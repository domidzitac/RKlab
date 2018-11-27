/* Match every k-character snippet of the query_doc document
 among a collection of documents doc1, doc2, ....
 
 ./rkmatch snippet_size query_doc doc1 [doc2...]
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <strings.h>
#include <assert.h>
#include <time.h>
#include "bloom.h"
#include <ctype.h>

enum algotype { SIMPLE = 0, RK, RKBATCH};

/* a large prime for RK hash (BIG_PRIME*256 does not overflow)*/
long long BIG_PRIME = 5003943032159437;

/* constants used for printing debug information */
const int PRINT_RK_HASH = 5;
const int PRINT_BLOOM_BITS = 160;

/* modulo addition */
long long
madd(long long a, long long b)
{
    return ((a+b)>BIG_PRIME?(a+b-BIG_PRIME):(a+b));
}

/* modulo substraction */
long long
mdel(long long a, long long b)
{
    return ((a>b)?(a-b):(a+BIG_PRIME-b));
}

/* modulo multiplication*/
long long
mmul(long long a, long long b)
{
    return ((a*b) % BIG_PRIME);
}

/* read the entire content of the file 'fname' into a
 character array allocated by this procedure.
 Upon return, *doc contains the address of the character array
 *doc_len contains the length of the array
 */
void
read_file(const char *fname, char **doc, int *doc_len)
{
    struct stat st;
    int fd;
    int n = 0;
    
    fd = open(fname, O_RDONLY);
    if (fd < 0) {
        perror("read_file: open ");
        exit(1);
    }
    
    if (fstat(fd, &st) != 0) {
        perror("read_file: fstat ");
        exit(1);
    }
    
    *doc = (char *)malloc(st.st_size);
    if (!(*doc)) {
        fprintf(stderr, " failed to allocate %d bytes. No memory\n", (int)st.st_size);
        exit(1);
    }
    
    n = read(fd, *doc, st.st_size);
    if (n < 0) {
        perror("read_file: read ");
        exit(1);
    }else if (n != st.st_size) {
        fprintf(stderr,"read_file: short read!\n");
        exit(1);
    }
    
    close(fd);
    *doc_len = n;
}


/* The normalize procedure examines a character array of size len
 in ONE PASS and does the following:
 1) turn all upper case letters into lower case ones
 2) turn any white-space character into a space character and,
 shrink any n>1 consecutive spaces into exactly 1 space only
 Hint: use C library function isspace()
 You must do the normalization IN PLACE so that when the procedure
 returns, the character array buf contains the normalized string and
 the return value is the length of the normalized string.
 */

int
normalize(char *buf,    /* The character array containing the string to be normalized*/
          int len            /* the size of the original character array */)
{
    int i;
    int j=0; //j is used for backtracking
    
    //white spaces at the beginning
    while (isspace(buf[i]))
        i++;
    int k=i;
    for (k=i; k < len; k++)
    {
        //if there is a whitespace, increase j
        if (isspace(buf[k]))
        {
            if (j > 0){
                if(!(isspace(buf[j-1])))
                {
                    buf[j++] = ' ';
                }
            }
        }
        //otherwise, if it's not a space, lowercase
        else
        {
            buf[j++] = tolower(buf[k]);
        }
    }
    
    if( j > 0){
        if (isspace(buf[j-1]))
            j--;
    }
    //NULL TERMINATOR
    buf[j] = '\0';
    return j;
}

/* check if a query string ps (of length k) appears
 in ts (of length n) as a substring
 If so, return 1. Else return 0
 You may want to use the library function strncmp
 */
int
simple_match(const char *ps,    /* the query string */
             int k,                     /* the length of the query string */
             const char *ts,    /* the document string (Y) */
             int n                        /* the length of the document Y */)
{
    
    int i;
    for (i=0; i < n+k-1; i++){
        if(strncmp(ts, ps, k)==0){
            return 1;
        }
        
        ts++;
    }
    return 0;
}

/* Check if a query string ps (of length k) appears
 in ts (of length n) as a substring using the rabin-karp algorithm
 If so, return 1. Else return 0
 In addition, print the first 'PRINT_RK_HASH' hash values of ts
 Example:
 $ ./rkmatch -t 1 -k 20 X Y
 605818861882592 812687061542252 1113263531943837 1168659952685767 4992125708617222
 0.01 matched: 1 out of 148
 */

//Helper function for the rabin karp algorithm that produces the initial hashing

long long initial_hashing(const char *s, int len) {
    
    int i;
    long long sum=0;
    
    for (i =0; i < len; i++)
    {
        sum = madd(s[i], mmul(sum, 256));
    }
    return sum;
}


int
rabin_karp_match(const char *ps,    /* the query string */
                 int k,                     /* the length of the query string */
                 const char *ts,    /* the document string (Y) */
                 int n                        /* the length of the document Y */ )
{
    
    long long hash1 = initial_hashing(ps,k); //hash value for query string
    long long hash2 = initial_hashing(ts,k); //hash value for document
    long long f = 1; //hash factor
    int i;
    int p = PRINT_RK_HASH;
    int newline=0;
    
    //if the document Y is shorter than the query, exit with 0
    if(n<k)
        return 0;
    
    for (i=0; i<k-1; i++){
        f= mmul(f,256);
    }
    
    
    for (i = 0; i < n - k + 1; i++){
        //prints
        if(i<p)
            printf("%lld ", hash2);
        if(i==p){
            printf("\n");
            newline = 1;
        }
        
        //check if there is a match, making sure it is an actual match => return true
        if(hash1==hash2 && strncmp((ts+i),ps, k)==0){
            if (newline == 0)
                printf("\n");
            return 1;
        }
        //compute the last part of the algorithm, move the hash
        hash2 = madd(mmul(mdel(hash2, mmul(*(ts + i), f)), 256), *(ts + k + i));
    }
    
    return 0;
}

/* Initialize the bitmap for the bloom filter using bloom_init().
 Insert all m/k RK hashes of qs into the bloom filter using bloom_add().
 Then, compute each of the n-k+1 RK hashes of ts and check if it's in the filter using bloom_query().
 Use the given procedure, hash_i(i, p), to compute the i-th bloom filter hash value for the RK value p.
 
 Return the number of matched chunks.
 Additionally, print out the first PRINT_BLOOM_BITS of the bloom filter using the given bloom_print
 after inserting m/k substrings from qs.
 */
int
rabin_karp_batchmatch(int bsz,        /* size of bitmap (in bits) to be used */
                      int k,          /* chunk length to be matched */
                      const char *qs, /* query docoument (X)*/
                      int m,          /* query document length */
                      const char *ts, /* to-be-matched document (Y) */
                      int n           /* to-be-matched document length*/)
{
    int count=0;
    int i;
    int j;
    long long f=1;
    long long hash1 = initial_hashing(qs,k); //hash value for query string
    long long hash2 = initial_hashing(ts,k); //hash value for document
    long long bloom_hash;
    
    //if the document Y is shorter than the query, exit with 0
    if(n<k)
        return 0;
    
    //Initialize the bitmap for the bloom filter using bloom_init().
    
    bloom_filter blFilter = bloom_init(bsz);
    
    for (i=0; i<k-1; i++){
        f= mmul(f,256);
    }
    
    //Insert all m/k RK hashes of qs into the bloom filter using bloom_add().
    
    for (i=0; i< m/k; i++){
        bloom_hash = initial_hashing(qs + i * k, k);
        bloom_add(blFilter, bloom_hash);
    }
    
    //Additionally, print out the first PRINT_BLOOM_BITS of the bloom filter using the given bloom_print
    //after inserting m/k substrings from qs.
    bloom_print(blFilter, PRINT_BLOOM_BITS);
    
    //Then, compute each of the n-k+1 RK hashes of ts and check if it's in the filter using bloom_query().
    
    for (i = 0 ; i < n - k + 1 ; i++)
    {
        if (bloom_query(blFilter, hash2))
        {
            //compare the stings on each part and increase the counter if there is a match
            for (j = 0 ; j < m / k ; j++)
            {
                
                if (strncmp(qs + j * k, ts + i , k)==0)
                {
                    count++;
                    break;
                }
            }
        }
        hash2 = madd(mmul(mdel(hash2, mmul(*(ts + i), f)), 256), *(ts + k + i)); //move the hash further
    }
    //Return the number of matched chunks.
    return count;
    
}

int
main(int argc, char **argv)
{
    int k = 100; /* default match size is 100*/
    int which_algo = SIMPLE; /* default match algorithm is simple */
    
    char *qdoc, *doc;
    int qdoc_len, doc_len;
    int i;
    int num_matched = 0;
    int to_be_matched;
    int c;
    
    /* Refuse to run on platform with a different size for long long*/
    assert(sizeof(long long) == 8);
    
    /*getopt is a C library function to parse command line options */
    while (( c = getopt(argc, argv, "t:k:q:")) != -1) {
        switch (c)
        {
            case 't':
                /*optarg is a global variable set by getopt()
                 it now points to the text following the '-t' */
                which_algo = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'q':
                BIG_PRIME = atoi(optarg);
                break;
            default:
                fprintf(stderr,
                        "Valid options are: -t <algo type> -k <match size> -q <prime modulus>\n");
                exit(1);
        }
    }
    
    /* optind is a global variable set by getopt()
     it now contains the index of the first argv-element
     that is not an option*/
    if (argc - optind < 1) {
        printf("Usage: ./rkmatch query_doc doc\n");
        exit(1);
    }
    
    /* argv[optind] contains the query_doc argument */
    read_file(argv[optind], &qdoc, &qdoc_len);
    qdoc_len = normalize(qdoc, qdoc_len);
    
    /* argv[optind+1] contains the doc argument */
    read_file(argv[optind+1], &doc, &doc_len);
    doc_len = normalize(doc, doc_len);
    
    switch (which_algo)
    {
        case SIMPLE:
            /* for each of the qdoc_len/k chunks of qdoc,
             check if it appears in doc as a substring*/
            for (i = 0; (i+k) <= qdoc_len; i += k) {
                if (simple_match(qdoc+i, k, doc, doc_len)) {
                    num_matched++;
                }
            }
            break;
        case RK:
            /* for each of the qdoc_len/k chunks of qdoc,
             check if it appears in doc as a substring using
             the rabin-karp substring matching algorithm */
            for (i = 0; (i+k) <= qdoc_len; i += k) {
                if (rabin_karp_match(qdoc+i, k, doc, doc_len)) {
                    num_matched++;
                }
            }
            break;
        case RKBATCH:
            /* match all qdoc_len/k chunks simultaneously (in batch) by using a bloom filter*/
            num_matched = rabin_karp_batchmatch(((qdoc_len*10/k)>>3)<<3, k, qdoc, qdoc_len, doc, doc_len);
            break;
        default :
            fprintf(stderr,"Wrong algorithm type, choose from 0 1 2\n");
            exit(1);
    }
    
    to_be_matched = qdoc_len / k;
    printf("%.2f matched: %d out of %d\n", (double)num_matched/to_be_matched,
           num_matched, to_be_matched);
    
    free(qdoc);
    free(doc);
    
    return 0;
}

