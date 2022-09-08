#include <stdbool.h>
#include <assert.h>

#include <math.h>		// for ceil, pow, log2
#include <stdio.h>		// for fprintf, snprintf, fflush, stderr, FILE
#include <stdlib.h>		// for free, malloc, exit, abort, realloc
//#include <string.h>		// for basic_string
//#include <type_traits>         // for remove_reference<>::type
#include <gmp.h>
//#include <time.h>
#include <sys/time.h>


/* structure to compute on-line a product tree, avoiding to first compute a
   list of mpz_t (which might take too much memory) */
typedef struct {
  mpz_t *l;     /* the value stored is l[0] * l[1] * ... * l[size-1],
                   where l[0] is the product of n[0] elements, l[1] is
                   the product of n[1] elements, ..., with n[0]=0 or 1,
                   n[1]=0 or 2, ..., n[k]=0 or 2^k */
  unsigned long *n;
  size_t size;
} mpz_product_tree_t;
typedef mpz_product_tree_t mpz_product_tree[1];


static unsigned long tree_height(unsigned long n)
{
	unsigned long h = 0;

	while (n > 1) {
		h++;
		n = (n + 1) / 2;
	}
	return h;
}

#define MAX_N 1024*256*60

mpz_t ** init_product_tree(int n)
{
        int h = tree_height(n);
        mpz_t **T = malloc((h + 1) * sizeof(mpz_t *));
        
        /* initialize tree */
        for (int i = 0; i <= h; i++) {
                int w = 1 + ((n - 1) >> i);
                T[i] = malloc(w * sizeof(mpz_t));
                for (size_t j = 0; j < w; j++)
                        mpz_init(T[i][j]);
        }
        return T;
}

/* return the product tree formed from R[0..n-1].
 * Put in w[i] the number of elements of level i:
 * w[0] = n, w[1] = ceil(n/2), ...
 *
 * returns the height
 */
int product_tree(mpz_t **T, mpz_t * R, int n, size_t *w)
{
	int h = tree_height(n);

	/* initialize tree */
	for (int i = 0; i <= h; i++)
		w[i] = 1 + ((n - 1) >> i);

	/* initialize T[0] to R */
        // CB: il y a moyen d'éviter ceci, en faisant directement pointer
        //     T[0] sur R...
	for (int j = 0; j < n; j++)
		mpz_set(T[0][j], R[j]);

	/* compute product tree */
	for (int i = 1; i <= h; i++) {
		for (int j = 0; j < w[i - 1] / 2; j++)
			mpz_mul(T[i][j], T[i - 1][2 * j], T[i - 1][2 * j + 1]);
		if (w[i - 1] & 1)
			mpz_set(T[i][w[i] - 1], T[i - 1][w[i - 1] - 1]);
	}
        return h;
}

/* Auxiliary routine: a node T[i][j] with left son T[i-1][2*j] and right
   son T[i-1][2*j+1] corresponds to a node say t with left son u and right
   son v in the original product tree, where t = u*v.
   Now the current value of T[i][j], say t', is an approximation of
   (P*2^k)/t mod 2^k, where k = nbits(t) + g (g is the guard).
   This routine puts in T[i-1][2*j] and T[i-1][2*j+1] an approximation
   u' of (P*2^ku)/u mod 2^ku and v' of (P*2^kv)/v mod 2^kb respectively, where
   ku = nbits(u) + g, and kv = nbits(v) + g.
   Assume at input we have |t' - (P*2^k)/t| < x mod 2^k, where "mod 2^k"
   means that all quantities are reduced modulo 2^k with remainder in
   [-2^k/2, 2^k/2]. Then we first compute v*t':
   |v*t' - (P*2^k)/u| < v*x mod 2^k   [since t = u*v]
   then we divide by 2^(k-ku):
   |v*t'/2^(k-ku) - (P*2^ku)/u| < v*x/2^(k-ku) mod 2^ku
   thus since v = t/u < 2^k/2^(ku-1) and |v*t'/2^(k-ku) - u'| < 1,
   |u' - (P*2^ku)/u| < 2*x+1 mod 2^ku
   This proves that the error goes from x to 2x+1 at each step
   of the tree, thus since it is at most 1 at the root of the tree, it is
   at most 2^(h+1)-1 at the leaves. */
static void remainder_tree_aux(mpz_t ** T, unsigned long **nbits, unsigned long i, unsigned long j, unsigned long guard)
{
	/* T[i][j]/2^(nbits[i][j] + guard) ~ P/T[i][j] */
	mpz_mul(T[i - 1][2 * j], T[i][j], T[i - 1][2 * j]);
	/* same for the right part */
	mpz_mul(T[i - 1][2 * j + 1], T[i][j], T[i - 1][2 * j + 1]);
	/* swap */
	mpz_swap(T[i - 1][2 * j], T[i - 1][2 * j + 1]);

	/* get the fractional part, i.e., the low nbits[i][j] + guard bits */
	mpz_tdiv_r_2exp(T[i - 1][2 * j], T[i - 1][2 * j], nbits[i][j] + guard);
	/* now keep only nbits[i-1][2*j] + guard significant bits */
	mpz_div_2exp(T[i - 1][2 * j], T[i - 1][2 * j], nbits[i][j] - nbits[i - 1][2 * j]);

	mpz_tdiv_r_2exp(T[i - 1][2 * j + 1], T[i - 1][2 * j + 1], nbits[i][j] + guard);
	mpz_div_2exp(T[i - 1][2 * j + 1], T[i - 1][2 * j + 1], nbits[i][j] - nbits[i - 1][2 * j + 1]);
}

/* Compute the remainder tree using the "scaled" variant
 * (http://cr.yp.to/arith/scaledmod-20040820.pdf).
 * At the root, we compute a floating-point
 * approximation of P/T[h][0] with m+guard bits, where m = nbits(T[h][0]).
 *
 * In the end, the new value of T[0][i] is P % T[0][i]. 
 */
void remainder_tree(mpz_t ** T, size_t *w, mpz_t P, mpz_t * R, size_t n)
{

	unsigned long h = tree_height(n), i, j, guard;
	unsigned long **nbits;
	mpz_t Q;

	guard = h + 2;		/* see error analysis above */
	nbits = (unsigned long **)malloc((h + 1) * sizeof(unsigned long *));
	for (i = 0; i <= h; i++) {
		nbits[i] = (unsigned long *)malloc(w[i] * sizeof(unsigned long));
		for (j = 0; j < w[i]; j++)
			nbits[i][j] = mpz_sizeinbase(T[i][j], 2);
	}

	mpz_init(Q);
	mpz_mod(Q, P, T[h][0]);	/* first reduce modulo T[h][0] in case P is huge */
	mpz_mul_2exp(Q, Q, nbits[h][0] + guard);
	mpz_tdiv_q(T[h][0], Q, T[h][0]);
	/* |T' - 2^k*P/T| < 1 mod 2^k, where T is the original value of T[h][0],
	   T' is the new value of T[h][0], k = nbits(T) + guard, and "mod 2^k"
	   means that all values are taken modulo 2^k, with remainder in
	   [-2^k/2,2^k/2]. */
	for (i = h; i > 0; i--) {
		for (j = 0; j < w[i - 1] / 2; j++)
			remainder_tree_aux(T, nbits, i, j, guard);
		if (w[i - 1] & 1)
			mpz_swap(T[i - 1][w[i - 1] - 1], T[i][w[i] - 1]);
	}

	/* now for all leaves, if R = R[j], and R' = T[0][j],
	   we have |R' - 2^k*P/R| < 2^h mod 2^k, where k = nbits(R) + g,
	   thus R' = 2^k*P/R + a*2^k + b, with a integer, and |b| < 2^(h+1),
	   thus R*R' = 2^k*P + a*R*2^k + b*R
	   thus P mod R = R*R'/2^k - b*R/2^k, with |b*R/2^k| < 2^(h+1-g).
	   Now it suffices to have g >= h+2 so that the 2^(h+1-g) term
	   is less than 1/2, and rounding R*R'/2^k to the nearest integer
	   gives P mod R. */

	/* from T[0][j] ~ P/R[j]*2^(nbits[0][j] + guard) mod 2^(nbits[0][j] + guard),
	   get T[0][j]*R[j]/2^(nbits[0][j] + guard) ~ P mod R[j] */

	//ASSERT_ALWAYS(R.size() == w[0]);

	for (size_t j = 0; j < n; ++j) {
		mpz_mul(T[0][j], T[0][j], R[j]);
		/* T[0][j] ~ P*2^(nbits[0][j] + guard) mod R[j]*2^(nbits[0][j]+guard) */
		mpz_div_2exp(T[0][j], T[0][j], nbits[0][j]);
		/* T[0][j] ~ P*2^guard mod R[j]*2^guard */
		mpz_add_ui(T[0][j], T[0][j], (1UL << guard) - 1UL);
		mpz_div_2exp(T[0][j], T[0][j], guard);
	}

	mpz_clear(Q);

	for (i = 0; i <= h; i++)
		free(nbits[i]);
	free(nbits);
}

/* Compute the remainder tree using the direct variant.
 *
 * In the end, the new value of T[0][i] is P % T[0][i]. 
 */
void remainder_tree_simple(mpz_t ** T, size_t *w, mpz_t P, int h)
{
        mpz_mod(T[h][0], P, T[h][0]);   // à priori, c'est cette opération la plus "dure"
        for (int i = h - 1; i >= 0; i--)
                for (int j = 0; j < w[i]; j++)
                        mpz_mod(T[i][j], T[i + 1][j / 2], T[i][j]);
        
}

/* Clear the product tree. */
void clear_product_tree(mpz_t ** T, unsigned long n, size_t *w)
{
	unsigned long i, j, h = tree_height(n);

	for (i = 0; i <= h; i++) {
		for (j = 0; j < w[i]; j++)
			mpz_clear(T[i][j]);
		free(T[i]);
	}
	free(T);
}

#define MAX_DEPTH 32

#define MAX_NORM_DIGITS 100

struct mpz_array {
        mpz_t *A;
        int len;
        int capacity;
};

void array_setup(struct mpz_array *T)
{
        T->len = 0;
        T->capacity = 0;
        T->A = NULL;
}

void array_pushback(struct mpz_array *T, mpz_t x)
{
        if (T->len == T->capacity) {
                T->capacity = 2 * T->capacity + 1;
                T->A = realloc(T->A, T->capacity * sizeof(mpz_t));
        }
        mpz_init_set(T->A[T->len], x);
        T->len += 1;
}

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

const char *PRIMES_FILE = "primes.txt";
const char *NORMS_FILE = "norms.txt";
const char *OUTPUT_FILE = "norms.txt";

void product(mpz_t *A,mpz_t x,int begin,int end)
{
	if(begin+1==end){
		mpz_set(x,A[begin]);
		return;
	}
	product(A,x,begin,(begin+end)/2);
	mpz_t y;
	mpz_init(y);
	product(A,y,(begin+end)/2,end);
	mpz_mul(x,x,y);
	mpz_clear(y);
}

#define LIM (pow(2,31))

int main(int argc, char *argv[])
{
        /* load norms */
	FILE * norm_file = fopen(argv[1], "r");
	struct mpz_array norms;
        array_setup(&norms);
        int Nsize = 0;
	while (true) {
                char norm_str[MAX_NORM_DIGITS];
		size_t read = fscanf(norm_file, "%s\n", norm_str);
		if (read == 0 || read == EOF)
                        break;
                mpz_t n;
                mpz_init(n);
		mpz_set_str(n, norm_str, 10);
		array_pushback(&norms, n);
                Nsize += mpz_size(n);
                mpz_clear(n);
	}
	fclose(norm_file);
        fprintf(stderr, "loaded %d norms (%d limbs)\n", norms.len, Nsize);

        /* load primes */
	FILE * prime_file = fopen(argv[2], "r");
        struct mpz_array primes;
        array_setup(&primes);
	//printf("%d\n",atoi(argv[3]));
        float max_prime = LIM / atof(argv[4]);
		fprintf(stderr, "max_prime = %.0f\n", max_prime);

	while (true) {
                long int prime;
		size_t read = fscanf(prime_file, "%lu\n", &prime);
                if (read == 0 || read == EOF || prime > max_prime)
                        break;
	       
		mpz_t p;
		mpz_init(p);
		mpz_set_ui(p, prime);
                array_pushback(&primes, p);
                mpz_clear(p);
	}
	
	fclose(prime_file);
        fprintf(stderr, "loaded %d primes\n", primes.len);

	
        /* compute product of primes */
        mpz_t P;
	mpz_init(P);
	product(primes.A,P,0,primes.len);
        int Psize = mpz_sizeinbase(P,2);
        fprintf(stderr, "prime product has %zd bits (%ld limbs)\n", mpz_sizeinbase(P, 2), mpz_size(P));

        double product_time = 0;
        double remainder_time = 0;
        mpz_t **T = init_product_tree(MAX_N);

	double start = wtime();

        int done = 0;     // norms up to done have been processed
        int batches = 0;  // number of batches processed
        mpz_t smooth;
        mpz_init(smooth);
        size_t w[MAX_DEPTH];
        
	while (done < norms.len) {
                /* delimit next batch of norms */
                int n = 0;
                int size = 0;
                while (done + n < norms.len && size < Psize/2) {
                        size += mpz_sizeinbase(norms.A[done + n],2);
                        n += 1;
                } 
                assert(n < MAX_N);
                assert(n > 0);

                double product_start = wtime();

		int h = product_tree(T, norms.A + done, n, w);
                
		double remainder_start = wtime();
		// remainder_tree(T, w, P, norms.A + done, n);
                remainder_tree_simple(T, w, P, h);

                remainder_time += wtime() - remainder_start;
                product_time += remainder_start - product_start;

		for (int j = 0; j < n; j++) {
			/* Divide out R by gcd(P, R) as much as we can. The first gcd may
			 * have some cost, the subsequent ones are supposedly cheap
			 * enough */
		        mpz_set_ui(smooth, 1);
                	for (;;) {
				mpz_gcd(T[0][j], T[0][j], norms.A[done + j]);
				if (mpz_cmp_ui(T[0][j], 1) == 0)
					break;
				mpz_divexact(norms.A[done + j], norms.A[done + j], T[0][j]);
                                mpz_mul(smooth, smooth, T[0][j]);
			}
                        mpz_swap(norms.A[done + j], smooth);
                        //version alternative: (les perf ont l'air rigoureusement équivalentes...)
                        //mpz_powm_ui(T[0][j], T[0][j], 17, norms.A[done + j]);
                        //mpz_gcd(norms.A[done + j], T[0][j], norms.A[done + j]);
		}
		// clear_product_tree(T, n, w);
                done += n;
                batches += 1;
	}
	double stop = wtime();

	fprintf(stderr, "%d batches processed in %.2fs\n", batches, stop - start);
        fprintf(stderr, "product time: %.2fs\n", product_time);
        fprintf(stderr, "remainder time: %.2fs\n", remainder_time);
        fprintf(stderr, "extra time: %.2fs\n", stop - start - product_time - remainder_time);

        /* open output file */
        FILE * output_file = fopen(argv[3], "w");
        for (int i = 0; i < norms.len; i++) {
                mpz_out_str(output_file, 10, norms.A[i]);
                fprintf(output_file, "\n");
        }
        fclose(output_file);


	return 0;
}
