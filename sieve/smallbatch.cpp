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

#include "smallbatch.hpp"







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

mpz_t ** sm_init_product_tree(int n)
{
        int h = tree_height(n);
        mpz_t **T = (mpz_t **)malloc((h + 1) * sizeof(mpz_t *));
        
        /* initialize tree */
        for (int i = 0; i <= h; i++) {
                int w = 1 + ((n - 1) >> i);
                T[i] = (mpz_t *)malloc(w * sizeof(mpz_t));
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
int sm_product_tree(mpz_t **T, std::vector<cofac_standalone> R, int n, size_t *w, int side)
{
	int h = tree_height(n);
	/* initialize tree */
	for (int i = 0; i <= h; i++)
		w[i] = 1 + ((n - 1) >> i);



	/* initialize T[0] to R */
        // CB: il y a moyen d'éviter ceci, en faisant directement pointer
        //     T[0] sur R...
	//for (int j = 0; j < n; j++)
	for (int j = 0; j < n; j++)
		mpz_set(T[0][j], R[j].norm[side].x);
	
// gmp_fprintf(stderr, "norm[%d] = %Zd (n=%d)\n",j ,R[j].norm[side].x, n);

	/* compute product tree */
	for (int i = 1; i <= h; i++) {
		for (int j = 0; j < w[i - 1] / 2; j++)
			mpz_mul(T[i][j], T[i - 1][2 * j], T[i - 1][2 * j + 1]);
		if (w[i - 1] & 1)
			mpz_set(T[i][w[i] - 1], T[i - 1][w[i - 1] - 1]);
	}
        return h;
}



/* Compute the remainder tree using the direct variant.
 *
 * In the end, the new value of T[0][i] is P % T[0][i]. 
 */
void remainder_tree_simple(mpz_t ** T, size_t *w, mpz_t P, int h, int side)
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


void clear_product_tree2(mpz_t ** T, unsigned long n, size_t *w)
{
	int h = tree_height(MAX_N);
	for (int i = 0; i <= h; i++) {
                int w = 1 + ((n - 1) >> i);
                
                for (size_t j = 0; j < w; j++)
                        mpz_clear(T[i][j]);
                free(T[i]);
        }
        free(T);
}

void clear_product_tree3(mpz_t ** T)
{
	int h = tree_height(MAX_N);
	for (int i = 0; i <= h; i++) {
                int w = 1 + ((MAX_N - 1) >> i);

                for (size_t j = 0; j < w; j++)
                        mpz_clear(T[i][j]);
                free(T[i]);
        }
        free(T);
}



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
                T->A = (mpz_t *)realloc(T->A, T->capacity * sizeof(mpz_t));
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

void product(mpz_t *A, mpz_t x, int begin, int end)
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



//- instead of a main
int sm_batch(std::vector<cofac_standalone> &surv, cxx_mpz fb_product, int side){


	// todo : change return (to actual facto)
	// todo : make it so that facto tree is not computed _every fucking time_

	/* compute product of primes */

        int Psize = mpz_sizeinbase(fb_product,2);
        fprintf(stderr, "prime product has %zd bits (%ld limbs)\n", mpz_sizeinbase(fb_product, 2), mpz_size(fb_product));

        double product_time = 0;
        double remainder_time = 0;

        
        mpz_t **T = sm_init_product_tree(MAX_N);
        //mpz_t **T2= sm_init_product_tree(MAX_N);

	double start = wtime();

        int done = 0;     // norms up to done have been processed
        int batches = 0;  // number of batches processed
        mpz_t smooth;   // smooth part of norm        
        mpz_init(smooth);
        size_t w[MAX_DEPTH];
        

        int n = surv.size();
        fprintf(stderr, "size is %d\n", n);
        //int size = mpz_sizeinbase(toutes les normes)
        //                n += 1; // ?????

        assert(n < MAX_N); //- n ???? todoo
        assert(n > 0);

        double product_start = wtime();

	int h = sm_product_tree(T, surv, n, w, side); // result is T // don't open seg fault inside

	double remainder_start = wtime();
	// remainder_tree(T, w, fb_product, norms.A + done, n);
        remainder_tree_simple(T, w, fb_product, h, side); // result is also T

        remainder_time += wtime() - remainder_start;
        product_time += remainder_start - product_start;

	for (int j = 0; j < n; j++) { // le résultat va dans smooth, smooth est réinjecté dans T -> ok cool mais il faut mettre les facteurs dans  (?)
		/* Divide out R by gcd(fb_product, R) as much as we can. The first gcd may
		 * have some cost, the subsequent ones are supposedly cheap
		 * enough */
	        mpz_set_ui(smooth, 1);
	        //-gmp_fprintf(stderr, "[PRE  BATCH] norm[%d] = %Zd           (side = %d)\n", j, surv[j].norm[side].x, side);
        	for (;;) {
        		//if (mpz_cmp_ui(T[0][j], 1) == 0)
			//	break;
			mpz_gcd(T[0][j], T[0][j], surv[j].norm[side]);
			if (mpz_cmp_ui(T[0][j], 1) == 0)
				break;
			mpz_divexact(surv[j].norm[side], surv[j].norm[side], T[0][j]);
                        mpz_mul(smooth, smooth, T[0][j]);
		}
                // mpz_swap(surv[j].norm[side], smooth); //- no swap, remainder/non-smooth part should be in surv[j].norm[side], instead, smooth part should go in surv[j].factors or smthg?
		
		mpz_init_set(surv[j].sm_smoothpart[side], smooth);

		//-gmp_fprintf(stderr, "[POST BATCH] norm[%d] = %Zd\n", j, surv[j].norm[side]);
		//-gmp_fprintf(stderr, "[POST BATCH] sm_s[%d] = %Zd\n", j, surv[j].sm_smoothpart[side]);
		//-fprintf(stderr, "-------------------\n");


                //version alternative: (les perf ont l'air rigoureusement équivalentes...)
                //mpz_powm_ui(T[0][j], T[0][j], 17, norms.A[done + j]);
                //mpz_gcd(norms.A[done + j], T[0][j], norms.A[done + j]);
	}
	//clear_product_tree2(T, n, w);
	clear_product_tree3(T);



	//clear_product_tree(T2, n, w);

        done += n;
        batches += 1;

        fprintf(stderr, "-- batch_fac done\n");

        mpz_clear(smooth);

        //delete &fb_product; //no
	
	double stop = wtime();

	fprintf(stderr, "%d batches processed in %.2fs\n", batches, stop - start);
        fprintf(stderr, "product time: %.2fs\n", product_time);
        fprintf(stderr, "remainder time: %.2fs\n", remainder_time);
        fprintf(stderr, "extra time: %.2fs\n", stop - start - product_time - remainder_time);


	return 0;
}


#if 0
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
        mpz_t **T = sm_init_product_tree(MAX_N);

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

		int h = sm_product_tree(T, norms.A + done, n, w);
                
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
#endif