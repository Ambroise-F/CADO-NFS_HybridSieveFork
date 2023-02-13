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

#include "las-cofac-standalone.hpp"       // for cofac_standalone
#include "cxx_mpz.hpp"




#define MAX_DEPTH 32

#define MAX_NORM_DIGITS 100

struct mpz_array {
        mpz_t *A;
        int len;
        int capacity;
};


/* structure to compute on-line a product tree, avoiding to first compute a
   list of mpz_t (which might take too much memory) */
typedef struct {
  mpz_t *l;     /* the value stored is l[0] * l[1] * ... * l[size-1],
                   where l[0] is the product of n[0] elements, l[1] is
                   the product of n[1] elements, ..., with n[0]=0 or 1,
                   n[1]=0 or 2, ..., n[k]=0 or 2^k */
  unsigned long *n;
  size_t size;
} sm_mpz_product_tree_t;



typedef sm_mpz_product_tree_t sm_mpz_product_tree[1];

mpz_t ** sm_init_product_tree(int);

int sm_product_tree(mpz_t **, std::vector<cofac_standalone>, int , size_t *);

int sm_batch(std::vector<cofac_standalone>&, cxx_mpz, int);
int sm_batch_initalized_tree(mpz_t **, std::vector<cofac_standalone>&, cxx_mpz, int);


