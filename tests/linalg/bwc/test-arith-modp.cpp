#include "cado.h"
#include <iostream>
#include <gmp.h>
#include "tests_common.h"
#include "arith-modp.hpp"


template<typename F, int summands, int cbound>
void do_tests(unsigned long iter)
{
    mp_size_t m = F::n*mp_bits_per_limb;
    // mp_size_t ell = F::extra*mp_bits_per_limb;

    typename F::elt x[summands];
    typename F::elt_ur y;
    typename F::elt r;
    typename F::elt p;

    mpz_t xz[summands];
    mpz_t yz;
    mpz_t rz;
    mpz_t pz;

    int coeffs[summands];

    mpz_init(pz);
    mpz_init(yz);
    mpz_init(rz);
    for(int i = 0 ; i < summands ; i++) {
        mpz_init(xz[i]);
    }

    for(unsigned long test = 0 ; test < iter ; test++) {
        for(int i = 0 ; i < summands ; i++) {
            coeffs[i] = gmp_urandomm_ui(state, 2 * cbound);
            coeffs[i] -= cbound - (coeffs[i] >= cbound);
        }
        /* Take a random modulus which is not a power of 2, and fits
         * within m bits */
        do {
            mpz_rrandomb(pz, state, m);
        } while (mpz_scan1(pz, 0) == mpz_sizeinbase(pz, 2) - 1);
        for(int i = 0 ; i < summands ; i++) {
            mpz_urandomm(xz[i], state, pz);
        }
        mpz_set_ui(yz, 0);

        MPN_SET_MPZ(p.x, F::n, pz);
        mpn_zero(y.x, F::n + F::extra);
        for(int i = 0 ; i < summands ; i++) {
            MPN_SET_MPZ(x[i].x, F::n, xz[i]);
        }

        typename F::preinv j;
        F::compute_preinv(j, p);

        /* Do the computation */
        for(int i = 0 ; i < summands ; i++) {
            if (coeffs[i] > 0) {
                mpz_addmul_ui(yz, xz[i], coeffs[i]);
            } else {
                mpz_submul_ui(yz, xz[i], -coeffs[i]);
            }
        }
        mpz_mod(yz, yz, pz);


        for(int i = 0 ; i < summands ; i++) {
            if (coeffs[i] > 0) {
                F::addmul(y, x[i], coeffs[i]);
            } else {
                F::submul(y, x[i], -coeffs[i]);
            }
        }
        F::reduce(r, y, p, j);

        MPZ_SET_MPN(rz, r.x, F::n);

        ASSERT_ALWAYS(mpz_cmp(yz, rz) == 0);

        if (tests_common_get_verbose()) {
            gmp_printf("ok [%d+%d] %Zd\n", F::n, F::extra, pz);
        }
    }
    printf("%lu tests ok [%d+%d]\n", iter, F::n, F::extra);

    mpz_clear(pz);
    mpz_clear(rz);
    mpz_clear(yz);
    for(int i = 0 ; i < summands ; i++) {
        mpz_clear(xz[i]);
    }
}

int main(int argc, const char * argv[])
{
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);

    unsigned long iter = 100;
    tests_common_get_iter(&iter);

#define SUMMANDS 200
#define CBOUND 1000
    do_tests< gfp<1>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<3>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<4>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<5>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<6>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<7>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<8>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<9>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<2, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<3, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<4, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<5, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<6, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<7, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<8, 2>, SUMMANDS, CBOUND>(iter);
    do_tests< gfp<9, 2>, SUMMANDS, CBOUND>(iter);
    tests_common_clear();
}

