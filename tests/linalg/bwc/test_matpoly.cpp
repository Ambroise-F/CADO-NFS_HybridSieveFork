#include "cado.h"
#include <gmp.h>
#include "lingen_matpoly_ft.hpp"
#include "utils.h"
#include "gmp_aux.h"

struct matpoly_checker_base {
    abfield ab;
    unsigned int m;
    unsigned int n;
    unsigned int len1;
    unsigned int len2;
    matpoly::memory_guard dummy;

    gmp_randstate_t rstate;
    /* tests are free to seed and re-seed the checker's private rstate, a
     * priori with this random seed which was taken from the initial
     * random state */
    unsigned long seed;

    matpoly_checker_base(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int len1, unsigned int len2, gmp_randstate_t rstate0)
        : m(m)
        , n(n)
        , len1(len1)
        , len2(len2)
        , dummy(SIZE_MAX)
        , seed(gmp_urandomm_ui(rstate0, ULONG_MAX))
    {
        abfield_init(ab);
        abfield_specify(ab, MPFQ_PRIME_MPZ, (mpz_srcptr) p);
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, seed);
    }
    matpoly_checker_base(matpoly_checker_base const & o)
        : m(o.m)
        , n(o.n)
        , len1(o.len1)
        , len2(o.len2)
        , dummy(SIZE_MAX)
        , seed(o.seed)
    {
        cxx_mpz p;
        abfield_init(ab);
        abfield_specify(ab, MPFQ_PRIME_MPZ, abfield_characteristic_srcptr(o.ab));
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, seed);
    }
    ~matpoly_checker_base() {
        gmp_randclear(rstate);
        abfield_clear(ab);
    }

    int ctor_and_pre_init() {
        matpoly A(ab, 0, 0, 0);
        if (!A.check_pre_init()) return 0;
        matpoly P;
        matpoly Q(ab, m, m+n, len1);
        Q.fill_random(len1, rstate);
        return P.check_pre_init() && !Q.check_pre_init();
    }

    int move_ctor() {
        /* when moving P to Q, we should de-init P. In fact, it's not
         * guaranteed, since we rely on the behaviour of the default move
         * ctor.  On the other hand, I'm making this assumption quite
         * often in the code, and it's good if I have an occasion to
         * check that it holds
         */
        matpoly P(ab, m, m+n, len1);
        P.fill_random(len1, rstate);
        matpoly Q(std::move(P));
        matpoly R(ab, m, n, len1);
        R.fill_random(len1, rstate);
        R = matpoly();
        return P.check_pre_init() && !Q.check_pre_init() && R.check_pre_init();
    }

    int copy_ctor() {
        matpoly P(ab, m, n, len1);
        P.fill_random(len1, rstate);
        matpoly Q;
        Q.set(P);
        return P.cmp(Q) == 0;
    }

    int fill_random_is_deterministic() {
        matpoly P0(ab, n, n, len1);
        matpoly P1(ab, n, n, len1 + len2);
        gmp_randseed_ui(rstate, seed); P0.fill_random(len1, rstate);
        gmp_randseed_ui(rstate, seed); P1.fill_random(len1, rstate);
        int ok = P0.capacity() >= len1 && P1.capacity() >= len1+len2 && P0.cmp(P1) == 0;
        return ok;
    }

    int realloc_does_what_it_says() {
        /* begin like the previous test. In particular, we  */
        matpoly P0(ab, n, n, len1);
        matpoly P1(ab, n, n, len1 + len2);
        gmp_randseed_ui(rstate, seed); P0.fill_random(len1, rstate);
        gmp_randseed_ui(rstate, seed); P1.fill_random(len1, rstate);
        int ok;
        /* If data is shrunk at or above the previous value of 'size',
         * then old data is kept, and the 'size' field is unchanged.
         */
        P1.realloc(len1);
        ok = P0.cmp(P1) == 0;
        /* If data is grown, old data is kept, and the 'size' field is
         * unchanged.  */
        P1.realloc(len1 + len2);
        ok = ok && P0.cmp(P1) == 0;
        /* If data is shrunk below the previous value of 'size', then
         * 'size' is set to zero.  */
        P1.realloc(len1-1);
        ok = ok && P1.size == 0;
        /* Note: The content of the data area above 'size' on return is
         * unspecified.
         */
        return ok;
    }

    int mulx_then_divx() {
        matpoly P(ab, m,   n, len1 + n);
        P.fill_random(len1, rstate);
        matpoly Q;
        Q.set(P);
        /* take some columns, do multiplies */
        std::vector<int> jlen(n, len1);
        gmp_randseed_ui(rstate, seed);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = gmp_urandomm_ui(rstate, n);
            P.multiply_column_by_x(j, jlen[j]++);
        }
        /* Arrange so that we pick the same list, and divide */
        gmp_randseed_ui(rstate, seed);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = gmp_urandomm_ui(rstate, n);
            P.divide_column_by_x(j, jlen[j]--);
        }
        return P.cmp(Q) == 0;
    }

    int truncate_is_like_mulx_then_divx_everywhere() {
        matpoly P(ab, m,   n, len1);
        unsigned int trmax = std::min(128u, len1 / 2);
        unsigned int tr = gmp_urandomm_ui(rstate, trmax + 1);
        P.fill_random(len1, rstate);
        matpoly Q;
        Q.set(P);
        for(unsigned int s = 0 ; s < tr ; s++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                P.multiply_column_by_x(k, len1 - 1);
            }
        }
        for(unsigned int s = 0 ; s < tr ; s++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                P.divide_column_by_x(k, len1);
            }
        }
        /* shouldn't we have an api call to chop off zero coefficients at
         * high degrees ? */
        if (!P.tail_is_zero(len1 - tr)) return false;
        P.truncate(len1 - tr);
        matpoly R;
        /* check truncate-self as well as truncate-foreign */
        R.truncate(Q, len1 - tr);
        Q.truncate(Q, len1 - tr);
        if (Q.cmp(R) != 0) return false;
        if (P.cmp(Q) != 0) return false;
        return true;
    }

    int rshift_is_like_divx_everywhere() {
        matpoly P(ab, m,   n, len1);
        unsigned int trmax = std::min(128u, len1 / 2);
        unsigned int tr = gmp_urandomm_ui(rstate, trmax + 1);
        P.fill_random(len1, rstate);
        matpoly Q;
        Q.set(P);
        for(unsigned int s = 0 ; s < tr ; s++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                P.divide_column_by_x(k, len1);
            }
        }
        /* shouldn't we have an api call to chop off zero coefficients at
         * high degrees ? */
        if (!P.tail_is_zero(len1 - tr)) return false;
        P.truncate(len1 - tr);
        matpoly R;
        /* check rshift-self as well as rshift-foreign */
        R.rshift(Q, tr);
        Q.rshift(Q, tr);
        if (Q.cmp(R) != 0) return false;
        if (P.cmp(Q) != 0) return false;
        return true;
    }

    int test_extract_column() {
        unsigned int s = n;
        matpoly P(ab, m,   n, s+1);
        P.fill_random(1, rstate);
        for(unsigned int k = 0 ; k < s ; k++)
            for(unsigned int j = 0 ; j < s ; j++)
                P.extract_column((j+1)%s, k+1, P, j, k);
        /* We've cycled the columns exactly s times, so the head matrix
         * should be equal to the very first one.
         */
        matpoly Q;
        Q.set(P);
        for(unsigned int k = 0 ; k < s ; k++)
            for(unsigned int j = 0 ; j < s ; j++)
                Q.divide_column_by_x(j, s+1-k);
        Q.truncate(Q, 1);
        P.truncate(P, 1);
        return P.cmp(Q) == 0;
    }

    int divx_then_mulx_is_like_zero_column() {
        /* This is a bit like doing the mulx_then_divx test, but in
         * reverse order */
        matpoly P(ab, m,   n, len1 + n);
        P.fill_random(len1, rstate);
        matpoly Q;
        Q.set(P);
        /* take some columns, divide */
        std::vector<int> jlen(n, len1);
        gmp_randseed_ui(rstate, seed);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = gmp_urandomm_ui(rstate, n);
            P.divide_column_by_x(j, jlen[j]);
            Q.zero_column(j, len1-jlen[j]);
            jlen[j]--;
        }
        /* Arrange so that we pick the same list, and divide */
        gmp_randseed_ui(rstate, seed);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = gmp_urandomm_ui(rstate, n);
            P.multiply_column_by_x(j, jlen[j]++);
        }
        return P.cmp(Q) == 0;
    }

    int add_and_sub() {
        matpoly P(ab, m,   n, len1);
        matpoly Q(ab, m,   n, len2);
        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);
        matpoly R;
        R.add(P, Q); P.add(Q);
        if (R.cmp(P) != 0) return 0;

        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);
        R.add(Q, P); P.add(Q, P);
        if (R.cmp(P) != 0) return 0;
        
        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);
        R.sub(P, Q); P.sub(Q);
        if (R.cmp(P) != 0) return 0;

        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);
        R.sub(Q, P); P.sub(Q, P);
        if (R.cmp(P) != 0) return 0;
        
        /* Also check that P + Q - Q == P, whether P is smaller or larger
         * than Q (hence we do Q + P - P == Q as well)
         */
        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);
        R.add(P, Q);
        R.sub(Q);
        R.truncate(P.size);
        if (R.cmp(P) != 0) return 0;

        R.add(Q, P);
        R.sub(P);
        R.truncate(Q.size);
        if (R.cmp(Q) != 0) return 0;

        return 1;
    }

    int mul_is_distributive()
    {
        /* Complexity of the other operations is usually m*n*(len1+len2)
         * at most. Here we'll have m*n*n*mlen1*mlen2, so let's arrange
         * so that n*mlen1*mlen2 is approximately the same as len1+len2.
         */
        unsigned int mlen1 = len1;
        unsigned int mlen2 = len2;
        for( ; mlen1 >= 4 && mlen2 >= 4 && n*mlen1*mlen2 >= len1+len2 ; ) {
            mlen1 /= 2;
            mlen2 /= 2;
        }
        matpoly P(ab, m, n, mlen1);
        matpoly Q(ab, m, n, mlen2);
        matpoly R(ab, n, n, mlen2);
        matpoly PQ, PR, QR, PQ_R, PR_QR;
        P.fill_random(mlen1, rstate);
        Q.fill_random(mlen2, rstate);
        PQ.add(P, Q);
        PR.mul(P, R);
        QR.mul(Q, R);
        PQ_R.mul(PQ, R);
        PR_QR.add(PR, QR);
        if (PQ_R.cmp(PR_QR) != 0) return 0;
        matpoly testz(ab, m, n, 0);
        testz.sub(PR_QR);
        testz.addmul(PQ, R);
        if (!testz.tail_is_zero(0)) return 0;
        return 1;
    }
    /* This does not really perform a check that we're doing a middle
     * product, of course.
     */
    int mp_is_distributive()
    {
        unsigned int mlen1 = len1;
        unsigned int mlen2 = len2;
        for( ; mlen1 >= 4 && mlen2 >= 4 && n*mlen1*mlen2 >= len1+len2 ; ) {
            mlen1 /= 2;
            mlen2 /= 2;
        }
        matpoly P(ab, m, n, mlen1);
        matpoly Q(ab, m, n, mlen2);
        matpoly R(ab, n, n, mlen2);
        matpoly PQ, PR, QR, PQ_R, PR_QR;
        P.fill_random(mlen1, rstate);
        Q.fill_random(mlen2, rstate);
        PQ.add(P, Q);
        PR.mp(P, R);
        QR.mp(Q, R);
        PQ_R.mp(PQ, R);
        PR_QR.add(PR, QR);
        if (PQ_R.cmp(PR_QR) != 0) return 0;
        matpoly testz(ab, m, n, 0);
        testz.sub(PR_QR);
        testz.addmp(PQ, R);
        if (!testz.tail_is_zero(0)) return 0;
        return 1;
    }
    int coeff_is_zero_and_zero_column_agree()
    {
        matpoly P(ab, m,   n, len1 + 2);
        P.fill_random(len1, rstate);
        unsigned int k = P.size / 2;
        for(unsigned int j = 0 ; j < n ; j++)
            P.zero_column(j, k);
        return P.coeff_is_zero(k);
    }
};

template<typename fft_type>
struct matpoly_checker_ft : public matpoly_checker_base {
    typename matpoly_ft<fft_type>::memory_guard dummy_ft;

    matpoly_checker_ft(matpoly_checker_base const & base)
        : matpoly_checker_base(base)
        , dummy_ft(SIZE_MAX)
    {}
    matpoly_checker_ft(cxx_mpz const & p, unsigned int m, unsigned int n, unsigned int len1, unsigned int len2, gmp_randstate_t rstate0)
        : matpoly_checker_base(p, m, n, len1, len2, rstate0)
        , dummy_ft(SIZE_MAX)
    {}
    int mul_and_mul_caching_are_consistent() {
        matpoly P(ab, n, n, len1);
        matpoly Q(ab, n, n, len2);
        matpoly R0;
        matpoly R1;

        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);

        R0.mul(P, Q);
        matpoly_ft<fft_type>::mul_caching(R1, P, Q, NULL);
        return (R0.cmp(R1) == 0);
    }

    int mp_and_mp_caching_are_consistent() {
        matpoly P(ab, m,   n, len1);
        matpoly Q(ab, n, n, len2);
        matpoly M0;
        matpoly M1;

        P.fill_random(len1, rstate);
        Q.fill_random(len2, rstate);

        M0.mp(P, Q);
        matpoly_ft<fft_type>::mp_caching(M1, P, Q, NULL);
        return M0.cmp(M1) == 0;
    }

};

void declare_usage(cxx_param_list & pl)
{
#ifndef SELECT_MPFQ_LAYER_u64k1
    param_list_decl_usage(pl, "prime", "(mandatory) prime defining the base field");
#else
    param_list_decl_usage(pl, "prime", "(unused) prime defining the base field -- we only use 2");
#endif
    param_list_decl_usage(pl, "m", "dimension m");
    param_list_decl_usage(pl, "n", "dimension n");
    param_list_decl_usage(pl, "len1", "length 1");
    param_list_decl_usage(pl, "len2", "length 2");
    param_list_decl_usage(pl, "seed", "random seed");
}

int main(int argc, char * argv[])
{
    cxx_mpz p;
    gmp_randstate_t rstate;

    unsigned int m = 4;
    unsigned int n = 2;
    unsigned int len1 = 1000;
    unsigned int len2 = 600;
    unsigned long seed = 0;

    cxx_param_list pl;

    const char * argv0 = argv[0];
    argv++,argc--;
    /* read all command-line parameters */
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unexpected argument %s\n", argv[0]);
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
#ifndef SELECT_MPFQ_LAYER_u64k1
    if (!param_list_parse_mpz(pl, "prime", (mpz_ptr) p)) {
        fprintf(stderr, "--prime is mandatory\n");
        param_list_print_command_line (stdout, pl);
        exit(EXIT_FAILURE);
    }
#else
    mpz_set_ui(p, 2);   /* unused anyway */
    param_list_parse_mpz(pl, "prime", (mpz_ptr) p);
#endif
    param_list_parse_uint(pl, "m", &m);
    param_list_parse_uint(pl, "n", &n);
    param_list_parse_uint(pl, "len1", &len1);
    param_list_parse_uint(pl, "len2", &len2);
    param_list_parse_ulong(pl, "seed", &seed);
    if (param_list_warn_unused(pl))
        exit(EXIT_FAILURE);
#ifdef SELECT_MPFQ_LAYER_u64k1
    if (m & 63) {
        unsigned int nm = 64 * iceildiv(m, 64);
        printf("Round m=%u to m=%u\n", m, nm);
        m = nm;
    }
    if (n & 63) {
        unsigned int nn = 64 * iceildiv(n, 64);
        printf("Round n=%u to n=%u\n", n, nn);
        n = nn;
    }
#endif

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);

    matpoly_checker_base checker(p, m, n, len1, len2, rstate);

    ASSERT_ALWAYS(checker.ctor_and_pre_init());
    ASSERT_ALWAYS(checker.move_ctor());
    ASSERT_ALWAYS(checker.copy_ctor());
    ASSERT_ALWAYS(checker.fill_random_is_deterministic());
    ASSERT_ALWAYS(checker.realloc_does_what_it_says());
    ASSERT_ALWAYS(checker.mulx_then_divx());
    ASSERT_ALWAYS(checker.truncate_is_like_mulx_then_divx_everywhere());
    ASSERT_ALWAYS(checker.rshift_is_like_divx_everywhere());
    ASSERT_ALWAYS(checker.test_extract_column());
    ASSERT_ALWAYS(checker.divx_then_mulx_is_like_zero_column());
    ASSERT_ALWAYS(checker.add_and_sub());
    ASSERT_ALWAYS(checker.mul_is_distributive());
    ASSERT_ALWAYS(checker.mp_is_distributive());
    ASSERT_ALWAYS(checker.coeff_is_zero_and_zero_column_agree());

#ifdef SELECT_MPFQ_LAYER_u64k1
    {
        matpoly_checker_ft<gf2x_fake_fft_info> checker_ft(checker);
        ASSERT_ALWAYS(checker_ft.mul_and_mul_caching_are_consistent());
        ASSERT_ALWAYS(checker_ft.mp_and_mp_caching_are_consistent());
    }
    {
        matpoly_checker_ft<gf2x_cantor_fft_info> checker_ft(checker);
        ASSERT_ALWAYS(checker_ft.mul_and_mul_caching_are_consistent());
        ASSERT_ALWAYS(checker_ft.mp_and_mp_caching_are_consistent());
    }
    {
        matpoly_checker_ft<gf2x_ternary_fft_info> checker_ft(checker);
        ASSERT_ALWAYS(checker_ft.mul_and_mul_caching_are_consistent());
        ASSERT_ALWAYS(checker_ft.mp_and_mp_caching_are_consistent());
    }
#else
    {
        matpoly_checker_ft<fft_transform_info> checker_ft(checker);
        ASSERT_ALWAYS(checker_ft.mul_and_mul_caching_are_consistent());
        ASSERT_ALWAYS(checker_ft.mp_and_mp_caching_are_consistent());
    }
#endif

    gmp_randclear(rstate);
}
