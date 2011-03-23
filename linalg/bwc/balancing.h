#ifndef BALANCING_H_
#define BALANCING_H_

#include "utils.h"

// #include "mf.h"

/* This is used only for sanity checking after dispatch. */
/* The primes are 2^64-59 and 2^63-25 */
#define DUMMY_VECTOR_COORD_VALUE(j)     \
                (UINT64_C(0xffffffffffffffc5) / (1 + (j))  \
              + UINT64_C(0x7fffffffffffffe7) * (1 + (j)))
#define DUMMY_VECTOR_COORD_VALUE2(j)     \
        (DUMMY_VECTOR_COORD_VALUE((j) ^ 0xbeef) ^ ((DUMMY_VECTOR_COORD_VALUE((j) ^ 0xcafe) << 32)))

/* balancing structure */

#define FLAG_COLPERM    1
#define FLAG_ROWPERM    2
#define FLAG_PADDING    4       /* pad to largest dimension */
#define FLAG_REPLICATE  8       /* only balancing in one dimension */
#define FLAG_SHUFFLED_MUL  16   /* corresponds to using shuffled product */

/* Shuffled (a.k.a. reordered) product is when matmul_top code dispatches data
 * with allgather + reduce-scatter instead of bcast + reduce, for a gain in
 * efficiency. This amounts to considering the matrix Sr*M*Sc^-1*P (for left
 * multiplication) instead of Sr*M*Sc^-1. Therefore the same permutation pair
 * cannot be used compatibly in the shuffled and non-shuffled case.
 */

struct balancing_header_s {
    // FIXME: add a magic number here ? This header is read directly in
    // binary format, so it might be a good idea.
    uint32_t nh;
    uint32_t nv;
    uint32_t nrows;
    uint32_t ncols;
    uint64_t ncoeffs;
    uint32_t checksum;
    uint32_t flags;
    // rshuf indicates two integers a,b such that the column i of the input
    // matrix is in fact mapped to column a*i+b mod n in the matrix we work
    // with. rshuf_inv indicates the inverse permutation. a and b do
    // **NOT** depend on nh and nv. Therefore all balancing are
    // compatible even though a permutation is being used.
    // a == rhuf[0] is the smallest integer greater than floor(sqrt(n)) which
    // is coprime to n. b is 42. rshuf_inv[0,1] are computed accordingly.
    uint32_t pshuf[2];
    uint32_t pshuf_inv[2];
};
typedef struct balancing_header_s balancing_header[1];

struct balancing_s {
    balancing_header h;
    uint32_t trows;     // target number of rows. ==nrows if no padding
    uint32_t tcols;     // target number of cols. ==ncols if no padding
    uint32_t * rowperm; // row index for new mat. --> row index for old mat.
    uint32_t * colperm; // might be equal to colperm.
};
typedef struct balancing_s balancing[1];
typedef struct balancing_s * balancing_ptr;

#ifdef __cplusplus
extern "C" {
#endif

/* Once the flags and perm[] fields have been provided, the caller must
 * call _finalize() in order to 1) update the trows and tcols fields 2)
 * compute the checksum of the balancing.
 */
extern void balancing_set_row_col_count(balancing_ptr bal);
extern void balancing_finalize(balancing_ptr bal);
extern void balancing_write_inner(balancing_ptr bal, const char *);
extern void balancing_write(balancing_ptr bal, const char * , const char *);
extern void balancing_read(balancing_ptr bal, const char *);
extern void balancing_read_header(balancing_ptr bal, const char * filename);
extern void balancing_clear(balancing_ptr bal);
extern void balancing_init(balancing_ptr bal);

/* helper for the functions below */
static inline unsigned long balancing_row_shuffle_common_(unsigned long r, unsigned long n, uint32_t * shuf)
{
    modulusul_t M;
    modul_initmod_ul(M, n);
    residueul_t x,a,b;
    modul_init_noset0(a, M);
    modul_init_noset0(b, M);
    modul_init_noset0(x, M);
    modul_set_ul(a, shuf[0], M);
    modul_set_ul(b, shuf[1], M);
    modul_set_ul(x, r, M);
    modul_mul(x, x, a, M);
    modul_add(x, x, b, M);
    unsigned long r1 = modul_get_ul(x, M);
    modul_clear(a, M);
    modul_clear(b, M);
    modul_clear(x, M);
    modul_clearmod(M);
    return r1;
}

/* These two relate to the global permutation represented by the rshuf /
 * rshuf_inv arrays */
static inline unsigned long balancing_pre_shuffle(balancing_ptr bal, unsigned long r)
{
    return balancing_row_shuffle_common_(r, bal->h->ncols, bal->h->pshuf);
}
static inline unsigned long balancing_pre_unshuffle(balancing_ptr bal, unsigned long r)
{
    return balancing_row_shuffle_common_(r, bal->h->ncols, bal->h->pshuf_inv);
}


static inline int balancing_progressive_dispatch_block(int n, int t, int k)
{
    int f = n / t;
    int r = n % t;
    int c = f + 1;
    if (k < r * c) {
        return k / c;
    } else {
        k -= r * c;
        return r + k / f;
    }
}
#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_H_ */
