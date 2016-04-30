#if !(defined(__OpenBSD__) || defined(__FreeBSD__))
#ifndef __STDC_FORMAT_MACROS
#error "Please define __STDC_FORMAT_MACROS before including lingen_mattypes.h"
#endif
#endif

#ifndef LINGEN_MAT_TYPES_HPP_
#define LINGEN_MAT_TYPES_HPP_

#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include <algorithm>
#ifdef  HAVE_OPENMP
#include <omp.h>
#define OMP_ROUND(k) ((k) % omp_get_num_threads() == omp_get_thread_num())
#else
#define OMP_ROUND(k) (1)
#endif

#include "bwc_config.h"
#include "alloc_proxy.h"
#include "utils.h"

/* Number of words holding B bits ; better naming sought. */
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))

/* Starting with gcc 4.3, -Wempty-body moans for loops like
 * for(;(x=x->next)!=NULL;y++);
 * It must shut up.
 *
 * Unfortunately the apple-bastardized gcc-4.2.1 backported this feature,
 * which gives rise to spurious warnings in that case as well. Since
 * there is no way to tell whether this pragma will be recognized or not,
 * we accept the change...
 */
// #if defined(__cplusplus) && GNUC_VERSION_ATLEAST(4,3,0)
// #pragma GCC diagnostic ignored "-Wempty-body"
// #endif

/* See the discussion in lingen_binary about the pros and cons of data
 * ordering schemes */

template<typename POD>/*{{{*/
inline void podswap(POD& a, POD& b) { POD c; c = a; a = b; b = c; } /*}}}*/

struct bcol;
struct bmat;
struct polmat;

struct bcol {/*{{{*/
    friend struct bmat;
    unsigned int nrows;
    inline unsigned int stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    private:
    unsigned long * x;
    public:
    bcol(unsigned int nrows)
        : nrows(nrows), x(mynew<unsigned long>(stride()))
    {}
    bcol() { nrows = 0; x = NULL; }
    void clear() { mydelete(x, stride()); nrows = 0; }
    private:
    bcol(bcol const& a) : nrows(a.nrows), x(a.x) {}
    bcol& operator=(bcol const& a) {
        nrows = a.nrows;
        x = a.x;
        return *this;
    }
    public:
    ~bcol() { mydelete(x, stride()); }
    void swap(bcol& a) {
        podswap(nrows, a.nrows);
        podswap(x, a.x);
    }
    void zero() {
        memset(x, 0, stride() * sizeof(unsigned long));
    }
    bool is_zero_likely() const {
        for(unsigned int i = 0 ; i < stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline unsigned long coeff(unsigned int i) const
    {
        ASSERT(i < nrows);
        unsigned int offset = i / ULONG_BITS;
        unsigned long shift = i % ULONG_BITS;
        return (x[offset] >> shift) & 1UL;
    }
    unsigned int ffs() const {
        unsigned int z = 0;
        const unsigned long * src = x;
        unsigned long w;
        unsigned int k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + cado_ctzl(w);
    }
    void add(bcol const& a, unsigned long mask = 1UL) {
        for(unsigned int l = 0 ; l < stride() ; l++) {
            x[l] ^= a.x[l] & -mask;
        }
    }
};/*}}}*/
struct bmat { /* This is for example small-e, an m times b matrix *//*{{{*/
    friend struct polmat;
    unsigned int nrows;
    unsigned int ncols;
    private:
    unsigned long * x;
    inline unsigned int stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    public:
    bmat(unsigned int nrows, unsigned int ncols)
        : nrows(nrows), ncols(ncols), x(mynew<unsigned long>(ncols*stride()))
        
    {
        memset(x, 0, ncols*stride()*sizeof(unsigned long));
    }
    bmat() {
        nrows = ncols = 0;
        x = NULL;
    }
    void clear() {
        mydelete(x, ncols*stride());
        nrows = ncols = 0;
    }
private:
    bmat(bmat const& a) : nrows(a.nrows), ncols(a.ncols), x(a.x) {}
    bmat& operator=(bmat const& a) {
        nrows = a.nrows;
        ncols = a.ncols;
        x = a.x;
        return *this;
    }
public:
    void swap(bmat& a) {
        podswap(nrows, a.nrows);
        podswap(ncols, a.ncols);
        podswap(x,a.x);
    }
    ~bmat() {
        mydelete(x, ncols*stride());
    };
#if 0
    bmat clone() const {
        bmat dst(nrows, ncols);
        memcpy(dst.x, x, ncols*stride()*sizeof(unsigned long));
        return dst;
    }
#endif
    /* zero column j */
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        memset(x + j * stride(), 0, stride() * sizeof(unsigned long));
    }
    bool col_is_zero(unsigned int j) const {
        ASSERT(j < ncols);
        const unsigned long * src = x + j * stride();
        for(unsigned int i = 0 ; i < stride() ; i++) {
            if (src[i] != 0) return false;
        }
        return true;
    }
    bool is_zero() const {
        for(unsigned int i = 0 ; i < ncols * stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline unsigned long coeff(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned int offset = i / ULONG_BITS;
        unsigned long shift = i % ULONG_BITS;
        return (x[j * stride() + offset] >> shift) & 1UL;
    }
    /* add column j0 to column j. Optionally, b indicates a number of rows
     * which are known to be zero in both columns.
     * The mask corresponds to a multiplicating coefficient.
     */
    void acol(unsigned int j, unsigned int j0, unsigned int b = 0, unsigned long mask = 1UL) {
        ASSERT(j < ncols);
        ASSERT(j0 < ncols);
        unsigned long * dst = x + j * stride();
        const unsigned long * src = x + j0 * stride();
        for(unsigned int l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    void acol(unsigned int j, bcol& a, unsigned int b = 0, unsigned long mask = 1UL) {
        ASSERT(j < ncols);
        unsigned long * dst = x + j * stride();
        const unsigned long * src = a.x;
        for(unsigned int l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    unsigned int ffs(unsigned int j) const {
        ASSERT(j < ncols);
        unsigned int z = 0;
        const unsigned long * src = x + j * stride();
        unsigned long w;
        unsigned int k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + cado_ctzl(w);
    }
    void extract_col(bcol& a, unsigned int j) const
    {
        ASSERT(j < ncols);
        bcol tmp_a(nrows);
        memcpy(tmp_a.x, x + j * stride(), stride() * sizeof(unsigned long));
        a.swap(tmp_a);
    }
    void addcoeff(unsigned int i, unsigned int j, unsigned long z)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned int offset = i / ULONG_BITS;
        x[j * stride() + offset] ^= z << (i % ULONG_BITS);
    }
};/*}}}*/
struct polmat { /* {{{ */
    /* polmat may represent:
     * - big-E, an m times b matrix of polynomials having nc max coeffs
     * - PI,     a b times b matrix of polynomials having nc max coeffs
     *
     * coeffs of big-E will eventually be stored in reversed bit order for
     * efficiency, but this is of little concern here.
     *
     * polynomials are stored packed, so this is not an array of bmat's.
     *
     * polmat may also represent other things, but for which speed is not
     * important.
     */
    unsigned int nrows;
    unsigned int ncols;
    unsigned long ncoef;
    static bool critical;
    private:
    unsigned long * x;
    long   * _deg;
    inline size_t stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }/*{{{*/
    inline size_t colstride() const { return nrows * stride(); }/*}}}*/
    public:
    long& deg(unsigned int j) { ASSERT(j < ncols); return _deg[j]; }/*{{{*/
    long deg(unsigned int j) const { ASSERT(j < ncols); return _deg[j]; }
    long maxdeg() const {
        long m = -1;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            if (_deg[j] > m) m = _deg[j];
        }
        return m;
    }/*}}}*/
    inline unsigned long maxlength() const { return 1 + maxdeg(); }
    void alloc() {
        /* we don't care about exceptions */
        x = mynew<unsigned long>(ncols*colstride());
        _deg = mynew<long>(ncols);
    }
    void clear() {
        mydelete(x, ncols*colstride());
        mydelete(_deg, ncols);
    }
    /* ctors dtors etc {{{ */
    polmat(unsigned int nrows, unsigned int ncols, unsigned long ncoef)
        : nrows(nrows), ncols(ncols), ncoef(ncoef)
    {
        alloc();
        memset(x, 0, ncols*colstride()*sizeof(unsigned long));
        for(unsigned int j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    polmat() {
        nrows = ncols = ncoef = 0;
        x = NULL;
        _deg = NULL;
    }
    private:
    polmat(polmat const& a MAYBE_UNUSED) { }
    polmat& operator=(polmat const&){ return *this;}
    public:
    ~polmat() {
        clear();
    }
    inline void swap(polmat & n) {
        podswap(nrows,n.nrows);
        podswap(ncols,n.ncols);
        podswap(ncoef,n.ncoef);
        podswap(x,n.x);
        podswap(_deg,n._deg);
    }
    inline void copy(polmat & n) {
        clear();
        nrows=n.nrows;
        ncols=n.ncols;
        ncoef=n.ncoef;
        alloc();
        memcpy(_deg,n._deg,ncols*sizeof(long));
        memcpy(x, n.x, ncols*colstride()*sizeof(unsigned long));
    }
    /* }}} */
    public:
#if 0
    polmat clone() {
        polmat dst(nrows, ncols, ncoef);
        memcpy(dst.x, x, ncols*colstride()*sizeof(unsigned long));
        memcpy(dst._deg, _deg, ncols*sizeof(long));
        return dst;
    }
#endif
    /* this handles expansion and shrinking */
    void resize(unsigned long ncoef2) {/*{{{*/
        size_t newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            ncoef = ncoef2;
            return;
        }
        size_t minstride = std::min(stride(),newstride);
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(nrows,ncols,ncoef2);
        unsigned long * dst = n.x;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            const unsigned long * src = x + j * colstride();
            for(unsigned int i = 0 ; i < nrows ; i++) {
                memcpy(dst, src, minstride * sizeof(unsigned long));
                dst += newstride;
                src += stride();
            }
        }
        swap(n);
    }/*}}}*/
    /* Divide by X^k, keep ncoef2 coefficients *//*{{{*/
    void xdiv_resize(unsigned long k, unsigned long ncoef2) {
        ASSERT(k + ncoef2 <= ncoef);
        ASSERT(k);
        size_t newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            /*
            mp_size_t sw = k / GMP_LIMB_BITS;
            mp_size_t sb = k & (GMP_LIMB_BITS - 1);
            mpn_rshift(x, x + sw, ncols * colstride() - sw, sb);
            */
            ASSERT(k && k < GMP_LIMB_BITS);
            mpn_rshift(x, x, ncols * colstride(), k);
            ncoef = ncoef2;
            return;
        }
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(nrows,ncols,ncoef2);

        ASSERT(GMP_LIMB_BITS == ULONG_BITS);
        mp_size_t sw = k / GMP_LIMB_BITS;
        mp_size_t sb = k & (GMP_LIMB_BITS - 1);
        size_t input_length = BITS_TO_WORDS(k + ncoef2, ULONG_BITS) - sw;

        /* sb might equal zero here, in which case sw >0. So k is
         * sw * ULONG_BITS, and input_length == newstride
         */

        /* We know that input_length <= stride */
        if (newstride < input_length) {
            ASSERT(sb && sb < GMP_LIMB_BITS);
            /* {{{ In this case we're encumbered by the fact that
             * mpn_rshift might land outside our main area for the last
             * cell. So we'll use an ugly hack */

            unsigned int i = 0 ;
            unsigned int j = 0 ;
            unsigned long * dst = n.poly(0,0);
            for(j = 0 ; j < ncols-1 ; j++) {
                const unsigned long * src = poly(0,j);
                for(i = 0 ; i < nrows ; i++) {
                    mpn_rshift(dst, src + sw, input_length, sb);
                    dst += newstride;
                    src += stride();
                }
            }

            {
                const unsigned long * src = poly(0,j);
                for(i = 0 ; i < nrows - 1 ; i++) {
                    mpn_rshift(dst, src + sw, input_length, sb);
                    dst += newstride;
                    src += stride();
                }

                unsigned long * tmp = mynew<unsigned long>(input_length);
                mpn_rshift(tmp, src + sw, input_length, sb);
                memcpy(dst, tmp, newstride * sizeof(unsigned long));
                mydelete(tmp, input_length);
            }
            /*}}}*/
        } else {
            unsigned long * dst = n.poly(0,0);
            if (sb == 0) {
                ASSERT(input_length == newstride);
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    const unsigned long * src = poly(0,j);
                    for(unsigned int i = 0 ; i < nrows ; i++) {
                        memcpy(dst, src + sw,
                                newstride * sizeof(unsigned long));
                        dst += newstride;
                        src += stride();
                    }
                }
            } else {
                ASSERT(sb && sb < GMP_LIMB_BITS);
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    const unsigned long * src = poly(0,j);
                    for(unsigned int i = 0 ; i < nrows ; i++) {
                        mpn_rshift(dst, src + sw, input_length, sb);
                        dst += newstride;
                        src += stride();
                    }
                }
            }
        }
        swap(n);
        for(unsigned int j = 0 ; j < ncols ; j++) {
            if (deg(j)) deg(j) -= k;
        }
    }
    /*}}}*/
    void xmul_col(unsigned int j, unsigned long s=1) {/*{{{*/
        ASSERT(j < ncols);
        mp_limb_t * dst = x + j * colstride();
        ASSERT(1ul <= s && s <= GMP_LIMB_BITS-1);
        mpn_lshift(dst, dst, colstride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        deg(j) += deg(j) >= 0;
        // BUG_ON(deg(j) >= (int) ncoef);
        // we do NOT consider this a bug. Instead, this means that te
        // corresponding column has gone live. 
        // Normally not a problem since it is supposed to occur only
        // when sufficiently many generators are known.
    }/*}}}*/
    void xmul_poly(unsigned int i, unsigned int j, unsigned long s=1) {/*{{{*/
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned long * dst = x + (j * nrows + i) * stride();
        ASSERT(1 <= s && s <= GMP_LIMB_BITS-1);
        mpn_lshift(dst, dst, stride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        // deg(j) += deg(j) >= 0;
    }/*}}}*/
    void addpoly(unsigned int i, unsigned int j, polmat const& y, unsigned int iy, unsigned int jy) {/*{{{*/
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(iy < y.nrows);
        ASSERT(jy < y.ncols);
        unsigned long * dst = x + (j * nrows + i) * stride();
        unsigned long * src = y.x + (jy * y.nrows + iy) * y.stride();
        mpn_xor_n(dst, dst, src, MIN(stride(), y.stride()));
    }/*}}}*/
    void xclean0_col(unsigned int j) {/*{{{*/
        ASSERT(j < ncols);
        unsigned long * dst = x + j * colstride();
        for(unsigned int i = 0 ; i < nrows ; i++, dst += stride())
            dst[0] &= ~1UL;
        deg(j) -= deg(j) == 0;
    }/*}}}*/
    /* {{{ add column j times mask to column i. */
    void acol(unsigned int i, unsigned int j, unsigned long mask = 1UL) {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned long * dst = x + i * colstride();
        const unsigned long * src = x + j * colstride();
        for(unsigned int l = 0 ; l < colstride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
        // ASSERT(_deg[i] >= _deg[j]);
        deg(i) = std::max(deg(j), deg(i));
    }
    /* }}} */
    unsigned int ffs(unsigned int j, unsigned int k) const {/*{{{*/
        ASSERT(k < nrows);
        ASSERT(j < ncols);
        unsigned int offset = k / ULONG_BITS;
        unsigned long mask = 1UL << (k % ULONG_BITS);
        const unsigned long * src = x + j * colstride() + offset;
        for(unsigned int z = 0 ; z < nrows ; z++) {
            if (*src & mask) return z;
            src += stride();
        }
        return UINT_MAX;
    }/*}}}*/

    unsigned long valuation() const {/*{{{*/
        /* stride() is the number of words it takes for one polynomial.
         * We'll just be OR'ing all polynomials into one, that's it.
         */
        /* It's fairly ugly, but not critical */
        unsigned long y[stride() + 1];
        for(size_t k = 0 ; k < stride() ; k++) {
            y[k] = 0;
        }
        y[stride()] = 1;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            const unsigned long * src = poly(0,j);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int k = 0 ; k < stride() ; k++) {
                    y[k] |= src[k];
                }
                src += stride();
            }
        }
        return mpn_scan1(y,0);
    }/*}}}*/
    void setdeg(unsigned int j) { /* {{{ */
        unsigned long y[stride()];
        for(unsigned int k = 0 ; k < stride() ; k++) {
            y[k] = 0;
        }
        const unsigned long * src = poly(0,j);
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int k = 0 ; k < stride() ; k++) {
                y[k] |= src[k];
            }
            src += stride();
        }
        long k;
        /* Keep only (ULONG_BITS-ncoef) mod ULONG BITS in the top word. */
        y[stride()-1] &= ~0UL >> ((-ncoef) & (ULONG_BITS-1));
        for(k = stride() - 1 ; k >= 0 && y[k] == 0 ; k--) ;
        if (k < 0) {
            deg(j) = -1;
        } else {
            deg(j) = (k+1) * ULONG_BITS - cado_clzl(y[k]) - 1;
        }
    } /* }}} */
    unsigned long * poly(unsigned int i, unsigned int j)/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (j * nrows + i) * stride();
    }
    unsigned long const * poly(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (j * nrows + i) * stride();
    }/*}}}*/
    unsigned long * col(unsigned int j) { return poly(0, j); }/*{{{*/
    unsigned long const * col(unsigned int j) const { return poly(0, j); }/*}}}*/
    /* zero column j *//*{{{*/
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        memset(col(j), 0, colstride() * sizeof(unsigned long));
        deg(j) = -1;
    }/*}}}*/
    /* is zero column j *//*{{{*/
    bool is_zcol(unsigned int j) const {
        ASSERT(j < ncols);
        unsigned long const * y = col(j);
        for(unsigned int x = colstride(); x--;y++)
            if (*y) return false;
        return true;
    }/*}}}*/
    /* shift is understood ``shift left'' (multiply by X) */
    void import_col_shift(unsigned int k, polmat const& a, unsigned int j, long s)/*{{{*/
    {
        ASSERT(k < ncols);
        ASSERT(j < a.ncols);
        ASSERT_ALWAYS(a.nrows == nrows);
        unsigned long const * src = a.col(j);
        unsigned long * dst = col(k);
        unsigned long as = s > 0 ? s : -s;
        mp_size_t sw = as / GMP_LIMB_BITS;
        mp_size_t sb = as & (GMP_LIMB_BITS - 1);
        /* Beware, trailing bits are lurking here and there */
        if (colstride() == a.colstride()) {
            /* Then this is fast *//*{{{*/
            // notice that since we shift the whole column, bits shifted
            // out from a row reappear as noise bit higher up.
            if (sb) {
                if (s > 0) {
                    mpn_lshift(dst + sw, src, colstride() - sw, sb);
                } else {
                    mpn_rshift(dst, src + sw, colstride() - sw, sb);
                }
            } else {
                if (s > 0) {
                    memcpy(dst + sw, src, (colstride()-sw)*sizeof(mp_limb_t));
                } else {
                    memcpy(dst, src + sw, (colstride()-sw)*sizeof(mp_limb_t));
                }
            }/*}}}*/
        } else if (colstride() < a.colstride()) {
            /* Otherwise, we have to resample... *//*{{{*/
            unsigned long tmp[a.colstride()];
            ASSERT_ALWAYS(!critical);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                if (sb) {
                    if (s > 0) {
                        mpn_lshift(tmp + sw, src, a.stride() - sw, sb);
                    } else {
                        mpn_rshift(tmp, src + sw, a.stride() - sw, sb);
                    }
                } else {
                    if (s > 0) {
                        memcpy(tmp + sw, src,(a.stride()-sw)*sizeof(mp_limb_t));
                    } else {
                        memcpy(tmp, src + sw,(a.stride()-sw)*sizeof(mp_limb_t));
                    }
                }
                memcpy(dst, tmp, stride() * sizeof(unsigned long));
                src += a.stride();
                dst += stride();
            }/*}}}*/
        } else {
            ASSERT_ALWAYS(0);/*{{{*/
            // I doubt we'll ever end up here, so I do not have a test
            // case. It should be ok though 
            ASSERT_ALWAYS(!critical);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                if (sb) {
                    if (s > 0) {
                        dst[stride()] = mpn_lshift(dst + sw, src, a.stride() -sw, sb);
                    } else {
                        mpn_rshift(dst, src + sw, a.stride() - sw, sb);
                    }
                } else {
                    if (s > 0) {
                        memcpy(dst+sw,src,(a.stride()-sw)*sizeof(mp_limb_t));
                    } else {
                        memcpy(dst,src+sw,(a.stride()-sw)*sizeof(mp_limb_t));
                    }
                }
                src += a.stride();
                dst += stride();
            }/*}}}*/
        }
    }/*}}}*/
    /* these accessors are not the preferred ones for performance */
    inline unsigned long coeff(unsigned int i, unsigned int j, unsigned long k) const/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
        ASSERT_ALWAYS(!critical);
        unsigned long offset = k / ULONG_BITS;
        unsigned long shift = k % ULONG_BITS;
        return (poly(i,j)[offset] >> shift) & 1UL;
    }/*}}}*/
    void extract_coeff(bmat& a, unsigned long k) const/*{{{*/
    {
        ASSERT(k < ncoef);
        ASSERT_ALWAYS(!critical);
        bmat tmp_a(nrows, ncols);
        for(unsigned int j = 0 ; j < ncols ; j++) {
            for(unsigned int i = 0 ; i < nrows ; i++) {
                unsigned int boffset = i / ULONG_BITS;
                unsigned long bmask = 1UL << (i % ULONG_BITS);
                unsigned long coeff = this->coeff(i,j,k);
                tmp_a.x[j*tmp_a.stride() + boffset] ^= bmask & -coeff;
            }
        }
        a.swap(tmp_a);
    }/*}}}*/
    void addcoeff(unsigned int i, unsigned int j, unsigned long k, unsigned long z)/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
        ASSERT_ALWAYS(!critical);
        // brev_warning();
        size_t offset = k / ULONG_BITS;
        poly(i,j)[offset] ^= z << (k % ULONG_BITS);
    }/*}}}*/
    void clear_highbits() {/*{{{*/
        size_t offset = ncoef / ULONG_BITS;
        unsigned long mask = (1UL << (ncoef % ULONG_BITS)) - 1UL;
        if (!mask) return;
        unsigned long * where = x + offset;
        for(unsigned int i = 0 ; i < nrows * ncols ; i++) {
            *where &= mask;
            where += stride();
        }
    }/*}}}*/
    uint32_t crc() const {
        if (stride() * sizeof(unsigned long) == BITS_TO_WORDS(ncoef, 64) * sizeof(uint64_t)) {
            return crc32(x, ncols * nrows * stride() * sizeof(unsigned long));
        } else {
            /* otherwise it's gonna be painful...  */
            ASSERT_ALWAYS(sizeof(unsigned long) == sizeof(uint32_t));
            ASSERT_ALWAYS(stride() & 1);
            cado_crc_lfsr l;
            cado_crc_lfsr_init(l);
            const uint32_t * xx = (const uint32_t *)(x);
            uint32_t z = 0;
            uint32_t w = 0;
            for(unsigned int count = nrows * ncols ; count-- ; xx += stride()) {
                w = cado_crc_lfsr_turn32_little(l, xx, stride() * sizeof(unsigned long));
                w = cado_crc_lfsr_turn32_little(l, &z, sizeof(unsigned long));
            }
            cado_crc_lfsr_clear(l);
            return w;
        }
    }
};


/*}}}*/
template<typename fft_type> struct tpolmat /* {{{ */
{
    /* tpolmat is essentially the same as polmat, except that it is
     * expected to hold transform data. Only a few operations are
     * allowed.
     */
    unsigned int nrows;
    unsigned int ncols;
    fft_type * po;
    private:
    typename fft_type::ptr x;
    int   * _deg;
    // inline unsigned int stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }
    // inline unsigned int colstride() const { return nrows * stride(); }
    public:
    int& deg(unsigned int j) { ASSERT(j < ncols); return _deg[j]; }
    int const & deg(unsigned int j) const { ASSERT(j < ncols); return _deg[j]; }
    tpolmat(unsigned int nrows, unsigned int ncols, fft_type & o)
        : nrows(nrows), ncols(ncols), po(&o),
        x(o.alloc(nrows * ncols)),
        _deg(mynew<int>(ncols))
    {
        for(unsigned int j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    void zero()
    {
        po->zero(x, nrows * ncols);
    }
    void clear() {
        po->free(x, nrows * ncols); x = NULL;
        mydelete(_deg, ncols);
        nrows = ncols = 0;
    }
    tpolmat() {
        nrows = ncols = 0;
        x = NULL;
        _deg = NULL;
        po = NULL;
    }
    private:
    tpolmat(tpolmat const&) { }
    tpolmat& operator=(tpolmat const&){ return *this;}
    public:
    ~tpolmat() {
        if (po)
            po->free(x, nrows * ncols);
        x = NULL;
        mydelete(_deg, ncols);
    }
    inline void swap(tpolmat & n) {
        podswap(nrows,n.nrows);
        podswap(ncols,n.ncols);
        podswap(po,n.po);
        podswap(x,n.x);
        podswap(_deg,n._deg);
    }
    public:
#if 0
    tpolmat clone() {
        tpolmat dst(nrows, ncols, o);
        memcpy(dst.x, x, ncols*colstride()*sizeof(unsigned long));
        memcpy(dst._deg, _deg, ncols*sizeof(unsigned int));
        return dst;
    }
#endif
    typename fft_type::ptr poly(unsigned int i, unsigned int j)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return po->get(x, j * nrows + i);
    }
    typename fft_type::srcptr poly(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return po->get(x, j * nrows + i);
    }
    typename fft_type::ptr col(unsigned int j) { return poly(0, j); }
    typename fft_type::srcptr col(unsigned int j) const { return poly(0, j); }
    /* zero column j */
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        po->zero(col(j), 1);
        deg(j) = -1;
    }
    uint32_t crc() const { return crc32((unsigned long *) x, nrows * ncols * po->size() * sizeof(*x)); }
};
/*}}}*/

template<typename fft_type>
void transform(tpolmat<fft_type>& dst, polmat& src, fft_type& o, int d)
{
    // clock_t t = clock();
    tpolmat<fft_type> tmp(src.nrows, src.ncols, o);
    tmp.zero();
    src.clear_highbits();
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
    {
        int k MAYBE_UNUSED = 0;
        for(unsigned int j = 0 ; j < src.ncols ; j++) {
            for(unsigned int i = 0 ; i < src.nrows ; i++) {
                if (OMP_ROUND(k++))
                    o.dft(tmp.poly(i,j), src.poly(i,j), d);
            }
        }
    }
    dst.swap(tmp);
    /*
    if (o.size() > 10) {
        printf("t(%u,%u,%u,%u): %.2fs\n",
                src.nrows,src.ncols,d,o.size(),
                (double)(clock()-t)/CLOCKS_PER_SEC);
    }
    */
}

template<typename fft_type>
void glue4(
        tpolmat<fft_type> & dst,
        tpolmat<fft_type> const & s00,
        tpolmat<fft_type> const & s01,
        tpolmat<fft_type> const & s10,
        tpolmat<fft_type> const & s11,
        fft_type& o)
{
    unsigned int nr2 = dst.nrows >> 1;
    unsigned int nc2 = dst.ncols >> 1;
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
    {
        int k MAYBE_UNUSED = 0;
        for(unsigned int i = 0; i < nr2; ++i) {
            for(unsigned int j = 0; j < nc2; ++j) {
                if (OMP_ROUND(k++))
                    o.cpy(dst.poly(i,j), s00.poly(i,j));
                if (OMP_ROUND(k++))
                    o.cpy(dst.poly(i,j+nc2), s01.poly(i,j));
                if (OMP_ROUND(k++))
                    o.cpy(dst.poly(i+nr2,j), s10.poly(i,j));
                if (OMP_ROUND(k++))
                    o.cpy(dst.poly(i+nr2,j+nc2), s11.poly(i,j));
            }
        }
    }
}

template<typename fft_type>
void splitin4(
        tpolmat<fft_type>& dst00,
        tpolmat<fft_type>& dst01,
        tpolmat<fft_type>& dst10,
        tpolmat<fft_type>& dst11,
        tpolmat<fft_type> const & s,
        fft_type& o)
{
    ASSERT((s.nrows & 1) == 0);
    ASSERT((s.ncols & 1) == 0);
    unsigned int nr2 = s.nrows >> 1;
    unsigned int nc2 = s.ncols >> 1;
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
    {
        int k MAYBE_UNUSED = 0;
        for(unsigned int i = 0; i < nr2; ++i) {
            for(unsigned int j = 0; j < nc2; ++j) {
                if (OMP_ROUND(k++))
                    o.cpy(dst00.poly(i,j), s.poly(i,j));
                if (OMP_ROUND(k++))
                    o.cpy(dst01.poly(i,j), s.poly(i,j+nc2));
                if (OMP_ROUND(k++))
                    o.cpy(dst10.poly(i,j), s.poly(i+nr2,j));
                if (OMP_ROUND(k++))
                    o.cpy(dst11.poly(i,j), s.poly(i+nr2,j+nc2));
            }
        }
    }
}

template<typename fft_type>
void add(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o)
{
    ASSERT(s1.nrows == s2.nrows);
    ASSERT(s1.ncols == s2.ncols);
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
    {
        int k MAYBE_UNUSED = 0;
        for(unsigned int i = 0 ; i < s1.nrows ; i++) {
            for(unsigned int j = 0 ; j < s1.ncols ; j++) {
                if (OMP_ROUND(k++))
                    o.add(dst.poly(i,j), s1.poly(i,j), s2.poly(i,j));
            }
        }
    }
}


template<typename fft_type, typename strassen_selector>
void compose_strassen(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o, strassen_selector const& s)
{
    ASSERT(s1.ncols == s2.nrows);
    // Build submatrices
    // We don't even try do to it in place
    ASSERT((s1.nrows & 1) == 0);
    ASSERT((s1.ncols & 1) == 0);
    unsigned int nr2 = s1.nrows >> 1;
    unsigned int nc2 = s1.ncols >> 1;
    tpolmat<fft_type> A11(nr2, nc2, o);
    tpolmat<fft_type> A12(nr2, nc2, o);
    tpolmat<fft_type> A21(nr2, nc2, o);
    tpolmat<fft_type> A22(nr2, nc2, o);
    tpolmat<fft_type> tmpA(nr2, nc2, o);
    splitin4(A11, A12, A21, A22, s1, o);

    ASSERT((s2.nrows & 1) == 0);
    ASSERT((s2.ncols & 1) == 0);
    nr2 = s2.nrows >> 1;
    nc2 = s2.ncols >> 1;
    tpolmat<fft_type> B11(nr2, nc2, o);
    tpolmat<fft_type> B12(nr2, nc2, o);
    tpolmat<fft_type> B21(nr2, nc2, o);
    tpolmat<fft_type> B22(nr2, nc2, o);
    tpolmat<fft_type> tmpB(nr2, nc2, o);
    splitin4(B11, B12, B21, B22, s2, o);

    // Build partial products
    unsigned int nr = s1.nrows >> 1;
    unsigned int nc = s2.ncols >> 1;
    // M1 = (A11 + A22)*(B11 +  B22)
    tpolmat<fft_type> M1(nr, nc, o);
    add(tmpA, A11, A22, o); add(tmpB, B11, B22, o);
    compose_inner(M1, tmpA, tmpB, o, s);
    // M2 = (A21 + A22)*B11
    tpolmat<fft_type> M2(nr, nc, o);
    add(tmpA, A21, A22, o); 
    compose_inner(M2, tmpA, B11, o, s);
    // M3 = A11*(B12-B22)
    tpolmat<fft_type> M3(nr, nc, o);
    add(tmpB, B12, B22, o); 
    compose_inner(M3, A11, tmpB, o, s);
    // M4 = A22*(B21-B11)
    tpolmat<fft_type> M4(nr, nc, o);
    add(tmpB, B21, B11, o); 
    compose_inner(M4, A22, tmpB, o, s);
    // M5 = (A11+A12)*B22
    tpolmat<fft_type> M5(nr, nc, o);
    add(tmpA, A11, A12, o);
    compose_inner(M5, tmpA, B22, o, s);
    // M6 = (A21-A11)*(B11+B12)
    tpolmat<fft_type> M6(nr, nc, o);
    add(tmpA, A21, A11, o); add(tmpB, B11, B12, o);
    compose_inner(M6, tmpA, tmpB, o, s);
    // M7 = (A12-A22)*(B21+B22)
    tpolmat<fft_type> M7(nr, nc, o);
    add(tmpA, A12, A22, o); add(tmpB, B21, B22, o);
    compose_inner(M7, tmpA, tmpB, o, s);

    // Reconstruct result
    // C11 = M1+M4-M5+M7  Store it in M7
    add(M7, M7, M1, o); add(M7, M7, M4, o); add(M7, M7, M5, o);
    // C12 = M3+ M5    Store it in M5
    add(M5, M5, M3, o);
    // C21 = M2+M4     Store it in M4
    add(M4, M4, M2, o);
    // C22 = M1 - M2 + M3 + M6  Store it in M6
    add(M6, M6, M1, o); add(M6, M6, M2, o); add(M6, M6, M3, o);

    glue4(dst, M7, M5, M4, M6, o);
}

struct strassen_default_selector {

// How to tune a threshold for Strassen's algo ????
// This is a 3-d problem...
// And this probably depends also on the degrees of the polynomials.

int operator()(unsigned int m, unsigned int n, unsigned int p, unsigned int nbits) const {
    static const unsigned int minsize = 4;
    static const unsigned int minbits = 12800;
    if (m < minsize || n < minsize || p < minsize || nbits < minbits) 
        return 0;
#if 0
    // not so much of a problem for the moment, and anyway it sould be
    // smarter to do one top-level recursion with simple addmuls,
    // avoiding extra temps, and then keep the possibility to descend to
    // strassen mults
    uint64_t temp_strassen = (n>>1)*(p>>1);
    temp_strassen *= nbits >> 3;
    if (temp_strassen >> 30) {
        /* Avoid strassen multiplication when temporaries will lead to
         * large extra memory requirements */
        fprintf(stderr, "%u*%u*%u %u avoids Strassen\n", m,n,p,nbits);
        return 0;
    }
#endif
    return 1;
}

};

template<typename T> struct remove_pointer {};
template<typename T> struct remove_pointer<T*> {typedef T t;};

template<typename fft_type, typename selector_type>
void compose_inner(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o, selector_type const& s)
{
    tpolmat<fft_type> tmp(s1.nrows, s2.ncols, o);
    ASSERT(s1.ncols == s2.nrows);
    unsigned int nbits;
    nbits = o.size() * sizeof(typename remove_pointer<typename fft_type::ptr>::t) * CHAR_BIT;
    if (s(s1.nrows, s1.ncols, s2.ncols, nbits)) {
        compose_strassen(tmp, s1, s2, o, s);
    } else {
        tmp.zero();
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
        {
            typename fft_type::ptr x = o.alloc(1);
	    /* This way of doing matrix multiplication is better for locality:
	       in the inner loop, the first element s1[i,k] is constant,
	       and in the second one s2[k,j], only the second index changes.
	       If the inner loop was on k, both s1[i,k] and s2[k,j] would
	       change, and a change on the first index is worse if matrices
	       are stored row by row. */
	    for(unsigned int i = 0 ; i < s1.nrows ; i++) {
	      for(unsigned int k = 0 ; k < s1.ncols ; k++) {
		for(unsigned int j = 0 ; j < s2.ncols ; j++) {
		  if (OMP_ROUND((int) (i * s2.ncols + j))) {
		    o.compose(x, s1.poly(i,k), s2.poly(k,j));
		    o.add(tmp.poly(i,j), tmp.poly(i,j), x);
		  }
		}
	      }
            }
            o.free(x,1);
            x = NULL;
        }
    }
    dst.swap(tmp);
}
template<typename fft_type>
inline void compose(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o)
{
    // clock_t t = clock();
    compose_inner(dst, s1, s2, o, strassen_default_selector());
    /*
    if (o.size() > 10) {
        printf("c(%u,%u,%u,%u): %.2fs [strassen]\n",
                s1.nrows,s1.ncols,s2.ncols,o.size(),
                (double)(clock()-t)/CLOCKS_PER_SEC);
    }
    */
}

template<typename fft_type, typename selector_type>
inline void compose2(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o, selector_type const& s)
{
    // clock_t t = clock();
    compose_inner(dst, s1, s2, o, s);
    /*
    if (o.size() > 10) {
        printf("c(%u,%u,%u,%u): %.2fs [strassen]\n",
                s1.nrows,s1.ncols,s2.ncols,o.size(),
                (double)(clock()-t)/CLOCKS_PER_SEC);
    }
    */
}


template<typename fft_type>
void itransform(polmat& dst, tpolmat<fft_type>& src, fft_type& o, int d)
{
    // clock_t t = clock();
    polmat tmp(src.nrows, src.ncols, d + 1);
#ifdef  HAVE_OPENMP
#pragma omp parallel
#endif  /* HAVE_OPENMP */
    {
        int k MAYBE_UNUSED = 0;
        for(unsigned int j = 0 ; j < src.ncols ; j++) {
            for(unsigned int i = 0 ; i < src.nrows ; i++) {
                if (OMP_ROUND(k++))
                    o.ift(tmp.poly(i,j), d, src.poly(i,j));
            }
        }
    }
    for(unsigned int j = 0 ; j < src.ncols ; j++) {
        tmp.setdeg(j);
    }
    dst.swap(tmp);
    /*
       if (o.size() > 10) {
       printf("i(%u,%u,%u,%u): %.2fs\n",
       src.nrows,src.ncols,d,o.size(),
       (double)(clock()-t)/CLOCKS_PER_SEC);
       }
       */
}


#endif	/* LINGEN_MAT_TYPES_HPP_ */
