#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#ifdef HAVE_SSE2
#include <emmintrin.h>                    // for __m128i, _mm_setzero_si128
#endif
#include <algorithm>    // for min
#include <cinttypes>                      // for PRId64, PRIu64
#include <cmath>                          // for log2
#include <cstdarg>             // IWYU pragma: keep
#include <cstring>                        // for size_t, NULL, memset
#include <sys/types.h>                    // for ssize_t
#include <array>                          // for array, array<>::value_type
#include <cstdint>                        // for uint32_t, uint8_t
#include <iterator>                       // for begin, end
#include <memory>                         // for allocator, __shared_ptr_access
#include <mutex>                          // for lock_guard, mutex
#include <utility>                        // for move
#include <vector>                         // for vector
#include <gmp.h>                          // for gmp_vfprintf, mpz_srcptr
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "las-process-bucket-region.hpp"  // for process_bucket_region_spawn
#include "bucket.hpp"                     // for bare_bucket_update_t<>::br_...
#include "fb.hpp"                         // for fb_factorbase::slicing, fb_...
#include "las-apply-buckets.hpp"          // for apply_one_bucket
#include "las-auxiliary-data.hpp"         // for nfs_aux, nfs_aux::thread_data
#include "las-cofac-standalone.hpp"       // for cofac_standalone
#include "las-cofactor.hpp"               // for check_leftover_norm
#include "las-config.h"                   // for LOG_BUCKET_REGION, BUCKET_R...
#include "las-coordinates.hpp"            // for NxToIJ, adjustIJsublat
#include "las-descent-trees.hpp"          // for descent_tree
#include "las-descent.hpp"                // for register_contending_relation
#include "las-detached-cofac.hpp"         // for cofac_standalone, detached_...
#include "las-divide-primes.hpp"          // for divide_known_primes, factor...
#include "las-dumpfile.hpp"               // for dumpfile_t
#include "las-globals.hpp"                // for exit_after_rel_found, globa...
#include "las-info.hpp"                   // for las_info, las_info::batch_p...
#include "las-multiobj-globals.hpp"       // for dlp_descent
#include "las-norms.hpp"                  // for lognorm_smart
#include "las-output.hpp"                 // for TRACE_CHANNEL, las_output
#include "las-qlattice.hpp"               // for qlattice_basis
#include "las-report-stats.hpp"           // for TIMER_CATEGORY, las_report
#include "las-siever-config.hpp"          // for siever_config, siever_confi...
#include "las-smallsieve-types.hpp"       // for small_sieve_data_t
#include "las-smallsieve.hpp"             // for resieve_small_bucket_region
#include "las-threads-work-data.hpp"      // for nfs_work, nfs_work::side_data
#include "las-todo-entry.hpp"             // for las_todo_entry
#include "las-unsieve.hpp"                // for search_survivors_in_line
#include "las-where-am-i-proxy.hpp"            // for where_am_I
#include "las-where-am-i.hpp"             // for where_am_I, WHERE_AM_I_UPDATE
#include "macros.h"                       // for ASSERT_ALWAYS, ASSERT, MAX
#include "tdict.hpp"                      // for slot, timetree_t, CHILD_TIMER
#include "threadpool.hpp"                 // for worker_thread, thread_pool
#include "verbose.h"

#include <cstdio>
#include <iostream>

#include "smallbatch.hpp"

MAYBE_UNUSED static inline void subusb(unsigned char *S1, unsigned char *S2, ssize_t offset)
{
    int ex = (unsigned int) S1[offset] - (unsigned int) S2[offset];
    if (UNLIKELY(ex < 0)) S1[offset] = 0; else S1[offset] = ex;	     
}

/* S1 = S1 - S2, with "-" in saturated arithmetic,
 * and memset(S2, 0, EndS1-S1).
 */
static void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2) {/*{{{*/
#ifndef HAVE_SSE2
    ssize_t mysize = EndS1 - S1;
    unsigned char *cS2 = S2;
    while (S1 < EndS1) {
        subusb(S1,S2,0);
        subusb(S1,S2,1);
        subusb(S1,S2,2);
        subusb(S1,S2,3);
        subusb(S1,S2,4);
        subusb(S1,S2,5);
        subusb(S1,S2,6);
        subusb(S1,S2,7);
        S1 += 8; S2 += 8;
    }
    memset(cS2, 0, mysize);
#else
    __m128i *S1i = (__m128i *) S1, *EndS1i = (__m128i *) EndS1, *S2i = (__m128i *) S2,
            z = _mm_setzero_si128();
    while (S1i < EndS1i) {
        __m128i x0, x1, x2, x3;
        __asm__ __volatile__
            ("prefetcht0 0x1000(%0)\n"
             "prefetcht0 0x1000(%1)\n"
             "movdqa (%0),%2\n"
             "movdqa 0x10(%0),%3\n"
             "movdqa 0x20(%0),%4\n"
             "movdqa 0x30(%0),%5\n"
             "psubusb (%1),%2\n"
             "psubusb 0x10(%1),%3\n"
             "psubusb 0x20(%1),%4\n"
             "psubusb 0x30(%1),%5\n"
             "movdqa %6,(%1)\n"
             "movdqa %6,0x10(%1)\n"
             "movdqa %6,0x20(%1)\n"
             "movdqa %6,0x30(%1)\n"
             "movdqa %2,(%0)\n"
             "movdqa %3,0x10(%0)\n"
             "movdqa %4,0x20(%0)\n"
             "movdqa %5,0x30(%0)\n"
             "add $0x40,%0\n"
             "add $0x40,%1\n"
             : "+&r"(S1i), "+&r"(S2i), "=&x"(x0), "=&x"(x1), "=&x"(x2), "=&x"(x3) : "x"(z));
        /* I prefer use ASM than intrinsics to be sure each 4
         * instructions which use exactly a cache line are together. I'm
         * 99% sure it's not useful...  but it's more beautiful :-)
         */
        /*
           __m128i x0, x1, x2, x3;
           _mm_prefetch(S1i + 16, _MM_HINT_T0); _mm_prefetch(S2i + 16, _MM_HINT_T0);
           x0 = _mm_load_si128(S1i + 0);         x1 = _mm_load_si128(S1i + 1);
           x2 = _mm_load_si128(S1i + 2);         x3 = _mm_load_si128(S1i + 3);
           x0 = _mm_subs_epu8(S2i[0], x0);       x1 = _mm_subs_epu8(S2i[1], x1);
           x2 = _mm_subs_epu8(S2i[2], x2);       x3 = _mm_subs_epu8(S2i[3], x3);
           _mm_store_si128(S2i + 0, z);          _mm_store_si128(S1i + 1, z);
           _mm_store_si128(S2i + 2, z);          _mm_store_si128(S1i + 3, z);
           _mm_store_si128(S1i + 0, x0);         _mm_store_si128(S1i + 1, x1);
           _mm_store_si128(S1i + 2, x2);         _mm_store_si128(S1i + 3, x3);
           S1i += 4; S2i += 4;
           */
    }
#endif 
}/*}}}*/

struct process_bucket_region_run : public process_bucket_region_spawn {/*{{{*/
    worker_thread * worker;
    nfs_aux::thread_data & taux;
    nfs_work::thread_data & tws;
    timetree_t & timer;
    int bucket_relative_index;
    las_report& rep;
    std::vector<unsigned char *> S;
    /* We will have this point to the thread's where_am_I data member.
     * (within nfs_aux::th). However it might be just as easy to let this
     * field be defined here, and drop the latter.
     */
    where_am_I & w;
    bool do_resieve;

    /* A note on SS versus S[side]
     *
     * SS is temp data. It's only used here, and it could well be defined
     * here only. We declare it at the thread_data level to avoid
     * constant malloc()/free().
     *
     * S[side] is where we compute the norm initialization. Some
     * tolerance is subtracted from these lognorms to account for
     * accepted cofactors.
     *
     * SS is the bucket region where we apply the buckets, and also later
     * where we do the small sieve.
     *
     * as long as SS[x] >= S[side][x], we are good.
     */

    unsigned char *SS;
    
    struct side_data {/*{{{*/
        bucket_array_complete purged;   /* for purge_buckets */
        bucket_primes_t primes;         /* for resieving */
    };/*}}}*/

    std::vector<side_data> sides;

    process_bucket_region_run(process_bucket_region_spawn const & p, timetree_t & timer, worker_thread * worker, int id);

    /* will be passed as results of functions
    std::vector<uint32_t> survivors;
    std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors2;
     * */

    /* most probably useless, I guess
    int N;
    int cpt;
    int copr;
    */

    void init_norms(int side);

    template<bool with_hints>
    void apply_buckets_inner(int side);

    void apply_buckets(int side);
    void small_sieve(int side);
    void SminusS(int side);
    typedef std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors_t;
    survivors_t search_survivors();
    void purge_buckets(int side);
    void resieve(int side);
    void cofactoring_sync (survivors_t & survivors2);
    void operator()();
};
/*}}}*/

/* process_bucket_region, split into pieces. */
process_bucket_region_run::process_bucket_region_run(process_bucket_region_spawn const & p, timetree_t & timer, worker_thread * worker, int id): /* {{{ */
    process_bucket_region_spawn(p),
    worker(worker),
    taux(aux_p->th[worker->rank()]),
    tws(ws.th[worker->rank()]),
    timer(timer),
    bucket_relative_index(id),
    rep(taux.rep),
    S(ws.las.cpoly->nb_polys),
    w(taux.w),
    sides(ws.las.cpoly->nb_polys)
{
    int nsides = sides.size();

    w = w_saved;
    WHERE_AM_I_UPDATE(w, N, first_region0_index + already_done + bucket_relative_index);

    /* This is local to this thread */
    for(int side = 0 ; side < nsides ; side++) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb()) {
            S[side] = NULL;
        } else {
            S[side] = tws.sides[side].bucket_region;
            ASSERT_ALWAYS(S[side]);
        }
    }

    SS = tws.SS;
    memset(SS, 0, BUCKET_REGION);

    /* see comment in process_bucket_region_run::operator()() */
    do_resieve = ws.conf.sides[0].lim && ws.conf.sides[1].lim;

    /* we're ready to go ! processing is in the operator() method.
    */
}/*}}}*/
void process_bucket_region_spawn::operator()(worker_thread * worker, int id) /*{{{{*/
{
    timetree_t & timer(aux_p->get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    /* create a temp object with more fields, and dispose it shortly
     * afterwards once we're done.  */
    process_bucket_region_run(*this, timer, worker, id)();
}/*}}}*/
void process_bucket_region_run::init_norms(int side)/*{{{*/
{
    CHILD_TIMER(timer, "init norms");
    TIMER_CATEGORY(timer, norms(side));

    int N = first_region0_index + already_done + bucket_relative_index;

    ws.sides[side].lognorms.fill(S[side], N);

#if defined(TRACE_K) 
    if (trace_on_spot_N(w->N))
        verbose_output_print(TRACE_CHANNEL, 0, "# After side %d init_norms_bucket_region, N=%u S[%u]=%u\n",
                side, w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
}/*}}}*/

template<bool with_hints> void process_bucket_region_run::apply_buckets_inner(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);

    typedef typename hints_proxy<with_hints>::l my_longhint_t;
    typedef typename hints_proxy<with_hints>::s my_shorthint_t;
    {
        CHILD_TIMER(timer, "apply buckets");
        TIMER_CATEGORY(timer, sieving(side));
        for (auto const & BA : wss.bucket_arrays<1, my_shorthint_t>())
            apply_one_bucket(SS, BA, already_done + bucket_relative_index, wss.fbs->get_part(1), w);
    }

    /* Apply downsorted buckets, if necessary. */
    if (ws.toplevel > 1) {
        CHILD_TIMER(timer, "apply downsorted buckets");
        TIMER_CATEGORY(timer, sieving(side));

        for (auto const & BAd : wss.bucket_arrays<1, my_longhint_t>()) {
            // FIXME: the updates could come from part 3 as well,
            // not only part 2.
            ASSERT_ALWAYS(ws.toplevel <= 2);
            apply_one_bucket(SS, BAd, already_done + bucket_relative_index, wss.fbs->get_part(2), w);
        }
    }
}/*}}}*/
void process_bucket_region_run::apply_buckets(int side)/*{{{*/
{
    if (do_resieve) {
        apply_buckets_inner<true>(side);
    } else {
        apply_buckets_inner<false>(side);
    }
}/*}}}*/

static void update_checksums(nfs_work::thread_data & tws, nfs_aux::thread_data & taux)
{
    for(int side = 0 ; side < 2 ; side++)
        taux.update_checksums(side, tws.sides[side].bucket_region, BUCKET_REGION);
}

void process_bucket_region_run::small_sieve(int side)/*{{{*/
{
    CHILD_TIMER(timer, "small sieve");
    TIMER_CATEGORY(timer, sieving(side));

    nfs_work::side_data & wss(ws.sides[side]);

    sieve_small_bucket_region(SS,
            first_region0_index + already_done + bucket_relative_index,
            wss.ssd,
            wss.ssd.ssdpos_many[bucket_relative_index],
            ws.conf.logI, ws.Q.sublat,
            w);
}/*}}}*/
void process_bucket_region_run::SminusS(int side)/*{{{*/
{
    /* compute S[side][x] = max(S[side][x] - SS[x], 0),
     * and clear SS.  */
    CHILD_TIMER(timer, "S minus S (2)");
    TIMER_CATEGORY(timer, cofactoring(side));

    ::SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
    if (trace_on_spot_N(w->N))
        verbose_output_print(TRACE_CHANNEL, 0,
                "# Final value on side %d, N=%u S[%u]=%u\n",
                side, w->N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
}/*}}}*/
process_bucket_region_run::survivors_t process_bucket_region_run::search_survivors() /*{{{*/
{
    typedef std::vector<uint32_t> surv1_t;

    surv1_t temp_sv;

    CHILD_TIMER(timer, __func__);
    TIMER_CATEGORY(timer, search_survivors());

    int N = first_region0_index + already_done + bucket_relative_index;

    /* change N, which is a bucket number, to
     * (i0, i1, j0, j1) */
    int logI = ws.conf.logI;
    /* This bit of code is replicated from las-smallsieve.cpp */
    const unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);
    const unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);
    const unsigned int regions_per_line = 1 << log_regions_per_line;           
    const unsigned int region_rank_in_line = N & (regions_per_line - 1);       
    const unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    
    const unsigned int j1 MAYBE_UNUSED = j0 + (1 << log_lines_per_region);    
    const int I = 1 << logI;                                            
    const int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          
    const int i1 = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     

    ASSERT(j1 > j0); /* even when we have a line fragment */


#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# When entering factor_survivors for bucket %u, "
                "S[0][%u]=%u, S[1][%u]=%u\n",
                trace_Nx.N, trace_Nx.x,
                S[0] ? S[0][trace_Nx.x] : ~0u,
                trace_Nx.x,
                S[1] ? S[1][trace_Nx.x] : ~0u);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n",
                (mpz_srcptr) traced_norms[0],
                (mpz_srcptr) traced_norms[1]);
    }
#endif  /* }}} */

#ifdef TRACE_K /* {{{ */
    nfs_work::side_data & side0(ws.sides[0]);
    nfs_work::side_data & side1(ws.sides[1]);
    for (int x = 0; x < 1 << LOG_BUCKET_REGION; x++) {
        if (trace_on_spot_Nx(N, x)) {
            verbose_output_print(TRACE_CHANNEL, 0,
                    "# side0.Bound[%u]=%u, side1.Bound[%u]=%u\n",
                    S[0] ? S[0][trace_Nx.x] : ~0u,
                    S[0] ? (S[0][x] <= side0.lognorms.bound ? 0 : side0.lognorms.bound) : ~0u,
                    S[1] ? S[1][trace_Nx.x] : ~0u,
                    S[1] ? (S[1][x] <= side1.lognorms.bound ? 0 : side1.lognorms.bound) : ~0u);
        }
    }
#endif /* }}} */

    rep.survivors.before_sieve += 1U << LOG_BUCKET_REGION;

    temp_sv.reserve(128);

    for (unsigned int j = j0; j < j1; j++)
    {
        int offset = (j-j0) << logI;

        unsigned char * const both_S[2] = {
            S[0] ? S[0] + offset : NULL,
            S[1] ? S[1] + offset : NULL,
        };
        /* TODO FIXME XXX that's weird. How come don't we merge that with
         * the lognorm computation that goes in the ws.sides[side]
         * regions before apply_buckets + small_sieve ?? Could it help
         * save a bit of time in search_survivors_in_line ?
         */
        const unsigned char both_bounds[2] = {
            ws.sides[0].lognorms.bound,
            ws.sides[1].lognorms.bound,
        };
        size_t old_size = temp_sv.size();

        ASSERT(j < ws.J);

        search_survivors_in_line(both_S, both_bounds,
                j,
                i0, i1,
                N,
                *ws.jd,
                ws.conf.unsieve_thresh,
                *ws.us,
                temp_sv,
                ws.Q.sublat);

        /* Survivors written by search_survivors_in_line() have index
         * relative to their j-line. We need to convert to index within
         * the bucket region by adding line offsets.
         */


        /* When several bucket regions are in a line, we have nothing to
         * change. Note in particular that we must not adjust with
         * respect to the starting i -- this info is already encoded with
         * the bucket number N.
         */

        if (!offset) continue;

        for (size_t i_surv = old_size; i_surv < temp_sv.size(); i_surv++)
            temp_sv[i_surv] += offset;
    }

    /* This used to be called convert_survivors */
    return survivors_t(begin(temp_sv), end(temp_sv));
}/*}}}*/
void process_bucket_region_run::purge_buckets(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);

    SIBLING_TIMER(timer, "purge buckets");
    TIMER_CATEGORY(timer, cofactoring(side));

    unsigned char * Sx = S[0] ? S[0] : S[1];

    for (auto & BA : wss.bucket_arrays<1, shorthint_t>()) {
#if defined(HAVE_SSE2) && defined(SMALLSET_PURGE)
        sides[side].purged.purge(BA, already_done + bucket_relative_index, Sx, survivors);
#else
        sides[side].purged.purge(BA, already_done + bucket_relative_index, Sx);
#endif
    }

    /* Add entries coming from downsorting, if any */
    for (auto const & BAd : wss.bucket_arrays<1, longhint_t>()) {
        sides[side].purged.purge(BAd, already_done + bucket_relative_index, Sx);
    }

    /* Sort the entries to avoid O(n^2) complexity when looking for
       primes during trial division */
    sides[side].purged.sort();
}/*}}}*/
void process_bucket_region_run::resieve(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);
    SIBLING_TIMER(timer, "resieve");
    TIMER_CATEGORY(timer, sieving(side));

    unsigned char * Sx = S[0] ? S[0] : S[1];

    /* Resieve small primes for this bucket region and store them 
       together with the primes recovered from the bucket updates */
    resieve_small_bucket_region (&sides[side].primes,
            Sx,
            first_region0_index + already_done + bucket_relative_index,
            wss.ssd,
            wss.ssd.ssdpos_many[bucket_relative_index],
            ws.conf.logI, ws.Q.sublat,
            w);

    /* same reason as above */
    sides[side].primes.sort();
}/*}}}*/


void process_bucket_region_run::cofactoring_sync (survivors_t & survivors)/*{{{*/
{
    int nsides = sides.size();
    int surv_pre_batch_bit_size[2] = {0, 0};

    /* by declaring this timer "fuzzy", we make the child timers use only
     * userspace calls, and not system calls. This makes it possible to
     * be really fine-grain, at only little expense.
     */
    CHILD_TIMER_FUZZY(this->timer, timer, __func__);
    // CHILD_TIMER(timer, __func__);
    TIMER_CATEGORY(timer, cofactoring_mixed());

    int N = first_region0_index + already_done + bucket_relative_index;
    unsigned char * Sx = S[0] ? S[0] : S[1];

    cofac_standalone cur;

    siever_config cf = ws.conf; //- obsolete?
    
    std::vector<cofac_standalone> surv_pre_batch; //-

    for (size_t i_surv = 0 ; i_surv < survivors.size(); i_surv++) {
        if (dlp_descent && ws.las.tree.must_take_decision())
            break;
        const size_t x = survivors[i_surv];
        ASSERT_ALWAYS (Sx[x] != 255);
        ASSERT(x < ((size_t) 1 << LOG_BUCKET_REGION));

        rep.survivors.after_sieve++;

        if (S[0] && S[1])
            rep.mark_survivor(S[0][x], S[1][x]);

        /* For factor_leftover_norm, we need to pass the information of the
         * sieve bound. If a cofactor is less than the square of the sieve
         * bound, it is necessarily prime. we implement this by keeping the
         * log to base 2 of the sieve limits on each side, and compare the
         * bitsize of the cofactor with their double.
         */

        SIBLING_TIMER(timer, "check_coprime");

        /* start building a new object. This is a swap operation */
        cur = cofac_standalone(nsides, N, x, ws.conf.logI, ws.Q);

        for(int side = 0 ; side < 2 ; side++) {
            if (ws.sides[side].no_fb()) continue;
            cur.S[side] = S[side][x];
        }

        if (cur.trace_on_spot())
            verbose_output_print(TRACE_CHANNEL, 0, "# about to start cofactorization for (%"
                    PRId64 ",%" PRIu64 ")  %zu %u\n", cur.a, cur.b, x, Sx[x]);

        /* since a,b both even were not sieved, either a or b should
         * be odd. However, exceptionally small norms, even without
         * sieving, may fall below the report bound (see tracker
         * issue #15437). Therefore it is safe to continue here. */
        // ASSERT((a | b) & 1);
        if (UNLIKELY(cur.both_even()))
            continue;

        rep.survivors.not_both_even++;

        if (!cur.gcd_coprime_with_q(ws.Q.doing))
            continue;

        rep.survivors.not_both_multiples_of_p++;

        BOOKKEEPING_TIMER(timer);
        int pass = 1;

        int i;
        unsigned int j;
        // Note that are, (i,j) must be true coordinates, not the
        // ones reduced to (-I/2, I/2) using sublattices.
        NxToIJ (i, j, N, x, ws.conf.logI);
        adjustIJsublat(i, j, ws.Q.sublat);

        //- ws.conf.sides[side].fbK.k_fb_product vs ws.sides[side].fbK.k_fb_product;
        
        // cxx_mpz fbp0; //- = ws.sides[0].fbK.k_fb_product;
        // cxx_mpz fbp1; //- = ws.sides[1].fbK.k_fb_product;

        
        // fbp0 = ws.sides[0].fbs->slicing_fb_product;
        // fbp1 = ws.sides[1].fbs->slicing_fb_product;

        // const cxx_mpz fbp0sd = ws.las.get_fb_product_shared_data(0);
        // const cxx_mpz fbp1sd = ws.las.get_fb_product_shared_data(1);

        // if (fbp0 > (cxx_mpz)0) 
        // {
        //     std::cout << "fb_prod a : " << fbp0 << "::" << fbp1 << "\n";
        //     std::cout << "fb_prod b : " << fbp0sd << "::" << fbp1sd << "\n";
            
        // }   
            

        //- ws.conf.sides[side].fbK.k_small_batch_max_prime

        /*
        //- auto truc = ws.sides[0].fbK.no_trial_div;
        
        unsigned int sbmp0, sbmp1;
        */
        if (do_resieve) {

            for(int pside = 0 ; pass && pside < 2 ; pside++) {
                int side = trialdiv_first_side ^ pside;
                nfs_work::side_data & wss(ws.sides[side]);






                //- auto truc = ws.sides[0].fbK.no_trial_div;
                unsigned int no_trial_div = wss.fbK.no_trial_div;
                
                /* 
                sbmp0 = ws.sides[0].fbK.k_small_batch_max_prime;
                sbmp1 = ws.sides[1].fbK.k_small_batch_max_prime;

                std::cout << "SBMP : " << sbmp0 << "::" << sbmp1 << "\n";
                std::cout << "FBPR : " << fbp0 << "::" << fbp1 << "\n";
                */



                CHILD_TIMER_PARAMETRIC(timer, "side ", side, " pre-cofactoring checks");
                TIMER_CATEGORY(timer, cofactoring(side));


                SIBLING_TIMER(timer, "recompute complete norm");

                // Trial divide norm on side 'side'
                /* Compute the norms using the polynomials transformed to 
                   i,j-coordinates. The transformed polynomial on the 
                   special-q side is already divided by q */
                wss.lognorms.norm(cur.norm[side], i, j);
                //- gmp_printf("## (%li, %li)(%Zd)[%d] \n", cur.a, cur.b, cur.norm[side], side);

                if (cur.trace_on_spot()) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0,
                            gmp_vfprintf, "# start trial division for norm=%Zd ", (mpz_srcptr) cur.norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            "on side %d for (%" PRId64 ",%" PRIu64 ")\n", side, cur.a, cur.b);
                }

                if (wss.no_fb()) {
                    /* This is a shortcut. We're probably replacing sieving
                     * by a product tree, there's no reason to bother doing
                     * trial division at this point (or maybe there is ?
                     * would that change the bit size significantly ?) */
                    rep.survivors.check_leftover_norm_on_side[side] ++; //- stats
                    continue;
                }

                SIBLING_TIMER(timer, "trial division");

                verbose_output_print(1, 2, "FIXME %s, line %d\n", __FILE__, __LINE__);
                const bool handle_2 = true; /* FIXME */
                rep.survivors.trial_divided_on_side[side]++;


                divide_known_primes (cur.factors[side], cur.norm[side], N, x,
                        handle_2,
                        &sides[side].primes,
                        &sides[side].purged,
                        *wss.td,
                        cur.a, cur.b,
                        *wss.fbs,
                        no_trial_div);

                /* if q is composite, its prime factors have not been sieved.
                 * Check if they divide. They probably don't, since we
                 * have computed the norm with the polynomials adapted to
                 * the (i,j) plane, and q divided out. But still,
                 * valuations are not desired here.
                 */
                if ((side == ws.Q.doing.side) && (!ws.Q.doing.is_prime())) {
                    for (const auto &x : ws.Q.doing.prime_factors) {
                        if (mpz_divisible_uint64_p(cur.norm[side], x)) {
                            mpz_divexact_uint64(cur.norm[side], cur.norm[side], x);
                            cur.factors[side].push_back(x);
                        }
                    }
                }

                SIBLING_TIMER(timer, "check_leftover_norm");

                siever_side_config & scs = ws.conf.sides[side]; //-
                

                pass = check_leftover_norm_pre_batch (cur.norm[side], scs); //-ws.conf.sides[side]);
                /*if (cur.b == 18){
                    fprintf(stderr, "\n\n---------------------\n(%ld, %lu), pass = %d\n---------------------\n\n", cur.a, cur.b, pass);
                }*/



                /*

                sbmp0 = ws.sides[0].fbK.k_small_batch_max_prime;
                sbmp1 = ws.sides[1].fbK.k_small_batch_max_prime;

                std::cout << "SBMPb : " << sbmp0 << "::" << sbmp1 << "\n";
                std::cout << "FBPRb : " << fbp0 << "::" << fbp1 << "\n";
                */

                if (cur.trace_on_spot()) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                            "# checked leftover norm=%Zd", (mpz_srcptr) cur.norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            " on side %d for (%" PRId64 ",%" PRIu64 "): %d\n",
                            side, cur.a, cur.b, pass);
                }
                rep.survivors.check_leftover_norm_on_side[side] += pass;
            }
        } else {
            ASSERT_ALWAYS(ws.las.batch || ws.las.batch_print_survivors.filename);

            /* no resieve, so no list of prime factors to divide. No
             * point in doing trial division anyway either.
             */
            for(int side = 0 ; side < 2 ; side++) {
                CHILD_TIMER_PARAMETRIC(timer, "side ", side, " pre-cofactoring checks");
                TIMER_CATEGORY(timer, cofactoring(side));

                SIBLING_TIMER(timer, "recompute complete norm");

                nfs_work::side_data & wss(ws.sides[side]);

                /* factor() in batch.cpp recomputes the complete norm, so
                 * there's no need to compute the norm right now for the
                 * side we've sieved with.
                 */
                if (wss.no_fb()) {
                    wss.lognorms.norm(cur.norm[side], i, j);
                } else {
                    /* This is recognized specially in the
                     * factor_simple_minded() code in batch.cpp
                     */
                    mpz_set_ui(cur.norm[side], 0);
                }

                /* We don't even bother with q and its prime factors.
                 * We're expecting to recover just everything after the
                 * game anyway */

                /* Note that we're *NOT* doing the equivalent of
                 * check_leftover_norm here. This is explained by two
                 * things:
                 *
                 *  - while the "red zone" of post-sieve values that we
                 *  know can't yield relations is quite wide (from L to
                 *  B^2), it's only a marginal fraction of the total
                 *  number of reports. Even more so if we take into
                 *  account the necessary tolerance near the boundaries
                 *  of the red zone.
                 *
                 *  - we don't have the complete norm (with factors taken
                 *  out) at this point, so there's no way we can do a
                 *  primality check -- which is, in fact, the most
                 *  stringent check because it applies to the bulk of the
                 *  candidates.
                 *
                 * Bottom line: we just hand over *everything* to the
                 * batch cofactorization.
                 */
            }
        }

        if (!pass) continue;

        

        if (do_resieve)
            for (int side_i = 0; side_i < 2; side_i++)
            {
                siever_side_config sc = cf.sides[side_i];
                unsigned long lim = sc.lim;
                std::sort(cur.factors[side_i].begin(), cur.factors[side_i].end());
                gmp_printf("(%li, %li)(%Zd)[%d]:", cur.a, cur.b, cur.norm[side_i], side_i);
                for (auto const& z : cur.factors[side_i]) {
                        assert(z < lim);
                        printf("%lx,", z);
                }
                printf("\n");
            }

                

#if 0
        if (ws.las.batch || ws.las.batch_print_survivors) {
            /* in these cases, we won't go through detached_cofac, hence
             * we won't check for potential (a,b) duplicates. Those can
             * and will happen if we encounter overflowing buckets at
             * level 2.
             *
             * However, in the long term we expect that this won't happen
             * for real, as overflowing buckets are unlikely to occur
             * unless very early in the process. Therefore we can omit
             * the check here, and do it only at relation printing time.
             */
            nfs_aux::rel_hash_t& rel_hash(aux_p->get_rel_hash());
            nfs_aux::abpair_t ab(cur.a, cur.b);
            std::lock_guard<std::mutex> foo(rel_hash.mutex());
            if (!rel_hash.insert(ab).second) continue;
        }
#endif

        rep.survivors.enter_cofactoring++; // seg fault si commenté ici et mis à la place pile avant l'ECM, mais tel quel on incrémente plus qu'on n'envoie à ECM ?

        //- bypass
        surv_pre_batch.push_back(cur);
        surv_pre_batch_bit_size[0] += mpz_sizeinbase(cur.norm[0], 2);
        surv_pre_batch_bit_size[1] += mpz_sizeinbase(cur.norm[1], 2);
/*
        gmp_fprintf(stderr, "%li,%li:", cur.a, cur.b);
        for (auto const& z : cur.factors[0]) gmp_fprintf(stderr, "%lx,", z);    
        gmp_fprintf(stderr, ":");
        for (auto const& z : cur.factors[1]) gmp_fprintf(stderr, "%lx,", z);
        gmp_fprintf(stderr, "\n");

*/        

        continue;
        assert(0);

        // we'll do the printing later.
        if (ws.las.batch || ws.las.batch_print_survivors.filename)
        {
            /* see above */
            rep.reports++;
            if (ws.conf.sublat_bound && !cur.ab_coprime()) continue;
            /* make sure threads don't write the cofactor list at the
             * same time !!! */
            cur.transfer_to_cofac_list(ws.cofac_candidates, aux_p->doing);
            std::cout << "transfer_to_cofac_list\n";
            continue; /* we deal with all cofactors at the end of subjob */
        }

        auto D = new detached_cofac_parameters(wc_p, aux_p, std::move(cur));

        if (!dlp_descent) {
            /* We must make sure that we join the async threads at some
             * point, otherwise we'll leak memory. It seems more appropriate
             * to batch-join only, so this is done at the las_subjob level */
            // worker->get_pool().get_result(1, false);
            worker->get_pool().add_task(detached_cofac, D, N, 1); /* id N, queue 1 */
        } else {
            /* We must proceed synchronously for the descent */
            auto res = dynamic_cast<detached_cofac_result*>(detached_cofac(worker, D, N));
            bool cc = false;
            if (res->rel_p) {
                cc = register_contending_relation(ws.las, ws.Q.doing, *res->rel_p);
            }
            delete res;
            if (cc)
                break;
        }
    }
    //- end for !

#if 0 //- old version + DEBUG
    //-std::cout << "end for\n";
    for (auto cur_last : surv_pre_batch) {
        if (ws.las.batch || ws.las.batch_print_survivors.filename)
        {
            /* see above */
            rep.reports++;
            if (ws.conf.sublat_bound && !cur_last.ab_coprime()) continue;
            /* make sure threads don't write the cofactor list at the
             * same time !!! */
            cur_last.transfer_to_cofac_list(ws.cofac_candidates, aux_p->doing);
            std::cout << "transfer_to_cofac_list\n";
            continue; /* we deal with all cofactors at the end of subjob */
        }

        auto D = new detached_cofac_parameters(wc_p, aux_p, std::move(cur_last));

        if (!dlp_descent) {
            /* We must make sure that we join the async threads at some
             * point, otherwise we'll leak memory. It seems more appropriate
             * to batch-join only, so this is done at the las_subjob level */
            // worker->get_pool().get_result(1, false);
            worker->get_pool().add_task(detached_cofac, D, N, 1); /* id N, queue 1 */
        } else {
            /* We must proceed synchronously for the descent */
            auto res = dynamic_cast<detached_cofac_result*>(detached_cofac(worker, D, N));
            bool cc = false;
            if (res->rel_p) {
                cc = register_contending_relation(ws.las, ws.Q.doing, *res->rel_p);
            }
            delete res;
            if (cc)
                break;
        }
    }
#endif
#if 1
    int res; //- todo
    //- cxx_mpz fbp0 = ws.sides[0].fbs->slicing_fb_product;
    //- cxx_mpz fbp1 = ws.sides[1].fbs->slicing_fb_product;

    if (surv_pre_batch.size() > 0){




#if 0
        long unsigned min0 = 0xffffffff;
        long unsigned min1 = 0xffffffff;

        for (auto const& z : cur.factors[0]) 
            min0 = min0 > z && z != 2 ? z : min0;
        for (auto const& z : cur.factors[1]) 
            min1 = min1 > z && z != 2 ? z : min1;

        gmp_fprintf(stderr, "facteurs min du batch : %lx : %lx\n", min0, min1);
        

#endif


        /* 
         * it batches on first_side then gets rid of survivors with bigger leftover norms
         * it then batches a second time on fewer survivors on first_side^1
         * not clear which side it is better to start with :
         * side 0 means batching with a smaller factor base
         * side 1 might mean getting rid of more survivors
         */

        int first_side = 1;

#if PRINTBATCH
        gmp_fprintf(stderr, "================ BATCHING (starting side %d  :  %d bits) ================\n", first_side, surv_pre_batch_bit_size[first_side]);
#endif
        for(int i = 0 ; i < 2; i++) {

            int side = first_side ^ i;


            auto conf_side = ws.conf.sides[side];
            //nfs_work::side_data & wss(ws.sides[side]);

            cxx_mpz fbp = ws.sides[side].fbs->slicing_fb_product;
            mpz_t **  tree = ws.sides[side].fbs->slicing_tree;

            //-gmp_fprintf(stderr, "1st surv pre batch   = %Zd\n", surv_pre_batch[0].norm[side]);
            //-std::cout << "1st surv pre batch   = " << surv_pre_batch[0].norm[side] << "\n"; 

#if PRINTBATCH
            gmp_fprintf(stderr, "* Batching %ld survivors on side %d\n", surv_pre_batch.size(), side);
#endif
            //-res = sm_batch(surv_pre_batch, fbp, side); //- primes ?
            res = sm_batch_initalized_tree_split(tree, surv_pre_batch, fbp, side);
            assert(res != -1);      

            /*
            for (auto &surv_tmp : surv_pre_batch)
            if (surv_tmp.b == 18){
                    gmp_fprintf(stderr, "\n\n---------------------\n(%ld, %lu), leftover norm = %Zd\n---------------------\n\n", surv_tmp.a, surv_tmp.b, surv_tmp.norm[side]);
            }
            */
/*
            std::vector<int> filtered;
            //- for (auto &surv : surv_pre_batch)
            for (int j = 0; j < surv.size(); j++)
            {
                int pass = check_leftover_norm_post_batch(surv.norm[side], ws.conf.sides[side]);
                if (!pass)
                    filtered.push_back(j)
            }
*/
            surv_pre_batch.erase(
                remove_if(
                    surv_pre_batch.begin(),
                    surv_pre_batch.end(),
                    [conf_side, side](cofac_standalone checked_surv) {
                        return (!check_leftover_norm_post_batch(checked_surv.norm[side], conf_side));}
                ), 
                surv_pre_batch.end());
#if PRINTBATCH
            gmp_fprintf(stderr, "  | filtering side %d, %ld survivors\n", side, surv_pre_batch.size());
#endif
/*
            for (int j : filtered) {
                surv[j].remove
            }
*/
            if (surv_pre_batch.size() == 0)
                break;
        }


        if (surv_pre_batch.size() > 0) {
            for(int side = 0 ; side < 2; side++) {
                for (auto &surv : surv_pre_batch) {
                    nfs_work::side_data & wss(ws.sides[side]);

#if 0  //- DEBUG
                    gmp_fprintf(stderr, "[BEFORE TRIALDIV], nb_fac = %lu\n", surv.factors[side].size());
                    for (auto fac : surv.factors[side])
                    {
                        gmp_fprintf(stderr, "%ld; ", fac);
                    }
                    gmp_fprintf(stderr, "\n");
                    
                    gmp_fprintf(stderr, "[BEFORE TRIALDIV], nb_lps = %lu\n", surv.lps[side].size());
                    for (auto fac : surv.lps[side])
                    {
                        gmp_fprintf(stderr, "%ld; ", fac);
                    }
                    gmp_fprintf(stderr, "\n");
#endif

#if 1
                    wss.td->trial_divide(surv.factors[side], surv.sm_smoothpart[side]);
                    // delete &surv.sm_smoothpart[side];
                    std::sort (surv.factors[side].begin(), surv.factors[side].end()); //- important ?? dunno
#endif

#if 0 //- DEBUG
                    // tmp test
                    cxx_mpz tmp_mpz;

                    for (auto fac : surv.factors[side]){
                        mpz_init_set_ui (tmp_mpz, fac);
                        surv.lps[side].push_back(tmp_mpz);
                    }

                    gmp_fprintf(stderr, "[AFTER  TRIALDIV], nb_fac = %lu\n", surv.factors[side].size());
                    gmp_fprintf(stderr, "[TEST *D0], side = %d, nb_fac = %lu, ptr = %p\n", side, surv.factors[side].size(), &surv);
                    
                    for (auto fac : surv.factors[side])
                    {
                        gmp_fprintf(stderr, "%ld; ", fac);
                    }
                    gmp_fprintf(stderr, "\n");

                    gmp_fprintf(stderr, "[TEST *D1], side = %d, nb_fac = %lu, ptr = %p\n", side, surv.factors[side].size(), &surv);

                    gmp_fprintf(stderr, "[AFTER  TRIALDIV], nb_lps = %lu\n", surv.lps[side].size());
                    std::sort (surv.lps[side].begin(), surv.lps[side].end()); //- important ??
                    for (auto fac : surv.lps[side])
                    {
                        gmp_fprintf(stderr, "%ld; ", fac);
                    }
                    gmp_fprintf(stderr, "\n");
#endif

                }
            }

            for (auto &cur_last : surv_pre_batch) {

#if 0
                if (ws.las.batch || ws.las.batch_print_survivors.filename)
                {
                    /* see above */
                    rep.reports++;
                    if (ws.conf.sublat_bound && !cur_last.ab_coprime()) continue;
                    /* make sure threads don't write the cofactor list at the
                     * same time !!! */
                    cur_last.transfer_to_cofac_list(ws.cofac_candidates, aux_p->doing);
                    std::cout << "transfer_to_cofac_list\n";
                    continue; /* we deal with all cofactors at the end of subjob */
                }
#endif

#if 1 // ECM
                //-rep.survivors.enter_cofactoring++;

                /*if (cur_last.b == 18){
                    gmp_fprintf(stderr, "\n\n---------------------\n(%ld, %lu), leftover norms = %Zd, %Zd\n---------------------\n\n", cur_last.a, cur_last.b, cur_last.norm[0], cur_last.norm[1]);
                }*/

                if (ws.las.batch || ws.las.batch_print_survivors.filename)
                    fprintf(stderr, "[[ERROR]] ws.las.batch || ws.las.batch_print_survivors.filename\n");
                if (dlp_descent)
                    fprintf(stderr, "[[ERROR]] dlp_descent\n");

#if 0
                gmp_fprintf(stderr, "[TEST *D3], side = 0, nb_fac = %lu, ptr = %p\n", cur_last.factors[0].size(), &cur_last);
                gmp_fprintf(stderr, "[TEST *D3], side = 1, nb_fac = %lu, ptr = %p\n", cur_last.factors[1].size(), &cur_last);
#endif
                auto D = new detached_cofac_parameters(wc_p, aux_p, std::move(cur_last)); //- ICI ?
                worker->get_pool().add_task(detached_cofac, D, N, 1); 

#endif

#if 0
                if (!dlp_descent) {
                    /* We must make sure that we join the async threads at some
                     * point, otherwise we'll leak memory. It seems more appropriate
                     * to batch-join only, so this is done at the las_subjob level */
                    // worker->get_pool().get_result(1, false);
                    worker->get_pool().add_task(detached_cofac, D, N, 1); /* id N, queue 1 */
                } else {
                    /* We must proceed synchronously for the descent */
                    auto res = dynamic_cast<detached_cofac_result*>(detached_cofac(worker, D, N));
                    bool cc = false;
                    if (res->rel_p) {
                        cc = register_contending_relation(ws.las, ws.Q.doing, *res->rel_p);
                    }
                    delete res;
                    if (cc)
                        break;
                }
#endif
            }
        }
    }
    
    // surv_pre_batch.clear();
    /// where leak??




#endif

//-fprintf(stderr, "OUT L-P-B-R\n");


}/*}}}*/


    
void process_bucket_region_run::operator()() {/*{{{*/

    int nsides = sides.size();

    // This is too verbose.
    // fprintf(stderr, "=== entering PBR for report id %lu\n", rep.id);

    /* first check some early abort conditions. */
    if (recursive_descent) {
        /* For the descent mode, we bail out as early as possible. We
         * need to do so in a multithread-compatible way, though.
         * Therefore the following access is mutex-protected within
         * las.tree. */
        if (ws.las.tree.must_take_decision())
            return;
    } else if (exit_after_rel_found) {
        if (rep.reports) {
            std::lock_guard<std::mutex> foo(protect_global_exit_semaphore);
            global_exit_semaphore=1;
            return;
        }
    }

    for (int side = 0; side < nsides; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb()) {
            ASSERT_ALWAYS(S[side] == NULL);
            continue;
        }

        MARK_TIMER_FOR_SIDE(timer, side);
        TIMER_CATEGORY(timer, sieving(side));

        /* Compute norms in S[side] */
        init_norms(side);

        /* Accumulate sieve contributions in SS */
        apply_buckets(side);
        small_sieve(side);

        /* compute S[side][x] = max(S[side][x] - SS[x], 0),
         * and clear SS.  */
        SminusS(side);

        ws.sides[side].dumpfile.write(S[side], BUCKET_REGION);

        BOOKKEEPING_TIMER(timer);
    }

    if (main_output.verbose >= 2)
        update_checksums(tws, taux);

    auto survivors = search_survivors();

    /* The "do_resieve" flag (set above in the ctor) checks if one of the
     * factor bases is empty. This means that we may have decided to
     * *not* sieve on that side, refusing to pay a per-area time a second
     * time. This is based on the rationale that we expect the *other*
     * side to be such a selective test that it isn't worth the trouble.
     * But then, it means that purge_buckets and resieving are not worth
     * the trouble either.
     */

    /* These two steps used to be called "prepare_cofactoring" */
    for(int side = 0 ; !survivors.empty() && do_resieve && side < nsides ; side++) {
        MARK_TIMER_FOR_SIDE(timer, side);
        sides[side].purged.allocate_memory(ws.local_memory, BUCKET_REGION); //- ICI ???
        purge_buckets(side);
        size_t ns = survivors.size();
        double maxnorm = ws.sides[side].lognorms.get_maxlog2();
        double logp_lb = log2(ws.sides[side].fbK.td_thresh);
        size_t nprimes_max = ns * maxnorm / logp_lb;
        sides[side].primes.allocate_memory(ws.local_memory, nprimes_max); //- OU ICI ???
        resieve(side);
    }

#ifdef TRACE_K
    int N = first_region0_index + already_done + bucket_relative_index;
    /* FIXME FIXME FIXME MNFS -- what do we want to do here? */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        unsigned char * Sx = S[0] ? S[0] : S[1];
        verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u has value %u\n",
                trace_Nx.x, trace_Nx.N, Sx[trace_Nx.x]);
    }
#endif

    cofactoring_sync(survivors);
}/*}}}*/


void process_many_bucket_regions(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool, int first_region0_index, where_am_I & w)/*{{{*/
{
    /* first_region0_index is always 0 when toplevel == 1, but the
     * present function is also called from within downsort_tree when
     * toplevel > 1, and then first_region0_index may be larger.
     */
    auto P = thread_pool::make_shared_task<process_bucket_region_spawn>(ws, wc_p, aux_p, w);

    /* Make sure we don't schedule too many tasks when J was truncated
     * anyway */

    int first_skipped_br = ws.J;

    if (ws.conf.logI >= LOG_BUCKET_REGION)
        first_skipped_br <<= ws.conf.logI - LOG_BUCKET_REGION;
    else
        first_skipped_br >>= LOG_BUCKET_REGION - ws.conf.logI;

    size_t small_sieve_regions_ready = std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]);

    for(int done = 0, ready = small_sieve_regions_ready ; done < ws.nb_buckets[1] ; ) {

        /* yes, it's a bit ugly */
        P->first_region0_index = first_region0_index;
        P->already_done = done;

        for(int i = 0 ; i < ready ; i++) {
            if (first_region0_index + done + i >= first_skipped_br) {
                /* Hmm, then we should also make sure that we truncated
                 * fill_in_buckets, right ? */
                break;
            }
            pool.add_shared_task(P, i, 0);
        }

        /* it's only really done when we do drain_queue(0), of course */
        done += ready;

        if (done < ws.nb_buckets[1]) {

            /* We need to compute more init positions */
            int more = std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1] - done);

            for(int side = 0 ; side < 2 ; side++) {
                nfs_work::side_data & wss(ws.sides[side]);
                if (wss.no_fb()) continue;
                pool.add_task_lambda([=,&ws](worker_thread * worker, int){
                        timetree_t & timer(aux_p->get_timer(worker));
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        SIBLING_TIMER(timer, "prepare small sieve");
                        nfs_work::side_data & wss(ws.sides[side]);
                        // if (wss.no_fb()) return;
                        SIBLING_TIMER(timer, "small sieve start positions");
                        /* When we're doing 2-level sieving, there is probably
                         * no real point in doing ssdpos initialization in
                         * several passes.
                         */
                        small_sieve_prepare_many_start_positions(
                                wss.ssd,
                                first_region0_index + done,
                                more,
                                ws.conf.logI, ws.Q.sublat);
                        },0);
            }

            pool.drain_queue(0);

            ready = more;

            /* Now these new start positions are ready to be used */
            for(int side = 0 ; side < 2 ; side++) {
                nfs_work::side_data & wss(ws.sides[side]);
                if (wss.no_fb()) continue;
                small_sieve_activate_many_start_positions(wss.ssd);
            }
        }
    }
}/*}}}*/
