#include "cado.h"

#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "async.h"
#include "filenames.h"
#include "xdotprod.h"
#include "rolling.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"

/*
 * Relatively common manipulation in fact. Move to matmul_top ?
static void xvec_to_vec(matmul_top_data_ptr mmt, uint32_t * gxvecs, int m, unsigned int nx, int d)
{
    mmt_comm_ptr mcol = mmt->wr[d];
    pi_comm_ptr picol = mmt->pi->wr[d];
    int shared = mcol->v->flags & THREAD_SHARED_VECTOR;
    if (!shared || picol->trank == 0) {
        abzero(mmt->abase, mcol->v->v, mcol->i1 - mcol->i0);
        for(int j = 0 ; j < m ; j++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                uint32_t i = gxvecs[j*nx+k];
                // set bit j of entry i to 1.
                if (i < mcol->i0 || i >= mcol->i1)
                    continue;
                abt * where;
                where = mcol->v->v + aboffset(mmt->abase, i-mcol->i0);
                abset_ui(mmt->abase, where, j, 1);
            }
        }
    }
    if (shared)
        serialize_threads(picol);
}
 */

void * krylov_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    matmul_top_data mmt;
    struct timing_data timing[1];

    int flags[2];
    flags[bw->dir] = THREAD_SHARED_VECTOR;
    flags[!bw->dir] = 0;

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        ASSERT_ALWAYS((bw->ys[1]-bw->ys[0]) % 2 == 0);
        ys[0] = bw->ys[0] + pi->interleaved->idx * (bw->ys[1]-bw->ys[0])/2;
        ys[1] = ys[0] + (bw->ys[1]-bw->ys[0])/2;
    }

    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpfq_vbase A;
    mpfq_vbase_oo_field_init_byfeatures(A, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, ys[1]-ys[0],
            MPFQ_DONE);

    pi_datatype_ptr A_pi = pi_alloc_mpfq_datatype(pi, A);

    /* Hmmm. This would deserve better thought. Surely we don't need 64
     * in the prime case. Anything which makes checks relevant will do.
     * For the binary case, we used to work with 64 as a constant, but
     * for the prime case we want to make this tunable (or maybe 1 ?)
     */
    mpfq_vbase Ac;
    mpfq_vbase_oo_field_init_byfeatures(Ac,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    pi_datatype_ptr Ac_pi = pi_alloc_mpfq_datatype(pi, Ac);

    block_control_signals();

    matmul_top_init(mmt, A, A_pi, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_comm_ptr mcol = mmt->wr[bw->dir];
    mmt_comm_ptr mrow = mmt->wr[!bw->dir];

    serialize(pi->m);
    
    A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);


    uint32_t * gxvecs = NULL;
    unsigned int nx = 0;
    if (!fake) {
        load_x(&gxvecs, bw->m, &nx, pi);
    } else {
        set_x_fake(&gxvecs, bw->m, &nx, pi);
    }

    indices_apply_S(mmt, gxvecs, nx * bw->m, bw->dir);
    if (bw->dir == 0)
        indices_apply_P(mmt, gxvecs, nx * bw->m, bw->dir);

    // xvec_to_vec(mmt, gxvecs, bw->m, nx, bw->dir);
    // matmul_top_save_vector(mmt, "Z", bw->dir, 0);

    /* let's be generous with interleaving protection. I don't want to be
     * bothered, really */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);
    matmul_top_comm_bench(mmt, bw->dir);
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    char * v_name = NULL;
    int rc = asprintf(&v_name, V_FILE_BASE_PATTERN, ys[0], ys[1]);
    ASSERT_ALWAYS(rc >= 0);
    if (!fake) {
        if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
        matmul_top_load_vector(mmt, v_name, bw->dir, bw->start, unpadded);
        if (tcan_print) { printf("done\n"); }
    } else {
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
#if 0
        /* This is for setting the source vector to something consistent
         * across mappings, so that given a fixed (fake, here) matrix, any
         * splitting will give the same source vector. This is just a
         * copy of the mechanism which exists in prep for doing exactly
         * this. Alas, in what we denote as a fake situation, there is
         * no chance of course that two different splittings lesad to
         * identical matrices ! Hence, we'd rather not bother with
         * generating something consistent.
         */
        if (pi->m->trank == 0 && !bw->seed) {
            bw->seed = time(NULL);
            MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
        }
        serialize_threads(pi->m);
        gmp_randseed_ui(rstate, bw->seed);
        if (tcan_print) {
            printf("// Random generator seeded with %d\n", bw->seed);
        }
        if (tcan_print) { printf("Creating fake %s...", v_name); fflush(stdout); }
        matmul_top_set_random_and_save_vector(mmt, v_name, bw->dir, bw->start, unpadded, rstate);
        if (tcan_print) { printf("done\n"); }
#else
        unsigned long g = pi->m->jrank * pi->m->ncores + pi->m->trank;
        gmp_randseed_ui(rstate, bw->seed + g);
        matmul_top_set_random_inconsistent(mmt, bw->dir, rstate);
#endif
        gmp_randclear(rstate);
    }

    mmt_vec check_vector;
    mmt_vec ahead;

    mpfq_vbase_tmpl AxAc;
    mpfq_vbase_oo_init_templates(AxAc, A, Ac);

    if (!bw->skip_online_checks) {
        matmul_top_vec_init_generic(mmt, Ac, Ac_pi,
                check_vector, !bw->dir, THREAD_SHARED_VECTOR);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        matmul_top_load_vector_generic(mmt,
                check_vector, CHECK_FILE_BASE, !bw->dir, bw->interval, unpadded);
        if (tcan_print) { printf("done\n"); }
    }

    if (!bw->skip_online_checks) {
        vec_init_generic(pi->m, A, A_pi, ahead, 0, nchecks);
    }

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    mmt_vec xymats;

    if (tcan_print) {
        printf("Each thread allocates %zd kb for the A matrices\n",
                A->vec_elt_stride(A, bw->m*bw->interval) >> 10);
    }
    vec_init_generic(pi->m, A, A_pi, xymats, 0, bw->m*bw->interval);
    
    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        /* The padded dimension is not the important one */
        length = MAX(mmt->n0[0], mmt->n0[1]);
        length = iceildiv(length, bw->m) + iceildiv(length, bw->n);
        length += 2 * iceildiv(bw->m, bw->n);
        length += 2 * iceildiv(bw->n, bw->m);
        length += 10;
        /* Because bw is a global variable, we protect its use */
        if (serialize_threads(pi->m)) {
            bw->end = length;
        }
        serialize(pi->m);
    }

    if (tcan_print) {
        fprintf(stderr, "Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }

#if 0
    /* FIXME -- that's temporary ! only for debugging */
    pi_log_init(pi->m);
    pi_log_init(pi->wr[0]);
    pi_log_init(pi->wr[1]);
#endif

    timing_init(timing, bw->start, bw->interval * iceildiv(bw->end, bw->interval));

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);


    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        /* Note that the check vector is always stored untwisted in
         * memory */

        if (!bw->skip_online_checks) {
            A->vec_set_zero(A, ahead->v, nchecks);
            unsigned int how_many;
            unsigned int offset_c;
            unsigned int offset_v;
            how_many = intersect_two_intervals(&offset_c, &offset_v,
                    mrow->i0, mrow->i1,
                    mcol->i0, mcol->i1);

            if (how_many) {
                AxAc->dotprod(A, Ac, ahead->v,
                        SUBVEC(check_vector, v, offset_c),
                        SUBVEC(mcol->v, v, offset_v),
                        how_many);
            }
        }

        /* Create an empty slot in program execution, so that we don't
         * impose strong constraints on twist/untwist_vector being free of
         * MPI calls.
         */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);
        matmul_top_twist_vector(mmt, bw->dir);

        A->vec_set_zero(A, xymats->v, bw->m*bw->interval);
        serialize(pi->m);
        for(int i = 0 ; i < bw->interval ; i++) {
            /* The first part of this loop must be guaranteed to be free
             * of any mpi calls */
            pi_interleaving_flip(pi);

            /* Compute the product by x */
            x_dotprod(mmt, gxvecs, nx, xymats, i * bw->m, bw->m, 1);
            matmul_top_mul_cpu(mmt, bw->dir);

#ifdef MEASURE_LINALG_JITTER_TIMINGS
            timing_next_timer(timing); /* now timer is [1] (cpu-wait) */
#endif
            pi_interleaving_flip(pi);
            /* from this point on in the loop, mpi calls are allowed */
#ifdef MEASURE_LINALG_JITTER_TIMINGS
            /* This is *only* in order to be able to measure the wait
             * times */
            serialize(pi->m);
#endif

            timing_next_timer(timing); /* now timer is [2] (COMM) */ 
            /* Now we can resume MPI communications. */
            matmul_top_mul_comm(mmt, bw->dir);

#ifdef MEASURE_LINALG_JITTER_TIMINGS
            timing_next_timer(timing); /* now timer is [3] (comm-wait) */
            /* This is *only* in order to be able to measure the wait
             * times */
            serialize(pi->m);
#endif

            timing_next_timer(timing);
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        /* See remark above. */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);

        if (!bw->skip_online_checks) {
            /* Last dot product. This must cancel ! */
            x_dotprod(mmt, gxvecs, nx, ahead, 0, nchecks, -1);

            pi_allreduce(NULL, ahead->v, nchecks, A_pi, BWC_PI_SUM, pi->m);
            if (!A->vec_is_zero(A, ahead->v, nchecks)) {
                printf("Failed check at iteration %d\n", s + bw->interval);
                exit(1);
            }
        }

        matmul_top_untwist_vector(mmt, bw->dir);

        /* Now (and only now) collect the xy matrices */
        pi_allreduce(NULL, xymats->v,
                bw->m * bw->interval,
                A_pi, BWC_PI_SUM, pi->m);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            int rc;
            rc = asprintf(&tmp, A_FILE_PATTERN, ys[0], ys[1], s, s+bw->interval);
            FILE * f = fopen(tmp, "wb");
            rc = fwrite(xymats->v, A->vec_elt_stride(A, 1), bw->m*bw->interval, f);
            if (rc != bw->m*bw->interval) {
                fprintf(stderr, "Ayee -- short write\n");
                // make sure our input data won't be deleted -- this
                // chunk will have to be redone later, maybe the disk
                // failure is temporary (?)
            }
            fclose(f);
            free(tmp);
        }

        matmul_top_save_vector(mmt, v_name, bw->dir, s + bw->interval, unpadded);

        if (pi->m->trank == 0 && pi->m->jrank == 0)
            keep_rolling_checkpoints(v_name, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, mmt->mm->ncoeffs, tcan_print, 0);
    }

    timing_final_tally(pi, timing, mmt->mm->ncoeffs, tcan_print, 0);

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (tcan_print) {
        printf("Done.\n");
    }
    serialize(pi->m);

    vec_clear_generic(pi->m, xymats, bw->m*bw->interval);

    if (!bw->skip_online_checks) {
        matmul_top_vec_clear_generic(mmt, check_vector, !bw->dir);
        vec_clear_generic(pi->m, ahead, nchecks);
    }

    free(gxvecs);
    free(v_name);

    matmul_top_clear(mmt);
    pi_free_mpfq_datatype(pi, A_pi);
    pi_free_mpfq_datatype(pi, Ac_pi);
    A->oo_field_clear(A);
    Ac->oo_field_clear(Ac);

    timing_clear(timing);

    return NULL;
}


int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init_new(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(krylov_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}

