#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <errno.h>
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

void * mksol_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
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

    int nsolvecs = bw->nsolvecs;
    /* In order to use only a subset of the number of solutions which can
     * be extracted with bw, follow this procedure.
     *
     * Split F0-nn into F0-1,  F1-2,  etc. The split program can do this.
     * This amounts to taking one column at a time, for whichever meaning
     * of column one could choose. In the binary case, we have in mind
     * splitting F0-512 into F0-64, F64-128, etc. This splitting still
     * has the complete set of solutions, but spanned across files. The
     * different files F we have just created are suitable for running on
     * distinct clusters: F0-64 will go to cluster number 0 which is
     * iterating with vector V0-64, etc.
     *
     * Now in order to work with nsolvecs < bw->n, an extraction must be
     * done within the F0-64 (or F0-1) files. This consists in picking
     * one set of rows, maybe one row only, for every bw->n rows.
     * This could be done, at least in theory, by the split program, as
     * it's really the same operation. However, it's not the case as it
     * is now.
     *
     * Anything with nsolvecs < bw->n requires the data files to be
     * modified ahead of time.
     */

    unsigned int multi = 1;
    unsigned int nsolvecs_pervec = nsolvecs;

    if (mpz_cmp_ui(bw->p,2) != 0 && nsolvecs > 1) {
        if (tcan_print) {
            fprintf(stderr,
"Note: the present code is a quick hack.\n"
"Some SIMD operations like\n"
"    w_i += c_i * v (with v a constant vector, c_i scalar)\n"
"are performed with a simple-minded loop on i, while there is an acknowledged\n"
"potential for improvement by using more SIMD-aware code.\n");
        }
        multi = nsolvecs;
        nsolvecs_pervec = 1;
    }

    ASSERT_ALWAYS(multi * nsolvecs_pervec == (unsigned int) nsolvecs);

    mpfq_vbase Ar;
    mpfq_vbase_oo_field_init_byfeatures(Ar,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nsolvecs_pervec,
            MPFQ_DONE);

    pi_datatype_ptr Ar_pi = pi_alloc_mpfq_datatype(pi, Ar);

    block_control_signals();

    matmul_top_init(mmt, A, A_pi, pi, flags, pl, bw->dir);
    unsigned int unpadded = MAX(mmt->n0[0], mmt->n0[1]);

    mmt_comm_ptr mcol = mmt->wr[bw->dir];
    mmt_comm_ptr mrow = mmt->wr[!bw->dir];
    pi_comm_ptr picol = mmt->pi->wr[bw->dir];

    serialize(pi->m);
    
    A->vec_set_zero(A, mrow->v->v, mrow->i1 - mrow->i0);

    uint32_t * gxvecs = NULL;
    unsigned int nx = 0;
    
    if (!bw->skip_online_checks) {
        if (!fake) {
            load_x(&gxvecs, bw->m, &nx, pi);
        } else {
            set_x_fake(&gxvecs, bw->m, &nx, pi);
        }
        /*
        indices_apply_S(mmt, gxvecs, nx * bw->m, bw->dir);
        if (bw->dir == 0)
            indices_apply_P(mmt, gxvecs, nx * bw->m, bw->dir);
            */
    }

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
        gmp_randclear(rstate);
    }


    mmt_vec check_vector;
    /* Our ``ahead zone'' will also be used for storing the data prior to
     * transposing */
    mmt_vec ahead;
    vec_init_generic(pi->m, A, A_pi, ahead, 0, bw->n);

    mpfq_vbase_tmpl AxAc;
    mpfq_vbase_oo_init_templates(AxAc, A, Ac);

    mpfq_vbase_tmpl ArxA;
    mpfq_vbase_oo_init_templates(ArxA, Ar, A);

    mpfq_vbase_tmpl AxAr;
    mpfq_vbase_oo_init_templates(AxAr, A, Ar);

    if (!bw->skip_online_checks) {
        matmul_top_vec_init_generic(mmt, Ac, Ac_pi,
                check_vector, !bw->dir, THREAD_SHARED_VECTOR);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        matmul_top_load_vector_generic(mmt,
                check_vector, CHECK_FILE_BASE, !bw->dir, bw->interval, unpadded);
        if (tcan_print) { printf("done\n"); }
    }

    /* F plays the role of a right-hand-side in mksol. Vector iterates
     * have to be multiplied on the right by a matrix corresponding to
     * some coefficient of F. Because F has been written to disk in the
     * proper order (=reversed order, leading coeff first) by lingen,
     * iterate i is to be multiplied by the i-th coefficient of F.
     *
     * The F-coefficients, from the overall description of the algorithm,
     * are square matrices having bw->n rows and columns.  Columns
     * correspond to candidate solutions. In some cases we are interested
     * in only a few solutions, maybe even one for the prime field case.
     * Having bw->n rows correspond to the fact that contributions from
     * all sequences have to be added in order to form a candidate
     * solution.
     *
     * However, for the purpose of splitting the work across sites, F is
     * stored on disk in the transposed order, so that the already
     * existing split program is able to do the split. This means that
     * the data we read here (which is presumed to come out of the split
     * program) consists of matrices having unconditionally
     * bw->n rows, and ys[1]-ys[0] columns.
     *
     * Because our multiply-by-smallish-matrix code does not like this
     * ordering, we have to transpose again in order to obtain a matrix
     * with ys[1]-ys[0] rows, and bw->n columns. This matrix is called
     * rhs in the text below.
     */

    ASSERT(A->vec_elt_stride(A,nsolvecs_pervec) == Ar->vec_elt_stride(Ar, A->groupsize(A)));

    /* We'll load all the F coefficient matrices before the main loop.
     * It's dual to the treatment of the xy matrices in krylov. Same
     * considerations apply w.r.t the memory footprint.
     */
    mmt_vec * fcoeffs;
    fcoeffs = malloc(multi * sizeof(mmt_vec));
    if (tcan_print) {
        printf("Each thread allocates %zd kb for the F matrices\n",
                (multi * Ar->vec_elt_stride(Ar, A->groupsize(A)*bw->interval)) >> 10);
    }
    for(unsigned int k = 0 ; k < multi ; k++) {
        vec_init_generic(pi->m, Ar, Ar_pi, fcoeffs[k], 0, A->groupsize(A)*bw->interval);
    }
    
    /* Our sum vector is only part of the story of course. All jobs, all
     * threads, will participate to the computation of the products by
     * tiny matrices. When dir==1, within a column, all threads have
     * access to the mcol->i0..mcol->i1 interval (which is roughly a
     * fraction equal to 1 / pirow->totalsize of the total
     * vector size. Now this interval will be split in picol->totalsize
     * parts.
     */

    unsigned int ii0, ii1;
    unsigned int eblock = (mcol->i1 - mcol->i0) / picol->totalsize;

    ii0 = mcol->i0 + eblock * (picol->jrank * picol->ncores + picol->trank);
    ii1 = ii0 + eblock;

    mmt_vec * sum;
    sum = malloc(multi * sizeof(mmt_vec));
    for(unsigned int k = 0 ; k < multi ; k++) {
        vec_init_generic(mmt->pi->m, Ar, Ar_pi, sum[k], 0, eblock);
    }

    if (bw->end == 0) {
        // for mksol we use EOF as an ending indication.
        serialize_threads(pi->m);
        if (pi->m->trank == 0)
            bw->end = INT_MAX;
        serialize_threads(pi->m);
        if (tcan_print) {
            fprintf(stderr, "Target iteration is unspecified ;"
                    " going to end of F file\n");
        }
    } else if (tcan_print) {
        fprintf(stderr, "Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }
    
    /* We can either look up the F files to get an idea of the number of
     * coefficients to be considered, or make our guess based on the
     * expected degree of the generator. The latter is obiviously less
     * accurate, but not by much, and anyway for the purpose of making up
     * an ETA, it's good enough.
     */

    int exp_end;
    exp_end = MAX(mmt->n[0], mmt->n[1]);
    exp_end = iceildiv(exp_end, bw->n);

    timing_init(timing, bw->start, bw->interval * iceildiv(exp_end, bw->interval));

    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        serialize(pi->m);
        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            for(unsigned int k = 0 ; k < multi ; k++) {
                Ar->vec_set_zero(Ar, fcoeffs[k]->v, A->groupsize(A)*bw->interval);
                char * tmp;
                rc = asprintf(&tmp, F_FILE_SLICE_PATTERN2, k, k + nsolvecs_pervec, ys[0], ys[1]);

                FILE * f = fopen(tmp, "rb");
                DIE_ERRNO_DIAG(f == NULL, "fopen", tmp);
                rc = fseek(f, A->vec_elt_stride(A, s * nsolvecs_pervec), SEEK_SET);
                if (rc < 0) {
                    bw->end = s;
                }

                for(int i = 0 ; i < bw->interval ; i++) {
                    rc = fread(ahead->v, A->vec_elt_stride(A,1), nsolvecs_pervec, f);
                    ASSERT_ALWAYS(rc == 0 || rc == (int) nsolvecs_pervec);
                    if (rc == 0) {
                        printf("Exhausted F data after %d coefficients\n", s+i);
                        bw->end = s+i;
                        break;
                    }

#if 0
                    {
                        void * test;
                        test = electric_alloc(rstride * abnbits(abase));
                        abvtranspose(mpfq_rhs, abase, test, ahead->v);
                        electric_free(test,rstride * abnbits(abase));
                    }
#endif
                    ArxA->transpose(Ar, A,
                            SUBVEC(fcoeffs[k], v, i * A->groupsize(A)),
                            ahead->v);
                }
                fclose(f);
                free(tmp);
            }
        }
        if (pi->m->trank == 0) {
            /* we're good friends with the global leader. Try to grasp
             * the information about the end value -- it might have
             * changed because of EOF. We're not so disturbed
             * by serialization points here. */
            int err = MPI_Bcast(&(bw->end), 1, MPI_INT, 0, pi->m->pals);
            ASSERT_ALWAYS(!err);
        }
        serialize_threads(pi->m);

        /* broadcast f */
        for(unsigned int k = 0 ; k < multi ; k++) {
            pi_bcast(fcoeffs[k]->v, A->groupsize(A) * bw->interval, Ar_pi, 0, 0, pi->m);
        }

        for(unsigned int k = 0 ; k < multi ; k++) {
            Ar->vec_set_zero(Ar, sum[k]->v, ii1 - ii0);
        }


        // Plan ahead. The check vector is here to predict the final A matrix.
        // Our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
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

        serialize(pi->m);
        /* Create an empty slot in program execution, so that we don't
         * impose strong constraints on twist/untwist_vector being free of
         * MPI calls.
         */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);
        matmul_top_twist_vector(mmt, bw->dir);
        serialize(pi->m);

        /* Despite the fact that the bw->end value might lead us
         * to stop earlier, we want to stop at multiples of the checking
         * interval. That's important, otherwise the last bunch of
         * computations won't be checked.
         */

        for(int i = 0 ; i < bw->interval ; i++) {
            /* The first part of this loop must be guaranteed to be free
             * of any mpi calls */

            /* first timer is [0] (CPU) */
            pi_interleaving_flip(pi);
            for(unsigned int k = 0 ; k < multi ; k++) {
                AxAr->addmul_tiny(A, Ar,
                        sum[k]->v,
                        SUBVEC(mcol->v, v, ii0 - mcol->i0),
                        SUBVEC(fcoeffs[k], v, i * A->groupsize(A)),
                        ii1 - ii0);
            }
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

        /* See remark above. */
        pi_interleaving_flip(pi);
        pi_interleaving_flip(pi);

        matmul_top_untwist_vector(mmt, bw->dir);

        // FIXME. serialize here ??
        // serialize(pi->m);

        if (!bw->skip_online_checks) {
            /* Last dot product. This must cancel ! Recall that x_dotprod
             * adds to the result. */
            x_dotprod(mmt, gxvecs, nx, ahead, 0, nchecks, -1);

            pi_allreduce(NULL, ahead->v, nchecks, A_pi, BWC_PI_SUM, pi->m);
            if (!A->vec_is_zero(A, ahead->v, nchecks)) {
                printf("Failed check at iteration %d\n", s + bw->interval);
                exit(1);
            }
        }

        matmul_top_save_vector(mmt, v_name, bw->dir, s + bw->interval, unpadded);

        /* We need to write our S files, untwisted. There are several
         * ways to do this untwisting. In order to keep the memory
         * footprint to a minimum, we reuse the vector area mcol->v. 
         */
        /* Since S is possibly wider than one single vector, we need
         * several passes. Each corresponds to one batch of potential
         * solutions in the end. There are two levels:
         *
         *  - npasses: this is made of nsolvecs_pervec vectors over the
         *  base field, which correspond to one of several of the vectors
         *  we work with. Those get saved *together* in S files.
         *
         *  - multi: these are saved to distinct files.
         *
         * in reality this corresponds to the same functionality. It's
         * just a design choice to force one behaviour in a situation and
         * not in the other. This mess ought to be fixed. (TODO)
         */
        serialize(pi->m);
        size_t  stride =  A->vec_elt_stride(A, 1);
        size_t rstride =  Ar->vec_elt_stride(Ar, 1);
        ASSERT_ALWAYS(rstride % stride == 0);
        int npasses = rstride / stride;
        for(unsigned int k = 0 ; k < multi ; k++) {
            for(int i = 0 ; i < npasses ; i++) {
                serialize(pi->m);
                matmul_top_zero_vec_area(mmt, bw->dir);
                serialize(pi->m);
                /* Each job/thread copies its data share to mcol->v */
                void * dst;
                const void * src;
                src = pointer_arith_const(sum[k]->v, i * stride);
                dst = pointer_arith(mcol->v->v, (ii0 - mcol->i0) * stride);
                for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
                    /* XXX If we dereference a function pointer for each such
                     * iteration here, we're going to pay a lot. And the mpfq
                     * API does not contain this ``strided copy''. So for the
                     * moment, since all types considered are flat anyway,
                     * this looks good enough.
                     */
                    memcpy(dst, src, stride);
                    dst = pointer_arith(dst, stride);
                    src = pointer_arith_const(src, rstride);
                }
                /* Now we collect the data in mcol->v. Note that each
                 * location is non-zero only at one location, so unreduced
                 * addition is unnecessary.
                 */
                allreduce_across(mmt, bw->dir);
                /* untwist the result */
                matmul_top_untwist_vector(mmt, bw->dir);

                serialize_threads(pi->m);
                dst = pointer_arith(sum[k]->v, i * stride);
                src = pointer_arith_const(mcol->v->v, (ii0 - mcol->i0) * stride);
                for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
                    memcpy(dst, src, stride);
                    dst = pointer_arith(dst, rstride);
                    src = pointer_arith_const(src, stride);
                }
            }
            // TODO: There is not much preventing us from using simply
            // something like:
            // matmul_top_save_vector(mmt, s_name, bw->dir, s +
            // bw->interval, unpadded);
            // since after the matmul_top_untwist_vector call, we really
            // have something which is ready to go trhough
            // matmul_top_save_vector.
            // the difficulty is that if we do so, we need to do so
            // within the loop over i. Which means that we'll create
            // (npasses) files instead of just one. When we provide the
            // preferred-io-width modification (see TODO), this whole
            // code branch should be simplified, and at least the
            // pi_save_file_2d call here should go away.
            char * s_name;
            rc = asprintf(&s_name, S_FILE_BASE_PATTERN ".%u", k, k + nsolvecs_pervec, ys[0], ys[1], s+bw->interval);
            pi_file_handle f;
            pi_file_open(f, pi, bw->dir, s_name, "wb", Ar->vec_elt_stride(Ar, unpadded));
            ssize_t s = pi_file_write(f, sum[k]->v, Ar->vec_elt_stride(Ar, ii1 - ii0));
            
            if(s < 0 || s < A->vec_elt_stride(A, unpadded)) {
                if (tcan_print) {
                    fprintf(stderr, "ERROR: could not save %s\n", s_name);
                    unlink(s_name);
                }
            }
            pi_file_close(f);
            free(s_name);
        }
        serialize(pi->m);

        /* recover the vector */
        matmul_top_load_vector(mmt, v_name, bw->dir, s + bw->interval, unpadded);

        if (pi->m->trank == 0 && pi->m->jrank == 0)
            keep_rolling_checkpoints(v_name, s + bw->interval);

        serialize(pi->m);

        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, mmt->mm->ncoeffs, tcan_print, 1);
    }
    timing_final_tally(pi, timing, mmt->mm->ncoeffs, tcan_print, 1);

    if (tcan_print) {
        printf("Done.\n");
    }
    serialize(pi->m);

    for(unsigned int k = 0 ; k < multi ; k++) {
        matmul_top_vec_clear_generic(mmt, sum[k], bw->dir);
    }
    if (!bw->skip_online_checks) {
        matmul_top_vec_clear_generic(mmt, check_vector, !bw->dir);
        vec_clear_generic(pi->m, ahead, nchecks);
        free(gxvecs);
    }
    for(unsigned int k = 0 ; k < multi ; k++) {
        vec_clear_generic(pi->m, fcoeffs[k], A->groupsize(A)*bw->interval);
    }
    free(v_name);
    free(fcoeffs);
    free(sum);

    matmul_top_clear(mmt);
    pi_free_mpfq_datatype(pi, A_pi);
    pi_free_mpfq_datatype(pi, Ar_pi);
    pi_free_mpfq_datatype(pi, Ac_pi);

    Ar->oo_field_clear(Ar);
    Ac->oo_field_clear(Ac);
    A->oo_field_clear(A);

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
    /* interpret our parameters: none here (so far). */
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(mksol_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}


