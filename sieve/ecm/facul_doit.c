#include "cado.h"
#include "facul_doit.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"

extern unsigned long stats_called[];
extern unsigned long stats_found_n[];

int
primetest (const modulus_t m)
{
  residue_t one, r;
  int isprime;
  
  isprime = mod_sprp2 (m);
  if (isprime)
    {
      mod_init_noset0 (one, m);
      mod_init_noset0 (r, m);
      mod_set1 (one, m);
      mod_add (r, one, one, m);
      mod_add (r, r, one, m);   /* r = 3 */
      isprime = mod_sprp (r, m);
      mod_clear (one, m);
      mod_clear (r, m);
    }
  
  return isprime;
}

int
facul_doit (unsigned long *factors, const modulus_t m, 
	    const facul_strategy_t *strategy, const int method_start)
{
  residue_t r;
  modint_t n, f;
  modulusredcul_t fm_ul, cfm_ul;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
  modulusredc15ul_t fm_15ul, cfm_15ul;
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
  modulusredc2ul2_t fm_2ul2, cfm_2ul2; /* Modulus for factor and cofactor */
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
  modulusmpz_t fm_mpz, cfm_mpz;
#endif
  int i, found = 0, bt, fprime, cfprime;
  /* The f_arith and cf_arith variables tell which modulus type has been 
     initiaised for arithmetic on the factor and the cofactor, resp. 
     The have a value of CHOOSE_NONE if no modulus for factor resp. 
     cofactor is currently initialised. */
  enum {
      CHOOSE_NONE,
      CHOOSE_UL,
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
      CHOOSE_15UL,
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
      CHOOSE_2UL2,
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
      CHOOSE_MPZ,
#endif
  } f_arith = CHOOSE_NONE, cf_arith = CHOOSE_NONE;
  
  mod_intinit (n);
  mod_intinit (f);
  mod_getmod_int (n, m);
  mod_intset_ul (f, 1UL);
  mod_init (r, m);
  
  for (i = method_start; strategy->methods[i].method != 0; i++)
    {
      /* Simple-minded early abort for large input.
         Note: before the test was "mod_intbits (n) > LONG_BIT" which was
         machine-dependent. However it would be better if the early abort
         test depends on the size of the number we are trying to factor,
         since for a large number we can invest more in cofactorization. */
#if 0
      if (i > 3 && mod_intbits (n) > 64)
        break;
#endif
      
      if (i < STATS_LEN)
	  stats_called[i]++;
      
      if (strategy->methods[i].method == PM1_METHOD)
	bt = pm1 (f, m, (pm1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == PP1_METHOD)
	bt = pp1 (f, m, (pp1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == EC_METHOD)
	bt = ecm (f, m, (ecm_plan_t *) (strategy->methods[i].plan));
      else 
	{
	  /* A method value we don't know about. Something's wrong, bail out */
	  found = -1;
	  break;
	}
      
      /* The following possibilities exist:
	 bt:   Factor:    Cofactor:   Action:
	 0           1    composite   Try next method
	 1           1    composite   Could try again with careful bt
	 0    prime<lp     prime<lp   Store both, exit successfully
	 0    prime>lp     prime<lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime>lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 1    prime<lp     prime<lp   Store both, exit successfully
	 1    prime>lp     prime<lp   Not smooth, exit NOT_SMOOTH
	 1    prime<lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 1    prime>lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp    composite   Store prime, continue with cofactor
	 0    prime>lp    composite   Not smooth, exit NOT_SMOOTH
	 1    prime<lp    composite   Store prime, try same method with cofactor
	 1    prime>lp    composite   Not smooth, exit NOT_SMOOTH
	 0   composite     prime<lp   Store prime, continue next method
	 0   composite     prime>lp   Not smooth, exit NOT_SMOOTH
	 1   composite     prime<lp   Store prime, retry this method
	 1   composite     prime>lp   Not smooth, exit NOT_SMOOTH
	 0   composite            1   Could try with lower bounds
	 1   composite            1   Could try again with careful bt

	 Simplified:

	 bt:   Factor:    Cofactor:   Action:
	 0           1    composite   Try next method
	 1           1    composite   Could try again with careful bt
	 ?    prime<lp     prime<lp   Store both, exit successfully
	 ?    prime>lp            ?   Not smooth, exit NOT_SMOOTH
	 ?           ?     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp    composite   Store prime, continue with cofactor
	 1    prime<lp    composite   Store prime, same method with cofactor
	 0   composite     prime<lp   Store prime, continue next method
	 1   composite     prime<lp   Store prime, retry this method
	 0   composite            1   Could try with lower bounds
	 1   composite            1   Could try again with careful bt
	 
      */
      
      
      if (mod_intequal_ul (f, 1UL))
	{
	  if (bt == 0)
	    {
	      /* No factor found, no backtracking... this was a simple miss. */
	      continue;
	    }
	  else
	    {
	      /* Backtracking was used, so there was a point where all
		 factors had been found simultaneously, but backing up
		 to the previous checkpoint resulted in no factors being
		 found. We could try to do some more clever backtracking 
		 to discover the factors yet. TODO. For now, just continue
		 to the next method. */
	      continue;
	    }
	}
      
      if (mod_intequal (f, n))
	{
	  if (i < STATS_LEN)
	    stats_found_n[i]++;
	  if (bt == 0)
	    {
	      /* Input number was found without any backtracking happening?
		 Find out when this can occur and how to get a chance of
		 finding the factors yet. TODO. */
	      continue;
	    }
	  else
	    {
	      /* Backtracking was used, but could not separate the factors,
	         e.g. if both factors are found in stage 1 without 
		 multiplying/exponentiating by 2 at all. Better backtracking
		 might recover the factors yet. TODO. */
	      continue;
	    }
	}
      
      /* So we found a non-trivial factor. See if it is prime, if the 
	 cofactor is prime, and if one of them is, whether they are too
	 large for our smoothness bounds */
      
      /* A quick test if the factor is <= fbb^2 and >lpb */
      /* FIXME: must always use same width for comparison */
      fprime = (mod_intcmp_uint64 (f, strategy->assume_prime_thresh) <= 0); 
      if (fprime && mod_intcmp_ul (f, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Compute the cofactor */
      mod_intdivexact (n, n, f);
      
      /* See if cofactor is <= fbb^2 and > lpb */
      cfprime = (mod_intcmp_uint64 (n, strategy->assume_prime_thresh) <= 0);
      if (cfprime && mod_intcmp_ul (n, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Determine for certain if the factor is prime */
      if (!fprime)
	{
	  const size_t bits = mod_intbits (f);
	  ASSERT(bits <= MOD_MAXBITS);
	  if (bits <= MODREDCUL_MAXBITS)
	    {
	      f_arith = CHOOSE_UL;
	      modredcul_initmod_ul (fm_ul, mod_intget_ul(f));
	      fprime = primetest_ul (fm_ul);
            }
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
          else if (bits <= MODREDC15UL_MAXBITS)
            {
              unsigned long t1[2];
              modintredc15ul_t t2;
              size_t nr_words = mod_intget_uls(t1, f);
              ASSERT_ALWAYS(nr_words <= 2);
              modredc15ul_intset_uls (t2, t1, nr_words);
              f_arith = CHOOSE_15UL;
	      modredc15ul_initmod_int (fm_15ul, t2);
	      fprime = primetest_15ul (fm_15ul);
            }
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
	  else if (bits <= MODREDC2UL2_MAXBITS)
            {
              unsigned long t1[2];
              modintredc2ul2_t t2;
              size_t nr_words = mod_intget_uls(t1, f);
              ASSERT_ALWAYS(nr_words <= 2);
              modredc2ul2_intset_uls (t2, t1, nr_words);
              f_arith = CHOOSE_2UL2;
	      modredc2ul2_initmod_int (fm_2ul2, t2);
	      fprime = primetest_2ul2 (fm_2ul2);
            }
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
	  else if (bits <= MODMPZ_MAXBITS)
            {
              /* We assume for now that f is a modintmpz_t */
              f_arith = CHOOSE_MPZ;
	      modmpz_initmod_int (fm_mpz, f);
	      fprime = primetest_mpz (fm_mpz);
            }
#endif
          else
              abort();

	  if (fprime && mod_intcmp_ul (f, strategy->lpb) > 0)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	      break;
	    }
	}
      
      /* Determine for certain if the cofactor is prime */
      if (!cfprime)
	{
	  const size_t bits = mod_intbits (n);
	  ASSERT (bits <= MOD_MAXBITS);
	  if (bits <= MODREDCUL_MAXBITS)
	    {
	      cf_arith = CHOOSE_UL;
	      modredcul_initmod_ul (cfm_ul, mod_intget_ul(n));
	      cfprime = primetest_ul (cfm_ul);
            }
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
	  else if (bits <= MODREDC15UL_MAXBITS)
	    {
              unsigned long t1[2];
              modintredc15ul_t t2;
              size_t nr_words = mod_intget_uls(t1, n);
              ASSERT_ALWAYS(nr_words <= 2);
              modredc15ul_intset_uls (t2, t1, nr_words);
	      cf_arith = CHOOSE_15UL;
	      modredc15ul_initmod_int (cfm_15ul, t2);
	      cfprime = primetest_15ul (cfm_15ul);
            }
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
	  else if (bits <= MODREDC2UL2_MAXBITS)
            {
              unsigned long t1[2];
              modintredc2ul2_t t2;
              size_t nr_words = mod_intget_uls(t1, n);
              ASSERT_ALWAYS(nr_words <= 2);
              modredc2ul2_intset_uls (t2, t1, nr_words);
              cf_arith = CHOOSE_2UL2;
	      modredc2ul2_initmod_int (cfm_2ul2, t2);
	      cfprime = primetest_2ul2 (cfm_2ul2);
            }
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
	  else if (bits <= MODMPZ_MAXBITS)
            {
              /* We assume for now that f is a modintmpz_t */
              cf_arith = CHOOSE_MPZ;
	      modmpz_initmod_int (cfm_mpz, n);
	      cfprime = primetest_mpz (cfm_mpz);
            }
#endif
          else
            abort ();

	  if (cfprime && mod_intcmp_ul (n, strategy->lpb) > 0)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	      break;
	    }
	}
      
      /* So each of factor and cofactor is either a prime < lpb, 
	 or is composite */

      if (fprime)
	factors[found++] = mod_intget_ul(f); /* f < lpb, so it fits in 1 unsigned long */
      else
	{
            int f2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite factor. Use the same method again so that
	     backtracking can separate the factors */
          switch (f_arith) {
              case CHOOSE_UL:
                  f2 = facul_doit_ul (factors + found, fm_ul, strategy, i);
                  modredcul_clearmod (fm_ul);
                  f_arith = CHOOSE_NONE;
                  break;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
              case CHOOSE_15UL:
                  f2 = facul_doit_15ul (factors + found, fm_15ul, strategy, i);
                  modredc15ul_clearmod (fm_15ul);
                  f_arith = CHOOSE_NONE;
                  break;
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
              case CHOOSE_2UL2:
                  f2 = facul_doit_2ul2 (factors + found, fm_2ul2, strategy, i);
                  modredc2ul2_clearmod (fm_2ul2);
                  f_arith = CHOOSE_NONE;
                  break;
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
              case CHOOSE_MPZ:
                  f2 = facul_doit_mpz (factors + found, fm_mpz, strategy, i);
                  modmpz_clearmod (fm_mpz);
                  f_arith = CHOOSE_NONE;
                  break;
#endif
              default: abort();
          }
          
	  if (f2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += f2;
	}

      if (cfprime)
	factors[found++] = mod_intget_ul(n); /* n < lp, so it fits in 1 unsigned long */
      else
	{
	  int f2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite cofactor */
          switch(cf_arith) {
              case CHOOSE_UL:
                  f2 = facul_doit_ul (factors + found, cfm_ul, strategy, i + 1);
                  modredcul_clearmod (cfm_ul);
                  cf_arith = CHOOSE_NONE;
                  break;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
              case CHOOSE_15UL:
                  f2 = facul_doit_15ul (factors + found, cfm_15ul, strategy, i + 1);
                  modredc15ul_clearmod (cfm_15ul);
                  cf_arith = CHOOSE_NONE;
                  break;
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
              case CHOOSE_2UL2:
                  f2 = facul_doit_2ul2 (factors + found, cfm_2ul2, strategy, i + 1);
                  modredc2ul2_clearmod (cfm_2ul2);
                  cf_arith = CHOOSE_NONE;
                  break;
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
              case CHOOSE_MPZ:
                  f2 = facul_doit_mpz (factors + found, cfm_mpz, strategy, i+1);
                  modmpz_clearmod (cfm_mpz);
                  cf_arith = CHOOSE_NONE;
                  break;
#endif
              default: abort();
          }
          
	  if (f2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += f2;
	}
      
      /* We found a non-trivial factorization and any composite 
	 factors/cofactors have been treated in recursive calls, 
	 so we can stop here */
      ASSERT_ALWAYS(f_arith == CHOOSE_NONE);
      ASSERT_ALWAYS(cf_arith == CHOOSE_NONE);
      break;
    }

  /* The breaks in the loop may leave initialised moduli behind */
  if (f_arith == CHOOSE_UL)
    modredcul_clearmod (fm_ul);
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
  else if (f_arith == CHOOSE_15UL)
    modredc15ul_clearmod (fm_15ul);
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
  else if (f_arith == CHOOSE_2UL2)
    modredc2ul2_clearmod (fm_2ul2);
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
  else if (f_arith == CHOOSE_MPZ)
    modmpz_clearmod (fm_mpz);
#endif

  if (cf_arith == CHOOSE_UL)
    modredcul_clearmod (cfm_ul);
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
  else if (cf_arith == CHOOSE_15UL)
    modredc15ul_clearmod (cfm_15ul);
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
  else if (cf_arith == CHOOSE_2UL2)
    modredc2ul2_clearmod (cfm_2ul2);
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
  else if (cf_arith == CHOOSE_MPZ)
    modmpz_clearmod (cfm_mpz);
#endif
  
  
  mod_clear (r, m);
  mod_intclear (n);
  mod_intclear (f);
  return found;
}
