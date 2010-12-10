/* Usage: cut_n_roots -poly xxx.poly -q0 nnn -q1 nnn -N nnn 
     split [q0, q1[ into maximal intervals of
     N special q's each. */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "modul_poly.h"
#include <gmp.h>
#include "utils.h"

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

void usage(const char *argv0)
{
    fprintf(stderr, "Usage: %s -poly xxx.poly -q0 nnn -q1 nnn -N nnn\n", argv0);
    fprintf(stderr, "  split [q0, q1[ into maximal intervals of N special q's each.\n");
    exit(1);
}

int
main (int argc0, char *argv0[])
{
  const char *polyfilename = NULL;
  cado_poly pol;
  param_list pl;
  uint64_t q0 = 0, q1 = 0, N = 0;
  mpz_t P, Q;
  int argc = argc0;
  char **argv = argv0;

  /* print the command line */
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (int k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");

  param_list_init(pl);
  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    /* Could also be a file */
    FILE * f;
    if ((f = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, f);
      fclose(f);
      argv++,argc--;
      continue;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage(argv0[0]);
  }

  param_list_parse_uint64(pl, "q0", &q0);
  param_list_parse_uint64(pl, "q1", &q1);
  param_list_parse_uint64(pl, "N", &N);
  if (q0 == 0 || q1 ==0 || N == 0) usage(argv0[0]);

  /* check that q1 fits into an unsigned long */
  if (q1 > (uint64_t) ULONG_MAX) {
    fprintf (stderr, "Error, q1=%" PRIu64 " exceeds ULONG_MAX\n", q1);
    exit (EXIT_FAILURE);
  }

  cado_poly_init (pol);
  ASSERT_ALWAYS((polyfilename = param_list_lookup_string(pl, "poly")) != NULL);
  if (polyfilename) param_list_read_file(pl, polyfilename);
  if (!cado_poly_set_plist (pol, pl)) {
    fprintf (stderr, "Error reading polynomial file\n");
    usage(argv[0]);
    exit (EXIT_FAILURE);
  }
  if (pol->skew <= 0.0) {
    fprintf (stderr, "Error, please provide a positive skewness\n");
    exit (EXIT_FAILURE);
  }

  mpz_init_set_ui (P, q0);
  mpz_init(Q);
  mpz_nextprime (P, P);
  mpz_nextprime(Q, P);
  unsigned long fence = 1000000 * (1 + (q0 / 1000000));
  unsigned long totnr = 0;
  unsigned long base = q0;
  while (mpz_cmp_ui (P, q1) < 0) {
    assert(mpz_cmp_ui(P, fence) < 0);
    assert(base < fence);
    unsigned long p = mpz_get_ui (P);
    unsigned long nr = modul_poly_roots (NULL, pol->f, pol->degree, &p);
    totnr+=nr;
    if (totnr > N) {
      unsigned long len = p + 1 - base;
      if (mpz_cmp_ui(Q, fence) >= 0) {
        /* swallow the tail */
        len = fence - base;
        fence += 1000000;
      }
      printf ("%lu,%lu\n", base, len);
      base += len;
      totnr = 0;
    } else if (mpz_cmp_ui(Q, fence) >= 0) {
      unsigned long len = fence - base;
      fence += 1000000;
      printf ("%lu,%lu\n", base, len);
      base += len;
      totnr = 0;
    }
    mpz_set(P, Q);
    mpz_nextprime(Q,Q);
  }
  if (totnr > 0)
    printf ("%lu,%" PRIu64 "\n", base, q1 - base);

  mpz_clear (P);
  mpz_clear (Q);
  cado_poly_clear (pol);
  param_list_clear(pl);
  return 0;
}
