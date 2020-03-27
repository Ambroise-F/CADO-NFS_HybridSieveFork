#include "cado.h"       /* HAVE_* macros ! */

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cinttypes>
#include <vector>
#include <string>
#include <algorithm>

#include "portability.h"
#include "macros.h"
#include "bblas.hpp"
#include "utils.h"
#include "params.h"

int test_accel = 0;

#include "test_bblas_level2.hpp"
#include "test_bblas_level3.hpp"
#include "test_bblas_level4.hpp"
#include "test_bblas_level5.hpp"

void print_features()/*{{{*/
{
    printf("## compile-time features\n");
#ifdef HAVE_M4RI
    printf("## HAVE_M4RI\n");
#endif /* HAVE_M4RI */
#ifdef HAVE_M4RIE
    printf("## HAVE_M4RIE\n");
#endif /* HAVE_M4RIE */
#ifdef HAVE_PCLMUL
    printf("## HAVE_PCLMUL\n");
#endif /* HAVE_PCLMUL */
#ifdef HAVE_SSE2
    printf("## HAVE_SSE2\n");
#endif /* HAVE_SSE2 */
#ifdef HAVE_SSE41
    printf("## HAVE_SSE41\n");
#endif /* HAVE_SSE41 */
#ifdef VALGRIND
    printf("## VALGRIND\n");
#endif /* VALGRIND */
    printf("## ULONG_BITS=%d\n", ULONG_BITS);
}/*}}}*/

void declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "seed", "seed for random data generation");
    param_list_decl_usage(pl, "n", "n value for some size-dependent tests");
    param_list_decl_usage(pl, "fast", "do quick tests only\n");
    param_list_decl_usage(pl, "tests", "list of tests to perform\n");
}/*}}}*/

int main(int argc, char * argv[])
{
    cxx_param_list pl;
    unsigned int n = 2 * 1000 * 1000;
    int seed = 0;
    param_list_configure_switch(pl, "-fast", &test_accel);
    const char * argv0 = argv[0];
    argc--,argv++;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_uint(pl, "n", &n);
    char ** tests_list;
    int ntests;
    param_list_parse_string_list_alloc(pl, "tests", &tests_list, &ntests, ",");
    std::vector<std::string> tests;
    for(int i = 0 ; i < ntests ; i++) {
        tests.emplace_back(tests_list[i]);
        free(tests_list[i]);
    }
    free(tests_list);
    tests_list = NULL;
    if (tests.empty())
        tests.emplace_back("all");

    if (!seed) seed = rand();

    print_features();

    printf("# seeding with seed %d\n", seed);

    test_bblas_level2 A2(n); A2.set_seed(seed);
    test_bblas_level3 A3(n); A3.set_seed(seed);
    test_bblas_level4 A4(n); A4.set_seed(seed);
    test_bblas_level5 A5(n); A5.set_seed(seed);

    for(auto const & s : tests) {
        bool match = false;
        if (A2(s)) match = true;
        if (A3(s)) match = true;
        if (A4(s)) match = true;
        if (A5(s)) match = true;
        if (!match) {
            fprintf(stderr, "## no test for key %s\n", s.c_str());
        }
    }

    return 0;
}
