#ifndef LAS_SIEVER_CONFIG_HPP_
#define LAS_SIEVER_CONFIG_HPP_

#include <cstring>      // for memcmp, NULL
#include <map>           // for operator==, map<>::const_iterator, _Rb_tree_...
#include <tuple>         // for tie, operator<, tuple
#include <utility>       // for pair
#include "fb.hpp"        // for fb_factorbase, fb_factorbase::key_type
#include "las-base.hpp"  // for _padded_pod
#include "params.h"     // param_list_ptr
#include "las-side-config.hpp"

struct las_todo_entry; // IWYU pragma: keep


/* siever_config */
 
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate. This is set by
     * siever_config_pool::get_config_for_q  */
    // unsigned int bitsize __attribute__((deprecated));  /* bitsize == 0 indicates end of table */
    // int side __attribute__((deprecated));              /* special-q side */

    int logA;
    int logI;   /* see below. logI is initialized late in the game */

    /* This does not really belong here. I'd rather have it at the
     * las_info level. However for obscure reasons,
     * las.get_strategies(siever_config&)
     * wants it.
     *
     * FIXME: er. now that get_strategies lives below las, we can do this
     * move, right ?
     */
    unsigned int sublat_bound;


    /* These four parameters are as they are provided in the command
     * line. In truth, the ones that really matter are the ones in the
     * fb_factorbase::key_type object that is stored within the
     * nfs_work structure (in sides[side].fbK).
     *
     * In particular, bucket_thresh and bucket_thresh1 below have a
     * default value that is dependent on:
     *  - logI -- which is not set here. Well, it exists here, but it is
     *    set late in the game, after all other fields
     *    (sieve_range_adjust does that).  And we don't want to go the
     *    route of making the default value for some of the fields in
     *    this struct be dynamic depending on a setter function for logI
     *  - the side, as well.
     */
    private:
    /* access should rather be
     * via ws.sides[side].fbK.thresholds[0,1]
     * and ws.sides[side].fbK.{td_thresh, skipped}
     */
    unsigned long bucket_thresh = 0;  // bucket sieve primes >= bucket_thresh
    unsigned long bucket_thresh1 = 0; // primes above are 2-level bucket-sieved
    unsigned int td_thresh = 1024;
    unsigned int skipped = 0;         // don't sieve below this
    unsigned int no_trial_div = 0; //- disable trial div
    unsigned int small_batch_max_prime = 3;

    public:
    /* the only way to access the four fields above */
    fb_factorbase::key_type instantiate_thresholds(int side) const;

    public:

    /* unsieve threshold is not related to the factor base. */
    unsigned int unsieve_thresh = 100;

    std::vector<siever_side_config> sides;

    void display(int side, unsigned int bitsize) const;

    static void declare_usage(cxx_param_list & pl);
    static bool parse_default(siever_config & sc, cxx_param_list & pl, int);

    /*{{{ has_same_config */
    bool operator==(siever_config const & o) const { return memcmp(this, &o, sizeof(*this)) == 0; }

    struct has_same_config {
        siever_config const & sc;
        has_same_config(siever_config const & sc) : sc(sc) {}
        bool operator()(siever_config const& o) const { return o == sc; }
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
    };
    has_same_config same_config() const { return has_same_config(*this); }
    /*}}}*/
#if 0
    /*{{{ has_same_config_q */
    struct has_same_config_q {
        siever_config const & sc;
        has_same_config_q(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return sc.side == o.side && sc.bitsize == o.bitsize;
        }
    };
    has_same_config_q same_config_q() const {
        return has_same_config_q(*this);
    }
    /*}}}*/
#endif
    /* {{{ has_same_fb_parameters */
    struct has_same_fb_parameters {
        siever_config const & sc;
        has_same_fb_parameters(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            bool ok = true;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lim == o.sides[side].lim;
                ok = ok && sc.sides[side].powlim == o.sides[side].powlim;
            }
            return ok;
        }
    };
    has_same_fb_parameters same_fb_parameters() const { return has_same_fb_parameters(*this); }
    /*}}}*/
#if 0
    /*{{{ has_same_sieving -- currently duplicates has_same_fb_parameters */
    struct has_same_sieving {
        siever_config const & sc;
        has_same_sieving(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return has_same_fb_parameters(sc)(o);
        }
    };
    has_same_sieving same_sieving() const { return has_same_sieving(*this); }
    /*}}}*/
#endif
    /*{{{ has_same_cofactoring */
    struct has_same_cofactoring {
        siever_config const & sc;
        has_same_cofactoring(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const { 
            bool ok = true;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lambda == o.sides[side].lambda;
                ok = ok && sc.sides[side].lpb == o.sides[side].lpb;
                ok = ok && sc.sides[side].mfb == o.sides[side].mfb;
                ok = ok && sc.sides[side].mfbb == o.sides[side].mfbb;
                ok = ok && sc.sides[side].sbmp == o.sides[side].sbmp;
                ok = ok && sc.sides[side].ncurves == o.sides[side].ncurves;
            }
            return ok;
        }
        struct comparison {
            bool operator()(siever_config const& a, siever_config const& b) const {
                return std::tie(
                        a.sides[0].lambda, a.sides[0].lpb, a.sides[0].mfb, a.sides[0].ncurves,
                        a.sides[1].lambda, a.sides[1].lpb, a.sides[1].mfb, a.sides[1].ncurves) <
                    std::tie(
                        b.sides[0].lambda, b.sides[0].lpb, b.sides[0].mfb, b.sides[0].ncurves,
                        b.sides[1].lambda, b.sides[1].lpb, b.sides[1].mfb, b.sides[1].ncurves);
            }
        };

    };
    has_same_cofactoring same_cofactoring() const { return has_same_cofactoring(*this); }
    /*}}}*/
};


/* {{{ descent_hint
 *
 * This is used for the descent. For each factor size, we provide a
 * reasonable siever_config value
 *
 * We also provide, based on experience, info relative to how long it
 * takes to finish the smoothing process for a prime factor of this size.
 */
/* }}} */

struct siever_config_pool {
    typedef std::pair<int, unsigned int> key_type;

    struct descent_hint : public siever_config {
        double expected_time;
        double expected_success;
    };

    typedef std::map<key_type, descent_hint> hint_table_t;
    hint_table_t hints;

    descent_hint const * get_hint(int side, unsigned int bitsize) const {
        hint_table_t::const_iterator it = hints.find(key_type(side, bitsize));
        if (it == hints.end())
            return NULL;
        else
            return &it->second;
    }

    /* The siever_config in [base] needs not be complete. The
     * default_config_ptr field points here if it is complete. If not,
     * the fields here are just used as a base for initializing the other
     * configurations.
     *
     * Note that while the "base for initializing" functionality is
     * likely to stay, the notion of a "default config" seems to be
     * screwed altogether, and we would rather like to see it disappear
     * someday. Currently it is used for displaying memory usage, setting
     * defaults for dupqmin, and getting the lim abd lpb parameters
     * before going to the batch step.
     */

    siever_config const * default_config_ptr;
    siever_config base;

    siever_config get_config_for_q(las_todo_entry const& doing) const;

    siever_config_pool(cxx_param_list& pl, int nb_polys);

    double hint_expected_time(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_time;
    }
    double hint_expected_success(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_success;
    }

    static void declare_usage(cxx_param_list & pl);
};

#endif	/* LAS_SIEVER_CONFIG_HPP_ */
