#pragma once

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
extern proj::Lot::SharedPtr rng;

namespace proj {

    struct G {
        // Program settings used in processCommandLineOptions
        static string                   _sim_filename_prefix;
        static unsigned                 _sim_ntaxa;
        static double                   _sim_lambda;
        static double                   _sim_mu;
        static double                   _sim_rho;
        static double                   _sim_root_age;
        static unsigned                 _rnseed;
        static unsigned                 _nthreads;
        static unsigned                 _verbosity;
        static string                   _start_mode;
        static bool                     _save_memory;
        static string                   _filename;
        static unsigned                 _nparticles;
        static string                   _model;
        static double                   _kappa;
        static string                   _string_base_frequencies;
        static string                   _string_relative_rates;
        static string                   _proposal;
        static bool                     _est_lambda;
        static bool                     _est_mu;
        static bool                     _est_root_age;
        static unsigned                 _ngroups;
        static bool                     _run_on_empty;
        
        // Candidates for setting status
        static double                   _asrv_shape;
        static double                   _comphet;
        static double                   _occupancy;
        
        // Dimensions
        static unsigned                 _nstates;
        static unsigned                 _ntaxa;
        static unsigned                 _nloci;
        static vector<unsigned>         _nsites_per_locus;
        
        // Names
        static vector<string>           _taxon_names;
        static vector<string>           _locus_names;

        // Other globals
        static bool                     _simulating;
        static bool                     _debugging;
        static map<unsigned,unsigned>   _nexus_taxon_map;
        static map<unsigned, double>    _relrate_for_locus;
        static double                   _lambda;
        static double                   _mu;
        static double                   _root_age;
        static vector<double>           _base_frequencies;
        static vector<double>           _double_relative_rates;
        static unsigned                 _save_every;
        
#if defined(FOSSILS)
        static vector<Fossil>           _fossils;
        static vector<TaxSet>           _taxsets;
        static vector<Fossil>           _sim_fossils;
#endif
        
        // validation
        static bool                     _ruv;

        // Useful constants
        static double                   _small_enough;
        static double                   _infinity;
        static double                   _negative_infinity;
        
        // Program info
        static string                   _program_name;
        static unsigned                 _major_version;
        static unsigned                 _minor_version;
        
        // Utility functions
        static string   inventName(unsigned k, bool lower_case);
        static double   calcLogSum(const vector<double> & log_values);
        static string   memoryAddressAsString(const void * ptr);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
    };
    
    inline string G::inventName(unsigned k, bool lower_case) {
        // If   0 <= k < 26, returns A, B, ..., Z,
        // If  26 <= k < 702, returns AA, AB, ..., ZZ,
        // If 702 <= k < 18278, returns AAA, AAB, ..., ZZZ, and so on.
        //
        // For example, k = 19009 yields ABCD:
        // ABCD 19009 = 26 + 26*26 + 26*26*26 + 0*26*26*26 + 1*26*26 + 2*26 + 3
        //              <------- base ------>   ^first       ^second   ^third ^fourth
        // base = (26^4 - 1)/25 - 1 = 18278
        //   26^1 + 26^2 + 26^3 = 26^0 + 26^1 + 26^2 + 26^3 - 1 = (q^n - 1)/(q - 1) - 1, where q = 26, n = 4
        //   n = 1 + floor(log(19009)/log(26))
        // fourth = ((19009 - 18278                           )/26^0) % 26 = 3
        // third  = ((19009 - 18278 - 3*26^0                  )/26^1) % 26 = 2
        // second = ((19009 - 18278 - 3*26^0 - 2*26^1         )/26^2) % 26 = 1
        // first  = ((19009 - 18278 - 3*26^0 - 2*26^1 - 1*26^2)/26^3) % 26 = 0
                
        // Find how long a species name string must be
        double logibase26 = (k > 0 ? log(k)/log(26) : 0);
        unsigned n = 1 + (unsigned)floor(logibase26);
        vector<char> letters;
        unsigned base = (unsigned)((pow(26,n) - 1)/25.0 - 1);
        unsigned cum = 0;
        int ordA = (unsigned)(lower_case ? 'a' : 'A');
        for (unsigned i = 0; i < n; ++i) {
            unsigned ordi = (unsigned)((k - base - cum)/pow(26,i)) % 26;
            letters.push_back(char(ordA + ordi));
            cum += (unsigned)(ordi*pow(26,i));
        }
        string species_name(letters.rbegin(), letters.rend());
        return species_name;
    }

    inline double G::calcLogSum(const vector<double> & log_values) {
        double max_logv = *max_element(log_values.begin(), log_values.end());
        
        double factored_sum = 0.0;
        for (auto & logv : log_values) {
            factored_sum += exp(logv - max_logv);
        }
        double log_sum_values = max_logv + log(factored_sum);
        return log_sum_values;
    }
    
    inline string G::memoryAddressAsString(const void * ptr) {
        ostringstream memory_address;
        memory_address << ptr;
        return memory_address.str();
    }
    
    inline unsigned G::multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cumprobs(probs.size());
        partial_sum(probs.begin(), probs.end(), cumprobs.begin());
        assert(fabs(*(cumprobs.rbegin()) - 1.0) < 0.0001);

        // Draw a Uniform(0,1) random deviate
        double u = lot->uniform();

        // Find first element in cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6 <-- u
        //
        // cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cumprobs.begin(), cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == cumprobs.end()) {
            double last_cumprob = *(cumprobs.rbegin());
            throw XProj(format("G::multinomialDraw failed: u = %.9f, last cumprob = %.9f") % u % last_cumprob);
        }

        auto d = distance(cumprobs.begin(), it);
        assert(d >= 0);
        assert(d < probs.size());
        return (unsigned)d;
    }
}

