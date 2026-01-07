#pragma once

#include <codecvt>
#include <random>

using boost::filesystem::current_path;
using boost::algorithm::join;
using boost::is_any_of;
using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::store;
using boost::program_options::parse_command_line;
using boost::program_options::parsed_options;
using boost::program_options::parse_config_file;
using boost::program_options::reading_file;
using boost::program_options::notify;

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
extern proj::StopWatch stopwatch;
extern proj::Lot::SharedPtr rng;

namespace proj {

    class Proj {
        public:
            Proj();
            ~Proj();
   
            void processCommandLineOptions(int argc, const char * argv[]);
            void run();
            void initializeParticles();
        
        private:
            void clear();
            
            void debugSaveParticleVectorInfo(string fn, unsigned step);
            
            void simulate();
#if defined (INCREMENT_COMPARISON_TEST)
            void simulateTree();
#endif
            void simulateData(Lot::SharedPtr lot, unsigned startat, unsigned locus_length);
            void simulateSave(string fnprefix);

            void smc();
            void summarizeData(Data::SharedPtr);
            void proposeParticles(unsigned step_number);
            void simProposeParticles(unsigned step_number);
            void proposeParticleRange(unsigned first, unsigned last, unsigned step_number);
            double filterParticles(unsigned step, unsigned group_number);
            void filterParticlesThreading(unsigned step);
            void filterParticlesRange(unsigned first, unsigned last, unsigned g);
            double computeEffectiveSampleSize(const vector<double> & probs) const;
            void writeTreeFile ();
            void writeLogFile();
            void writeLogMarginalLikelihoodFile() const;
            void handleBaseFrequencies();
            void handleRelativeRates();
            void writePartialCount();
            void removeUnecessaryTaxsets();
            void buildNonzeroMap(vector<Particle> & particles, map<const void *, list<unsigned> > & nonzero_map, const vector<unsigned> & nonzeros, vector<unsigned> particle_indices, unsigned start, unsigned end);
        
            Fossil parseFossilDefinition(string & fossil_def);
            TaxSet parseTaxsetDefinition(string & taxset_def);

            vector<Particle>            _particle_vec;
        
            Forest::SharedPtr           _sim_tree;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood;
            vector<Lot::SharedPtr>      _group_rng;
            vector<unsigned>            _indices_to_keep; // indices of particles to write to output files
            vector<pair<double, double>> _hpd_values;
    };

    inline Proj::Proj() {
        clear();
    }

    inline Proj::~Proj() {
    }
    
    inline void Proj::clear() {
        _log_marginal_likelihood = 0.0;
        _data = nullptr;
        _partition.reset(new Partition());
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        variables_map vm;
        
        vector<string> taxsets;
        vector<string> fossils;
        vector<string> sim_fossils;
        
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("lambda",  value(&G::_lambda)->default_value(10.0), "per lineage speciation rate assumed for the Yule model")
        ("root_age",  value(&G::_root_age)->default_value(1.0), "root age for birth death model")
        ("mu",  value(&G::_mu)->default_value(1.0), "per lineage extinction rate assumed for the birth death model")
        ("rnseed",  value(&G::_rnseed)->default_value(1), "pseudorandom number seed")
        ("nthreads",  value(&G::_nthreads)->default_value(1), "number of threads")
        ("startmode",  value(&G::_start_mode)->default_value("smc"), "smc or sim")
        ("filename",  value(&G::_filename), "name of file containing data")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("clock_rate",  value(&G::_clock_rate)->default_value(1.0), "clock rate if fixed, mean clock rate if estimating clock rate")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles)")
        ("savememory",  value(&G::_save_memory)->default_value(false), "save memory by recalculating partials each round")
        ("model",  value(&G::_model)->default_value("JC"), "model to use for likelihood calculations")
        ("verbosity",  value(&G::_verbosity)->default_value(1), "0 (nothing but essential output), 1, 2, ...")
        ("simfnprefix",  value(&G::_sim_filename_prefix), "prefix of files in which to save simulated trees and data if startmode is 'sim' (e.g. specifying 'sim' results in files named 'sim.tre' and 'sim.nex')")
        ("simntaxa",  value(&G::_sim_ntaxa)->default_value(4), "number of taxa to simulate if startmode is 'sim'")
        ("simlambda",  value(&G::_sim_lambda)->default_value(1.0), "true speciation rate for simulating tree under the Yule model if startmode is 'sim'")
        ("simmu",  value(&G::_sim_mu)->default_value(0.0), "true extinction rate for simulating tree under the constant-rates Birth-Death model if startmode is 'sim'")
        ("simrho",  value(&G::_sim_rho)->default_value(1.0), "true extant taxon sampling rate for simulating tree under the constant-rates Birth-Death model if startmode is 'sim'")
        ("simrootage",  value(&G::_sim_root_age)->default_value(1.0), "true root age for simulating tree under the constant-rates Birth-Death model if startmode is 'sim'")
        ("simclockrate",  value(&G::_sim_clock_rate)->default_value(1.0), "true clock rate for simulating tree under constant-rates Birth-Death model if startmode is 'sim'")
        ("kappa",  boost::program_options::value(&G::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&G::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("relative_rates", boost::program_options::value(&G::_string_relative_rates)->default_value("null"), "relative rates by locus")
        ("proposal", boost::program_options::value(&G::_proposal)->default_value("prior-prior"), "prior-prior or prior-post")
        ("estimate_lambda", boost::program_options::value(&G::_est_lambda)->default_value(false), "estimate birth rate")
        ("estimate_mu", boost::program_options::value(&G::_est_mu)->default_value(false), "estimate death rate")
        ("estimate_root_age", boost::program_options::value(&G::_est_root_age)->default_value(false), "estimate root age")
        ("estimate_clock_rate", boost::program_options::value(&G::_est_clock_rate)->default_value(false), "estimate clock rate")
        ("ngroups", boost::program_options::value(&G::_ngroups)->default_value(1.0), "number of subgroups")
        ("save_every", boost::program_options::value(&G::_save_every)->default_value(1.0), "save one out of this number of trees in params and tree files")
        ("run_on_empty", boost::program_options::value(&G::_run_on_empty)->default_value(false), "run without likelihood")
        ("sim_dir", boost::program_options::value(&G::_sim_dir)->default_value("."), "run without likelihood")
        ("fossil",  value(&fossils), "a string defining a fossil, e.g. 'Ursus_abstrusus         1.8–5.3 4.3' (4.3 is time, 1.8-5.3 is prior range)")
        ("taxset",  value(&taxsets), "a string defining a taxon set, e.g. 'Ursinae: Helarctos_malayanus Melursus_ursinus Ursus_abstrusus Ursus_americanus Ursus_arctos Ursus_maritimus Ursus_spelaeus Ursus_thibetanus'")
        ("simfossil",  value(&sim_fossils), "a string defining a fossil, e.g. 'Ursus_abstrusus         1.8–5.3 4.3' (4.3 is time, 1.8-5.3 is prior range)")
        // the following variables relate to validation analyses
        ("ruv", boost::program_options::value(&G::_ruv)->default_value(false), "run ruv analysis")
        ("coverage", boost::program_options::value(&G::_coverage)->default_value(false), "run coverage analysis")
        ;
        
        store(parse_command_line(argc, argv, desc), vm);
        try {
            const parsed_options & parsed = parse_config_file< char >("smctree.conf", desc, false);
            store(parsed, vm);
        }
        catch(reading_file & x) {
            throw XProj("Configuration file (smctree.conf) not found\n");
        }
        notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            output(format("%s\n") % desc, 2);
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            output(format("This is %s version %d.%d\n") % G::_program_name % G::_major_version % G::_minor_version, 1);
            exit(1);
        }
        
        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            _partition.reset(new Partition());
            for (auto s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }
        
        // If user specified "base_frequencies" in conf file, convert them to a vector<double>
        if (vm.count("base_frequencies") > 0) {
            handleBaseFrequencies();
        }
        
        // If user specified "relative_rates" in conf file, convert them to a vector<double>
        if (vm.count("relative_rates") > 0) {
            handleRelativeRates();
        }
        
        // If user specified a start mode other than "prior-prior" or "prior-post", throw an exception
        if (G::_proposal != "prior-prior" && G::_proposal != "prior-post") {
            throw XProj(format("must specify a proposal of either prior-prior or prior-post but %s was specified")%G::_proposal);
        }
        
        // If user specified lambda <= 0, throw an exception
        if (G::_lambda <= 0.0) {
            throw XProj(format("must specify a value of lambda > 0 but %d was specified")%G::_lambda);
        }
        
        // If user specified mu < 0, throw an exception
        if (G::_mu < 0.0) {
            throw XProj(format("must specify a value of mu >= 0 but %d was specified")%G::_mu);
        }
        
        // If user specified estimate_mu but not estimate_lambda or vice versa, throw an exception
        if ((G::_est_mu && !G::_est_lambda) || (G::_est_lambda && !G::_est_mu && G::_mu != 0.0)) {
            throw XProj(format("must estimate both lambda and mu or neither under birth-death model"));
        }
        
        // If user specified root age <= 0, throw an exception
        if (G::_root_age <= 0.0) {
            throw XProj(format("must specify a value of root age >= 0 but %d was specified")%G::_root_age);
        }
        
        // If user specified mu = 0 and estimate_mu, print a warning that Yule model will be chosen and mu will not be estimated
        if (G::_mu == 0.0 && G::_est_mu) {
            output(format("warning: extinction rate specified to be estimated but mean rate set to %d; extinction rate will be set to 0 and Yule model will be used\n")%G::_mu, 1);
        }
        
        // If user specified --fossil on command line, break specified
        // fossil definition into species name, taxset name, and age
        if (vm.count("fossil") > 0) {
            G::_fossils.clear();
            for (auto fdef : fossils) {
                G::_fossils.push_back(parseFossilDefinition(fdef));
            }
        }
        
        // If user specified --simfossil on command line, break specified
        // fossil definition into species name, taxset name, and age
        if (vm.count("simfossil") > 0) {
            G::_fossils.clear();
            for (auto fdef : sim_fossils) {
                G::_fossils.push_back(parseFossilDefinition(fdef));
            }
        }
        
        // If user specified --taxset on command line, break specified
        // taxset definition into name and species included
        if (vm.count("taxset") > 0) {
            G::_taxsets.clear();
            for (auto tdef : taxsets) {
                G::_taxsets.push_back(parseTaxsetDefinition(tdef));
            }
        }
        
        if (G::_model == "JC") {
            G::_model_type = G::ModelType::MODEL_TYPE_JC;
        }
        else if (G::_model == "HKY") {
           G:: _model_type = G::ModelType::MODEL_TYPE_HKY;
        }
    }

    // Note: This code is unduly complex because of the fact that there are
    // many different kinds of dashes that can be used in a utf8-encoded text file.
    // One solution is to eliminate the dash separating the lower and upper
    // bounds for the age range, but I felt that using a dash makes it clear that
    // it is an age range and I didn't want to put the onus on the user to use one
    // particular kind of dash (especially since I could not figure out how to
    // use the right kind of dash myself in a text file created using BBEdit).
    
    // The ws_to_utf8 and utf8_to_ws functions below are
    // slightly modified from Galik's answer at
    // stackoverflow.com/questions/43302279
    //   /any-good-solutions-for-c-string-code-point-and-code-unit
    //   /43302460#43302460
    string ws_to_utf8(wstring const & s) {
        wstring_convert<codecvt_utf8<wchar_t>, wchar_t> cnv;
        string utf8 = cnv.to_bytes(s);
        if(cnv.converted() < s.size())
            throw XProj("incomplete conversion to utf8");
        return utf8;
    }

    wstring utf8_to_ws(string const & utf8) {
        wstring_convert<codecvt_utf8<wchar_t>, wchar_t> cnv;
        wstring s = cnv.from_bytes(utf8);
        if(cnv.converted() < utf8.size())
            throw XProj("incomplete conversion to wstring");
        return s;
    }

    inline Fossil Proj::parseFossilDefinition(string & fossil_def) {
        // Examples showing range of possible inputs:
        //   fossil_def = "Ursus_abstrusus 1.8–5.3 4.3"
        //   fossil_def = "Parictis_montanus      33.90  – 37.20  36.6"

        // Vector to hold strings obtained by parsing fossil_def
        vector<string> v;
        
        // Separate fossil_def into 4 strings:
        //  v[0] = fossil species name
        //  v[1] = fossil age lower bound
        //  v[2] = fossil age upper bound
        //  v[3] = fossil age

        // regex_pattern specifies string of characters that do not
        // represent a tab (\u0009), space (\u0020), or dash of any kind:
        //  \u002D Hyphen-minus
        //  \u058A Armenian Hyphen
        //  \u05BE Hebrew Punctuation Maqaf
        //  \u2010 Hyphen
        //  \u2011 Non-Breaking Hyphen
        //  \u2012 Figure Dash
        //  \u2013 En dash
        //  \u2014 Em dash
        //  \u2015 Horizontal bar
        //  \u2E3A Two-Em Dash
        //  \u2E3B Three-Em Dash
        //  \uFE58 Small Em Dash
        //  \uFE63 Small Hyphen-Minus
        //  \uFF0D Fullwidth Hyphen-Minus
        wregex regex_pattern(utf8_to_ws("[^\\u0020\\u0009\\u002D\\u058A\\u05BE\\u2010\\u2011\\u2012\\u2013\\u2014\\u2015\\u2E3A\\u2E3B\\uFE58\\uFE63\\uFF0D]+"));

        // 0 means keep whole match:
        //     1 would mean keep first match only
        //     2 would mean keep second match only
        //     {1,2} would mean keep first and secod matches only
        //    -1 means use regex_pattern as delimiter
        wstring ws_fossil_def = utf8_to_ws(fossil_def);
        wsregex_token_iterator regex_iter(ws_fossil_def.begin(), ws_fossil_def.end(), regex_pattern, 0);

        // regex_end is, by default, signifies the end of the input
        wsregex_token_iterator regex_end;

        // iterate to get matches
        for ( ; regex_iter != regex_end; ++regex_iter) {
            wstring s = *regex_iter;
            v.push_back(ws_to_utf8(s));
        }

        string fossil_species_name = v[0];
        string fossil_age_lower_str  = v[1];
        string fossil_age_upper_str  = v[2];
        string fossil_age_str   = v[3];
        
        double fossil_age_lower;
        try {
            fossil_age_lower = stof(fossil_age_lower_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age lower bound \"%s\" to a floating point value") % fossil_age_lower_str);
        }

        double fossil_age_upper;
        try {
            fossil_age_upper = stof(fossil_age_upper_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age upper bound \"%s\" to a floating point value") % fossil_age_upper_str);
        }
        
        double fossil_age;
        try {
            fossil_age = stof(fossil_age_str);
        }
        catch (const std::invalid_argument& ia) {
            throw XProj(format("Could not convert fossil age \"%s\" to a floating point value") % fossil_age_str);
        }
        
        return Fossil(fossil_species_name, fossil_age_lower, fossil_age_upper, fossil_age);
    }
    
    inline TaxSet Proj::parseTaxsetDefinition(string & taxset_def) {
        // Example:
        //   taxset_def = "Ursinae = Helarctos_malayanus Melursus_ursinus Ursus_abstrusus Ursus_americanus Ursus_arctos Ursus_maritimus Ursus_spelaeus Ursus_thibetanus;"
        vector<string> v;
        
        // Separate taxset_def into 2 strings at equal sign:
        //  v[0] = taxset name
        //  v[1] = list of species in taxset
        split(v, taxset_def, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj(format("Expecting exactly 2 items separated by colon (:) in taxset definition but instead found %d") % v.size());
        string taxset_name = v[0];
        string taxset_species_list  = v[1];
        
        // Trim whitespace from both ends
        trim(taxset_name);
        trim(taxset_species_list);
        
        // Separate taxset_species_list into strings separated by spaces
        
        // regex_pattern specifies string of characters not a hyphen or whitespace
        regex regex_pattern("[^\\s]+");
        
        // 0 means keep whole match:
        //     1 would mean keep first match only
        //     2 would mean keep second match only
        //     {1,2} would mean keep first and secod matches only
        //    -1 means use regex_pattern as delimiter
        sregex_token_iterator regex_iter(taxset_species_list.begin(), taxset_species_list.end(), regex_pattern, 0);
        
        // regex_end is, by default, signifies the end of the input
        sregex_token_iterator regex_end;
        
        v.clear();
        for ( ; regex_iter != regex_end; ++regex_iter) {
            v.push_back(*regex_iter);
        }
        return TaxSet(taxset_name, v);
    }

    inline void Proj::handleBaseFrequencies() {
        vector <string> temp;
        split(temp, G::_string_base_frequencies, is_any_of(","));
        double sum = 0.0;
        // iterate through temp
        for (auto &i:temp) {
            double f = stof(i);
            G::_base_frequencies.push_back(f);
            sum +=f;
        }
        if (fabs(sum-1)>0.000001) {
            throw XProj(format("base frequencies (%s) don't add to 1")%G::_string_base_frequencies);
        }
        assert (fabs(sum-1) < 0.000001);
    }

    inline void Proj::handleRelativeRates() {
        vector <string> temp;
        assert (G::_double_relative_rates.size() == 0);
        unsigned nloci = _partition->getNumSubsets();
        
        split(temp, G::_string_relative_rates, is_any_of(","));
        
        if (G::_string_relative_rates == "null") {
            for (unsigned i=0; i<nloci; i++) {
                G::_double_relative_rates.push_back(1.0); // if no relative rates provided, set them all to 1
            }
        }
        else {
        // iterate through temp
            for (auto &i:temp) {
                double f = stof(i);
                G::_double_relative_rates.push_back(f);
            }
        }
        
        if (G::_double_relative_rates.size() != nloci) {
            throw XProj(format("%d relative rates were provided, but there are %d loci")% G::_double_relative_rates.size() % nloci);
        }
        
        // relative rates must average to 1.0
        double average = 0.0;
        for (auto &i:G::_double_relative_rates) {
            average += i;
        }
        average /= G::_double_relative_rates.size();
        
        if (fabs(average-1)>0.000001) {
            throw XProj(format("relative rates (%s) don't average to 1")%G::_string_relative_rates);
        }
        assert (fabs(average-1) < 0.000001);
    }

    inline void Proj::run() {
        rng->setSeed(G::_rnseed);

        if (G::_start_mode == "sim") {
            simulate();
        }
        else {
            smc();
        }
    }

    inline void Proj::simulate() {
        // Sanity checks
        assert(G::_sim_ntaxa > 0);
        
        // Simulate the tree
#if defined (INCREMENT_COMPARISON_TEST)
        simulateTree();
#endif
        
        // create vector of particles
        G::_nparticles = 100;
        G::_ngroups = 1;
        G::_est_lambda = false;
        G::_est_root_age = false;
        G::_est_mu = false;
        unsigned nleaves = G::_sim_ntaxa;
        G::_ntaxa = nleaves;
        G::_root_age = G::_sim_root_age;
        
        // make up taxon names
        G::_taxon_names.resize(nleaves);
        for (unsigned i = 0; i < nleaves; i++) {
            G::_taxon_names[i] = G::inventName(i, /*lower_case*/false);
        }
        
        if (G::_fossils.size() > 0) {
            // check that no fossil is older than the mean root age
            for (auto &f:G::_fossils) {
                if (f._upper >= G::_root_age) {
                    throw XProj("fossil range cannot exceed or equal root age");
                }
            }
            
            // sort fossils from youngest to oldest for use in later proposal
            sort(G::_fossils.begin(), G::_fossils.end(), [](Fossil & left, Fossil & right) {
                return left._age < right._age;
            });
            
            for (auto &f:G::_fossils) {
                if (G::_root_age < f._age) {
                    throw XProj(format("Root age set to %d but oldest fossil has age %d; root age must be older than fossils")% G::_root_age % f._age);
                }
            }
        }
        
        // simulate many particles, then choose among the survivors
        // must do this to simulate fossils
                
        // set group rng
        G::_ngroups = 1.0;
        _group_rng.resize(G::_ngroups);
        unsigned psuffix = 1;
        for (auto &g:_group_rng) {
            g.reset(new Lot());
            g->setSeed(rng->randint(1,9999)+psuffix);
            psuffix += 2;
        }
        
        _particle_vec.resize(G::_nparticles);
        initializeParticles();
        
        unsigned nsteps = (G::_ntaxa-1);
        
        for (unsigned n=0; n<nsteps; n++) {
            simProposeParticles(n);
            filterParticles(n, 0);
            G::_step++;
        }
        
        // save only one particle that survived
        if (_particle_vec.size() > 1) {
            _particle_vec.erase(_particle_vec.begin() + 1, _particle_vec.end());
        }
        
        // Interrogate _partition to determine number of genes, gene names, and
        // number of sites in each gene
        G::_nloci = _partition->getNumSubsets();
        G::_nsites_per_locus.resize(G::_nloci);
        G::_locus_names.resize(G::_nloci);
        unsigned total_nsites = 0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            unsigned locusnsites = _partition->numSitesInSubset(g);
            total_nsites += locusnsites;
            string locusname = _partition->getSubsetName(g);
            G::_nsites_per_locus[g] = locusnsites;
            G::_locus_names[g] = locusname;
        }

        // Create data object
        assert(!_data);
        _data = Data::SharedPtr(new Data());
        _data->setTaxonNames(G::_taxon_names);
        _data->setPartition(_partition);
        
        // Simulate data for each locus given the tree
        unsigned starting_site = 0;
        for (unsigned g = 0; g < G::_nloci; ++g) {
            _particle_vec[0].simulateData(::rng, _data, starting_site, G::_nsites_per_locus[g]);
            starting_site += G::_nsites_per_locus[g];
        }
        
        // Save simulated data to file
        simulateSave(G::_sim_filename_prefix);
    }

#if defined (INCREMENT_COMPARISON_TEST)
    inline void Proj::simulateTree() {
        // this function is only for the increment comparison test
        
        // Set global _ntaxa
        unsigned nleaves = G::_sim_ntaxa;
        G::_ntaxa = nleaves;
        
        // Make up taxon names
        G::_taxon_names.resize(nleaves);
        for (unsigned i = 0; i < nleaves; i++)
            G::_taxon_names[i] = G::inventName(i, /*lower_case*/false);
        
        _sim_tree = Forest::SharedPtr(new Forest);
        
        _sim_tree->incrementComparisonTest();
    }
#endif
        
    inline void Proj::simulateSave(string fnprefix) {
        // Show simulation settings to user
        string locus_or_loci = (G::_nloci != 1 ? "loci" : "locus");
        output(format("Simulating data for %d taxa and %d %s\n") % G::_sim_ntaxa % G::_nloci % locus_or_loci, 0);
        output("  Locus names and lengths:\n", 0);
        for (unsigned i = 0; i < G::_nloci; ++i) {
            output(format("    %s (%d sites)\n") % G::_locus_names[i] % G::_nsites_per_locus[i], 0);
        }
        output(format("  True speciation rate: %g\n") % G::_sim_lambda, 0);
        output(format("  True extinction rate: %g\n") % G::_sim_mu, 0);
        output(format("  True sampling rate: %g\n") % G::_sim_rho, 0);
        output(format("  True root age: %g\n") % G::_sim_root_age, 0);

        string tree_file_name = fnprefix + ".tre";
        ofstream treef(tree_file_name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree true = [&R] " << _particle_vec[0].makeNewick(/*precision*/9, /*use_names*/true) << ";\n";
        treef << "end;\n\n";
        treef.close();

        // Output data to file
        string data_file_name = fnprefix + ".nex";
        _data->compressPatterns();
        _data->writeDataToFile(data_file_name);
        output(format("  Sequence data saved in file \"%s\"\n") % data_file_name, 0);
        
        string height_first_splitf = "height_of_first_split.txt";
        ofstream splitf(height_first_splitf);
        double first_split_height = _particle_vec[0].getHeightFirstSplit();
        splitf << first_split_height << endl;

        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        string paup_command_file_name = fnprefix + "-paup.nex";
        output(format("  PAUP* commands saved in file \"%s\"\n") % paup_command_file_name, 1);
        ofstream paupf(paup_command_file_name);
        paupf << "#NEXUS\n\n";
        paupf << "begin paup;\n";
        paupf << "  log start file=pauplog.txt replace;\n";
        paupf << "  exe " << data_file_name << ";\n";
        paupf << "  set crit= like;\n";
        paupf << "  lset nst=1 basefreq=equal rates=equal pinvar=0 clock userbrlen;\n";
        paupf << "  gettrees file=" << tree_file_name << ";\n";
        paupf << "[!\n";
        paupf << "*********** true tree ***********]\n";
        paupf << "  describe 1 / plot=none brlens=sumonly;\n";
        paupf << "  lset nouserbrlen;\n";
        paupf << "  hsearch;\n";
        paupf << "  savetrees file=mltree.tre brlen format=altnexus;\n";
        paupf << "[!\n";
        paupf << "*********** ML tree ***********]\n";
        paupf << "  describe 1 / plot=none brlens=sumonly;\n";
        paupf << "[!\n";
        paupf << "*********** Tree distances ***********]\n";
        paupf << "  lset userbrlen;\n";
        paupf << "  gettrees file=sim.tre mode=3;\n";
        paupf << "  gettrees file=mltree.tre mode=7;\n";
        paupf << "  [gettrees file=beast/summary.tree mode=7;]\n";
        paupf << "  treedist reftree=1 measure=pathlength file=pldist.txt replace;\n";
        paupf << "  treedist reftree=1 measure=rfSymDiff file=rfdist.txt replace;\n";
        paupf << "  log stop;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
    }
    
    inline void Proj::debugSaveParticleVectorInfo(string fn, unsigned step) {
        ofstream outf(fn, ios::out | ios::app);
        outf << str(format("\n********** particles step %d **********\n\n") % step);
        unsigned i = 1;
        for (auto & p : _particle_vec) {
            outf << p.debugSaveParticleInfo(i++);
        }
        outf.close();
    }

    inline void Proj::smc() {
        output("\nStarting....\n",1);
        output(format("Current working directory: %s\n") % boost::filesystem::current_path(), 1);
        output(format("Random seed: %d\n") % G::_rnseed, 1);
        output(format("Number of threads: %d\n") % G::_nthreads, 1);
        
        // save indices of particles to keep when thinning output
        for (unsigned p=0; p<G::_ngroups * G::_nparticles; p++) {
            _indices_to_keep.push_back(p);
        }
        
        std::shuffle(_indices_to_keep.begin(), _indices_to_keep.end(), std::default_random_engine(G::_rnseed)); // shuffle particle indices
        
        // save first sample_size of indices
        unsigned sample_size = round(double (G::_nparticles * G::_ngroups) / double(G::_save_every) );
        if (sample_size == 0) {
            sample_size = G::_nparticles * G::_ngroups;
        }
        
        unsigned n_indices_to_delete = (unsigned) _indices_to_keep.size() - sample_size;
        
        // delete last n_indices_to_delete elements of vector
        _indices_to_keep.erase(_indices_to_keep.end() - n_indices_to_delete, _indices_to_keep.end());
        assert (_indices_to_keep.size() == sample_size);
        
        try {
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(G::_filename);
            
            G::_nloci = _data->getNumSubsets();
            assert(G::_nloci > 0);
            
            // Copy taxon names to global variable in G
            G::_ntaxa = _data->getNumTaxa();
            _data->copyTaxonNames(G::_taxon_names);
            
            if (G::_verbosity > 0) {
                summarizeData(_data);
            }
            
            assert (G::_nloci > 0);
            ps.setNLoci(G::_nloci);
            ps.setNElements(G::_nstates * _data->getNumPatterns());
//            for (unsigned locus = 1; locus < G::_nloci+1; locus++) {
//                // Set length of partials for gene g
//                ps.setNElements(G::_nstates*_data->getNumPatternsInSubset(locus-1));
//            }
                        
            // create vector of particles
            _particle_vec.resize(G::_nparticles * G::_ngroups);

            for (unsigned i=0; i<G::_nparticles * G::_ngroups; i++) {
                _particle_vec[i] = Particle();
            }
            
            // set group rng
            _group_rng.resize(G::_ngroups);
            unsigned psuffix = 1;
            for (auto &g:_group_rng) {
                g.reset(new Lot());
                g->setSeed(rng->randint(1,9999)+psuffix);
                psuffix += 2;
            }
            

            if (G::_fossils.size() > 0) {
                // check that no fossil is older than the mean root age
                for (auto &f:G::_fossils) {
                    if (f._upper >= G::_root_age) {
                        throw XProj("fossil range cannot exceed or equal root age");
                    }
                }
                
                // sort fossils from youngest to oldest for use in later proposal
                sort(G::_fossils.begin(), G::_fossils.end(), [](Fossil & left, Fossil & right) {
                    return left._age < right._age;
                });
                
                for (auto &f:G::_fossils) {
                    if (G::_root_age < f._age) {
                        throw XProj(format("Root age set to %d but oldest fossil has age %d; root age must be older than fossils")% G::_root_age % f._age);
                    }
                }
            }

            G::_step = 0;
            initializeParticles(); // TODO: make one template particle and copy it
            
            //debugSaveParticleVectorInfo("debug-initialized.txt", 0);

            // reset marginal likelihood
            _log_marginal_likelihood = 0.0;
            vector<double> starting_log_likelihoods = _particle_vec[0].calcGeneTreeLogLikelihoods();
            
            _log_marginal_likelihood = 0.0;
            for (auto &l:starting_log_likelihoods) {
                _log_marginal_likelihood += l;
            }
            
            _log_marginal_likelihood *= G::_ngroups;
            
            // initialize starting log likelihoods for all other particles
            // necessary for calculating the first weight
            
            for (unsigned p=1; p<_particle_vec.size(); p++) {
                _particle_vec[p].setStartingLogLikelihoods(starting_log_likelihoods);
            }
            
            if (G::_save_memory) {
                for (auto &p:_particle_vec) {
                    p.clearPartials();
                }
            }
            
            unsigned nsteps = (G::_ntaxa-1);
            
            for (unsigned g=0; g<nsteps; g++){
                G::_step = g;
                
                // set random number seeds
                unsigned psuffix = 1;
                for (auto &p:_particle_vec) {
                    p.setSeed(rng->randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                unsigned step_plus_one = g+1;
                output(format("Step %d of %d.\n") % step_plus_one % nsteps, 1);
                
                proposeParticles(g);
                
# if defined (INCREMENT_COMPARISON_TEST)
                
                if (G::_step == 2) {
                    ofstream outFile("test3.log", std::ios::app);
                    unsigned count = 0;
                    outFile << "sample" << "\t" << "increment";
                    outFile << endl;
                    for (auto &p:_particle_vec) {
//                        outFile << count << "\t" << p.getHeightFirstSplit() << endl;
//                        outFile << count << "\t" << p.getHeightSecondIncr() << endl;
                        outFile << count << "\t" << p.getHeightThirdIncr() << endl;
                        count++;
                    }
                }
#endif
                
//                debugSaveParticleVectorInfo("debug-proposed.txt", g+1);
                
                if (G::_nthreads == 1) {
                    for (unsigned a=0; a<G::_ngroups; a++) {
                        double ess = filterParticles(g, a);
                        output(format("     ESS = %d\n") % ess, 2);
                    }
                    
                    //debugSaveParticleVectorInfo("debug-filtered.txt", g+1);
                }
                else {
                    filterParticlesThreading(g);
                }
            }
            output(format("     log marginal likelihood = %d\n") % _log_marginal_likelihood, 2);
            
            writeTreeFile();
            writePartialCount();
            writeLogFile(); // TODO: fix this
            writeLogMarginalLikelihoodFile();
            
            if (G::_ruv) {
                // get height of first split from sim directory (read in as command line option)
                string splitfname;
                if (G::_sim_dir != ".") {
                    splitfname = G::_sim_dir + "height_of_first_split.txt";
                }
                else {
                    splitfname = "height_of_first_split.txt";
                }
                ifstream inputf(splitfname);
                double true_split_height = 0.0;
                string line;
                while (getline(inputf, line)) {
                    true_split_height = stod(line);
                    break;
                }
                
                vector<pair<double, bool>> first_split_heights;
                for (auto &p:_particle_vec) {
                    double first_split_height = p.getHeightFirstSplit();
                    first_split_heights.push_back(make_pair(first_split_height, false));
                }
                
                first_split_heights.push_back(make_pair(true_split_height, true)); // true first split height
                
                // sort split heights
                sort(first_split_heights.begin(), first_split_heights.end());

                // find rank of truth
                auto it = std::find_if(first_split_heights.begin(), first_split_heights.end(), [&](const pair<double, bool>& p) { return p.second == true;});
                unsigned index_value = (unsigned) std::distance(first_split_heights.begin(), it);
                
                // write rank value to file
                
                ofstream rankf("rank_first_split.txt");
                rankf << "rank: " << index_value << endl;
            }
            
            if (G::_coverage) {
                assert (_hpd_values.size() > 0);
                
                // get height of first split from sim directory (read in as command line option)
                string splitfname;
                if (G::_sim_dir != ".") {
                    splitfname = G::_sim_dir + "height_of_first_split.txt";
                }
                else {
                    splitfname = "height_of_first_split.txt";
                }
                ifstream inputf(splitfname);
                double true_split_height = 0.0;
                string line;
                while (getline(inputf, line)) {
                    true_split_height = stod(line);
                    break;
                }
                
                // get average first split height
                double total_first_split_heights = 0.0;
                for (auto &p:_particle_vec) {
                    total_first_split_heights += p.getHeightFirstSplit();
                }

                double observed_mean = total_first_split_heights /= _particle_vec.size();
                
                ofstream hpdf("hpd.txt");
                hpdf << "min    " << "max   " << "true    " << "observed mean   " << endl;

                // sort hpd values largest to smallest
                std::sort(_hpd_values.begin(), _hpd_values.end());
                std::reverse(_hpd_values.begin(), _hpd_values.end());

                // take first 95% of values (round down to nearest integer)
                double total = boost::size(_hpd_values);
                double ninety_five_index = floor(0.95*total);

                if (ninety_five_index == 0) {
                    ninety_five_index = 1;
                }

                vector<double> hpd_values_in_range;

                for (unsigned h=0; h<ninety_five_index; h++) {
                    hpd_values_in_range.push_back(_hpd_values[h].first);
                }

                auto max = *std::max_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                auto min = *std::min_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                assert (min < max || min == max);

                // write min and max to file
                hpdf << min << "\t" << max << "\t" << true_split_height << "\t" << observed_mean << endl;
            }
            
        }
        
        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }
    }

    inline double Proj::filterParticles(unsigned step, unsigned group_number) {
        // Copy log weights for all particles to probs vector
        vector<double> probs(G::_nparticles, 0.0);
        unsigned start = group_number * G::_nparticles;
        unsigned end = start + (G::_nparticles) - 1;
        
        vector<unsigned> particle_indices;
        for (unsigned a = 0; a < G::_nparticles * G::_ngroups; a++) {
            particle_indices.push_back(a);
        }
//        for (unsigned a=start; a <end+1; a++) {
//            particle_indices.push_back(a); // TODO: check this
//        }
        
        unsigned prob_count = 0;
        for (unsigned p=start; p < end + 1; p++) {
            probs[prob_count] = _particle_vec[p].getLogWeight();
            prob_count++;
        }
        
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(probs);
        
        // if all weights are -inf, exit
        if (log_sum_weights != log_sum_weights) {
            throw XProj("All particles have violated fossil constraints; try again with more particles");
        }
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood due to this step
        _log_marginal_likelihood += log_sum_weights - log(G::_nparticles);
        
        double ess = 0.0;
        if (G::_verbosity > 1) {
            // Compute effective sample size
            ess = computeEffectiveSampleSize(probs);
        }
        
#if defined (SYSTEMATIC_FILTERING)
        vector<unsigned> zeros;
        zeros.reserve(G::_nparticles);
        vector<unsigned> nonzeros;
        nonzeros.reserve(G::_nparticles);
                  
        // Zero vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (G::_nparticles, 0);
        
        double cump = probs[0];
        double delta = _group_rng[group_number]->uniform() / G::_nparticles;
        unsigned c = (unsigned)(floor(1.0 + G::_nparticles*(cump - delta)));
        if (c > 0) {
            nonzeros.push_back(0);
        }
        else {
            zeros.push_back(0);
        }
        counts[0] = c;
        unsigned prev_cum_count = c;
        for (unsigned i = 1; i < G::_nparticles; ++i) {
            cump += probs[i];
            double cum_count = floor(1.0 + G::_nparticles*(cump - delta));
            if (cum_count > G::_nparticles) {
                cum_count = G::_nparticles;
            }
            unsigned c = (unsigned)cum_count - prev_cum_count;
            if (c > 0) {
                nonzeros.push_back(i);
            }
            else {
                zeros.push_back(i);
            }
            counts[i] = c;
            prev_cum_count = cum_count;
        }
        
//        unsigned locus = _particle_vec[particle_indices[start]].getNextGene() - 1; // subtract 1 because vector of gene forests starts at 0
//        assert (locus == _particle_vec[particle_indices[end]].getNextGene() - 1);
        // Create map (nonzero_map) in which the key for an element
        // is the memory address of a gene forest and
        // the value is a vector of indices of non-zero counts.
        // This map is used to determine which of the nonzeros
        // that need to be copied (last nonzero count for any
        // memory address does not need to be copied and can be
        // modified in place).
        map<const void *, list<unsigned> > nonzero_map;
        buildNonzeroMap(_particle_vec, nonzero_map, nonzeros, particle_indices, start, end);

        
        // Example of following code that replaces dead
        // particles with copies of surviving particles:
        //             0  1  2  3  4  5  6  7  8  9
        // _counts  = {0, 2, 0, 0, 0, 8, 0, 0, 0, 0}  size = 10
        // zeros    = {0, 2, 3, 4, 6, 7, 8, 9}        size =  8
        // nonzeros = {1, 5}                          size =  2
        //
        //  next_zero   next_nonzero   k   copy action taken
        //  --------------------------------------------------------------
        //      0             0        0   _particles[1] --> _particles[0]
        //  --------------------------------------------------------------
        //      1             1        0   _particles[5] --> _particles[2]
        //      2             1        1   _particles[5] --> _particles[3]
        //      3             1        2   _particles[5] --> _particles[4]
        //      4             1        3   _particles[5] --> _particles[6]
        //      5             1        4   _particles[5] --> _particles[7]
        //      6             1        5   _particles[5] --> _particles[8]
        //      7             1        6   _particles[5] --> _particles[9]
        //  --------------------------------------------------------------
        unsigned next_zero = 0;
        unsigned next_nonzero = 0;
        while (next_nonzero < nonzeros.size()) {
            double index_survivor = nonzeros[next_nonzero];
            
            unsigned index_survivor_in_particles = particle_indices[index_survivor + start];
            
            // TODO: check this works with multiple groups
            _particle_vec[index_survivor_in_particles].finalizeLatestJoin(index_survivor_in_particles, nonzero_map);

            
            unsigned ncopies = counts[index_survivor] - 1;
            for (unsigned k = 0; k < ncopies; k++) {
                double index_nonsurvivor = zeros[next_zero++];

                // Replace non-survivor with copy of survivor
//                unsigned survivor_index_in_particles = index_survivor+start;
//                unsigned non_survivor_index_in_particles = index_nonsurvivor+start;

                unsigned survivor_index_in_particles = particle_indices[index_survivor+start];
                unsigned non_survivor_index_in_particles = particle_indices[index_nonsurvivor+start];
                
                _particle_vec[non_survivor_index_in_particles] = _particle_vec[survivor_index_in_particles];
            }
            ++next_nonzero;
        }
        
        return ess;
#else

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts(G::_nparticles, 0);

        // Throw _nparticles darts
        for (unsigned i=0; i<G::_nparticles; i++) {
            double u = _group_rng[group_number]->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)std::distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Copy particles

        int donor = -1;
        for (unsigned i = 0; i < counts.size(); i++) {
            if (counts[i] > 1) {
                donor = i;
                break;
            }
        }
        bool copying_needed = (donor >= 0);
 
        if (copying_needed) {
            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
            }

            // Count number of cells with zero count that can serve as copy recipients
            unsigned nzeros = 0;
            for (unsigned i = 0; i < G::_nparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }

            while (nzeros > 0) {
                assert(donor < G::_nparticles);
                assert(recipient < G::_nparticles);
                
                assert(donor < G::_nparticles * group_number + end + 1);
                assert(recipient < G::_nparticles * group_number + end + 1);

                // Copy donor to recipient
                _particle_vec[recipient + start] = _particle_vec[donor + start];

                counts[donor]--;
                counts[recipient]++;
                nzeros--;

                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < G::_nparticles && counts[donor] < 2) {
                        donor++;
                    }
                }

                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < G::_nparticles && counts[recipient] > 0) {
                    recipient++;
                }
            } // while nzeros > 0
        } // if copying_needed
        return ess;
#endif
    }

    inline void Proj::filterParticlesThreading(unsigned g) {
      // divide up the groups as evenly as possible across threads
        unsigned first = 0;
        unsigned last = 0;
        unsigned stride = (G::_ngroups) / G::_nthreads; // divisor
        unsigned r = (G::_ngroups) % G::_nthreads; // remainder
        
        // need a vector of threads because we have to wait for each one to finish
        vector<thread> threads;

        for (unsigned i=0; i<G::_nthreads; i++) {
            first = last;
            last = first + stride;
            
            if (r > 0) {
                last += 1;
                r -= 1;
            }
            
            if (last > (G::_ngroups)) {
                last = (G::_ngroups);
            }
            
            threads.push_back(thread(&Proj::filterParticlesRange, this, first, last, g));
        }


      // the join function causes this loop to pause until the ith thread finishes
      for (unsigned i = 0; i < threads.size(); i++) {
        threads[i].join();
      }
    }

    inline void Proj::filterParticlesRange(unsigned first, unsigned last, unsigned g) {
        for (unsigned i=first; i<last; i++){
            // i is the group number
            filterParticles(g, i);
        }
    }

    inline double Proj::computeEffectiveSampleSize(const vector<double> & probs) const {
        double ss = 0.0;
        for_each(probs.begin(), probs.end(), [&ss](double w){ss += w*w;});
        double ess = 1.0/ss;
        return ess;
    }

    inline void Proj::simProposeParticles(unsigned step_number) {
        assert(G::_nthreads > 0);
        // don't bother threading simulations
        for (auto & p : _particle_vec) {
            p.simProposal(step_number);
        }
    }

    inline void Proj::proposeParticles(unsigned step_number) {
        assert(G::_nthreads > 0);
        if (G::_nthreads == 1) {
            for (auto & p : _particle_vec) {
                p.proposal(step_number);
            }
        }
        else {
            // run each group separately
            // divide up the groups as evenly as possible across threads
            unsigned first = 0;
            unsigned last = 0;
            unsigned stride = (G::_nparticles*G::_ngroups) / G::_nthreads; // divisor
            unsigned r = (G::_nparticles*G::_ngroups) % G::_nthreads; // remainder

            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;
                
                if (r > 0) {
                    last += 1;
                    r -= 1;
                }
                
                if (last > G::_ngroups * G::_nparticles) {
                    last = G::_ngroups * G::_nparticles;
                }
                
                threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, step_number));
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, unsigned step_number) {
        for (unsigned i=first; i<last; i++){
            _particle_vec[i].proposal(step_number);
        }
    }

    inline void Proj::summarizeData(Data::SharedPtr) {
        // Report information about data partition subsets

        std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
        std::cout << "Number of partition subsets: " << G::_nloci << std::endl;
        std::cout << "Number of particles: " << G::_nparticles << std::endl;

        for (unsigned subset = 0; subset < G::_nloci; subset++) {
            DataType dt = _partition->getDataTypeForSubset(subset);
            std::cout << "  Subset " << (subset+1) << " (" << _data->getSubsetName(subset) << ")" << std::endl;
            std::cout << "    data type: " << dt.getDataTypeAsString() << std::endl;
            std::cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << std::endl;
            std::cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << std::endl;
        }
    }

    inline void Proj::initializeParticles() { // TODO: make one template particle and copy it
        // set partials for first particle under save_memory setting for initial marginal likelihood calculation
         assert (G::_nthreads > 0);

        unsigned psuffix = 1;
        
         bool partials = true;

         for (auto & p:_particle_vec ) {
             if (G::_start_mode != "sim") {
                 p.setParticleData(_data, partials);
                 partials = false;
             }
             else {
                 p.createTrivialForest();
             }
             
             // set particle seed for drawing new values
             p.setSeed(rng->randint(1,9999) + psuffix);
             if (G::_start_mode != "sim") {
                 if (G::_est_clock_rate) {
                     p.drawClockRate();
                 }
                 else {
                     p.setClockRate(G::_clock_rate);
                 }
             }
             else {
                 p.setSimClockRate();
             }
             
             psuffix += 2;
             
             if (G::_est_lambda) {
                 assert (G::_est_mu);
                 if (G::_mu > 0.0) {
                     p.drawBirthDiff();
                     p.drawTurnover();
                     p.calculateLambdaAndMu();
                 }
                 else {
                     // Yule model
                     p.drawLambda();
                 }
             }
         }
        
        if (G::_taxsets.size() > 0) {
            removeUnecessaryTaxsets(); // remove any taxsets with just one lineage and one fossil because they will not provide any constraints
        }
        
        for (auto &p:_particle_vec) {
            p.setParticleTaxSets();
            p.setOverlappingTaxSets();
            p.setTaxSetsNoFossils();
        }
        
        if (G::_fossils.size() > 0) {
            for (auto &p:_particle_vec) {
                p.setFossils();
                p.drawFossilAges(); // draw a separate fossil age for each particle
            }
        }
        
    if (G::_est_root_age) {
        for (auto &p:_particle_vec) {
                p.drawRootAge();
            }
        }
    }

    inline void Proj::removeUnecessaryTaxsets() {
        // remove any taxsets with just one real taxon because there is no constraint in that case
        vector<unsigned> taxsets_to_remove;
        
        for (int i= (int) G::_taxsets.size() - 1; i >=0; i--) {
            if (G::_taxsets[i]._species_included.size() == 2) {
                string name1 = G::_taxsets[i]._species_included[0];
                string name2 = G::_taxsets[i]._species_included[1];
                
                if (boost::ends_with(name1, "FOSSIL") || (boost::ends_with(name2, "FOSSIL") )) {
                    taxsets_to_remove.push_back(i);
                }
            }
        }
        
        if (taxsets_to_remove.size() > 0) {
            for (int a = 0; a < taxsets_to_remove.size(); a++) {
                G::_taxsets.erase(G::_taxsets.begin() + taxsets_to_remove[a]);
            }
            
            // renumber taxsets
            int name_num = 1;
            for (auto &g:G::_taxsets) {
                g._name = name_num;
                name_num++;
            }
        }
        
    }

    inline void Proj::writePartialCount() {
        // save all partials to a .txt file
        ofstream partialf("partials.txt");
        double total_partial_count = 0;
        for (auto &p:_particle_vec) {
            total_partial_count += p.getPartialCount();
        }
        partialf << total_partial_count << endl;
        partialf.close();
    }

    inline void Proj::writeTreeFile() {
        // save all trees in nexus files
        ofstream treef("trees.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        
        for (unsigned n = 0; n <_indices_to_keep.size(); n++) {
            unsigned index = _indices_to_keep[n];
            treef << "  tree sample = [&R] " << _particle_vec[index].saveForestNewick() << ";\n";
        }
        
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::writeLogFile() {
        // this function creates a params file that is comparable to output from beast
        ofstream logf("params.log");
        logf << "Replicate_ID ";
        logf << "\t" << "posterior ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "prior ";

        for (unsigned i=0; i<G::_nloci; i++) {
            logf << "\t" << "Tree" + to_string(i) + "Likelihood";
        }
        
        logf << "\t" << "branch_rates";
        logf << "\t" << "diversificationRate";
        logf << "\t" << "extant_mrca";
        logf << "\t" << "extinction_rate";
        logf << "\t" << "speciation_rate";
        logf << "\t" << "turnover";
        
        for (auto &t:G::_taxsets) {
            string taxset_name = t._name;
            logf << "\t" << "mrca.date-backward(" + taxset_name + ")";
        }

        logf << endl;
        
        int iter = 0;
        for (unsigned n = 0; n <_indices_to_keep.size(); n++) {
            unsigned index = _indices_to_keep[n];
            Particle p = _particle_vec[index];
            logf << iter;
            iter++;
            
            vector<double> tree_log_likelihoods = p.calcGeneTreeLogLikelihoods();
            double log_likelihood = 0.0;
            for (auto &t:tree_log_likelihoods) {
                log_likelihood += t;
            }
            
//            double log_likelihood = p.getLogLikelihood();
            double yule = 0.0;
            
            double total_log_prior = p.getAllPriors();
            
            double log_posterior = total_log_prior + log_likelihood;
            
            if (G::_coverage) {
                _hpd_values.push_back(make_pair(p.getHeightFirstSplit(), log_posterior));
            }
            
//            vector<double> tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            
//            double height = p.getTreeHeight();
//            double length = p.getTreeLength();
            
            double clock_rate = p.getClockRate();
//            double FBD = p.getFBDModel();
//            double sampling_proportion = 1.0;
            
            double birth_death = 0.0;
            if (G::_mu > 0.0) {
                birth_death = p.getBirthDeathModel();
            }
            else {
                yule = p.getYuleModel();
            }
            
            logf << "\t" << log_posterior;
            logf << "\t" << log_likelihood;
            logf << "\t" << total_log_prior;
            for (auto &l:tree_log_likelihoods) {
                logf << "\t" << l;
            }
//            logf << "\t" << height;
//            logf << "\t" << length;
            logf << "\t" << clock_rate;
//            logf << "\t" << FBD;
            
            double lambda = 0.0;
            double mu = 0.0;
            
            if (G::_est_lambda) {
                lambda = p.getEstLambda();
            }
            else {
                lambda = G::_lambda;
            }
            
            if (G::_est_mu) {
                    mu = p.getEstMu();
            }
            else {
                mu = G::_mu;
            }
            
            double root_age = 0.0;
            if (G::_est_root_age) {
                root_age = p.getEstRootAge();
            }
            else {
                root_age = G::_root_age;
            }
            
//            to be consistent with beast output, birth rate is "effective birth rate": lambda - mu
//            death rate is relative to birth rate (mu / lambda)
//            https://groups.google.com/g/beast-users/c/HtpPKHGYNYg/m/pC1rIS5iCgAJ
//            https://beast2-dev.github.io/hmc/hmc//Priors/DeathRatePrior/
//            https://beast2-dev.github.io/hmc/hmc//Priors/BirthRatePrior/
            
            double diversification_rate = lambda - mu;
            double turnover = mu / lambda;
//            double SACountFBD = 0;
            
            map<string, double> particle_taxset_map = p.getTaxsetAges();
            
            logf << "\t" << diversification_rate;
            
//            if (G::_mu > 0.0) {
//                logf << "\t" << turnover;
//            }
            
//            logf << "\t" << sampling_proportion;
            
            logf << "\t" << root_age;
            
            logf << "\t" << mu;
            logf << "\t" << lambda;
            
//            logf << "\t" << SACountFBD;
            
            logf << "\t" << turnover;
            
            for (auto &t:G::_taxsets) {
                string name = t._name;
                double age = particle_taxset_map[name];
                logf << "\t" << age;
            }

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::writeLogMarginalLikelihoodFile() const {
        // this function writes the log marginal likelihood to a file
        ofstream logf ("marginal_likelihood.txt");
        logf << "log marginal likelihood: " << _log_marginal_likelihood << "\n";
    }

    inline void Proj::buildNonzeroMap(vector<Particle> &particles, map<const void *, list<unsigned> > & nonzero_map, const vector<unsigned> & nonzeros, vector<unsigned> particle_indices, unsigned start, unsigned end) {
        
        for (auto i : nonzeros) {
            unsigned particle_index = particle_indices[start + i];
            void * ptr = particles[particle_index].getForestPtr().get();

    //            void * ptr = particles[i].getGeneForestPtr(locus).get();
            if (nonzero_map.count(ptr) > 0) {
    //                nonzero_map[ptr].push_back(i);
                nonzero_map[ptr].push_back(particle_index);
            }
            else {
    //                nonzero_map[ptr] = {i};
                nonzero_map[ptr] = {particle_index};
            }
        }
    }

}
