#pragma once

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
        
        private:
            void clear();
            
            void simulate();
            void simulateTree();
            void simulateData(string locus_name, unsigned locus_length);
            void simulateSave(string fnprefix);
            
            void smc();
            
            void parseLocusSpec(string s) const;
            
            Forest _sim_tree;
    };

    inline Proj::Proj() {
        clear();
    }

    inline Proj::~Proj() {
    }
    
    inline void Proj::clear() {
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> locus_specs;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("lambda",  value(&G::_lambda)->default_value(10.0), "per lineage speciation rate assumed for the Yule model")
        ("rnseed",  value(&G::_rnseed)->default_value(1), "pseudorandom number seed")
        ("nthreads",  value(&G::_nthreads)->default_value(1), "number of threads")
        ("startmode",  value(&G::_start_mode)->default_value("smc"), "smc or sim")
        ("simfnprefix",  value(&G::_sim_filename_prefix), "prefix of files in which to save simulated trees and data if startmode is 'sim' (e.g. specifying 'sim' results in files named 'sim.tre' and 'sim.nex')")
        ("simntaxa",  value(&G::_sim_ntaxa)->default_value(4), "number of taxa to simulate if startmode is 'sim'")
        ("simlocus",  value(&locus_specs), "string providing the name and number of sites for a locus if startmode is 'sim' (e.g. 'rbcl:1314')")
        ("simlambda",  value(&G::_sim_lambda)->default_value(1.0), "true speciation rate for simulating tree under the Yule model if startmode is 'sim'")
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
        
        // If user specified simlocus, break locus specification
        // into a locus name and number of sites
        if (vm.count("simlocus") > 0) {
            for (auto s : locus_specs) {
                parseLocusSpec(s);
            }
        }
        
    }
    
    inline void Proj::parseLocusSpec(string s) const {
        vector<string> v;
        
        // Separate part before colon from part after colon
        split(v, s, boost::is_any_of(":"));
        if (v.size() != 2) {
            throw XProj(format("Expecting exactly one colon in locus definition (\"%s\")") % s);
        }

        string locus_name = v[0];
        trim(locus_name);
        
        if (!all_of(v[1].begin(), v[1].end(), ::isdigit)) {
            throw XProj(format("Locus length specification (\"%s\") must contain only digits") % v[1]);
        }
        
        unsigned locus_length = (unsigned)stoi(v[1]);
        G::_sim_locus_name.push_back(locus_name);
        G::_sim_locus_length.push_back(locus_length);
    }

    inline void Proj::run() {
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
        assert(G::_sim_locus_name.size() > 0);
        assert(G::_sim_locus_name.size() == G::_sim_locus_length.size());
        
        // Show simulation settings to user
        output(format("Simulating data for %d taxa and %d loci\n") % G::_sim_ntaxa % G::_sim_locus_name.size(), 0);
        output("  Locus names and lengths:\n", 0);
        for (unsigned i = 0; i < G::_sim_locus_name.size(); ++i) {
            output(format("    %s (%d sites)\n") % G::_sim_locus_name[i] % G::_sim_locus_length[i], 0);
        }
        output(format("  True speciation rate: %g\n") % G::_sim_lambda, 0);
        
        // Simulate the tree
        simulateTree();
        
        // Simulate data for each locus given the tree
        unsigned nloci = (unsigned)G::_sim_locus_name.size();
        for (unsigned i = 0; i < nloci; i++) {
            simulateData(G::_sim_locus_name[i], G::_sim_locus_length[i]);
        }
        
        // Save simulated data to file
        simulateSave(G::_sim_filename_prefix);
    }

    inline void Proj::simulateTree() {
        // Set global _ntaxa
        unsigned nleaves = G::_sim_ntaxa;
        G::_ntaxa = nleaves;
        
        // Make up taxon names
        G::_taxon_names.resize(nleaves);
        for (unsigned i = 0; i < nleaves; i++)
            G::_taxon_names[i] = G::inventName(i, /*lower_case*/false);
        
        unsigned nsteps = nleaves - 1;
        _sim_tree.createTrivialForest();
        for (unsigned i = 0; i < nsteps; i++) {
            // Determine number of lineages remaining
            unsigned n = _sim_tree.getNumLineages();
            assert(n > 1);
            
            // Waiting time to speciation event is Exponential(rate = n*lambda)
            // u = 1 - exp(-r*t) ==> t = -log(1-u)/r
            double r = G::_sim_lambda*n;
            double u = rng->uniform();
            double t = -log(1.0 - u)/r;
            _sim_tree.advanceAllLineagesBy(t);
            
            // Join two random lineages
            _sim_tree.joinRandomLineagePair(rng);
        }
        assert(_sim_tree.getNumLineages() == 1);
        _sim_tree.refreshAllPreorders();
        _sim_tree.renumberInternals();
    }
    
    inline void Proj::simulateData(string locus_name, unsigned locus_length) {
        output(format("*** not simulating data for locus %s (length %d) ***\n") % locus_name % locus_length, 0);
    }
    
    inline void Proj::simulateSave(string fnprefix) {
        string tree_file_name = fnprefix + ".tre";
        ofstream treef(tree_file_name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree true = [&R] " << _sim_tree.makeNewick(/*precision*/9, /*use_names*/true) << ";\n";
        treef << "end;\n\n";
        treef.close();
    }

    inline void Proj::smc() {
        
    }
        
}
