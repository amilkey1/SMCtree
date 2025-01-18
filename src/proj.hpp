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
            void initializeParticles(vector<Particle> &particles);
        
        private:
            void clear();

            void simulate();
            void simulateTree();
            void simulateData(Lot::SharedPtr lot, unsigned startat, unsigned locus_length);
            void simulateSave(string fnprefix);

            void smc();
            void summarizeData(Data::SharedPtr);
            void proposeParticles(vector<Particle> &particles);
            void proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            double filterParticles(unsigned step, vector<Particle> & particles);
            double computeEffectiveSampleSize(const vector<double> & probs) const;
            void writeTreeFile (vector<Particle> &v) const;
            void writeLogFile(vector<Particle> &v) const;
            void handleBaseFrequencies();
            void handleRelativeRates();

        
            Forest::SharedPtr    _sim_tree;
            Partition::SharedPtr _partition;
            Data::SharedPtr      _data;
            double               _log_marginal_likelihood;
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
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("lambda",  value(&G::_lambda)->default_value(10.0), "per lineage speciation rate assumed for the Yule model")
        ("mu",  value(&G::_mu)->default_value(1.0), "per lineage extinction rate assumed for the birth death model")
        ("rnseed",  value(&G::_rnseed)->default_value(1), "pseudorandom number seed")
        ("nthreads",  value(&G::_nthreads)->default_value(1), "number of threads")
        ("startmode",  value(&G::_start_mode)->default_value("smc"), "smc or sim")
        ("filename",  value(&G::_filename), "name of file containing data")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles)")
        ("savememory",  value(&G::_save_memory)->default_value(false), "save memory by recalculating partials each round")
        ("model",  value(&G::_model)->default_value("JC"), "model to use for likelihood calculations")
        ("verbosity",  value(&G::_verbosity)->default_value(1), "0 (nothing but essential output), 1, 2, ...")
        ("simfnprefix",  value(&G::_sim_filename_prefix), "prefix of files in which to save simulated trees and data if startmode is 'sim' (e.g. specifying 'sim' results in files named 'sim.tre' and 'sim.nex')")
        ("simntaxa",  value(&G::_sim_ntaxa)->default_value(4), "number of taxa to simulate if startmode is 'sim'")
        ("simlambda",  value(&G::_sim_lambda)->default_value(1.0), "true speciation rate for simulating tree under the Yule model if startmode is 'sim'")
        ("kappa",  boost::program_options::value(&G::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&G::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("relative_rates", boost::program_options::value(&G::_string_relative_rates)->default_value("null"))
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
        simulateTree();
        
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
            _sim_tree->simulateData(::rng, _data, starting_site, G::_nsites_per_locus[g]);
            starting_site += G::_nsites_per_locus[g];
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
        
        _sim_tree = Forest::SharedPtr(new Forest);
        _sim_tree->buildYuleTree();
    }
        
    inline void Proj::simulateSave(string fnprefix) {
        // Show simulation settings to user
        string locus_or_loci = (G::_nloci != 1 ? "loci" : "locus");
        output(format("Simulating data for %d taxa and %d %s\n") % G::_sim_ntaxa % G::_nloci % locus_or_loci, 0);
        output("  Locus names and lengths:\n", 0);
        for (unsigned i = 0; i < G::_nloci; ++i) {
            output(format("    %s (%d sites)\n") % G::_locus_names[i] % G::_nsites_per_locus[i], 0);
        }
        output(format("  True speciation rate: %g\n") % G::_sim_lambda, 0);

        string tree_file_name = fnprefix + ".tre";
        ofstream treef(tree_file_name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree true = [&R] " << _sim_tree->makeNewick(/*precision*/9, /*use_names*/true) << ";\n";
        treef << "end;\n\n";
        treef.close();

        // Output data to file
        string data_file_name = fnprefix + ".nex";
        _data->compressPatterns();
        _data->writeDataToFile(data_file_name);
        output(format("  Sequence data saved in file \"%s\"\n") % data_file_name, 0);

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

    inline void Proj::smc() {
        output("\nStarting....\n",1);
        output(format("Current working directory: %s\n") % boost::filesystem::current_path(), 1);
        output(format("Random seed: %d\n") % G::_rnseed, 1);
        output(format("Number of threads: %d\n") % G::_nthreads, 1);
        
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
                        
            // create vector of particles
            vector<Particle> particle_vec;
            particle_vec.resize(G::_nparticles);

            for (unsigned i=0; i<G::_nparticles; i++) {
                particle_vec[i] = Particle();
            }
            
            initializeParticles(particle_vec); // initialize in parallel with multithreading
            
            // reset marginal likelihood
            _log_marginal_likelihood = 0.0;
            vector<double> starting_log_likelihoods = particle_vec[0].calcGeneTreeLogLikelihoods();
            
            _log_marginal_likelihood = 0.0;
            for (auto &l:starting_log_likelihoods) {
                _log_marginal_likelihood += l;
            }
            
            // initialize starting log likelihoods for all other particles
            // necessary for calculating the first weight
            
            for (unsigned p=1; p<particle_vec.size(); p++) {
                particle_vec[p].setStartingLogLikelihoods(starting_log_likelihoods);
            }
            
            if (G::_save_memory) {
                for (auto &p:particle_vec) {
                    p.clearPartials();
                }
            }
            
            unsigned nsteps = (G::_ntaxa-1);
            
            for (unsigned g=0; g<nsteps; g++){
                // set random number seeds
                unsigned psuffix = 1;
                for (auto &p:particle_vec) {
                    p.setSeed(rng->randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                unsigned step_plus_one = g+1;
                output(format("Step %d of %d.\n") % step_plus_one % nsteps, 1);
                proposeParticles(particle_vec);
                
                double ess = filterParticles(g, particle_vec);
                output(format("     ESS = %d\n") % ess, 2);
                
            }
            
            writeTreeFile(particle_vec);
            writeLogFile(particle_vec);
            
        }
        
        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }
    }

    inline double Proj::filterParticles(unsigned step, vector<Particle> & particles) {
        unsigned nparticles = (unsigned) particles.size();
        // Copy log weights for all bundles to prob vector
        vector<double> probs(nparticles, 0.0);
        
        for (unsigned p=0; p < nparticles; p++) {
            probs[p] = particles[p].getLogWeight();
        }
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = G::calcLogSum(probs);
        
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood due to this step
        _log_marginal_likelihood += log_sum_weights - log(nparticles);
        
        double ess = 0.0;
        if (G::_verbosity > 1) {
            // Compute effective sample size
            ess = computeEffectiveSampleSize(probs);
        }

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (nparticles, 0);

        // Throw _nparticles darts
        for (unsigned i=0; i<nparticles; i++) {
            double u = rng->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)std::distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Copy particles

      bool copying_needed = true;
      
        // Locate first donor
        unsigned donor = 0;
        while (counts[donor] < 2) {
            donor++;
            if (donor >= counts.size()) {
                copying_needed = false; // all the particle counts are 1
                break;
            }
        }

      if (copying_needed) {
            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
            }

            // Count number of cells with zero count that can serve as copy recipients
            unsigned nzeros = 0;
            for (unsigned i = 0; i < nparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }

            while (nzeros > 0) {
                assert(donor < nparticles);
                assert(recipient < nparticles);

                // Copy donor to recipient
                particles[recipient] = particles[donor];

                counts[donor]--;
                counts[recipient]++;
                nzeros--;

                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < nparticles && counts[donor] < 2) {
                        donor++;
                    }
                }

                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < nparticles && counts[recipient] > 0) {
                    recipient++;
                }
            }
      }
        return ess;
    }

    inline double Proj::computeEffectiveSampleSize(const vector<double> & probs) const {
        double ss = 0.0;
        for_each(probs.begin(), probs.end(), [&ss](double w){ss += w*w;});
        double ess = 1.0/ss;
        return ess;
    }

    inline void Proj::proposeParticles(vector<Particle> & particles) {
        assert(G::_nthreads > 0);
        if (G::_nthreads == 1) {
          for (auto & p : particles) {
              p.proposal();
          }
        }
        else {
          // divide up the particles as evenly as possible across threads
            unsigned first = 0;
            unsigned last = 0;
            unsigned stride = G::_nparticles / G::_nthreads; // divisor
            unsigned r = G::_nparticles % G::_nthreads; // remainder
            
            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;
                
                if (r > 0) {
                    last += 1;
                    r -= 1;
                }
                
                if (last > G::_nparticles) {
                    last = G::_nparticles;
                }
                
                threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles)));
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i].proposal();
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

    inline void Proj::initializeParticles(vector<Particle> &particles) { // TODO: can thread this
        // set partials for first particle under save_memory setting for initial marginal likelihood calculation
         assert (G::_nthreads > 0);

         bool partials = true;

         for (auto & p:particles ) {
             p.setParticleData(_data, partials);
             partials = false;
         }
    }

    inline void Proj::writeTreeFile(vector<Particle> &particles) const {
        // save all trees in nexus files
        ofstream treef("trees.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        
        for (auto &p:particles) {
            treef << "  tree sample = [&R] " << p.saveForestNewick() << ";\n";
        }
        
        treef << "end;\n";
        treef.close();
        
        // save unique trees only
        ofstream uniquetreef("unique_trees.trees");
        uniquetreef << "#nexus\n\n";
        uniquetreef << "begin trees;\n";
        
        vector<vector<pair<double, double>>> unique_increments_and_priors;
        for (auto &p:particles) {
            vector<pair<double, double>> increments_and_priors = p.getSpeciesTreeIncrementPriors();
            bool found = false;
            if(std::find(unique_increments_and_priors.begin(), unique_increments_and_priors.end(), increments_and_priors) != unique_increments_and_priors.end()) {
                found = true;
            }
            if (!found) {
                unique_increments_and_priors.push_back(increments_and_priors);
                uniquetreef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
            }
        }
        uniquetreef << "end;\n";
        uniquetreef.close();
    }

    inline void Proj::writeLogFile(vector<Particle> &v) const {
        // this function creates a params file that is comparable to output from beast
        ofstream logf("params.log");
        logf << "iter ";
        logf << "\t" << "posterior ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "prior ";

        for (unsigned i=0; i<G::_nloci; i++) {
            logf << "\t" << "Tree" + to_string(i) + "Likelihood";
        }
        
        logf << "\t" << "Tree.height";
        logf << "\t" << "Tree.treeLength";
        logf << "\t" << "YuleModel";
        logf << "\t" << "birthRate";


        logf << endl;

        int iter = 0;
        for (auto &p:v) {
            logf << iter;
            iter++;
            
            double log_likelihood = p.getLogLikelihood();
            double yule = p.getYuleModel();
            double log_prior = yule;
            
            double log_posterior = log_prior + log_likelihood;
            vector<double> tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            double height = p.getTreeHeight();
            double length = p.getTreeLength();
            


            logf << "\t" << log_posterior;
            logf << "\t" << log_likelihood;
            logf << "\t" << log_prior;
            for (auto &l:tree_log_likelihoods) {
                logf << "\t" << l;
            }
            logf << "\t" << height;
            logf << "\t" << length;
            logf << "\t" << yule;
            logf << "\t" << G::_lambda;

            logf << endl;
        }

        logf.close();
    }

}
