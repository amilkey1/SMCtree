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
   
            void                            processCommandLineOptions(int argc, const char * argv[]);
            void                            run();
            void                            initializeParticles(vector<Particle> &particles);
        
        private:
            void                            clear();
            void                            simulate();
            void                            smc();
            void                            summarizeData(Data::SharedPtr);
            void                            proposeParticles(vector<Particle> &particles);
            void                            proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            double                          filterParticles(unsigned step, vector<Particle> & particles);
            double                          computeEffectiveSampleSize(const vector<double> & probs) const;
        
            Partition::SharedPtr            _partition;
            Data::SharedPtr                 _data;
            double                          _log_marginal_likelihood;
                                    
    };

    inline Proj::Proj() {
        clear();
    }

    inline Proj::~Proj() {
    }
    
    inline void Proj::clear() {
        _log_marginal_likelihood = 0.0;
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        vector<string> partition_subsets;
        variables_map vm;
        options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("lambda",  value(&G::_lambda)->default_value(10.0), "per lineage speciation rate assumed for the Yule model")
        ("theta",  value(&G::_theta)->default_value(0.1), "theta for entire tree")
        ("rnseed",  value(&G::_rnseed)->default_value(1), "pseudorandom number seed")
        ("nthreads",  value(&G::_nthreads)->default_value(1), "number of threads")
        ("startmode",  value(&G::_start_mode)->default_value("smc"), "smc or sim")
        ("filename",  value(&G::_filename), "name of file containing data")
        ("verbosity",  value(&G::_verbosity), "level of output")
        ("subset",  value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("nparticles",  value(&G::_nparticles)->default_value(500), "number of particles)")
        ("savememory",  value(&G::_save_memory)->default_value(false), "save memory by recalculating partials each round)")
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
            
            rng->setSeed(G::_rnseed);
            
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
            
            unsigned nsteps = (G::_ntaxa-1);
            
            for (unsigned g=0; g<nsteps; g++){
                // set random number seeds
                unsigned psuffix = 1;
                for (auto &p:particle_vec) {
                    p.setSeed(rng->randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                unsigned step_plus_one = g+1;
                output(format("Step %d of %d steps.\n") % step_plus_one % nsteps, 1);
                proposeParticles(particle_vec);
                double ess = filterParticles(g, particle_vec);
                output(format("     ESS = %d\n") % ess, 2);
                
                for (auto &p:particle_vec) {
                    p.showParticle();
                }
            }
            
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
        
}
