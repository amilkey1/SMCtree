#pragma once

namespace proj {

class Particle {
    public:
        Particle();
        Particle(const Particle & other);
        ~Particle();

        typedef std::shared_ptr<Particle>               SharedPtr;
    
        void                                            setParticleData(Data::SharedPtr d, bool partials);
        void                                    operator=(const Particle & other);
        vector<double>                          calcGeneTreeLogLikelihoods();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    proposal(unsigned step_number);
        void                                    showParticle();
        vector<double>                          getGeneTreeLogLikelihoods();
        double                                  getTreeHeight();
        double                                  getTreeLength();
        double                                  getEstLambda() {return _forest._estimated_lambda;}
        double                                  getEstMu() {return _forest._estimated_mu;}
        double                                  getLogLikelihood();
        double                                  getYuleModel();
        double                                  getAllPriors();
        double                                  getBirthDeathModel();
        vector<pair<double, double>>            getTreeIncrementPriors();
        void                                    setStartingLogLikelihoods(vector<double> starting_log_likelihoods);
        void                                    clearPartials();
        string                                  saveForestNewick() {
                                                        return _forest.makeNewick(8, true);}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}

        string debugSaveParticleInfo(unsigned i) const;
        void drawBirthDiff();
        void drawTurnover();
        void drawRootAge();
        void calculateLambdaAndMu();
        void drawLambda();
    
        void calcStartingUPGMAMatrix();
        vector<vector<double>> getStartingUPGMAMatrix();
        void setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene);
        vector<map<Node*,  unsigned>> getStartingRowCount();
        void setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene);
        void calcStartingRowCount();
#if defined (FOSSILS)
        void setFossils() {_particle_fossils = G::_fossils;}
        void drawFossilAges();
#endif
    
    private:
        mutable                                 Lot::SharedPtr _lot;
        void                                    clear();
    
        Forest                                  _forest;
        double                                  _log_weight;
#if defined (FOSSILS)
        unsigned                                _fossil_number;
        vector<Fossil>                          _particle_fossils; // each particle needs it own set of fossils with their own ages
#endif
};

    inline string Particle::debugSaveParticleInfo(unsigned i) const {
        string s;
        s += str(format("Particle %d:\n") % i);
        s += str(format("  _log_weight = %.9f:\n") % _log_weight);
        s += _forest.debugSaveForestInfo();
        s += "\n";
        return s;
    }

    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
    }

    inline Particle::Particle(const Particle & other) {
        *this = other;
    }

    inline Particle::~Particle() {
        clear();
    }

    inline void Particle::clear() {
        _log_weight = 0.0;
#if defined (FOSSILS)
        _fossil_number = 0;
#endif
      }

    inline void Particle::setParticleData(Data::SharedPtr d, bool partials) {
        // one forest contains information about all loci since all loci share a tree
        _forest.setForestData(d, partials);
    }

    inline vector<double> Particle::calcGeneTreeLogLikelihoods() {
        vector<double> gene_forest_likelihoods;
        gene_forest_likelihoods.resize(G::_nloci);
        
        //calculate likelihood for each gene tree
        for (unsigned i=0; i<G::_nloci; i++) {
            double gene_tree_log_likelihood  = 0.0;
            gene_tree_log_likelihood = _forest.calcSubsetLogLikelihood(i);
            assert(!isnan (gene_tree_log_likelihood));
            gene_forest_likelihoods[i] = gene_tree_log_likelihood;
        }

        return gene_forest_likelihoods;
    }

    inline void Particle::proposal(unsigned step_number) {
        double prev_log_likelihood = _forest.getLogLikelihood();
        
        bool done = false;
        
        while (!done) {
            // TODO: need to draw an age for each fossil and sort the fossils by age
            // TODO: for now, just use existing ages
#if defined (FOSSILS)
                bool fossil_added = false;
            bool done_adding_increment = false;
            
            while (!done_adding_increment) {
                // loop through until you have ended on a non-fossil increment added
                if ((_fossil_number < G::_fossils.size()) || (_fossil_number == 0 && G::_fossils.size() == 1)) {
//                    fossil_added = _forest.addIncrementFossil(_lot, G::_fossils[_fossil_number]._age, G::_fossils[_fossil_number]._name);
                    fossil_added = _forest.addIncrementFossil(_lot, _particle_fossils[_fossil_number]._age, _particle_fossils[_fossil_number]._name);
                    if (fossil_added) {
                        _fossil_number++;
                        done_adding_increment = false;
                    }
                    else {
                        done_adding_increment = true;
                    }
                }
                else {
                    fossil_added = _forest.addIncrementFossil(_lot, -1, "placeholder");
                    assert (!fossil_added);
                    done_adding_increment = true;
                }
            }
#else
            _forest.addIncrement(_lot);
#endif
            pair<double, bool> output = _forest.joinTaxa(prev_log_likelihood, _lot);
            _log_weight = output.first;
            bool filter = output.second;
            // step is done when log weight is not 0 or when all the lineages are joined and all the fossils have been added
#if defined (FOSSILS)
            // TODO: if the branch is really long, joining nodes will not change the likelihood
//            if (_log_weight != 0.0 || (_forest._lineages.size() == 1 && _fossil_number == (unsigned) (_particle_fossils.size()))) {
            if (filter || (_forest._lineages.size() == 1 && _fossil_number == (unsigned) (_particle_fossils.size()))) {
                done = true;
            }
#else
            if (_log_weight != 0.0 || (_forest._lineages.size() == 1)) {
                done = true;
            }
#endif
            
            // if you are on the last step and still have fossils to join, keep going
            if (step_number == (G::_ntaxa-2) && _forest._lineages.size() > 1) {
                done = false;
            }
            
            if (G::_upgma_completion) {
                _log_weight = _forest.buildRestOfTreeUPGMA();
            }
        }
//        _forest.showForest();
    }

    inline void Particle::showParticle() {
        //print out weight of each particle
        output("\nParticle:\n", 1);
        output(format("log weight: %d") % _log_weight, 1);
        output("\nForest:\n", 1);
        _forest.showForest();
    }

    inline vector<double> Particle::getGeneTreeLogLikelihoods() {
        return _forest._gene_tree_log_likelihoods;
    }

    inline double Particle::getLogLikelihood() {
        double log_likelihood = 0.0;
        for (auto &l:_forest._gene_tree_log_likelihoods) {
            log_likelihood += l;
        }
        
        return log_likelihood;
    }

    inline double Particle::getTreeHeight() {
        return _forest.getTreeHeight();
    }

    inline double Particle::getTreeLength() {
        return _forest.getTreeLength();
    }

    inline double Particle::getYuleModel() {
        return _forest.getTreePrior();
    }

    inline double Particle::getAllPriors() {
        double tree_prior = _forest.getTreePrior();
        double param_prior = 0.0;
        
        // TODO: for now, these params are all drawn from exponential distributions - need to change priors if distribution is changed
        // TODO: this also assumes the mean of the exponential distribution is the user-specified param
        if (G::_est_mu && G::_mu > 0.0) {
            param_prior += log(G::_mu) - (_forest._estimated_mu * G::_mu);
        }
        
        if (G::_est_lambda) {
            param_prior += log(G::_lambda) - (_forest._estimated_lambda * G::_lambda);
        }
        
        if (G::_est_root_age) {
            param_prior += log(G::_root_age) - (_forest._estimated_root_age * G::_root_age);
        }
                
        double total_prior = tree_prior + param_prior;
        return total_prior;
    }

    inline double Particle::getBirthDeathModel() {
        // TODO: fix this for birth death model
        return _forest.getTreePrior();
    }

    inline void Particle::setStartingLogLikelihoods(vector<double> starting_log_likelihoods) {
        _forest._gene_tree_log_likelihoods = starting_log_likelihoods;
    }

    inline vector<pair<double, double>> Particle::getTreeIncrementPriors() {
        return _forest._increments_and_priors;
    }

    inline void Particle::clearPartials() {
        _forest.clearPartials();
    }

    inline void Particle::calcStartingUPGMAMatrix() {
        _forest.buildStartingUPGMAMatrix();
    }

    inline vector<vector<double>> Particle::getStartingUPGMAMatrix() {
        vector<vector<double>> starting_upgma_matrices_by_gene;
        starting_upgma_matrices_by_gene.push_back(_forest._starting_dij);
        return starting_upgma_matrices_by_gene;
    }

    inline void Particle::setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene) {
        _forest._starting_dij = starting_upgma_matrices_by_gene[0];
    }

    inline vector<map<Node*,  unsigned>> Particle::getStartingRowCount() {
        vector<map<Node*,  unsigned>> starting_row_count_by_gene;
        starting_row_count_by_gene.push_back(_forest._starting_row);
        return starting_row_count_by_gene;
    }

    inline void Particle::setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene) {
        _forest._starting_row = starting_row_count_by_gene[0];
    }

    inline void Particle::calcStartingRowCount() {
        _forest.buildStartingRow();
    }

    inline void Particle::drawBirthDiff() {
        assert (G::_est_lambda);
        assert (G::_est_mu);
        _forest.drawBirthDiff(_lot);
    }

    inline void Particle::drawTurnover() {
        assert (G::_est_mu);
        assert (G::_est_lambda);
        _forest.drawTurnover(_lot);
    }

    inline void Particle::drawRootAge() {
        assert (G::_est_root_age > 0.0);
        _forest.drawRootAge(_lot);
    }
    
    inline void Particle::calculateLambdaAndMu() {
        _forest.calculateLambdaAndMu();
    }

    inline void Particle::drawLambda() {
        _forest.drawLambda(_lot);
    }

    inline void Particle::drawFossilAges() {
        // draw an age for each fossil, using the upper and lower bounds as specified
         for (auto &f:_particle_fossils) {
            f._age = _lot->uniformConstrained(f._lower, f._upper); // TODO: unsure if this lot function is working correctly
        }
        
        sort(_particle_fossils.begin(), _particle_fossils.end(), [](Fossil & left, Fossil & right) {
            return left._age < right._age;
        });
    }

    inline void Particle::operator=(const Particle & other) {
        _forest = other._forest;
        _log_weight = other._log_weight;
#if defined (FOSSILS)
        _fossil_number = other._fossil_number;
        _particle_fossils = other._particle_fossils;
#endif
    }
    
}
