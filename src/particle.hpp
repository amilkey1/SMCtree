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
        void                                    proposal();
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
        void drawLambda();
        void drawMu();
        void drawRootAge();
    
#if defined (UPGMA_COMPLETION)
        void calcStartingUPGMAMatrix();
        vector<vector<double>> getStartingUPGMAMatrix();
        void setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene);
        vector<map<Node*,  unsigned>> getStartingRowCount();
        void setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene);
        void calcStartingRowCount();
#endif
    
    private:
        mutable                                 Lot::SharedPtr _lot;
        void                                    clear();
    
        Forest                                  _forest;
        double                                  _log_weight;
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

    inline void Particle::proposal() {
       _forest.addIncrement(_lot);
        _log_weight = _forest.joinTaxa(_lot);
#if defined (UPGMA_COMPLETION)
        _log_weight = _forest.buildRestOfTreeUPGMA();
#endif
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
        if (G::_est_mu) {
            param_prior += log(G::_mu) - (_forest._estimated_mu * G::_mu);
        }
        
        if (G::_est_lambda) {
            param_prior += log(G::_lambda) - (_forest._estimated_lambda * G::_lambda);
        }
        
        if (G::_est_root_age) {
            param_prior += log(G::_root_age) - (_forest._estimated_root_age * G::_root_age);
        }
        
//        total_prior += log(Forest::_theta_prior_mean) - (_forests[1]._theta * Forest::_theta_prior_mean);
        
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

#if defined (UPGMA_COMPLETION)
    inline void Particle::calcStartingUPGMAMatrix() {
            _forest.buildStartingUPGMAMatrix();
    }
#endif

#if defined (UPGMA_COMPLETION)
    inline vector<vector<double>> Particle::getStartingUPGMAMatrix() {
        vector<vector<double>> starting_upgma_matrices_by_gene;
        starting_upgma_matrices_by_gene.push_back(_forest._starting_dij);
        return starting_upgma_matrices_by_gene;
    }
#endif

#if defined (UPGMA_COMPLETION)
    inline void Particle::setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene) {
        _forest._starting_dij = starting_upgma_matrices_by_gene[0];
    }
#endif

#if defined (UPGMA_COMPLETION)
    inline vector<map<Node*,  unsigned>> Particle::getStartingRowCount() {
        vector<map<Node*,  unsigned>> starting_row_count_by_gene;
        starting_row_count_by_gene.push_back(_forest._starting_row);
        return starting_row_count_by_gene;
    }
#endif

#if defined (UPGMA_COMPLETION)
    inline void Particle::setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene) {
        _forest._starting_row = starting_row_count_by_gene[0];
    }
#endif

#if defined (UPGMA_COMPLETION)
    inline void Particle::calcStartingRowCount() {
        _forest.buildStartingRow();
    }
#endif

    inline void Particle::drawLambda() {
        assert (G::_est_lambda);
        _forest.drawLambda(_lot);
    }

    inline void Particle::drawMu() {
        assert (G::_est_mu);
        _forest.drawMu(_lot);
    }

    inline void Particle::drawRootAge() {
        assert (G::_est_root_age > 0.0);
        _forest.drawRootAge(_lot);
    }

    inline void Particle::operator=(const Particle & other) {
        _forest = other._forest;
        _log_weight = other._log_weight;
    }
    
}
