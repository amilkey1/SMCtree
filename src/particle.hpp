#pragma once

namespace proj {

class Particle {
    public:
        Particle();
//        Particle(const Particle & other); // TODO: why does this cause compilation to fail?
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
        double                                  getLogLikelihood();
        double                                  getYuleModel();
        vector<pair<double, double>>            getSpeciesTreeIncrementPriors();
        void                                    setStartingLogLikelihoods(vector<double> starting_log_likelihoods);
        string                                  saveForestNewick() {
                                                        return _forest.makeNewick(8, true);}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}

    
    private:
        mutable                                 Lot::SharedPtr _lot;
        void                                    clear();
    
        Forest                                  _forest;
        double                                  _log_weight;

    
};

    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
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
    }

    inline void Particle::showParticle() {
        //print out weight of each particle
        output("\nParticle:\n", 1);
        output(format("log weight: %d") % _log_weight, 1);
//        output(format("log likelihood: %d") % _log_likelihood, 1);
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
        return _forest.getSpeciesTreePrior();
    }

    inline void Particle::setStartingLogLikelihoods(vector<double> starting_log_likelihoods) {
        _forest._gene_tree_log_likelihoods = starting_log_likelihoods;
    }

    inline vector<pair<double, double>> Particle::getSpeciesTreeIncrementPriors() {
        return _forest._increments_and_priors;
    }

    inline void Particle::operator=(const Particle & other) {
        _forest = other._forest;
        _log_weight = other._log_weight;
    }
}
