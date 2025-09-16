#pragma once

namespace proj {

class Particle {
    public:
        Particle();
        Particle(const Particle & other);
        ~Particle();

        typedef std::shared_ptr<Particle>               SharedPtr;
    
        void                                    setParticleData(Data::SharedPtr d, bool partials);
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
        double                                  getEstRootAge() {return _forest._estimated_root_age;}
        double                                  getLogLikelihood();
        double                                  getYuleModel();
        double                                  getAllPriors();
        double                                  getFBDModel();
        map<string, double>                     getTaxsetAges();
        double                                  getBirthDeathModel();
        void                                    setStartingLogLikelihoods(vector<double> starting_log_likelihoods);
        void                                    clearPartials();
        double                                  getPartialCount();
        void                                    drawClockRate();
        void                                    createTrivialForest();
        double                                  getClockRate() {return _forest._clock_rate;}
        void                                    simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites);
        string                                  makeNewick(unsigned precision, bool use_names);
        string                                  saveForestNewick() {
                                                        return _forest.makeNewick(8, true);}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}

        string debugSaveParticleInfo(unsigned i) const;
        void drawBirthDiff();
        void drawTurnover();
        void drawRootAge();
        void calculateLambdaAndMu();
        void drawLambda();

#if defined (FOSSILS)
        void setFossils() {_particle_fossils = G::_fossils;}
        void drawFossilAges();
        void setParticleTaxSets();
        void setOverlappingTaxSets();
#endif
    
        // validation stuff
        double getLowestFossilAge();
    
    private:
        mutable                                 Lot::SharedPtr _lot;
        void                                    clear();
    
        Forest                                  _forest;
        double                                  _log_weight;
#if defined (FOSSILS)
        unsigned                                _fossil_number;
        vector<Fossil>                          _particle_fossils; // each particle needs it own set of fossils with their own ages
        vector<TaxSet>                          _particle_taxsets; // update this as nodes are joined
        vector<TaxSet>                          _unused_particle_taxsets; // if there are overlapping taxa in taxsets, put the largest groups here until they can be used
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
        _lot.reset(new Lot());
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

    inline void Particle::createTrivialForest() {
        _forest.createTrivialForest();
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
        bool last_step = false;
        if (step_number == G::_ntaxa - 2) {
            last_step = true;
        }
        double prev_log_likelihood = _forest.getLogLikelihood();
        
        bool done = false;
        
        while (!done) {
#if defined (FOSSILS)
            bool fossil_added = false;
            bool done_adding_increment = false;
            bool valid = false;
            
            while (!done_adding_increment || !valid) {
                _forest._valid_taxsets.clear();
                // loop through until you have ended on a non-fossil increment added or forest is finished
                if ((_fossil_number < G::_fossils.size()) || (_fossil_number == 0 && G::_fossils.size() == 1)) {
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
                
                // check that at least one taxon set is valid
                // if there is no valid taxon set, keep adding increments until a valid set has been reached
                
                valid = _forest.checkForValidTaxonSet(_particle_taxsets, _unused_particle_taxsets);
                
                if (_forest._lineages.size() == 1 && _fossil_number == G::_fossils.size() - 1) { // if forest is finished, break out of step even if we have ended on a fossil
                    done_adding_increment = true;
                    assert (valid == true);
                }
            }
#else
            _forest.addIncrement(_lot);
#endif
            pair<double, bool> output = _forest.joinTaxa(prev_log_likelihood, _lot, _particle_taxsets, _unused_particle_taxsets);

            _log_weight = output.first; // TODO: how to decide when the step is done if there is no data and no likelihood change?
            bool filter = output.second;
            
            // step is done when log weight is not 0 or when all the lineages are joined and all the fossils have been added
#if defined (FOSSILS)
            // if the branch is really long, joining nodes will not change the likelihood
            if (filter || (_forest._lineages.size() == 1 && _fossil_number == (unsigned) (_particle_fossils.size()))) {
                done = true; // if the forest isn't finished and we're on the last step, keep going
                if (last_step && (_fossil_number != (unsigned) _particle_fossils.size() || _forest._lineages.size() != 1)) {
                    done = false;
                }
            }

#else
            if (_log_weight != 0.0 || (_forest._lineages.size() == 1)) {
                done = true;
            }
#endif
        }
        if (step_number == G::_ntaxa - 2) {
            assert (_fossil_number == G::_fossils.size());
            assert (_forest._lineages.size() == 1);
        }
        
        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }
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

    inline double Particle::getFBDModel() {
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

    inline void Particle::clearPartials() {
        _forest.clearPartials();
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
        double max_fossil_age = -1;
#if defined (FOSSILS)
        max_fossil_age = _particle_fossils.back()._age;
#endif
        _forest.drawRootAge(_lot, max_fossil_age);
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
             if (f._lower == f._upper) {
                 f._lower -= 0.001;
                 f._upper += 0.001; // make upper and lower slightly different if they are equal
             }
            f._age = _lot->uniformConstrained(f._lower, f._upper); // TODO: unsure if this lot function is working correctly
        }
        
        sort(_particle_fossils.begin(), _particle_fossils.end(), [](Fossil & left, Fossil & right) {
            return left._age < right._age;
        });
    }

    inline void Particle::setParticleTaxSets() {
        _particle_taxsets = G::_taxsets;
    }

    inline void Particle::setOverlappingTaxSets() {
        for (auto &t:_particle_taxsets) {
            sort(t._species_included.begin(), t._species_included.end());
        }
        
        vector<unsigned> unused_taxsets_index;
        // go through taxsets one by one and check for duplicate taxa
        for (unsigned count = 0; count < _particle_taxsets.size(); count++) {
            for (unsigned comparison = count + 1; comparison < _particle_taxsets.size(); comparison++) {
                vector<string> common_elements;
                // Find the intersection of the two vectors

                std::set_intersection(_particle_taxsets[count]._species_included.begin(), _particle_taxsets[count]._species_included.end(),
                                  _particle_taxsets[comparison]._species_included.begin(), _particle_taxsets[comparison]._species_included.end(),
                                  std::back_inserter(common_elements));
                if (common_elements.size() > 0) { // the larger taxset  goes in unused taxsets
                    if (_particle_taxsets[count]._species_included.size() > _particle_taxsets[comparison]._species_included.size()) {
                        if (std::find(unused_taxsets_index.begin(), unused_taxsets_index.end(), count) == unused_taxsets_index.end()) {
                            unused_taxsets_index.push_back(count); // don't include duplicates
                        }
                    }
                    else {
                        if (std::find(unused_taxsets_index.begin(), unused_taxsets_index.end(), comparison) == unused_taxsets_index.end()) {
                            unused_taxsets_index.push_back(comparison); // don't include duplicates
                        }
                    }
                }
            }
        }
        
        if (unused_taxsets_index.size() > 0) {
            for (unsigned i=0; i<unused_taxsets_index.size(); i++) {
                _unused_particle_taxsets.push_back(_particle_taxsets[unused_taxsets_index[i]]);
            }
            for (int i = (int) unused_taxsets_index.size() - 1; i>=0; i--) {
                _particle_taxsets.erase(_particle_taxsets.begin() + unused_taxsets_index[i]);
            }
        }
    }

    inline map<string, double> Particle::getTaxsetAges() {
        return _forest._taxset_ages;
    }

    inline double Particle::getPartialCount() {
        return _forest._partial_count;
    }

    inline void Particle::drawClockRate() {
        double clock_rate = _lot->uniform(); // TODO: maybe not the best prior
        _forest._clock_rate = clock_rate;
    }

    inline void Particle::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        _forest.simulateData(lot, data, starting_site, nsites);
    }

    inline string Particle::makeNewick(unsigned precision, bool use_names) {
        return _forest.makeNewick(precision, use_names);
    }

    inline double Particle::getLowestFossilAge() {
        return _particle_fossils[0]._age;
    }

    inline void Particle::operator=(const Particle & other) {
        _forest = other._forest;
        _log_weight = other._log_weight;
#if defined (FOSSILS)
        _fossil_number = other._fossil_number;
        _particle_fossils = other._particle_fossils;
        _particle_taxsets = other._particle_taxsets;
        _unused_particle_taxsets = other._unused_particle_taxsets;
#endif
    }
    
}

