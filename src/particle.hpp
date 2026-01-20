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
        void                                    simProposal(unsigned step_number);
        void                                    showParticle();
        double                                  getTreeHeight();
        double                                  getTreeLength();
        double                                  getEstLambda() {return _estimated_lambda;}
        double                                  getEstMu() {return _estimated_mu;}
        double                                  getEstRootAge() {return _estimated_root_age;}
        double                                  getLogLikelihood();
        double                                  getYuleModel();
        double                                  getAllPriors();
        map<string, double>                     getTaxsetAges();
        double                                  getBirthDeathModel();
        void                                    setStartingLogLikelihoods(vector<double> starting_log_likelihoods);
        void                                    clearPartials();
        void                                    drawClockRate();
        double                                  getClockRate() {return _clock_rate;}
        void                                    setSimClockRate();
        void                                    createTrivialForest();
        void                                    setClockRate(double rate) {_clock_rate = rate;}
        void                                    simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites);
        string                                  makeNewick(unsigned precision, bool use_names);
        string                                  saveForestNewick() {
                                                        return _forest_ptr->makeNewick(8, true);}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}
        void                                    showTaxaJoined();

        string debugSaveParticleInfo(unsigned i) const;
        double getHeightFirstSplit(){return _forest_ptr->getHeightFirstSplit();}
        double getHeightSecondIncr() {return _forest_ptr->getHeightSecondIncr();}
        double getHeightThirdIncr() {return _forest_ptr->getHeightThirdIncr();}

        void setFossils() {_particle_fossils = G::_fossils;}
        void setTaxSetsNoFossils();
        void drawFossilAges();
        void setParticleTaxSets();
        void setOverlappingTaxSets();
    
        void    finalizeLatestJoin(unsigned index, map<const void *, list<unsigned> > & nonzero_map);
        void    finalizeThisParticle();
        Forest::SharedPtr getForestPtr() {return _forest_ptr;}
    
        double  getPartialCount() {return _total_particle_partials;}
        void    resetSubgroupPointers(Forest::SharedPtr gene_forest_copy);
    
        void    drawBirthDiff();
        void    drawTurnover();
        void    calculateLambdaAndMu();
        void    drawLambda();
        void    drawRootAge();

    
        // validation stuff
        double  getRootAge();
    
        vector<double>  getHeightNodesWithFossilCalibrations(){return _node_heights_with_fossil_calibrations;}
    
    
    private:
        mutable                                 Lot::SharedPtr _lot;
        void                                    clear();
        
        mutable ForestExtension                 _forest_extension;
        Forest::SharedPtr                       _forest_ptr;
    
        double                                  _log_weight;
        unsigned                                _fossil_number;
        vector<Fossil>                          _particle_fossils; // each particle needs it own set of fossils with their own ages
    
        vector<TaxSet>                          _particle_taxsets; // update this as nodes are joined - for node ages only
        vector<TaxSet>                          _unused_particle_taxsets; // if there are overlapping taxa in taxsets, put the largest groups here until they can be used - for node ages only
    
        vector<TaxSet>                          _particle_taxsets_no_fossils; // update this as nodes are joined
        vector<TaxSet>                          _unused_particle_taxsets_no_fossils; // if there are overlapping taxa in taxsets, put the largest groups here until they can be used
        vector<bool>                            _valid_taxsets;
        map<string, double>                     _taxset_ages;
        double                                  _clock_rate;
        double                                  _estimated_lambda;
        double                                  _estimated_mu;
        double                                  _estimated_root_age;
        double                                  _estimated_birth_difference;
        double                                  _turnover;
    
        double                                  _prev_log_likelihood;
        unsigned                                _total_particle_partials = 0;
        vector<double>                          _node_heights_with_fossil_calibrations;
};

    inline string Particle::debugSaveParticleInfo(unsigned i) const {
        string s;
        s += str(format("Particle %d:\n") % i);
        s += str(format("  _log_weight = %.9f:\n") % _log_weight);
        s += _forest_ptr->debugSaveForestInfo();
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
        _fossil_number = 0;
        _prev_log_likelihood = 0.0;
        _clock_rate = 0.0;
        _forest_ptr = nullptr;
        _total_particle_partials = 0.0;
        _node_heights_with_fossil_calibrations.clear();
        _estimated_lambda = G::_lambda;
        _estimated_mu = G::_mu;
        _estimated_root_age = G::_root_age;
        _estimated_birth_difference = 0.0;
        _turnover = 0.0;
      }

    inline void Particle::setParticleData(Data::SharedPtr d, bool partials) {
        // one forest contains information about all loci since all loci share a tree
        _prev_log_likelihood = 0.0;
        _forest_ptr = Forest::SharedPtr(new Forest());
        Forest::SharedPtr gfp = _forest_ptr;
        gfp->setForestData(d, partials);
    }

    inline void Particle::createTrivialForest() {
        _forest_ptr = Forest::SharedPtr(new Forest());
        _forest_ptr->createTrivialForest();
    }

    inline void Particle::setTaxSetsNoFossils() {
        if (_particle_taxsets.size() > 0) {
            string search_string = "FOSSIL"; // only add taxa that are not fossils to these taxsets
            _particle_taxsets_no_fossils = _particle_taxsets;
            _unused_particle_taxsets_no_fossils = _unused_particle_taxsets;
            
            for (auto &p:_particle_taxsets_no_fossils) {
                for (auto &t:p._species_included) {
                    // remove elements containing the word FOSSIL
                    p._species_included.erase(std::remove_if(p._species_included.begin(), p._species_included.end(),
                                                      [&](const string& s) {
                                                          return s.find(search_string) != string::npos;
                                                      }),
                                              p._species_included.end());
                }
            }
            for (auto &p:_unused_particle_taxsets_no_fossils) {
                for (auto &t:p._species_included) {
                    // remove elements containing the word FOSSIL
                    p._species_included.erase(std::remove_if(p._species_included.begin(), p._species_included.end(),
                                                      [&](const string& s) {
                                                          return s.find(search_string) != string::npos;
                                                      }),
                                              p._species_included.end());
                }
            }
        }
    }

    inline vector<double> Particle::calcGeneTreeLogLikelihoods() {
        vector<double> gene_forest_likelihoods;
        gene_forest_likelihoods.resize(G::_nloci);
        
        //calculate likelihood for each gene tree
        for (unsigned i=0; i<G::_nloci; i++) {
            double gene_tree_log_likelihood  = 0.0;
            gene_tree_log_likelihood = _forest_ptr->calcSubsetLogLikelihood(i);
            assert(!isnan (gene_tree_log_likelihood));
            gene_forest_likelihoods[i] = gene_tree_log_likelihood;
        }

        return gene_forest_likelihoods;
    }

    inline void Particle::simProposal(unsigned step_number) {
        proposal(step_number);
    }

    inline void Particle::proposal(unsigned step_number) {
        if (step_number == 0) {
            // set parameters and taxon sets
            if (G::_est_clock_rate) {
                drawClockRate();
            }
            else {
                setClockRate(G::_clock_rate);
            }
            
            if (G::_est_lambda) {
                if (G::_mu > 0) {
                    assert (G::_est_mu);
                    drawBirthDiff();
                    drawTurnover();
                    calculateLambdaAndMu();
                }
                else {
                    // Yule model
                    drawLambda();
                }
            }
            
            setParticleTaxSets();
            setOverlappingTaxSets();
            setTaxSetsNoFossils();
            
            if (G::_fossils.size() > 0) {
                setFossils();
                drawFossilAges();
            }
            
            if (G::_est_root_age) {
                drawRootAge();
            }
            
            pair<bool, vector<bool>> valid_output = _forest_ptr->checkForValidTaxonSet(_particle_taxsets_no_fossils, _unused_particle_taxsets_no_fossils);
            assert (valid_output.first); // there should always be a valid taxon set if things have been merged correctly
            _valid_taxsets = valid_output.second;
            
        }
        
        if (G::_start_mode != "sim") {
            _forest_extension.dock(_forest_ptr, _forest_ptr->pullPartial(), _lot);
        }
        else {
            _forest_extension.dockSim(_forest_ptr, _lot);
        }
        
        double increment = _forest_ptr->drawBirthDeathIncrement(_lot, -1, _estimated_lambda, _estimated_mu, _estimated_root_age);
        _forest_extension.addIncrement(increment);

        // save heights of each clade with a fossil calibration
        unsigned size_before = (unsigned) _particle_taxsets.size();
        
        _forest_extension.joinPriorPrior(_particle_taxsets, _unused_particle_taxsets, _particle_taxsets_no_fossils, _unused_particle_taxsets_no_fossils, _particle_fossils, _valid_taxsets, _taxset_ages, _clock_rate);
        _total_particle_partials++;
        
        unsigned size_after = (unsigned) _particle_taxsets.size();
        
        if (size_after != size_before) {
            assert (size_after == size_before - 1);
            double proposed_delta = _forest_extension.getProposedDelta();
            _node_heights_with_fossil_calibrations.push_back(_forest_ptr->_tree_height + proposed_delta);
        }
        
        if (step_number == G::_ntaxa - 2) {
            // if we are on the last step, check that the forest is down to 2 lineages (because last two lineags will be joined in the finalizing step in filtering)
            assert (_forest_ptr->_lineages.size() == 2);
        }
        
        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }
        else {
            _log_weight = _forest_extension.getLogWeight();
        }
    }

    inline void Particle::showParticle() {
        //print out weight of each particle
        output("\nParticle:\n", 1);
        output(format("log weight: %d") % _log_weight, 1);
        output("\nForest:\n", 1);
        _forest_ptr->showForest();
    }

    inline double Particle::getLogLikelihood() {
        double log_likelihood = 0.0;
        for (auto &l:_forest_ptr->_gene_tree_log_likelihoods) {
            log_likelihood += l;
        }
        
        return log_likelihood;
    }

    inline double Particle::getTreeHeight() {
        return _forest_ptr->getTreeHeight();
    }

    inline double Particle::getTreeLength() {
        return _forest_ptr->getTreeLength();
    }

    inline double Particle::getYuleModel() {
        return _forest_ptr->getTreePrior();
    }

    inline double Particle::getAllPriors() {
        double tree_prior = _forest_ptr->getTreePrior();
        double param_prior = 0.0;
        
        // these params are all drawn from exponential distributions
        // this also assumes the mean of the exponential distribution is the user-specified param
        if (G::_est_mu && G::_mu > 0.0) {
            param_prior += log(G::_mu) - (_estimated_mu * G::_mu);
        }
        
        if (G::_est_lambda) {
            param_prior += log(G::_lambda) - (_estimated_lambda * G::_lambda);
        }
        
        if (G::_est_root_age) {
            param_prior += log(G::_root_age) - (_estimated_root_age * G::_root_age);
        }
                
        double total_prior = tree_prior + param_prior;
        return total_prior;
    }

    inline double Particle::getBirthDeathModel() {
        return _forest_ptr->getTreePrior();
    }

    inline void Particle::setStartingLogLikelihoods(vector<double> starting_log_likelihoods) {
        _forest_ptr->_gene_tree_log_likelihoods = starting_log_likelihoods;
        
        double total = 0;
        for (auto &s:starting_log_likelihoods) {
            total += s;
        }
        
        _prev_log_likelihood = total;
    }

    inline void Particle::clearPartials() {
        _forest_ptr->clearPartials();
    }

    inline void Particle::drawBirthDiff() {
        assert (G::_est_lambda);
        assert (G::_est_mu);
            
        // birth diff = lambda - mu
        // Gamma(1, n) = Exp(1/n)
        // mean = n
        // for now, n = (G::_lambda - G::_mu) set by user
        double mean = G::_lambda - G::_mu;
        _estimated_birth_difference = _lot->gamma(1, mean);
    }

    inline void Particle::drawLambda() {
        // Yule model; don't estimate mu since it is fixed at 0.0
        // Gamma(1, n) = Exp(1/n)
        // mean = n
        // for now, n = G::_lambda set by user
        _estimated_lambda = _lot->gamma(1, G::_lambda);
    }

    inline void Particle::drawFossilAges() {
        // draw an age for each fossil, using the upper and lower bounds as specified
         for (auto &f:_particle_fossils) {
             if ((f._lower - f._upper) < 0.001) { // if upper and lower are equal, set fossil age
                 f._age = f._lower;
             }
             else {
                f._age = _lot->uniformConstrained(f._lower, f._upper); // TODO: unsure if this lot function is working correctly
             }
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
                // TODO: if two taxsets are equal except for the fossil, combine them and use the older fossil
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
        return _taxset_ages;
    }

    inline void Particle::drawClockRate() {
        _clock_rate = _lot->gamma(1, G::_clock_rate);
    }

    inline void Particle::setSimClockRate() {
        _clock_rate = G::_sim_clock_rate;
    }

    inline void Particle::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        _forest_ptr->simulateData(lot, data, starting_site, nsites, _clock_rate);
    }

    inline string Particle::makeNewick(unsigned precision, bool use_names) {
        return _forest_ptr->makeNewick(precision, use_names);
    }

    inline double Particle::getRootAge() {
        return _estimated_root_age;
    }

    inline void Particle::finalizeThisParticle() {
        // Makes join closest to leaf-level in _forest_extension
        // permanent, then undocks _forest_extension
        
        // Get reference to gene forest extension for this locus
        ForestExtension & gfx = _forest_extension;
        
        // Get pointer to gene forest for this locus
        Forest::SharedPtr gfp = _forest_ptr;
        
        // Copy log likelihood
        gfp->setLogLikelihood(_prev_log_likelihood + gfx.getLogWeight());
                        
        // Get splits for children of _proposed_anc
        const Node * anc = gfx.getProposedAnc();
        assert(anc);
        const Node * lchild = gfx.getProposedLChild();
        assert(lchild);
        const Node * rchild = gfx.getProposedRChild();
        assert(rchild);
        Split lsplit = lchild->_split;
        Split rsplit = rchild->_split;
        
        assert(anc->_split.isEquivalent(lsplit + rsplit));
        
        // Recreate extension's join in the actual gene forest
        double incr = gfx.getProposedDelta();
        assert(incr > 0.0);
        
        gfp->addIncrAndJoin(incr, lsplit, rsplit, gfx);
        
        // reset valid taxon sets
        if (gfp->_lineages.size() > 1) {
            // if not on the last step, reset the valid taxon sets
            pair<bool, vector<bool>> valid_output = _forest_ptr->checkForValidTaxonSet(_particle_taxsets_no_fossils, _unused_particle_taxsets_no_fossils);
            assert (valid_output.first); // there should always be a valid taxon set if things have been merged correctly
            _valid_taxsets = valid_output.second;
        }

        // Can now get rid of extension
        _forest_extension.undock();
    }

    inline void Particle::finalizeLatestJoin(unsigned index, map<const void *, list<unsigned> > & nonzero_map) {
        // Makes join closest to leaf-level in _forest_extension
        // permanent, then undocks _forest_extension
        
        // Get reference to gene forest extension for this locus
        ForestExtension & gfx = _forest_extension;
        
        // Get pointer to gene forest for this locus
        Forest::SharedPtr gfp = _forest_ptr;
        
        // If we are not finalizing the last particle for this
        // gene forest object, make a copy that can be modified
        // without affecting other surviving particles
        unsigned nz = (unsigned)nonzero_map[gfp.get()].size();
        if (nz > 1) {
            // Remove the element corresponding to index
            list<unsigned> & v = nonzero_map[gfp.get()];
            auto it = find(v.begin(), v.end(), index);
            assert(it != v.end());
            v.erase(it);
            
            // Make a copy of the object pointed to by gfp
            Forest::SharedPtr gfcpy = Forest::SharedPtr(new Forest());
            *gfcpy = *gfp;
            _forest_ptr = gfcpy;
            
            // Let gpf point to the copy
            gfp = gfcpy;
        }
        
        // Copy log likelihood
        gfp->setLogLikelihood(_prev_log_likelihood + gfx.getLogWeight());
                        
        // Get splits for children of _proposed_anc
        const Node * anc = gfx.getProposedAnc();
        assert(anc);
        const Node * lchild = gfx.getProposedLChild();
        assert(lchild);
        const Node * rchild = gfx.getProposedRChild();
        assert(rchild);
        Split lsplit = lchild->_split;
        Split rsplit = rchild->_split;
        
        assert(anc->_split.isEquivalent(lsplit + rsplit));
        
        // Recreate extension's join in the actual gene forest
        double incr = gfx.getProposedDelta();
        assert(incr > 0.0);
        
        gfp->addIncrAndJoin(incr, lsplit, rsplit, gfx);
        
        // reset valid taxon sets
        if (gfp->_lineages.size() > 1) {
            // if not on the last step, reset the valid taxon sets
            pair<bool, vector<bool>> valid_output = _forest_ptr->checkForValidTaxonSet(_particle_taxsets_no_fossils, _unused_particle_taxsets_no_fossils);
            assert (valid_output.first); // there should always be a valid taxon set if things have been merged correctly
            _valid_taxsets = valid_output.second;
        }

        // Can now get rid of extension
        _forest_extension.undock();
    }

    inline void Particle::resetSubgroupPointers(Forest::SharedPtr forest_copy) {
        Forest::SharedPtr gfp = _forest_ptr;
        Forest::SharedPtr gfcpy = forest_copy;
        
        *gfcpy = *gfp;
        _forest_ptr = gfcpy;
        
        // let gpf point to the copy
        gfp = gfcpy;
    }

    inline void Particle::showTaxaJoined() {
        cout << _forest_extension.getProposedLChild()->_name << "\t" << "\t" << _forest_extension.getProposedRChild()->_name << endl;
    }

    inline void Particle::drawTurnover() {
        assert (G::_est_mu);
        assert (G::_est_lambda);
        // turnover = mu / lambda
        // Uniform distribution - must be (0, 1)
        
        _turnover = _lot->uniform();
    }

    inline void Particle::calculateLambdaAndMu() {
        _estimated_lambda = _estimated_birth_difference / (1 - _turnover);
        _estimated_mu = (_turnover * _estimated_birth_difference) / (1 - _turnover);
    }
    inline void Particle::drawRootAge() {
        assert (G::_est_root_age > 0.0);
        double max_fossil_age = -1;
        if (_particle_fossils.size() > 0) {
            max_fossil_age = _particle_fossils.back()._age;
        }
            
        bool done = false;
        while (!done) {
            // Gamma(1, n) = Exp(1/n)
            // mean = n
            // for now, n = G::_root_age set by user
//            _estimated_root_age = _lot->gamma(1, G::_root_age);
            // TODO: trying uniform root age to match rev bayes
            _estimated_root_age = _lot->uniformConstrained(G::_root_age_min, G::_root_age_max);
            assert (_estimated_root_age >= G::_root_age_min);
            assert (_estimated_root_age <= G::_root_age_max);
            if (max_fossil_age != -1) {
                if (_estimated_root_age > max_fossil_age) {
                    done = true;
                }
            }
            else {
                done = true;
            }
        }
    }


    inline void Particle::operator=(const Particle & other) {
        _log_weight = other._log_weight;
        _clock_rate = other._clock_rate;
        _fossil_number = other._fossil_number;
        _particle_fossils = other._particle_fossils;
        _particle_taxsets = other._particle_taxsets;
        _unused_particle_taxsets = other._unused_particle_taxsets;
        _particle_taxsets_no_fossils = other._particle_taxsets_no_fossils;
        _unused_particle_taxsets_no_fossils = other._unused_particle_taxsets_no_fossils;
        _valid_taxsets = other._valid_taxsets;
        _taxset_ages = other._taxset_ages;
        _prev_log_likelihood = other._prev_log_likelihood;
        _total_particle_partials = other._total_particle_partials;
        _node_heights_with_fossil_calibrations = other._node_heights_with_fossil_calibrations;
        _estimated_lambda = other._estimated_lambda;
        _estimated_mu = other._estimated_mu;
        _estimated_root_age = other._estimated_root_age;
        _estimated_birth_difference = other._estimated_birth_difference;
        _turnover = other._turnover;
        
        // undock forest extension
        _forest_extension.undock();
        
        _forest_ptr = other._forest_ptr;
    }
    
}
