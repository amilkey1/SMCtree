#pragma once

namespace proj {

    class ForestExtension {
        public:
            
                            ForestExtension();
                                
            void            dock(const Forest::SharedPtr gf, PartialStore::partial_t partial, Lot::SharedPtr lot);
            void            dockSim(const Forest::SharedPtr gf, Lot::SharedPtr lot);
            void            undock();

            tuple<double, pair<unsigned, unsigned>, pair<bool, double>>  chooseNodesToJoin(vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils, vector<bool> &valid_taxsets);
            double          getProposedDelta() const;
            double          getHeight() const;
            double          getLogWeight() const;
            const Node *    getProposedAnc() const;
            const Node *    getProposedLChild() const;
            const Node *    getProposedRChild() const;
                        
            void            addIncrement(double increment);
            void            joinPriorPrior(vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils, vector<bool> &valid_taxsets, map<string, double> &taxset_ages, double clock_rate);
            pair<unsigned, unsigned> chooseTaxaToJoin(double s);
                    
            PartialStore::partial_t getExtensionPartial();
            double getLineageHeight(const Node* nd) ;
            pair<bool, vector<bool>>    checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets);

        private:

            Forest::ConstSharedPtr          _docked_gene_forest;
            double                          _log_weight;
            double                          _proposed_delta;
            Node                            _proposed_anc;
            const Node *                    _proposed_lchild;
            const Node *                    _proposed_rchild;
            Lot::SharedPtr                  _lot;
        
    };
    
    inline ForestExtension::ForestExtension() {
        undock();
    }
    
    inline void ForestExtension::dock(const Forest::SharedPtr gf, PartialStore::partial_t partial, Lot::SharedPtr lot) {
        // Check to make sure this extension was previously undocked
        assert(gf);
        assert(_docked_gene_forest == nullptr);
        
        // Reset the random number generator
        _lot = lot;

        // Attach Forest
        _docked_gene_forest = gf;
        
        _proposed_delta = 0.0;
        _proposed_anc._height = gf->getForestHeight();
        _proposed_anc._edge_length = 0.0;
        _proposed_anc._partials = partial;
    }

    inline void ForestExtension::dockSim(const Forest::SharedPtr gf, Lot::SharedPtr lot) {
        // Check to make sure this extension was previously undocked
        assert(gf);
        assert(_docked_gene_forest == nullptr);
        
        // Reset the random number generator
        _lot = lot;

        // Attach Forest
        _docked_gene_forest = gf;
        
        _proposed_delta = 0.0;
        _proposed_anc._height = gf->getForestHeight();
        _proposed_anc._edge_length = 0.0;
    }
    
    inline void ForestExtension::undock() {
        _lot.reset();
        _docked_gene_forest.reset();
        _log_weight = 0.0;
        _proposed_delta = 0.0;
        _proposed_anc.clearPointers();
        _proposed_anc._number = -2;
        _proposed_anc._name = "fake";
        _proposed_anc._height = 0.0;
        _proposed_anc._partials.reset();
        _proposed_anc._edge_length = 0.0;
//        _proposed_anc._flags = 0;
        _proposed_anc._split.clear();
        _proposed_lchild = nullptr;
        _proposed_rchild = nullptr;
    }

    inline const Node * ForestExtension::getProposedAnc() const {
        return &_proposed_anc;
    }
    
    inline const Node * ForestExtension::getProposedLChild() const {
        return _proposed_lchild;
    }
    
    inline const Node * ForestExtension::getProposedRChild() const {
        return _proposed_rchild;
    }
    
    inline double ForestExtension::getProposedDelta() const {
        return _proposed_delta;
    }
    
    inline double ForestExtension::getHeight() const {
        return _proposed_anc._height;
    }
    
    inline double ForestExtension::getLogWeight() const {
        return _log_weight;
    }

    inline void ForestExtension::addIncrement(double increment) {
        assert (increment > 0.0);
        _proposed_delta += increment;
        _proposed_anc._height += increment;
    }
    
    inline tuple<double, pair<unsigned, unsigned>, pair<bool, double>> ForestExtension::chooseNodesToJoin(vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils, vector<bool> &valid_taxsets) {
        
        double weight_correction = 0.0;

        vector<unsigned> set_counts;
        
        int chosen_taxset = -1;
        
        bool fossil_constraint = false;
        string match_string = "FOSSIL";
        double fossil_age = -1;
        
        pair <unsigned, unsigned> t = make_pair (0, 1); // if there is only choice, 0 and 1 will be chosen
        if (taxset.size() > 0) {
            // if no taxsets, can choose any remaining nodes
            vector<double> set_sizes;
            bool taxa_not_included_in_sets = false;
            unsigned total_nodes_in_sets = 0;
            
            vector<string> names_of_nodes_in_sets;
            for (auto &t:taxset_no_fossils) {
                for (auto &n:t._species_included) {
                    names_of_nodes_in_sets.push_back(n);
                }
            }
            
            // reset unused taxsets
            vector<vector<unsigned>> node_choices;
            
            // figure out which taxon set choices are valid (exist in _lineages)
            for (unsigned t=0; t<taxset_no_fossils.size(); t++) {
                if (valid_taxsets[t]) {
                    node_choices.resize(node_choices.size() + 1);
                    unsigned real_count = 0;
                    for (auto &n:taxset_no_fossils[t]._species_included) {
                        for (unsigned l=0; l<_docked_gene_forest->_lineages.size(); l++) {
                            if (_docked_gene_forest->_lineages[l]->_name == n) {
                                node_choices[node_choices.size()-1].push_back(_docked_gene_forest->_lineages[l]->_position_in_lineages);
                                real_count++;
                                break;
                            }
                        }
                    }
                    // need to use node_choices not taxset[t]._species_included because species_included may include a fossil that hasn't yet been added
                    set_sizes.push_back((unsigned) node_choices.back().size());
                    set_counts.push_back(t);
                    total_nodes_in_sets += (unsigned) node_choices.back().size();
                }
            }
                
            unsigned nnodes = (unsigned) _docked_gene_forest->_lineages.size();
            
            // add unused taxset taxa to names_of_nodes_in_sets
            for (auto &t:unused_taxset_no_fossils) {
                for (auto &name:t._species_included) {
                    if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), name) == names_of_nodes_in_sets.end()) {
                        names_of_nodes_in_sets.push_back(name); // don't include duplicates
                    }
                }
            }
                
            vector<unsigned> nodes_not_in_taxon_sets;
            for (auto &nd:_docked_gene_forest->_lineages) {
                if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), nd->_name) != names_of_nodes_in_sets.end()) {
                    nnodes--;
                }
                else {
                    // node is not in the taxon sets
                    nodes_not_in_taxon_sets.push_back(nd->_position_in_lineages);
                }
            }
                
            if (nnodes > 1) {
                assert (valid_taxsets.back()); // last taxset must be valid
                taxa_not_included_in_sets = true;
                set_counts.push_back((unsigned) taxset.size());
                set_sizes.push_back(nnodes);
                node_choices.push_back(nodes_not_in_taxon_sets);
                total_nodes_in_sets += nnodes;
            }
            
            vector<double> set_probabilities = set_sizes;
            unsigned total_choices = 0;
            
            for (unsigned s=0; s<set_sizes.size(); s++) {
                // node_choices includes only nodes that exist
                unsigned num_choices = (unsigned) (node_choices[s].size() * (node_choices[s].size() - 1)) / 2;
                set_probabilities[s] = num_choices;
                total_choices += num_choices;
            }
            
            double weight_posterior = 0;
            
            for (auto &s:set_probabilities) { //s * s  because you have ex (3/4) probability of choosing a taxset * 3 taxa within that taxset
                weight_posterior += (s * s / total_choices);
                s /= total_choices;
            }
            
            double nlineages = (double) _docked_gene_forest->_lineages.size();
            double weight_prior = 1 / ((nlineages * (nlineages - 1)) / 2);
            weight_posterior = 1 / weight_posterior;
            
            weight_correction = weight_prior / weight_posterior;
            
            assert (set_probabilities.size() > 0);
            
            chosen_taxset = G::multinomialDraw(_lot, set_probabilities);
            unsigned node_choice_index = chosen_taxset;
            
            int taxset_count = -1;
            
            for (unsigned v=0; v<valid_taxsets.size(); v++) {
                if (valid_taxsets[v]) {
                    taxset_count++;
                }
                if (taxset_count == chosen_taxset) {
                    chosen_taxset = v;
                    break;
                }
            }
            
            // node choice index corresponds to the taxset chosen
            unsigned size1 = (unsigned) taxset.size();
            if (node_choice_index < size1) { // only existing taxon sets will have associated fossils; there is no branch length constraint if the node chosen is not in a taxon set
                vector<string> species_included = taxset[node_choice_index]._species_included;
                if (species_included.size() == 3) {
                    // fossil constraint only applies when the taxon set is down to last 2 lineages + fossil (= 3 entries in species_included)
                    // TODO: double check this
                    for (auto &s:species_included) {
                        if (s.find(match_string) != std::string::npos) {
                            fossil_constraint = true;
                            string fossil_name = s;
                            for (auto &f:particle_fossils) {
                                if (f._name + "_FOSSIL" == fossil_name) {
                                    fossil_age = f._age;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            
            t = chooseTaxaToJoin(set_sizes[node_choice_index]);
            
            t.first = _docked_gene_forest->_lineages[node_choices[node_choice_index][t.first]]->_position_in_lineages;
            t.second = _docked_gene_forest->_lineages[node_choices[node_choice_index][t.second]]->_position_in_lineages;
        }
        else {
            unsigned nlineages = (unsigned) _docked_gene_forest->_lineages.size();
            if (nlineages > 2) {
                t = chooseTaxaToJoin(nlineages);
            }
        }
        
        pair<bool, double> fossil_info = make_pair(fossil_constraint, fossil_age);
        return make_tuple(weight_correction, t, fossil_info);
    }

    inline pair<unsigned, unsigned> ForestExtension::chooseTaxaToJoin(double s){
        assert (s>1);
        double nsubtrees = s;
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        
        if (nsubtrees > 2) {
            t1 = _lot->randint(0, nsubtrees-1);
            t2 = _lot->randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = _lot->randint(0, nsubtrees-1);
            }
        }
        assert(t1 < nsubtrees);
        assert (t2 < nsubtrees);

        return make_pair(t1, t2);
    }

    inline void ForestExtension::joinPriorPrior(vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils, vector<bool> &valid_taxsets, map<string, double> &taxset_ages, double clock_rate) {
        vector<unsigned> set_counts;
        
        string match_string = "FOSSIL";
        
        // Choose the two nodes to join
         tuple<double, pair<unsigned, unsigned>, pair<bool, double>> chosen = chooseNodesToJoin(taxset, unused_taxset, taxset_no_fossils, unused_taxset_no_fossils, particle_fossils, valid_taxsets);
        
        bool update_unused = false;
        double weight_correction = get<0>(chosen);
        bool fossil_constraint = get<2>(chosen).first;
        double fossil_age = get<2>(chosen).second;
        
        // Get pointers to the two nodes to join
        _proposed_lchild = _docked_gene_forest->_lineages[get<1>(chosen).first];
        _proposed_rchild = _docked_gene_forest->_lineages[get<1>(chosen).second];
        
        // Set _proposed_anc's split to union of the two child splits
        _proposed_anc._split.resize(G::_ntaxa);
        _proposed_anc._split += _proposed_lchild->_split;
        _proposed_anc._split += _proposed_rchild->_split;
        
        int proposed_anc_number = _docked_gene_forest->_lineages.back()->_number + 1;
        string proposed_anc_name = "node-" + to_string(proposed_anc_number);
        
        string name1 = _proposed_lchild->_name;
        string name2 = _proposed_rchild->_name;
        
        vector<bool> update_these_taxsets;
         vector<bool> update_these_unused_taxsets;
         
         for (auto &t:taxset_no_fossils) {
             bool update = false;
             unsigned count = 0;
             for (auto &n:t._species_included) {
                 if (n == name1 || n == name2) {
                     update_these_taxsets.push_back(true);
                     break;
                 }
                 count++;
                 if (count == t._species_included.size()) {
                     update_these_taxsets.push_back(update);
                 }
             }
         }
         
         for (auto &t:unused_taxset_no_fossils) {
             bool update = false;
             unsigned count = 0;
             for (auto &n:t._species_included) {
                 if (n == name1 || n == name2) {
                     update_these_unused_taxsets.push_back(true);
                     break;
                 }
                 count++;
                 if (count == t._species_included.size()) {
                     update_these_unused_taxsets.push_back(update);
                 }
             }
         }
         // update taxsets as needed
         for (unsigned i=0; i<update_these_taxsets.size(); i++) {
             if (update_these_taxsets[i] == true) {
                 taxset_no_fossils[i]._species_included.erase(remove(taxset_no_fossils[i]._species_included.begin(), taxset_no_fossils[i]._species_included.end(), _proposed_lchild->_name));
                 taxset_no_fossils[i]._species_included.erase(remove(taxset_no_fossils[i]._species_included.begin(), taxset_no_fossils[i]._species_included.end(), _proposed_rchild->_name));
                 taxset_no_fossils[i]._species_included.push_back(proposed_anc_name);
                 
                 if (taxset_no_fossils[i]._species_included.size() == 1) {
                     update_unused = true;
                     string name = taxset_no_fossils[i]._name;
                     double tree_height = getLineageHeight(_docked_gene_forest->_lineages.back());
                     taxset_ages[name] = tree_height + _proposed_delta;
                     taxset_no_fossils.erase(taxset_no_fossils.begin() + i);
                 }
             }
         }
         
         for (unsigned i=0; i<update_these_unused_taxsets.size(); i++) {
             if (update_these_unused_taxsets[i] == true) {
                 unused_taxset_no_fossils[i]._species_included.erase(remove(unused_taxset_no_fossils[i]._species_included.begin(), unused_taxset_no_fossils[i]._species_included.end(), _proposed_lchild->_name));
                 unused_taxset_no_fossils[i]._species_included.erase(remove(unused_taxset_no_fossils[i]._species_included.begin(), unused_taxset_no_fossils[i]._species_included.end(), _proposed_rchild->_name));
                 unused_taxset_no_fossils[i]._species_included.push_back(proposed_anc_name);
             }
         }
         
         // if a taxset is down to one taxon, find corresponding unused taxset and replace
         if (update_unused) {
             vector<unsigned> updateable_unused;
             vector<unsigned> updateable_unused_sizes;
             string new_name = proposed_anc_name;
             for (unsigned count=0; count < unused_taxset_no_fossils.size(); count++) {
                 for (auto &n:unused_taxset_no_fossils[count]._species_included) {
                     if (n == new_name) {
                         updateable_unused.push_back(count);
                         updateable_unused_sizes.push_back((unsigned) unused_taxset_no_fossils[count]._species_included.size());
                         break;
                     }
                 }
             }
             
             if (updateable_unused.size() > 0) {
                 // check for overlapping taxa with existing taxsets before adding anything in
                 auto min_it = min_element(updateable_unused_sizes.begin(), updateable_unused_sizes.end());
                 unsigned min_index = (unsigned) std::distance(updateable_unused_sizes.begin(), min_it);
                 vector<string> common_elements;
                 // Find the intersection of the two vectors

                 for (unsigned count = 0; count < taxset_no_fossils.size(); count++) {
                     std::set_intersection(taxset_no_fossils[count]._species_included.begin(), taxset_no_fossils[count]._species_included.end(),
                                           unused_taxset_no_fossils[min_index]._species_included.begin(), unused_taxset_no_fossils[min_index]._species_included.end(),
                                       std::back_inserter(common_elements));
                 }
                     
                 if (common_elements.size() == 0) {
                     // add unused taxset into taxsets
                     taxset_no_fossils.push_back(unused_taxset_no_fossils[updateable_unused[min_index]]);
                     unused_taxset_no_fossils.erase(unused_taxset_no_fossils.begin() + min_index);
                 }
             }
         }
         
         // update non fossil taxsets too
         update_these_taxsets.clear();
         update_these_unused_taxsets.clear();
         update_unused = false;
         
         for (auto &t:taxset) {
             bool update = false;
             unsigned count = 0;
             for (auto &n:t._species_included) {
                 if (n == name1 || n == name2) {
                     update_these_taxsets.push_back(true);
                     break;
                 }
                 count++;
                 if (count == t._species_included.size()) {
                     update_these_taxsets.push_back(update);
                 }
             }
         }
         
         for (auto &t:unused_taxset) {
             bool update = false;
             unsigned count = 0;
             for (auto &n:t._species_included) {
                 if (n == name1 || n == name2) {
                     update_these_unused_taxsets.push_back(true);
                     break;
                 }
                 count++;
                 if (count == t._species_included.size()) {
                     update_these_unused_taxsets.push_back(update);
                 }
             }
         }
         // update taxsets as needed
         for (unsigned i=0; i<update_these_taxsets.size(); i++) {
             if (update_these_taxsets[i] == true) {
                 taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), _proposed_lchild->_name));
                 taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), _proposed_rchild->_name));
                 taxset[i]._species_included.push_back(proposed_anc_name);
                 
                 unsigned n_non_fossil_lineages = 0;
                 for (auto &t:taxset[i]._species_included) {
                     string match_string = "FOSSIL";
                     if (t.find(match_string) == std::string::npos) {
                         n_non_fossil_lineages++;
                     }
                 }
                     
                 if (n_non_fossil_lineages == 1) {
                     update_unused = true;
                     string name = taxset[i]._name;
                     double tree_height = getLineageHeight(_docked_gene_forest->_lineages.back());
                     taxset_ages[name] = tree_height + _proposed_delta;
                     taxset.erase(taxset.begin() + i);
                 }
             }
         }
         
         for (unsigned i=0; i<update_these_unused_taxsets.size(); i++) {
             if (update_these_unused_taxsets[i] == true) {
                 unused_taxset[i]._species_included.erase(remove(unused_taxset[i]._species_included.begin(), unused_taxset[i]._species_included.end(), _proposed_lchild->_name));
                 unused_taxset[i]._species_included.erase(remove(unused_taxset[i]._species_included.begin(), unused_taxset[i]._species_included.end(), _proposed_rchild->_name));
                 unused_taxset[i]._species_included.push_back(proposed_anc_name);
             }
         }
         
         // if a taxset is down to one taxon, find corresponding unused taxset and replace
         if (update_unused) {
             vector<unsigned> updateable_unused;
             vector<unsigned> updateable_unused_sizes;
             string new_name = proposed_anc_name;
             for (unsigned count=0; count < unused_taxset.size(); count++) {
                 for (auto &n:unused_taxset[count]._species_included) {
                     if (n == new_name) {
                         updateable_unused.push_back(count);
                         updateable_unused_sizes.push_back((unsigned) unused_taxset[count]._species_included.size());
                         break;
                     }
                 }
             }
             
             if (updateable_unused.size() > 0) {
                 // check for overlapping taxa with existing taxsets before adding anything in
                 auto min_it = min_element(updateable_unused_sizes.begin(), updateable_unused_sizes.end());
                 unsigned min_index = (unsigned) std::distance(updateable_unused_sizes.begin(), min_it);
                 vector<string> common_elements;
                 // Find the intersection of the two vectors

                 for (unsigned count = 0; count < taxset.size(); count++) {
                     std::set_intersection(taxset[count]._species_included.begin(), taxset[count]._species_included.end(),
                                           unused_taxset[min_index]._species_included.begin(), unused_taxset[min_index]._species_included.end(),
                                       std::back_inserter(common_elements));
                 }
                     
                 if (common_elements.size() == 0) {
                     // add unused taxset into taxsets
                     taxset.push_back(unused_taxset[updateable_unused[min_index]]);
                     unused_taxset.erase(unused_taxset.begin() + min_index);
                 }
             }
         }


        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }
        else {
            if (fossil_constraint) {
                // TODO: I don't think this works if > 2 taxa in clade
                assert (fossil_age != -1);
                // figure out if associated branch length has violated the fossil constraint
                // node must be at least as deep as the fossil, which sets the minimum age for the group
                if ((getLineageHeight(_proposed_lchild) + _proposed_delta) < fossil_age) {
//                if ((getLineageHeight(_docked_gene_forest->_lineages.back()->_left_child) + _proposed_delta )< fossil_age) {
                    // fossil age is violated
                    _log_weight = -1 * G::_infinity;
                }
                else if (G::_start_mode != "sim") {
                    if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                        _log_weight = _docked_gene_forest->calcPartialArrayLazyJC(&_proposed_anc, _proposed_lchild, _proposed_rchild, clock_rate);
                    }
                    else {
                        _log_weight = _docked_gene_forest->calcPartialArrayLazyHKY(&_proposed_anc, _proposed_lchild, _proposed_rchild, clock_rate);
                    }
                    _log_weight += weight_correction;
                }
                else {
                    _log_weight = 0.0; // weight is 0 if simulation
                }
            }
            else if (G::_start_mode != "sim") {
                // Compute partial likelihood array of ancestral node
                if (G::_model_type == G::ModelType::MODEL_TYPE_JC) {
                    _log_weight = _docked_gene_forest->calcPartialArrayLazyJC(&_proposed_anc, _proposed_lchild, _proposed_rchild, clock_rate);
                }
                else {
                    _log_weight = _docked_gene_forest->calcPartialArrayLazyHKY(&_proposed_anc, _proposed_lchild, _proposed_rchild, clock_rate);
                }
                _log_weight += weight_correction;
            }
            else {
                _log_weight = 0.0; // weight is 0 if simulation
            }
        }
    }
    
    inline PartialStore::partial_t ForestExtension::getExtensionPartial() {
        return _proposed_anc._partials;
    }

    inline double ForestExtension::getLineageHeight(const Node* nd) {
        if (nd != nullptr) {
            double sum_height = 0.0;
            
            sum_height += nd->getEdgeLength();
            if (nd->_left_child) {
                for (Node* child = nd->_left_child; child; child=child->_left_child) {
                    sum_height += child->getEdgeLength();
                }
            }
            return sum_height;
        }
        else {
            return 0.0;
        }
    }

 }


