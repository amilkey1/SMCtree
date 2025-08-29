#pragma once

extern proj::PartialStore ps;

namespace proj {

class Forest {
    
    friend class Particle;
    
    public:
        Forest();
        ~Forest();
        Forest(const Forest & other);
        string makeNewick(unsigned precision, bool use_names);
        void operator=(const Forest & other);
        
        //POL added below
        void createTrivialForest();
        void simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites);
        void buildYuleTree();
        void buildBirthDeathTree();
        double getLogLikelihood();
    
#if defined (INCREMENT_COMPARISON_TEST)
    void buildBirthDeathTreeTest();
    void addBirthDeathTreeIncrementTest(Lot::SharedPtr lot);
#endif
    
        typedef shared_ptr<Forest> SharedPtr;
        string debugSaveForestInfo() const;
        //POL added above
    
    private:

        void clear();
        void setForestData(Data::SharedPtr d, bool partials);
        Node * findNextPreorder(Node * nd);
        void refreshPreorder();
        double calcSubsetLogLikelihood(unsigned i);
        void addIncrement(Lot::SharedPtr lot);
        pair<double, bool> joinTaxa(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset);
        pair<double, bool> joinPriorPrior(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset);
        double joinPriorPost(Lot::SharedPtr lot);
        pair<pair<Node*, Node*>, double> chooseAllPairs(Lot::SharedPtr lot);
        void calcPartialArray(Node* new_nd);
        double calcTransitionProbability(Node* child, double s, double s_child, unsigned locus);
        
        //POL added below
        Node * pullNode();
        void joinRandomLineagePair(Lot::SharedPtr lot);
        void advanceAllLineagesBy(double dt);
        void scaleAllEdgeLengthsBy(double scaling_factor);
        void renumberInternals();
        unsigned getNumLineages() const;
        unsigned getNumNodes() const;
        double calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
        //POL added above

        pair<unsigned, unsigned> chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        string makePartialNewick(unsigned precision, bool use_names);
        void showForest();
        void updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void revertNodeVector(vector<Node *> & node_vector, Node * addnode1, Node * addnode2, Node * delnode1);
        tuple<Node*, Node*, Node*> createNewSubtree(pair<unsigned, unsigned> t);
        double getRunningSumChoices(vector<double> &log_weight_choices);
        vector<double> reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood);
        int selectPair(vector<double> weight_vec, Lot::SharedPtr lot);
        void drawBirthDiff(Lot::SharedPtr lot);
        void drawTurnover(Lot::SharedPtr lot);
        void drawRootAge(Lot::SharedPtr lot);
        void calculateLambdaAndMu();
        void drawLambda(Lot::SharedPtr lot);

        double getTreeHeight();
        double getTreeLength();
        double getTreePrior();
        double calcTopologyPrior(unsigned nlineages);
        void clearPartials();
        double getLineageHeight(Node* nd);
        void addYuleTreeIncrement(Lot::SharedPtr lot);
        void addBirthDeathTreeIncrement(Lot::SharedPtr lot);
    
#if defined (FOSSILS)
        bool addIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name);
        bool addBirthDeathIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name);
        double calcTransitionProbabilityFossil(Node* child, double s, double s_child, unsigned locus);
        bool checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets);
#endif
    
        Data::SharedPtr _data;
        vector<Node *> _lineages;
        vector<Node> _nodes;
        vector<Node*> _preorder;
        unsigned _first_pattern;
        unsigned _npatterns;
        vector<double> _gene_tree_log_likelihoods;
        unsigned _ninternals;
        unsigned _nleaves;
        double _log_joining_prob;
        vector<double> _increments;
        vector<pair<Node*, Node*>> _node_choices;
        double _estimated_lambda;
        double _estimated_mu;
        double _estimated_root_age;
        double _estimated_birth_difference;
        double _turnover;
        double _partial_count;
    
#if defined (FOSSILS)
        double _tree_height;
        vector<bool> _valid_taxsets;
#endif
};

    inline string Forest::debugSaveForestInfo() const {
        string s = "  Forest:\n";

        s += str(format("    _nleaves    = %d\n") % _nleaves);
        s += str(format("    _ninternals = %d\n") % _ninternals);
            
        s += "    _lineages = ";
        for (auto nd : _lineages)
            s += str(format(" %d") % nd->_number);
        s += "\n";

        s += "    _preorder = ";
        for (auto nd : _preorder)
            s += str(format(" %d") % nd->_number);
        s += "\n";

        s += "    _nodes:\n";
        for (auto nd : _nodes) {
            s += nd.saveNodeInfo("      ");
            s += "\n";
        }
        return s;
    }

    inline Forest::Forest() {
        clear();
    }

    inline Forest::~Forest() {
    }

    inline void Forest::clear() {
        _data = nullptr;
        _lineages.clear();
        _nodes.clear();
        _first_pattern = 0;
        _npatterns = 0;
        _gene_tree_log_likelihoods.clear();
        _ninternals = 0;
        _log_joining_prob = 0.0;
        _increments.clear();
        _nleaves = 0;
        _node_choices.clear();
        _estimated_lambda = G::_lambda;
        _estimated_mu = G::_mu;
        _estimated_root_age = G::_root_age;
        _estimated_birth_difference = 0.0;
        _turnover = 0.0;
        _partial_count = 0;
#if defined (FOSSILS)
        _tree_height = 0.0;
        _valid_taxsets.clear();
#endif
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::createTrivialForest() {
        assert(G::_ntaxa > 0);
        assert(G::_ntaxa == G::_taxon_names.size());
        clear();
        unsigned nnodes = 2*G::_ntaxa - 1;
#if defined (FOSSILS)
        nnodes = (unsigned) 2*(G::_ntaxa + G::_sim_fossils.size()) - 1;
#endif
        _nodes.reserve(nnodes);
        _nodes.resize(G::_ntaxa);
        _nleaves = G::_ntaxa;
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            string taxon_name = G::_taxon_names[i];
            _nodes[i]._number = (int)i;
            _nodes[i]._name = taxon_name;
            _nodes[i].setEdgeLength(0.0);
            _nodes[i]._height = 0.0;
            _nodes[i]._split.resize(G::_ntaxa);
            _nodes[i]._split.setBitAt(i);
            _lineages.push_back(&_nodes[i]);
        }
        refreshPreorder();
        _partial_count = 0;
    }
    
    inline void Forest::setForestData(Data::SharedPtr d, bool partials) {
        _data = d;
        
        _gene_tree_log_likelihoods.resize(G::_nloci, 0.0);
        
        auto &data_matrix=_data->getDataMatrix();

        unsigned nnodes = 2*G::_ntaxa - 1;
#if defined (FOSSILS)
//        nnodes += G::_fossils.size();
        nnodes = (unsigned) 2*(G::_ntaxa + G::_fossils.size()) - 1;
#endif
        _nodes.reserve(nnodes);
        _nodes.resize(G::_ntaxa);

        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            _nodes[i]._right_sib=0;
            _nodes[i]._name=" ";
            _nodes[i]._left_child=0;
            _nodes[i]._right_sib=0;
            _nodes[i]._parent=0;
            _nodes[i]._number=i;
            _nodes[i]._edge_length=0.0;
            _nodes[i]._accumulated_height = 0.0;
            _nodes[i]._position_in_lineages=i;
            _nodes[i]._partials = nullptr;
            _nodes[i]._name = G::_taxon_names[i];
            _nodes[i]._set_partials = true;
            _lineages.push_back(&_nodes[i]);
        }
        
        _nleaves = G::_ntaxa;

        for (unsigned index = 0; index < G::_nloci; index ++) {
            for (auto &nd:_lineages) {
                if (index == 0) {
                    double npatterns_total = _data->getNumPatterns();
                    nd->_partials=ps.getPartial(npatterns_total*G::_nstates);
                }
                
                Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index);
                _first_pattern = gene_begin_end.first;
                _npatterns = _data->getNumPatternsInSubset(index);
                
                
                if (!nd->_left_child) {

                    if (!G::_save_memory || (G::_save_memory && partials)) { // if save memory setting, don't set tip partials yet
                            for (unsigned p=_first_pattern; p<_npatterns + _first_pattern; p++) {
                                unsigned pp = p;
                                for (unsigned s=0; s<G::_nstates; s++) {
                                    Data::state_t state = (Data::state_t)1 << s;
                                    Data::state_t d = data_matrix[nd->_number][pp];
                                    double result = state & d;
                                    (*nd->_partials)[p*G::_nstates + s] = (result == 0.0 ? 0.0:1.0);
                                }
                            }
                    }
                }
            }
        }
    }

    inline double Forest::calcSubsetLogLikelihood(unsigned i) {
        _gene_tree_log_likelihoods[i] = 0.0;
        
        auto &counts = _data->getPatternCounts();
        _npatterns = _data->getNumPatternsInSubset(i);
        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
        _first_pattern = gene_begin_end.first;
        
        for (auto &nd:_nodes) {
            if (nd._use_in_likelihood) {
                assert (nd._partials != nullptr); // ignore fossils and fake nodes in likelihood calculations
                double log_like = 0.0;
                for (unsigned p=_first_pattern; p<_npatterns + _first_pattern; p++) {
                    double site_like = 0.0;
                    for (unsigned s=0; s<G::_nstates; s++) {
                        double partial = (*nd._partials)[p*G::_nstates+s];
                        site_like += 0.25*partial;
                    }
                    if (site_like == 0) {
                        showForest();
                    }
                    assert(site_like>0);
                    log_like += log(site_like)*counts[p];
                }
                _gene_tree_log_likelihoods[i] += log_like;
                }
        }
        return _gene_tree_log_likelihoods[i];
    }

    inline void Forest::addYuleTreeIncrement(Lot::SharedPtr lot) {
        // Yule tree
        unsigned nlineages = getNumLineages();
        
        double rate = 0.0;
        
        if (G::_est_lambda) {
            rate = nlineages * _estimated_lambda;
        }
        else {
            rate = nlineages * G::_lambda;
        }
        
        double increment = -log(1.0 - lot->uniform())/rate;
        
        for (auto &nd:_lineages) {
            nd->_edge_length += increment;
            nd->_accumulated_height += increment;
        }
        
        // lorad only works if all topologies the same - then don't include the prior on joins because it is fixed
        double increment_prior = (log(rate)-increment*rate);
                    
        _increments.push_back(increment);
    }

    inline void Forest::addBirthDeathTreeIncrement(Lot::SharedPtr lot) {
        // birth death
        double cum_height = getLineageHeight(_lineages.back());
        // Determine number of lineages remaining
        unsigned n = getNumLineages();
        assert(n > 1);
        
        // Draw n-1 internal node heights and store in vector heights
        vector<double> heights(n - 1, 0.0);
        
        double rho = 1.0; // TODO: for now, assume rho = 1.0
        
        double birth_rate = _estimated_lambda;
        
        double death_rate = _estimated_mu;

        double exp_death_minus_birth = exp(death_rate - birth_rate);
        double phi = 0.0;
        phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
        phi += (death_rate - birth_rate)*exp_death_minus_birth;
        phi /= (exp_death_minus_birth - 1.0);
        for (unsigned i = 0; i < n - 2; i++) {
            double u = lot->uniform();
            double y = u/(1.0 + birth_rate*rho*(1.0 - u));
            if (birth_rate > death_rate) {
                y = log(phi - u*rho*birth_rate);
                y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                y /= (death_rate - birth_rate);
            }
            heights[i] = y;
        }
        heights[n-2] = 1.0;
        sort(heights.begin(), heights.end());
        
        // Waiting time to next speciation event is first height
        // scaled so that max height is mu - cum_height
        double t = heights[0]*(_estimated_root_age - cum_height); // TODO: not sure this is right
        
        assert (t > 0.0);
        
        for (auto &nd:_lineages) {
            nd->_edge_length += t;
            nd->_accumulated_height += t;
        }
        
        // lorad only works if all topologies the same - then don't include the prior on joins because it is fixed
        double rate = 0.0; // TODO: need to modify this for birth-death
        rate = getNumLineages() * birth_rate;
        double increment_prior = (log(rate)-t*rate);
        
        _increments.push_back(increment_prior);
    }

#if defined (INCREMENT_COMPARISON_TEST)
    inline void Forest::addBirthDeathTreeIncrementTest(Lot::SharedPtr lot) {
        // birth death
        
            ofstream logf("smc5.log");
        logf << "sample" << "\t" << "increment" << endl;
        double cum_height = 0.0;
        unsigned n = getNumLineages();
        for (unsigned a = 0; a < 1000000; a++) {
            for (unsigned b = 0; b < 5; b++) {
//        double cum_height = getLineageHeight(_lineages.back());
        // Determine number of lineages remaining
//        unsigned n = getNumLineages();
//        assert(n > 1);
//
//                n -= b;
        
        // Draw n-1 internal node heights and store in vector heights
        vector<double> heights(n - 1, 0.0);
        
        double rho = 1.0;
        
        double birth_rate = G::_lambda;
        if (_estimated_lambda > 0.0) {
            birth_rate = _estimated_lambda;
        }
        
        double death_rate = G::_mu;
        if (_estimated_mu > 0.0) {
            death_rate = _estimated_mu;
        }
            _estimated_root_age = 1.0;
        double exp_death_minus_birth = exp(death_rate - birth_rate);
        double phi = 0.0;
        phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
        phi += (death_rate - birth_rate)*exp_death_minus_birth;
        phi /= (exp_death_minus_birth - 1.0);
        for (unsigned i = 0; i < n - 2; i++) {
            double u = lot->uniform();
            double y = u/(1.0 + birth_rate*rho*(1.0 - u));
            if (birth_rate > death_rate) {
                y = log(phi - u*rho*birth_rate);
                y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                y /= (death_rate - birth_rate);
            }
            heights[i] = y;
        }
        heights[n-2] = 1.0;
        sort(heights.begin(), heights.end());
        
        // Waiting time to next speciation event is first height
        // scaled so that max height is mu - cum_height
        double t = heights[0]*(_estimated_root_age - cum_height); // TODO: not sure this is right
        
        assert (t > 0.0);
            
                if (b < 4) {
                    cum_height += t;
                    n -= 1;
                }
                else {
                    cum_height = 0.0;
                    logf << a << "\t";
                    logf << t << endl;
                    n = getNumLineages();
                }
        }
        }
    }
#endif

#if defined (FOSSILS)
    bool Forest::addIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name) {
        bool fossil_added = addBirthDeathIncrementFossil(lot, age, fossil_name);
        return fossil_added;
    }
#endif

#if defined (FOSSILS)
    bool Forest::addBirthDeathIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name) {
        bool fossil_added = false;
        
        // if there is only one lineage left, all extant taxa have been joined and remaining fossils must be added
        
        // birth death
        double cum_height = _tree_height;
        // Determine number of lineages remaining
        unsigned n = getNumLineages();
        
        double birth_rate = _estimated_lambda;
        
        double death_rate = _estimated_mu;
        
        if (n > 1) {
            // Draw n-1 internal node heights and store in vector heights
            vector<double> heights(n - 1, 0.0);
            
            double rho = 1.0; // TODO: for now, assume rho = 1.0
            
            double exp_death_minus_birth = exp(death_rate - birth_rate);
            double phi = 0.0;
            phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
            phi += (death_rate - birth_rate)*exp_death_minus_birth;
            phi /= (exp_death_minus_birth - 1.0);
            for (unsigned i = 0; i < n - 2; i++) {
                double u = lot->uniform();
                double y = u/(1.0 + birth_rate*rho*(1.0 - u));
                if (birth_rate > death_rate) {
                    y = log(phi - u*rho*birth_rate);
                    y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                    y /= (death_rate - birth_rate);
                }
                heights[i] = y;
            }
            heights[n-2] = 1.0;
            sort(heights.begin(), heights.end());
            
            // Waiting time to next speciation event is first height
            // scaled so that max height is mu - cum_height
            double t = heights[0]*(_estimated_root_age - cum_height); // TODO: not sure this is right
            
            assert (t > 0.0);
                    
            if ((t + _tree_height < age) || (age == -1)) { // don't add fossil because next fossil placement is deeper than current tree
                for (auto &nd:_lineages) {
                    nd->_edge_length += t;
                    nd->_accumulated_height += t;
                }
                _tree_height += t;
                _increments.push_back(t);
            }
            
            else {
                // add fossil
                double edge_len = age - _tree_height;
                _tree_height += edge_len;
                _increments.push_back(edge_len);
                
                for (auto &nd:_lineages) {
                    nd->_edge_length += edge_len;
                    nd->_accumulated_height += edge_len;
                }

                //new node is always needed
                Node* new_nd = pullNode();

                new_nd->_name = fossil_name + "_FOSSIL";
                new_nd->_set_partials = false; // do not include this node in likelihood calculation
                // don't add anything to fossil edge length; go back and draw another increment
                new_nd->_position_in_lineages = (unsigned) _lineages.size();
                new_nd->_use_in_likelihood = false;
                _lineages.push_back(new_nd);
                
                fossil_added = true;
            }
        }
            
        else {
            // add fossil because tree is done and remaining fossils must go at the base
            assert (age != -1.0);
            double edge_len = age - _tree_height;
            _tree_height += edge_len;
            _increments.push_back(edge_len);
            
            for (auto &nd:_lineages) {
                nd->_edge_length += edge_len;
                nd->_accumulated_height += edge_len;
            }

            //new node is always needed
            Node* new_nd = pullNode();

            new_nd->_name = fossil_name + "_FOSSIL";
            new_nd->_set_partials = false; // do not include this node in likelihood calculation
            new_nd->_edge_length = edge_len;
            new_nd->_accumulated_height = edge_len;
            new_nd->_position_in_lineages = (unsigned) _lineages.size();
            _lineages.push_back(new_nd);
            
            fossil_added = true;
        }
            
        // lorad only works if all topologies the same - then don't include the prior on joins because it is fixed
//        double rate = 0.0; // TODO: need to modify this for birth-death
//        rate = getNumLineages() * birth_rate;
//        double increment_prior = (log(rate)-t*rate);
//
//        _increments_and_priors.push_back(make_pair(t, increment_prior));
                
        return fossil_added;
    }
#endif

    inline void Forest::addIncrement(Lot::SharedPtr lot) {
#if defined (INCREMENT_COMPARISON_TEST)
        addBirthDeathTreeIncrementTest(lot);
#else
        addBirthDeathTreeIncrement(lot);
#endif
    }

    inline pair<double, bool> Forest::joinTaxa(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset) {
        double log_weight = 0.0;
        pair<double, bool> output;
        
        if (G::_proposal == "prior-prior") {
            output = joinPriorPrior(prev_log_likelihood, lot, taxset, unused_taxset);
        }
        else {
            log_weight = joinPriorPost(lot);
        }

        return output;
    }

    inline double Forest::joinPriorPost(Lot::SharedPtr lot) {
        if (G::_save_memory) {
            double npatterns_total = _data->getNumPatterns();
            for (auto &nd:_lineages) {
                if (nd->_partials == nullptr) {
                    nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                    calcPartialArray(nd);
                }
            }
        }
        
        pair<pair<Node*, Node*>, double> node_pair = chooseAllPairs(lot);
        
        Node* subtree1 = node_pair.first.first;
        Node* subtree2 = node_pair.first.second;
        
        double log_weight = node_pair.second;

        assert (subtree1 != subtree2);

        //new node is always needed
        Node* new_nd = pullNode();

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;

        // calculate new partials
        assert (new_nd->_partials == nullptr);
        double npatterns_total = _data->getNumPatterns();
        new_nd->_partials = ps.getPartial(G::_nstates*npatterns_total);
        assert(new_nd->_left_child->_right_sib);

        if (G::_save_memory) {
            double npatterns_total = _data->getNumPatterns();
            new_nd->_partials = ps.getPartial(npatterns_total*G::_nstates);
            
            for (auto &nd:_lineages) {
                if (nd->_partials == nullptr) {
                    nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                    calcPartialArray(nd);
                }
            }
        }
                
        calcPartialArray(new_nd);

        for (unsigned index = 0; index < G::_nloci; index++) {
            subtree1->_partials=nullptr; // throw away subtree partials now, no longer needed
            subtree2->_partials=nullptr;
        }
        
        //update node lists
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);

        for (unsigned index = 0; index<G::_nloci; index++) {
            calcSubsetLogLikelihood(index);
        }
        
       if (G::_save_memory) {
           for (auto &nd:_nodes) {
               nd._partials = nullptr;
           }
       }
        
        calcTopologyPrior((unsigned) getNumLineages() +1);

        return log_weight;
    }

    inline pair<pair<Node*, Node*>, double> Forest::chooseAllPairs(Lot::SharedPtr lot) {
        // this function tries all possible node pairs for "prior-post" proposal
        double prev_log_likelihood = 0.0;
        double log_weight = 0.0;
        assert (_node_choices.size() == 0);
        
        for (auto &g:_gene_tree_log_likelihoods) {
            prev_log_likelihood += g;
        }
                
        vector<double> log_likelihood_choices;
             
          // choose pair of nodes to try
        for (unsigned i = 0; i < _lineages.size()-1; i++) {
            for (unsigned j = i+1; j < _lineages.size(); j++) {

                // createNewSubtree returns subtree1, subtree2, new_nd
              
                tuple<Node*, Node*, Node*> t = createNewSubtree(make_pair(i,j));
              
                double log_likelihood_of_pair = 0.0;
              
                  for (unsigned index = 0; index<G::_nloci; index++) {
                      log_likelihood_of_pair += calcSubsetLogLikelihood(index);
                  }
              
                log_likelihood_choices.push_back(log_likelihood_of_pair);
              
                // revert _lineages if > 1 choice
                revertNodeVector(_lineages, get<0>(t), get<1>(t), get<2>(t));

                //reset siblings and parents of original nodes back to 0
                get<0>(t)->resetNode(); //subtree1
                get<1>(t)->resetNode(); //subtree2

                // delete new node -- TODO: just reset it and keep using it?
                _nodes.pop_back();
                // clear new node from _nodes
                //clear new node that was just created
//                get<2>(t)->clear(); //new_nd
                _ninternals--;
            }
        }
              
        assert (_node_choices.size() == log_likelihood_choices.size());
        // reweight each choice of pairs
        vector<double> log_weight_choices = reweightChoices(log_likelihood_choices, prev_log_likelihood);

        // sum unnormalized weights before choosing the pair
        // must include the likelihoods of all pairs in the final particle weight
        double log_weight_choices_sum = getRunningSumChoices(log_weight_choices);
        log_weight = log_weight_choices_sum;
        for (unsigned b=0; b < log_weight_choices.size(); b++) {
            log_weight_choices[b] -= log_weight_choices_sum;
        }
             
        // randomly select a pair
        unsigned index_of_choice = selectPair(log_weight_choices, lot);

        // find nodes to join in node_list
        Node* subtree1 = _node_choices[index_of_choice].first;
        Node* subtree2 = _node_choices[index_of_choice].second;
        
        _node_choices.clear();
             
        return make_pair(make_pair(subtree1, subtree2), log_weight);
      }

    inline tuple<Node*, Node*, Node*> Forest::createNewSubtree(pair<unsigned, unsigned> t) {
        Node* subtree1 = _lineages[t.first];
        Node* subtree2 = _lineages[t.second];

        Node* new_nd = pullNode();

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;

        // calculate partials for new node
        double npatterns_total = _data->getNumPatterns();
        new_nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
        calcPartialArray(new_nd);
        
         assert(new_nd->_left_child->_right_sib);
         calcPartialArray(new_nd);

         // update lineages vector
         updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        _node_choices.push_back(make_pair(subtree1, subtree2));
                  
         return make_tuple(subtree1, subtree2, new_nd);
     }

    inline pair<double, bool> Forest::joinPriorPrior(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset) {
        bool filter = false;
        // find the new_nd from the previous step and accumulate height if needed
        
        // if lineages.back() is a fossil, go to the previous node
        unsigned end_node = (unsigned) _lineages.size() - 1;
        vector<unsigned> set_counts;
        
        Node *node_to_check = _lineages[end_node];
        
        // check for the new node that was added in the previous step
        bool done = false;
        while (!done) {
            if (boost::ends_with(node_to_check->_name, "FOSSIL")) {
                end_node --;
                node_to_check = _lineages[end_node];
            }
            else {
                done = true;
            }
        }
        
        int chosen_taxset = -1;

        // TODO: can speed this up? for now, at every step, go through every node and reset accumulated height to 0
        for (auto &nd:_nodes) {
            nd._accumulated_height = nd._edge_length;
            if (!nd._set_partials && !boost::ends_with(nd._name, "FOSSIL")) {
                int next_node = nd._next_real_node;
                if (next_node != -1) {
                    _nodes[next_node]._accumulated_height += nd._edge_length;
                }
            }
        }
        
        unsigned nlineages = getNumLineages();
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        
        if (G::_save_memory) {
            double npatterns_total = _data->getNumPatterns();
            for (auto &nd:_lineages) {
                if (nd->_set_partials) {
                    if (nd->_partials == nullptr) {
                        nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                        calcPartialArray(nd);
                    }
                }
            }
        }
        
        pair <unsigned, unsigned> t = make_pair (0, 1); // if there is only choice, 0 and 1 will be chosen
        if (taxset.size() > 0) {
            // if no taxsets, can choose any remaining nodes
            vector<double> set_sizes;
            bool taxa_not_included_in_sets = false;
            unsigned total_nodes_in_sets = 0;
            
            vector<string> names_of_nodes_in_sets;
            for (auto &t:taxset) {
                for (auto &n:t._species_included) {
                    names_of_nodes_in_sets.push_back(n);
                }
            }
            
            // TODO: reset unused taxsets
            vector<vector<unsigned>> node_choices;
            
            // figure out which taxon set choices are valid (exist in _lineages)
            for (unsigned t=0; t<taxset.size(); t++) {
                if (_valid_taxsets[t]) {
                    node_choices.resize(node_choices.size() + 1);
                    unsigned real_count = 0;
                    for (auto &n:taxset[t]._species_included) {
                        for (unsigned l=0; l<_lineages.size(); l++) {
                            if (_lineages[l]->_name == n) {
                                node_choices[node_choices.size()-1].push_back(_lineages[l]->_position_in_lineages);
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
                
            unsigned nnodes = (unsigned) _lineages.size();
            
            // add unused taxset taxa to names_of_nodes_in_sets
            for (auto &t:unused_taxset) {
                for (auto &name:t._species_included) {
                    if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), name) == names_of_nodes_in_sets.end()) {
                        names_of_nodes_in_sets.push_back(name); // don't include duplicates
                    }
                }
            }
                
            vector<unsigned> nodes_not_in_taxon_sets;
            for (auto &nd:_lineages) {
                if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), nd->_name) != names_of_nodes_in_sets.end()) {
                    nnodes--;
                }
                else {
                    // node is not in the taxon sets
                    nodes_not_in_taxon_sets.push_back(nd->_position_in_lineages);
                }
            }
                
            if (nnodes > 1) {
                assert (_valid_taxsets.back()); // last taxset must be valid
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
            
            for (auto &s:set_probabilities) {
                s /= total_choices;
            }
            
            assert (set_probabilities.size() > 0);
            
            chosen_taxset = G::multinomialDraw(lot, set_probabilities);
            unsigned node_choice_index = chosen_taxset;
            
            int taxset_count = -1;
            for (unsigned v=0; v<_valid_taxsets.size(); v++) {
                if (_valid_taxsets[v]) {
                    taxset_count++;
                }
                if (taxset_count == chosen_taxset) {
                    chosen_taxset = v;
                    break;
                }
            }
            
            t = chooseTaxaToJoin(set_sizes[node_choice_index], lot);
            
            
            subtree1 = _lineages[node_choices[node_choice_index][t.first]];
            subtree2 = _lineages[node_choices[node_choice_index][t.second]];
        }
        else {
            if (nlineages > 2) {
                t = chooseTaxaToJoin(nlineages, lot);
            }
            subtree1 = _lineages[t.first];
            subtree2 = _lineages[t.second];
        }
        
        assert (subtree1 != subtree2);

        //new node is always needed
        Node* new_nd = pullNode();

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;
        
        // if new node has any child that is not a real node (set partials = false) and has no next node (next node != -1), the new node is not a real node either
        bool fake_node = false;
        if (!subtree1->_set_partials) {
            int next_node = subtree1->_next_real_node;
            if (next_node == -1) {
                fake_node = true;
            }
            else {
                if (next_node == -1 && !_nodes[next_node]._set_partials) {
                    fake_node = true;
                }
            }
        }
        
        if (!subtree2->_set_partials) {
            int next_node = subtree2->_next_real_node;
            if (next_node == -1) {
                fake_node = true;
            }
            else {
                if (next_node == -1 && !_nodes[next_node]._set_partials) {
                    fake_node = true;
                }
            }
        }

        if (fake_node) {
            new_nd->_set_partials = false;
            new_nd->_use_in_likelihood = false;
        }

        if (new_nd->_set_partials) {
            // calculate new partials
            assert (new_nd->_partials == nullptr);
            double npatterns_total = _data->getNumPatterns();
            new_nd->_partials = ps.getPartial(G::_nstates*npatterns_total);
            assert(new_nd->_left_child->_right_sib);

            if (G::_save_memory) {
                double npatterns_total = _data->getNumPatterns();
                new_nd->_partials = ps.getPartial(npatterns_total*G::_nstates);
                
                for (auto &nd:_lineages) {
                    if (nd->_partials == nullptr) {
                        nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                        calcPartialArray(nd);
                    }
                }
            }
            calcPartialArray(new_nd);
            
            filter = true; // must filter if a real node has been added
            
            subtree1->_use_in_likelihood = false;
            subtree2->_use_in_likelihood = false;
        }
        
        else {
            if (!subtree1->_set_partials) {
                subtree1->_use_in_likelihood = false;
    //                subtree1->_partials = nullptr;
            }
            if (!subtree2->_set_partials) {
                subtree2->_use_in_likelihood = false;
    //                subtree2->_partials = nullptr;
            }
        }
        
        // if new_nd is a real node, none of its children should be used in the likelihood calculation
        // TODO: faster way to do this?
        
        if (new_nd->_set_partials) {
            for (auto nd:_nodes){
                unsigned node_number = nd._number;
                bool done = false;
                while (!done) {
                    if (nd._parent) {
                        if (nd._parent == new_nd) {
                            _nodes[node_number]._use_in_likelihood = false;
                            done = true;
                        }
                        else {
                            if (nd._parent) {
                                nd = *nd._parent;
                            }
                            else {
                                done = true;
                            }
                        }
                    }
                    else {
                        done = true;
                    }
                }
            }
        }
        
        //update node lists
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);
        
        bool update_unused = false;
        // update taxset if needed
        // TODO: go through all taxsets in existence and look for chosen taxa, then update all of them
        // TODO: because there might be multiple taxsets with the same taxon names
        
        string name1 = subtree1->_name;
        string name2 = subtree2->_name;
        
        vector<bool> update_these_taxsets;
        vector<bool> update_these_unused_taxsets;
        
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
        // TODO: update taxsets as needed
        for (unsigned i=0; i<update_these_taxsets.size(); i++) {
            if (update_these_taxsets[i] == true) {
                taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), subtree1->_name));
                taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), subtree2->_name));
                taxset[i]._species_included.push_back(new_nd->_name);
                
                if (taxset[i]._species_included.size() == 1) {
                    taxset.erase(taxset.begin() + i);
                    update_unused = true;
                }
            }
        }
        
        for (unsigned i=0; i<update_these_unused_taxsets.size(); i++) {
            if (update_these_unused_taxsets[i] == true) {
                unused_taxset[i]._species_included.erase(remove(unused_taxset[i]._species_included.begin(), unused_taxset[i]._species_included.end(), subtree1->_name));
                unused_taxset[i]._species_included.erase(remove(unused_taxset[i]._species_included.begin(), unused_taxset[i]._species_included.end(), subtree2->_name));
                unused_taxset[i]._species_included.push_back(new_nd->_name);
            }
        }
        
//        if (chosen_taxset != -1 && chosen_taxset != _valid_taxsets.size()-1) { // nothing to update if this was not a real taxon set
//            taxset[chosen_taxset]._species_included.erase(remove(taxset[chosen_taxset]._species_included.begin(), taxset[chosen_taxset]._species_included.end(), subtree1->_name));
//            taxset[chosen_taxset]._species_included.erase(remove(taxset[chosen_taxset]._species_included.begin(), taxset[chosen_taxset]._species_included.end(), subtree2->_name));
//            taxset[chosen_taxset]._species_included.push_back(new_nd->_name);
//
//            if (taxset[chosen_taxset]._species_included.size() == 1) {
//                taxset.erase(taxset.begin() + chosen_taxset);
//                update_unused = true;
//            }
//        }
        
        // update unused taxset with new names
//        if (unused_taxset.size() > 0) {
//            string name1 = subtree1->_name;
//            string name2 = subtree2->_name;
//            string new_name = new_nd->_name;
//
//            for (auto &t:unused_taxset) {
//                bool new_name_added = false;
//                for (auto &name:t._species_included) {
//                    if (!new_name_added) {
//                        if (name == name1) {
//                            name = new_name;
//                            new_name_added = true;
//                        }
//                        else if (name == name2) {
//                            name = new_name;
//                            new_name_added = true;
//                        }
//                    }
//                    else {
//                        if (name == name1) {
//                            t._species_included.erase(remove(t._species_included.begin(), t._species_included.end(), name1));
//                        }
//                        else if (name == name2) {
//                            t._species_included.erase(remove(t._species_included.begin(), t._species_included.end(), name2));
//                        }
//                    }
//                }
//            }
//        }
        
        // TODO: find corresponding unused taxset and replace
        // TODO: find unused taxsets with new_name in them
        if (update_unused) {
            vector<unsigned> updateable_unused;
            vector<unsigned> updateable_unused_sizes;
            string new_name = new_nd->_name;
            for (unsigned count=0; count < unused_taxset.size(); count++) {
                for (auto &n:unused_taxset[count]._species_included) {
                    if (n == new_name) { // TODO: this doesn't work if there are multiple overlapping taxsets; need to check for that before adding in the new one
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

        for (unsigned index = 0; index<G::_nloci; index++) {
            calcSubsetLogLikelihood(index);
        }
        
        // after a new node is created, check if it's fake
        // if it is, set the next real node
        if (!new_nd->_set_partials) {
            if (new_nd->_left_child) {
                if (!boost::ends_with(new_nd->_left_child->_name, "FOSSIL")) {
                    // if new node's left child is not a fossil, use it to find the next real node
                    int next_node = new_nd->_left_child->_next_real_node;
                    if (next_node > -1) {
                        new_nd->_next_real_node = new_nd->_left_child->_number;
                    }
                    else {
                        if (new_nd->_left_child->_set_partials) {
                            new_nd->_next_real_node = new_nd->_left_child->_number;
                        }
                    }
                }
                if (!boost::ends_with(new_nd->_left_child->_right_sib->_name, "FOSSIL")) {
                    int next_node = new_nd->_left_child->_right_sib->_next_real_node;
                    if (next_node > -1) {
                        new_nd->_next_real_node = new_nd->_left_child->_right_sib->_number;
                    }
                    else {
                        if (new_nd->_left_child->_right_sib->_set_partials) {
                            new_nd->_next_real_node = new_nd->_left_child->_right_sib->_number;
                        }
                    }
                }
            }
        }
        
        double new_log_likelihood = 0.0;
        for (auto &g:_gene_tree_log_likelihoods) {
            new_log_likelihood += g;
        }
        
        double log_weight = new_log_likelihood - prev_log_likelihood;
               
       if (G::_save_memory) {
           for (auto &nd:_nodes) {
               nd._partials = nullptr;
           }
       }
        
        calcTopologyPrior(getNumLineages() +1);
        
        _valid_taxsets.clear();
        
        return make_pair(log_weight, filter);
    }

    inline void Forest::updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode) {
        // Delete delnode1 from node_vector
        auto it1 = find(node_vector.begin(), node_vector.end(), delnode1);
        assert(it1 != node_vector.end());
        node_vector.erase(it1);

        // Delete delnode2 from node_vector
        auto it2 = find(node_vector.begin(), node_vector.end(), delnode2);
        assert(it2 != node_vector.end());
        node_vector.erase(it2);

        // Add addnode to node_vector
        node_vector.push_back(addnode);

        // reset _position_in_lineages
        for (unsigned i=0; i < _lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
        }
    }

    inline vector<double> Forest::reweightChoices(vector<double> & likelihood_vec, double prev_log_likelihood) {
        vector<double> weight_vec;
        for (int a = 0; a < (int) likelihood_vec.size(); a++) {
            weight_vec.push_back(likelihood_vec[a]-prev_log_likelihood);
        }
        return weight_vec;
    }

    inline int Forest::selectPair(vector<double> weight_vec, Lot::SharedPtr lot) {
         // choose a random number [0,1]
         assert (lot != nullptr);
         double u = lot->uniform();
         
         double cum_prob = 0.0;
         int index = 0.0;
         for (int i=0; i < (int) weight_vec.size(); i++) {
             cum_prob += exp(weight_vec[i]);
             if (u <= cum_prob) {
                 index = i;
                 break;
             }
         }
         // return index of choice
         return index;
     }

    inline double Forest::getRunningSumChoices(vector<double> &log_weight_choices) {
        double running_sum = 0.0;
        double log_weight_choices_sum = 0.0;
        double log_max_weight = *max_element(log_weight_choices.begin(), log_weight_choices.end());
        for (auto & i:log_weight_choices) {
            running_sum += exp(i - log_max_weight);
        }
        log_weight_choices_sum = log(running_sum) + log_max_weight;
        return log_weight_choices_sum;
    }

    inline void Forest::revertNodeVector(vector<Node *> &node_vector, Node *addnode1, Node *addnode2, Node *delnode1) {
        // Delete delnode1 from node_vector
        auto it = find(node_vector.begin(), node_vector.end(), delnode1);
        assert (it != node_vector.end());
        node_vector.erase(it);

        // find positions of nodes to insert
        auto position1 = addnode1->_position_in_lineages;
        auto iter1 = addnode1;

        auto position2 = addnode2->_position_in_lineages;
        auto iter2 = addnode2;

        // lower position must be inserted first
        if (position1 < position2) {
            node_vector.insert(node_vector.begin()+position1, iter1);
            node_vector.insert(node_vector.begin()+position2, iter2);
        }
        else {
            node_vector.insert(node_vector.begin()+position2, iter2);
            node_vector.insert(node_vector.begin()+position1, iter1);
        }

        assert(_lineages[addnode1->_position_in_lineages] == addnode1);
        assert(_lineages[addnode2->_position_in_lineages] == addnode2);

        // reset _position_in_lineages
        for (unsigned i=0; i < _lineages.size(); i++) {
            _lineages[i] -> _position_in_lineages=i;
        }
    }

    inline void Forest::updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode) {
        // Delete delnode1 from node_list
        auto it1 = find(node_list.begin(), node_list.end(), delnode1);
        assert(it1 != node_list.end());
        node_list.erase(it1);

        // Delete delnode2 from node_list
        auto it2 = find(node_list.begin(), node_list.end(), delnode2);
        assert(it2 != node_list.end());
        node_list.erase(it2);

        // Add addnode to node_list
        node_list.push_back(addnode);
    }

    inline double Forest::getTreeHeight() {
        return _tree_height;
    }

    inline double Forest::getTreeLength() {
        // sum of all edge lengths in tree
        double sum_height = 0.0;

        for (auto &nd:_nodes) {
            // sum edge lengths from all nodes
            sum_height += nd._edge_length;
        }
        return sum_height;
    }

    inline double Forest::getTreePrior() {
        double birth_death_prior = 0.0;
        
        // https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon)/10%3A_Introduction_to_Birth-Death_Models/10.02%3A_The_Birth-Death_Model
        // Raup 1985
        
        // starting from the root, go increment by increment and calculate the following parameters:
        // TODO: do fossil increments count?
        
        double E = _estimated_mu / _estimated_lambda;
        double r = _estimated_lambda - _estimated_mu;
        double t = 0.0;
        unsigned count = (unsigned) _increments.size() - 1;
        
        for (unsigned i=1; i<G::_fossils.size() + G::_ntaxa + 1; i++) {
            // i is number of species in existence
            // density increment is the probability of having exactly i species at time t
            // TODO: is t the increment or the accumulated height?
            
            t = _increments[count]; // TODO: = or += ?
            double a = (E*(exp(r)*t - 1)) / (exp(r)*t - E);
//            double B = (exp(r)*t - 1) / (exp(r)*t - E);
            double B = a / E;
            
            if (birth_death_prior == 0) {
                double density_increment = a; // starting with a single lineage
                birth_death_prior += density_increment;
            }
            else {
                double density_increment = (1-a)*(1-B)*pow(B, i-1);
                birth_death_prior += density_increment;
            }
            
//            double log_prob_density_increment = log((1-a)*(1-B)*pow(B, i-1));
            
//            birth_death_prior += log_prob_density_increment;
            count--;
        }
        
        birth_death_prior = log(birth_death_prior);
        
//        assert(birth_death_prior == birth_death_prior); // TODO: why does this turn into Nan?
        // check for NaN
        return birth_death_prior;
    }

    inline double Forest::calcTopologyPrior(unsigned nlineages) {
        _log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        assert (!isinf(_log_joining_prob));
        return _log_joining_prob;
    }

    inline void Forest::refreshPreorder() {
       // Create vector of node pointers in preorder sequence
        Node *nd = &_nodes.back();
       _preorder.clear();
       _preorder.reserve(_nodes.size()); // _preorder must include root node

        _preorder.push_back(nd);
        
       while (true) {
           nd = findNextPreorder(nd);
           if (nd)
               _preorder.push_back(nd);
           else
               break;
       }   // end while loop
    }

    inline Node * Forest::findNextPreorder(Node * nd) {
        assert(nd);
        Node * next = 0;
        if (!nd->_left_child && !nd->_right_sib) {
            // nd has no children and no siblings, so next preorder is the right sibling of
            // the first ancestral node that has a right sibling.
            Node * anc = nd->_parent;
            while (anc && !anc->_right_sib)
                anc = anc->_parent;
            if (anc) {
                // We found an ancestor with a right sibling
                next = anc->_right_sib;
            }
            else {
                // nd is last preorder node in the tree
                next = 0;
            }
        }
        else if (nd->_right_sib && !nd->_left_child) {
            // nd has no children (it is a tip), but does have a sibling on its right
            next = nd->_right_sib;
        }
        else if (nd->_left_child && !nd->_right_sib) {
            // nd has children (it is an internal node) but no siblings on its right
            next = nd->_left_child;
        }
        else {
            // nd has both children and siblings on its right
            next = nd->_left_child;
        }
        return next;
    }

    inline pair<unsigned, unsigned> Forest::chooseTaxaToJoin(double s, Lot::SharedPtr lot){
        assert (s>1);
        double nsubtrees = s;
        unsigned t1=0;
        unsigned t2=1;
        //don't use this when there's only one choice (2 subtrees)
        
        if (nsubtrees > 2) {
            t1 = lot->randint(0, nsubtrees-1);
            t2 = lot->randint(0, nsubtrees-1);

            //keep calling t2 until it doesn't equal t1
            while (t2 == t1) {
                t2 = lot->randint(0, nsubtrees-1);
            }
        }
        assert(t1 < nsubtrees);
        assert (t2 < nsubtrees);

        return make_pair(t1, t2);
    }

    inline void Forest::calcPartialArray(Node* new_nd) {
        _partial_count++;
        auto &data_matrix=_data->getDataMatrix();

        for (unsigned i=0; i<G::_nloci; i++) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(i);
            
            if (new_nd->_set_partials) {
                if (!new_nd->_left_child) {
//                    assert (G::_save_memory || G::_start_mode == "sim");
                    if (!new_nd->_left_child) {
                        for (unsigned p=_first_pattern; p<_npatterns + _first_pattern; p++) {
                            unsigned pp = p;
                            for (unsigned s=0; s<G::_nstates; s++) {
                                Data::state_t state = (Data::state_t)1 << s;
                                Data::state_t d = data_matrix[new_nd->_number][pp];
                                double result = state & d;
                                (*new_nd->_partials)[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
                            }
                        }
                    }
                }
            }
        }
        
        double npatterns_total = _data->getNumPatterns();
//        new_nd->_partials = ps.getPartial(G::_nstates*npatterns_total);
        
        auto & parent_partial_array = *(new_nd->_partials);
        for (unsigned i=0; i<G::_nloci; i++) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(i);
            
            Node* first_child = new_nd->_left_child;
            Node* second_child = new_nd->_left_child->_right_sib;
            
//            for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
            for (unsigned a = 0; a < 2; a++) {
                Node* child;
                if (a == 0) {
                    child = first_child;
                }
                else {
                    child = second_child;
                }
                int next_real_node = -1;
                
                bool done = false;
                // find the next real node to use in partial calculation
                // even if you change the child, the right sib must stay the right sib from the original node
                while (!done) {
                    if (child->_next_real_node != -1) {
                        next_real_node = child->_next_real_node;
                        child = &_nodes[child->_next_real_node];
                    }
                    if (child->_next_real_node == -1) {
                        done = true;
                    }
                }
                
                if (child->_set_partials) {
                    if (child->_partials == nullptr) {
                        child->_partials = ps.getPartial(npatterns_total*G::_nstates);
                        calcPartialArray(child);
                    }
                }
                
                if (child->_set_partials) {
//                if (child->_set_partials && child->_use_in_likelihood) {
                    assert (child->_partials != nullptr);
                    auto & child_partial_array = *(child->_partials);

                    for (unsigned p = 0; p < _npatterns + _first_pattern; p++) {
                        for (unsigned s = 0; s <G::_nstates; s++) {
                            double sum_over_child_states = 0.0;
                            for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                                double child_transition_prob = 0.0;
# if defined (FOSSILS)
                                child_transition_prob = calcTransitionProbabilityFossil(child, s, s_child, i);
#else
                                child_transition_prob = calcTransitionProbability(child, s, s_child, i);
#endif
                                double child_partial = child_partial_array[p*G::_nstates + s_child];
                                sum_over_child_states += child_transition_prob * child_partial;
                            }   // child state loop
//                            if (!skip) {
                            if (next_real_node == -1) {
                                if (child == new_nd->_left_child) {
                                    parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                                }
                                else {
                                    parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
                                }
                            }
                            else {
                                if (a == 0) {
                                    parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                                }
                                else {
                                    assert (a == 1);
                                    parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
                                }
                            }
                        }   // parent state loop
                    }   // pattern loop
                }
            }   // child loop
        }
    }

    inline double Forest::calcTransitionProbability(Node* child, double s, double s_child, unsigned locus) {
        double relative_rate = G::_double_relative_rates[locus];
        assert (relative_rate > 0.0);
        
        double child_transition_prob = 0.0;

        if (G::_model == "JC" ) {
            double expterm = exp(-4.0*(child->_edge_length * relative_rate)/3.0);
            double prsame = 0.25+0.75*expterm;
            double prdif = 0.25 - 0.25*expterm;

            child_transition_prob = (s == s_child ? prsame : prdif);
            assert (child_transition_prob > 0.0);
            return child_transition_prob;
        }

        if (G::_model == "HKY") {
            double pi_A = G::_base_frequencies[0];
            double pi_C = G::_base_frequencies[1];
            double pi_G = G::_base_frequencies[2];
            double pi_T = G::_base_frequencies[3];

            double pi_j = 0.0;
            double PI_J = 0.0;

            double phi = (pi_A+pi_G)*(pi_C+pi_T)+G::_kappa*(pi_A*pi_G+pi_C*pi_T);
            double beta_t = 0.5*(child->_edge_length * relative_rate)/phi;

            // transition prob depends only on ending state
            if (s_child == 0) {
                // purine
                pi_j = pi_A;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 1) {
                // pyrimidine
                pi_j = pi_C;
                PI_J = pi_C + pi_T;
            }
            else if (s_child == 2) {
                // purine
                pi_j = pi_G;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 3) {
                // pyrimidine
                pi_j = pi_T;
                PI_J = pi_C + pi_T;
            }

            while (true) {
                if (s == s_child) {
                    // no transition or transversion
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (PI_J-pi_j)/PI_J*exp(-beta_t*(PI_J*G::_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*first_term+second_term;
                    break;
                }

                else if ((s == 0 && s_child == 2) || (s == 2 && s_child == 0) || (s == 1 && s_child == 3) || (s == 3 && s_child==1)) {
                    // transition
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (1/PI_J)*exp(-beta_t*(PI_J*G::_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*(first_term-second_term);
                    break;
                }

                else {
                    // transversion
                    child_transition_prob = pi_j*(1-exp(-beta_t));
                    break;
                }
            }
        }
        assert (child_transition_prob > 0.0);
        return child_transition_prob;
    }

#if defined (FOSSILS)
    inline bool Forest::checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets) {
        bool valid = false;
                
        vector<string> names_of_nodes_in_sets;
        for (auto &t:taxset) {
            for (auto &n:t._species_included) {
                names_of_nodes_in_sets.push_back(n);
            }
        }
        
        for (auto &t:taxset) {
            // if there are at least two nodes in the set that existing in _lineages, it's valid
            unsigned real_count = 0;
                for (auto &nd:_lineages) {
                    if (!valid) { // if already valid, can stop looking
                        for (auto &name:t._species_included) {
                            if (nd->_name == name) {
                                real_count++;
                                if (real_count == 2) {
                                    valid = true;
                                }
                                break;
                            }
                            if (real_count == 2) {
                                valid = true;
                                break;
                            }
                        }
                    }
                }
            
            if (real_count == t._species_included.size()) {
                valid = true;
            }
            _valid_taxsets.push_back(valid);
            valid = false;
            }
        
        // check if there are any taxa that are not in the taxset
        unsigned nnodes = (unsigned) _lineages.size();
        
        for (auto &t:unused_taxsets) {
            for (auto &name:t._species_included) {
                if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), name) == names_of_nodes_in_sets.end()) {
                    names_of_nodes_in_sets.push_back(name); // don't include duplicates
                }
            }
        }
        
        for (auto &nd:_lineages) {
            if (std::find(names_of_nodes_in_sets.begin(), names_of_nodes_in_sets.end(), nd->_name) != names_of_nodes_in_sets.end()) {
                nnodes--;
            }
        }
        
        if (nnodes > 1) { // TODO: if you have overlapping taxsets, this may not be true
            valid = true;
        }
        
        _valid_taxsets.push_back(valid);
        
        bool at_least_one_valid = false;
        for (unsigned v=0; v<_valid_taxsets.size(); v++) {
            if (_valid_taxsets[v]) {
                at_least_one_valid = true;
                break;
            }
        }
        
        return at_least_one_valid;
    }
#endif

#if defined (FOSSILS)
    inline double Forest::calcTransitionProbabilityFossil(Node* child, double s, double s_child, unsigned locus) {
        // this function uses accumulated branch lengths for the likelihood calculations
        double relative_rate = G::_double_relative_rates[locus];
        assert (relative_rate > 0.0);
        
        double child_transition_prob = 0.0;

        if (G::_model == "JC" ) {
    //            double expterm = exp(-4.0*(child->_edge_length * relative_rate)/3.0);
            double expterm = exp(-4.0*(child->_accumulated_height * relative_rate)/3.0);
            double prsame = 0.25+0.75*expterm;
            double prdif = 0.25 - 0.25*expterm;

            child_transition_prob = (s == s_child ? prsame : prdif);
            assert (child_transition_prob > 0.0);
            return child_transition_prob;
        }

        if (G::_model == "HKY") {
            double pi_A = G::_base_frequencies[0];
            double pi_C = G::_base_frequencies[1];
            double pi_G = G::_base_frequencies[2];
            double pi_T = G::_base_frequencies[3];

            double pi_j = 0.0;
            double PI_J = 0.0;

            double phi = (pi_A+pi_G)*(pi_C+pi_T)+G::_kappa*(pi_A*pi_G+pi_C*pi_T);
    //            double beta_t = 0.5*(child->_edge_length * relative_rate)/phi;
            double beta_t = 0.5*(child->_accumulated_height * relative_rate)/phi;


            // transition prob depends only on ending state
            if (s_child == 0) {
                // purine
                pi_j = pi_A;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 1) {
                // pyrimidine
                pi_j = pi_C;
                PI_J = pi_C + pi_T;
            }
            else if (s_child == 2) {
                // purine
                pi_j = pi_G;
                PI_J = pi_A + pi_G;
            }
            else if (s_child == 3) {
                // pyrimidine
                pi_j = pi_T;
                PI_J = pi_C + pi_T;
            }

            while (true) {
                if (s == s_child) {
                    // no transition or transversion
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (PI_J-pi_j)/PI_J*exp(-beta_t*(PI_J*G::_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*first_term+second_term;
                    break;
                }

                else if ((s == 0 && s_child == 2) || (s == 2 && s_child == 0) || (s == 1 && s_child == 3) || (s == 3 && s_child==1)) {
                    // transition
                    double first_term = 1+(1-PI_J)/PI_J*exp(-beta_t);
                    double second_term = (1/PI_J)*exp(-beta_t*(PI_J*G::_kappa+(1-PI_J)));
                    child_transition_prob = pi_j*(first_term-second_term);
                    break;
                }

                else {
                    // transversion
                    child_transition_prob = pi_j*(1-exp(-beta_t));
                    break;
                }
            }
        }
        assert (child_transition_prob > 0.0);
        return child_transition_prob;
    }
#endif

    inline void Forest::showForest() {
        for (unsigned g=0; g<_gene_tree_log_likelihoods.size(); g++) {
            output(format("log likelihood for gene tree %d: %d \n") % g % _gene_tree_log_likelihoods[g], 1);

        }
        
        output(makeNewick(15, true), 1);
        output("\n", 1);
    }

    inline string Forest::makeNewick(unsigned precision, bool use_names) {
        if (_lineages.size() > 1) {
            return makePartialNewick(precision, use_names);
        }
        else {
            string newick = "";

            const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
            const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
            const format internal_node_format( str(format("):%%.%df") % precision) );

            stack<Node *> node_stack;

            Node * nd = _lineages[0];
            while (nd) {
                if (nd->_left_child) {
                    // internal node
                    newick += "(";
                    node_stack.push(nd);
                }
                else {
                    // leaf node
                    if (use_names) {
                        newick += str(format(tip_node_name_format)
                            % nd->_name
                            % nd->_edge_length);
                    } else {
                        newick += str(format(tip_node_number_format)
                            % (nd->_number + 1)
                            % nd->_edge_length);
                    }
                    
                    if (nd->_right_sib) {
                        // Going to right sibling
                        newick += ",";
                    }
                    else if (nd->_parent) {
                        // Go down until we find an ancestor with a right sibling
                        // or until we reach the root node
                        Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                        while (popped && !popped->_right_sib) {
                            node_stack.pop();
                            if (node_stack.empty()) {
                                // We're back at the root node where we started,
                                // so finish by adding the final right parenthesis
                                newick += ")";
                                
                                // Set popped to nullptr to halt algorithm
                                popped = nullptr;
                            }
                            else {
                                // Save internal edge length and continue going down
                                newick += str(format(internal_node_format)
                                    % popped->_edge_length);
                                popped = node_stack.top();
                            }
                        }
                        
                        // Assuming we're not at the root node, move to right sibling
                        if (popped && popped->_right_sib) {
                            node_stack.pop();
                            newick += str(format(internal_node_format)
                                % popped->_edge_length);
                            newick += ",";
                        }
                    }
                }   // leaf node
                
                nd = findNextPreorder(nd);
                
            }   // while (nd)...

            return newick;
        }
    }

    inline string Forest::makePartialNewick(unsigned precision, bool use_names) {
            // this function makes a newick string for a partially constructed tree
            string newick = "(";
            const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
            const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
            const format internal_node_format( str(format("):%%.%df") % precision) );
            stack<Node *> node_stack;

            unsigned i = 0;
            unsigned a = 0;
            for (auto &lineage : _lineages) {
                Node * nd = lineage;
                while (nd) {
                    if (nd->_left_child) {
                        a++;
                        // internal node
                        newick += "(";
                        node_stack.push(nd);
                    }
                    else {
                        a++;
                        // leaf node
                        if (use_names) {
                            newick += str(format(tip_node_name_format)
                                % nd->_name
                                % nd->_edge_length);
                        } else {
                            newick += str(format(tip_node_number_format)
                                % (nd->_number + 1)
                                % nd->_edge_length);
                        }
                        if (nd->_right_sib)
                            newick += ",";
                        else {
                            Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                            while (popped && !popped->_right_sib) {
                                node_stack.pop();
                                if (node_stack.empty()) {
                                    //newick += ")";
                                    newick += str(format(internal_node_format) % lineage->_edge_length);
                                    popped = 0;
                                }
                                else {
                                    newick += str(format(internal_node_format) % popped->_edge_length);
                                    popped = node_stack.top();
                                }
                            }
                            if (popped && popped->_right_sib) {
                                node_stack.pop();
                                newick += str(format(internal_node_format) % popped->_edge_length);
                                newick += ",";
                            }
                        }
                    }   // leaf node
                    nd = findNextPreorder(nd);
                }   // while (subnd)...

                if (i < _lineages.size() - 1)
                    newick += ",";
                ++i;
            }
            newick += ")";

            return newick;
        }

    inline void Forest::clearPartials() {
        for (auto &nd:_lineages) {
            nd->_partials = nullptr;
        }
    }

    inline void Forest::operator=(const Forest & other) {
        _data                      = other._data;
        _first_pattern             = other._first_pattern;
        _npatterns                 = other._npatterns;
        _gene_tree_log_likelihoods = other._gene_tree_log_likelihoods;
        _ninternals                = other._ninternals;
        _nleaves                   = other._nleaves;
        _log_joining_prob          = other._log_joining_prob;
        _increments                = other._increments;
        _node_choices              = other._node_choices;
        _estimated_lambda          = other._estimated_lambda;
        _estimated_mu              = other._estimated_mu;
        _estimated_root_age        = other._estimated_root_age;
        _estimated_birth_difference = other._estimated_birth_difference;
        _turnover = other._turnover;
        _partial_count = other._partial_count;
#if defined (FOSSILS)
        _tree_height = other._tree_height;
        _valid_taxsets = other._valid_taxsets;
#endif

        // Copy _nodes
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        unsigned i = 0;
        for (const Node & othernd : other._nodes) {
            Node & nd = _nodes[i++];

            // copy parent
            if (othernd._parent) {
                unsigned parent_number = othernd._parent->_number;
                Node * parent = &*next(_nodes.begin(), parent_number);
                nd._parent = parent;
            }

            // copy left child
            if (othernd._left_child) {
                unsigned left_child_number = othernd._left_child->_number;
                Node* left_child = &*next(_nodes.begin(), left_child_number);
                nd._left_child = left_child;
            }
            else {
                nd._left_child = 0;
            }

            // copy right sibling
            if (othernd._right_sib) {
                unsigned right_sib_number = othernd._right_sib->_number;
                Node* right_sib = &*next(_nodes.begin(), right_sib_number);
                nd._right_sib = right_sib;
            }
            else
                nd._right_sib = 0;

            nd._number               = othernd._number;
            nd._name                 = othernd._name;
            nd._edge_length          = othernd._edge_length;
            nd._accumulated_height   = othernd._accumulated_height;
            nd._position_in_lineages = othernd._position_in_lineages;
            nd._partials             = othernd._partials;
            nd._split                = othernd._split;
            nd._height               = othernd._height;
            nd._set_partials         = othernd._set_partials;
            nd._use_in_likelihood    = othernd._use_in_likelihood;
            nd._next_real_node       = othernd._next_real_node;
            nd._renumber             = othernd._renumber;
            
            // Sanity check
            assert(nd._number >= 0);
            assert(nd._number < _nodes.size());
        }

        // Copy _lineages
        _lineages.resize(other._lineages.size());
        unsigned j = 0;
        for (auto othernd : other._lineages) {
            unsigned k = othernd->_number;
            Node * nd = &*next(_nodes.begin(), k);
            _lineages[j] = nd;
            j++;
        }
        
        // Copy _preorder
        _preorder.resize(other._preorder.size());
        unsigned m = 0;
        for (auto othernd : other._preorder) {
            unsigned n = othernd->_number;
            Node * nd = &*next(_nodes.begin(), n);
            _preorder[m] = nd;
            m++;
        }
    }

    inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) {
        assert(pi.size() == 4);
        assert(fabs(accumulate(pi.begin(), pi.end(), 0.0) - 1.0) < G::_small_enough);
        //assert(_relrate > 0.0);
        double _relrate = 1.0;
        double transition_prob = 0.0;
        
        // F81 transition probabilities
        double Pi[] = {pi[0] + pi[2], pi[1] + pi[3], pi[0] + pi[2], pi[1] + pi[3]};
        bool is_transition = (from == 0 && to == 2) || (from == 1 && to == 3) || (from == 2 && to == 0) || (from == 3 && to == 1);
        bool is_same = (from == 0 && to == 0) || (from == 1 && to == 1) | (from == 2 && to == 2) | (from == 3 && to == 3);
        bool is_transversion = !(is_same || is_transition);

        // HKY expected number of substitutions per site
        //  v = betat*(AC + AT + CA + CG + GC + GT + TA + TG) + kappa*betat*(AG + CT + GA + TC)
        //    = 2*betat*(AC + AT + CG + GT + kappa(AG + CT))
        //    = 2*betat*((A + G)*(C + T) + kappa(AG + CT))
        //  betat = v/[2*( (A + G)(C + T) + kappa*(AG + CT) )]
        double kappa = 1.0;
        double betat = 0.5*_relrate*edge_length/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
        if (is_transition) {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j - exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j);
        }
        else if (is_transversion) {
            double pi_j = pi[to];
            transition_prob = pi_j*(1.0 - exp(-betat));
        }
        else {
            double pi_j = pi[to];
            double Pi_j = Pi[to];
            transition_prob = pi_j*(1.0 + (1.0 - Pi_j)*exp(-betat)/Pi_j) + (Pi_j - pi_j)*exp(-betat*(kappa*Pi_j + 1.0 - Pi_j))/Pi_j;
        }
        return transition_prob;
    }
        
    inline void Forest::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        
        // Create vector of states for each node in the tree
        unsigned nnodes = (unsigned)_nodes.size();
#if defined (FOSSILS)
        nnodes = 0;
        for (auto &nd:_nodes) {
            if (nd._set_partials) {
                nnodes++;
            }
        }
#endif
        vector< vector<unsigned> > sequences(nnodes);
        for (unsigned i = 0; i < nnodes; i++) {
            sequences[i].resize(nsites, 4);
        }
        
        // Walk through tree in preorder sequence, simulating all sites as we go
        //    DNA   state      state
        //         (binary)  (decimal)
        //    A      0001        1
        //    C      0010        2
        //    G      0100        4
        //    T      1000        8
        //    ?      1111       15
        //    R      0101        5
        //    Y      1010       10
        
        // Draw equilibrium base frequencies from Dirichlet
        // having parameter G::_comphet
        vector<double> basefreq = {0.25, 0.25, 0.25, 0.25};
        if (G::_comphet != G::_infinity) {
            // Draw 4 Gamma(G::_comphet, 1) variates
            double A = lot->gamma(G::_comphet, 1.0);
            double C = lot->gamma(G::_comphet, 1.0);
            double G = lot->gamma(G::_comphet, 1.0);
            double T = lot->gamma(G::_comphet, 1.0);
            double total = A + C + G + T;
            basefreq[0] = A/total;
            basefreq[1] = C/total;
            basefreq[2] = G/total;
            basefreq[3] = T/total;
        }
        
#if defined (FOSSILS)
        // skip some nodes by renumbering everything
        unsigned new_node_number = 0;
        for (auto &nd:_nodes) {
            if (nd._set_partials) {
                nd._renumber = new_node_number;
                new_node_number++;
            }
            else {
                nd._renumber = -1;
            }
        }
        // Simulate starting sequence at the root node
        Node * nd;
        for (auto &node:_nodes) {
            if (node._renumber == nnodes-1) {
                nd = &node;
                break;
            }
        }
        int ndnum = nd->_renumber;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = G::multinomialDraw(lot, basefreq);
        }
#else
        // Simulate starting sequence at the root node
        Node * nd = *(_lineages.begin());
        int ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = G::multinomialDraw(lot, basefreq);
        }
#endif
        
# if defined (FOSSILS)
        nd = findNextPreorder(nd);
                
        while (nd) {
            if (nd->_set_partials) {
                ndnum = nd->_number;
                
                // Get reference to parent sequence
                assert(nd->_parent);
                
                bool parent_found = false;
                bool skip_node = false;
                unsigned parnum = nd->_parent->_number;
                Node parent = _nodes[parnum];
                
                while (!parent_found) {
                    if (!parent._set_partials) { // if parent is not a real node, use its parent
                        if (parent._parent) {
                            parnum = parent._parent->_number;
                            parent = _nodes[parnum];
                        }
                        else {
                            skip_node = true;
                            // skip node if the parent is not a real node and neither is its parent
                            parent_found = true;
                        }
                    }
                    else {
                        parent_found = true;
                    }
                }
                
                if (!skip_node) {
                    assert (_nodes[parnum]._set_partials);
                        
                        // Evolve nd's sequence given parent's sequence and edge length
                    for (unsigned i = 0; i < nsites; i++) {
                        // Choose relative rate for this site
                        double site_relrate = 1.0;
                        if (G::_asrv_shape != G::_infinity)
                            site_relrate = lot->gamma(G::_asrv_shape, 1.0/G::_asrv_shape);
                        unsigned from_state = sequences[_nodes[parnum]._renumber][i];
                        double cum_prob = 0.0;
                        double u = lot->uniform();
                        for (unsigned to_state = 0; to_state < 4; to_state++) {
                            cum_prob += calcSimTransitionProbability(from_state, to_state, basefreq, site_relrate*nd->_accumulated_height);
                            if (u < cum_prob) {
                                sequences[_nodes[ndnum]._renumber][i] = to_state;
                                break;
                            }
                        }
                        assert(sequences[_nodes[ndnum]._renumber][i] < 4);
                    }
                }
            }
            // else, node is not a real node and should not be included in the data
            // Move to next node in preorder sequence
            nd = findNextPreorder(nd);
        }

#else
        nd = findNextPreorder(nd);
        while (nd) {
            ndnum = nd->_number;
            assert(ndnum < nnodes);

            // Get reference to parent sequence
            assert(nd->_parent);
            unsigned parnum = nd->_parent->_number;
            assert(parnum < nnodes);
            
            // Evolve nd's sequence given parent's sequence and edge length
            for (unsigned i = 0; i < nsites; i++) {
                // Choose relative rate for this site
                double site_relrate = 1.0;
                if (G::_asrv_shape != G::_infinity)
                    site_relrate = lot->gamma(G::_asrv_shape, 1.0/G::_asrv_shape);
                unsigned from_state = sequences[parnum][i];
                double cum_prob = 0.0;
                double u = lot->uniform();
                for (unsigned to_state = 0; to_state < 4; to_state++) {
                    cum_prob += calcSimTransitionProbability(from_state, to_state, basefreq, site_relrate*nd->_edge_length);
                    if (u < cum_prob) {
                        sequences[ndnum][i] = to_state;
                        break;
                    }
                }
                assert(sequences[ndnum][i] < 4);
            }
            
            // Move to next node in preorder sequence
            nd = findNextPreorder(nd);
        }
#endif

        assert(data);
        Data::data_matrix_t & dm = data->getDataMatrixNonConst();
        
        // Copy sequences to data object
        for (unsigned t = 0; t < G::_ntaxa; t++) {
            // Allocate row t of _data's _data_matrix data member
            dm[t].resize(starting_site + nsites);
            
            // Get reference to nd's sequence
            unsigned ndnum = _nodes[t]._number;
            
            // Translate to state codes and copy
            for (unsigned i = 0; i < nsites; i++) {
                dm[t][starting_site + i] = (Data::state_t)1 << sequences[ndnum][i];
            }
        }
    }

    inline unsigned Forest::getNumLineages() const {
        return (unsigned)_lineages.size();
    }
    
    inline unsigned Forest::getNumNodes() const {
        return (unsigned)_preorder.size();
    }

    inline Node * Forest::pullNode() {
        // Add one node to the end of _nodes vector
        _nodes.resize(_nodes.size() + 1);
        
        // Get pointer to the new node
        Node * new_nd = &(*_nodes.rbegin());
        
        // Set up the new node
        new_nd->clear();
        new_nd->_number = _nleaves + _ninternals;
        new_nd->_split.resize(G::_ntaxa);
        new_nd->_name=boost::str(boost::format("node-%d")%new_nd->_number);

        // Increment number of internals
        _ninternals++;

        return new_nd;
    }
    
    inline void Forest::joinRandomLineagePair(Lot::SharedPtr lot) {
        unsigned n = getNumLineages();
        auto lineage_pair = lot->nchoose2(n);
        unsigned i = lineage_pair.first;
        unsigned j = lineage_pair.second;
        Node * ancnd = pullNode();
        Node * lchild = _lineages[i];
        Node * rchild = _lineages[j];
        
        ancnd->_left_child = lchild;
        ancnd->_right_sib = nullptr;
        ancnd->_parent = nullptr;
        
        lchild->_right_sib = rchild;
        lchild->_parent = ancnd;
        
        rchild->_right_sib = nullptr;
        rchild->_parent = ancnd;
        
        updateNodeVector(_lineages, lchild, rchild, ancnd);
        //removeTwoAddOne(_lineages, lchild, rchild, ancnd);
        
#if defined (FOSSILS)
        // if subtree1 or subtree2 is a fossil, set partials to false for the new node
        if (boost::ends_with(lchild->_name, "FOSSIL") || boost::ends_with(rchild->_name, "FOSSIL")) {
            ancnd->_set_partials = false;
            ancnd->_use_in_likelihood = false;
        }
        // check if both of child's children are fossils
        if (lchild->_left_child) {
            if (boost::ends_with(lchild->_left_child->_name, "FOSSIL") && boost::ends_with(lchild->_left_child->_right_sib->_name, "FOSSIL")) {
                ancnd->_set_partials = false;
                ancnd->_use_in_likelihood = false;
            }
        }
        if (rchild->_left_child) {
            if (boost::ends_with(rchild->_left_child->_name, "FOSSIL") && boost::ends_with(rchild->_left_child->_right_sib->_name, "FOSSIL")) {
                ancnd->_set_partials = false;
                ancnd->_use_in_likelihood = false;
            }
        }
        
        unsigned end_node = (unsigned) _lineages.size() - 1;
        Node *node_to_check = _lineages[end_node];
        
        // check for the new node that was added in the previous step
        bool done = false;
        while (!done) {
            if (boost::ends_with(node_to_check->_name, "FOSSIL")) {
                end_node --;
                node_to_check = _lineages[end_node];
            }
            else {
                done = true;
            }
        }
        
        if (!node_to_check->_set_partials &&node_to_check->_left_child) {
            if (!boost::ends_with(node_to_check->_left_child->_name, "FOSSIL")) {
                // add extra branch length to the next node
                int next_node = node_to_check->_left_child->_next_real_node;
                if (next_node > -1) {
                    _nodes[next_node]._accumulated_height += node_to_check->_edge_length;
                    node_to_check->_next_real_node = _nodes[next_node]._number;
                }
                else {
                    node_to_check->_left_child->_accumulated_height += node_to_check->_edge_length;
                    node_to_check->_next_real_node = node_to_check->_left_child->_number;
                }
            }
            else if (!boost::ends_with(node_to_check->_left_child->_right_sib->_name, "FOSSIL")) {
                // add extra branch length to the next node
                // can't accumulate to a fake node - if the node is fake, find its next real node
                int next_node = node_to_check->_left_child->_right_sib->_next_real_node;
                if (next_node > -1) {
                    _nodes[next_node]._accumulated_height += node_to_check->_edge_length;
                    node_to_check->_next_real_node = _nodes[next_node]._number;
                }
                else {
                    node_to_check->_left_child->_right_sib->_accumulated_height += node_to_check->_edge_length;
                    node_to_check->_next_real_node = node_to_check->_left_child->_right_sib->_number;
                }
            }
        }
#endif
    }
    
    inline void Forest::renumberInternals() {
        // First internal node number is the number of leaves
        int next_node_number = (unsigned)G::_taxon_names.size();

        // Renumber internal nodes in postorder sequence for each lineage in turn
        for (auto nd : boost::adaptors::reverse(_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                assert(nd->_height != G::_infinity);
                nd->_number = next_node_number++;
                assert(nd->_left_child->_right_sib);
                assert(nd->_left_child->_right_sib->_right_sib == nullptr);
            }
            else {
                // nd is a leaf node
                assert(nd->_number > -1);
                nd->_height = 0.0;
            }
                            
            if (nd->_parent) {
                // Set parent's height if nd is right-most child of its parent
                bool is_rightmost_child = !nd->_right_sib;
                double parent_height = nd->_height + nd->_edge_length;
                if (is_rightmost_child) {
                    nd->_parent->_height = parent_height;
                }
                
                // If nd is not its parent's rightmost child, check ultrametric assumption
                assert(!is_rightmost_child || fabs(nd->_parent->_height - parent_height) < G::_small_enough);
            }
        }
    }
    

    inline void Forest::buildBirthDeathTree() {
        int fossil_number = 0;
    # if defined (FOSSILS)
    //    fossil_number = (unsigned) G::_sim_fossils.size() - 1;
        // sort fossils from youngest to oldest for use in later proposal
        sort(G::_sim_fossils.begin(), G::_sim_fossils.end(), [](Fossil & left, Fossil & right) {
            return left._age < right._age;
        });
    #endif
        
        // Algorithm from Yang and Rannala. 1997. MBE 14(7):717-724.
        createTrivialForest();
    #if 1
        unsigned nsteps = G::_ntaxa - 1;
    # if defined (FOSSILS)
        nsteps += G::_sim_fossils.size();
    #endif
        double cum_height = 0.0;
        for (unsigned i = 0; i < nsteps; i++) {
            // Determine number of lineages remaining
            unsigned n = getNumLineages();
            assert(n > 1);
            
            // Draw n-1 internal node heights and store in vector heights
            vector<double> heights(n - 1, 0.0);
            
            double rho = G::_sim_rho;
            double birth_rate = G::_sim_lambda;
            double death_rate = G::_sim_mu;
            double exp_death_minus_birth = exp(death_rate - birth_rate);
            double phi = 0.0;
            phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
            phi += (death_rate - birth_rate)*exp_death_minus_birth;
            phi /= (exp_death_minus_birth - 1.0);
            for (unsigned i = 0; i < n - 2; i++) {
                double u = rng->uniform();
                double y = u/(1.0 + birth_rate*rho*(1.0 - u));
                if (birth_rate > death_rate) {
                    y = log(phi - u*rho*birth_rate);
                    y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                    y /= (death_rate - birth_rate);
                }
                heights[i] = y;
            }
            heights[n-2] = 1.0;
            sort(heights.begin(), heights.end());
            
            // Waiting time to next speciation event is first height
            // scaled so that max height is 1 - cum_height
            double t = heights[0]*(1.0 - cum_height);
            
    # if defined (FOSSILS)
            if (fossil_number < G::_sim_fossils.size()) {
                if (cum_height + t > G::_sim_fossils[fossil_number]._age) {
                    // add fossil
                    t = G::_sim_fossils[fossil_number]._age - cum_height;
                    
                    for (auto &nd:_lineages) {
                        nd->_edge_length += t;
                        nd->_accumulated_height += t;
                    }
                    
                    cum_height += t;
                    
                    //new node is always needed
                    Node* new_nd = pullNode();

                    new_nd->_name = G::_sim_fossils[fossil_number]._name + "_FOSSIL";
                    new_nd->_set_partials = false; // do not include this node in likelihood calculation
                    new_nd->_position_in_lineages = (unsigned) _lineages.size();
                    new_nd->_use_in_likelihood = false;
                    _lineages.push_back(new_nd);
                    
                    fossil_number++;
                }
                else {
                    for (auto &nd:_lineages) {
                        nd->_edge_length += t;
                        nd->_accumulated_height +=t;
                    }
                    cum_height += t;
                }
            }
            else {
                for (auto &nd:_lineages) {
                    nd->_edge_length += t;
                    nd->_accumulated_height +=t;
                }
                cum_height += t;
            }
    #endif
            
    #if !defined (FOSSILS)
            cum_height += t;
            advanceAllLineagesBy(t);
    #endif
            
            joinRandomLineagePair(rng);
        }
    #else
        // Draw n-1 internal node heights and store in vector heights
        unsigned n = getNumLineages();
        vector<double> heights(n - 1, 0.0);
        
        double rho = G::_sim_rho;
        double birth_rate = G::_sim_lambda;
        double death_rate = G::_sim_mu;
        double exp_death_minus_birth = exp(death_rate - birth_rate);
        double phi = 0.0;
        phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
        phi += (death_rate - birth_rate)*exp_death_minus_birth;
        phi /= (exp_death_minus_birth - 1.0);
        for (unsigned i = 0; i < n - 2; i++) {
            double u = rng->uniform();
            double y = u/(1.0 + birth_rate*rho*(1.0 - u));
            if (birth_rate > death_rate) {
                y = log(phi - u*rho*birth_rate);
                y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                y /= (death_rate - birth_rate);
            }
            heights[i] = y;
        }
        heights[n-2] = 1.0;
        sort(heights.begin(), heights.end());
        
        // Now that we have the increments, perform the joins
        //
        // 1   2   3  4  5 n = 5, n - 1 = 4
        // |   |   |  |  |
        // +-+-+   |  |  | i = 0
        //   |     |  |  |
        //   +--+--+  |  | i = 1
        //      |     |  |
        //      +--+--+  | i = 2
        //         |     |
        //         +--+--+ i = 3
        
        double t0 = 0.0;
        for (unsigned i = 0; i < n - 1; i++) {
            double t = heights[i];
            double dt = t - t0;
            advanceAllLineagesBy(dt);
            joinRandomLineagePair(rng);
            t0 = t;
        }
    #endif
        assert(getNumLineages() == 1);

        // Scale all edge lengths by G::_sim_root_age
        scaleAllEdgeLengthsBy(G::_sim_root_age);
        
        refreshPreorder();
    }

#if defined (INCREMENT_COMPARISON_TEST)
    inline void Forest::buildBirthDeathTreeTest() {
        // Algorithm from Yang and Rannala. 1997. MBE 14(7):717-724.
        createTrivialForest();
#if 0
        ofstream logf("sim8.log");
        for (unsigned a = 0; a < 100000; a++) {
        
        unsigned nsteps = G::_ntaxa - 1;
        double cum_height = 0.0;
//        for (unsigned i = 0; i < nsteps; i++) {
            for (unsigned i = 0; i < 2; i++) {
            // Determine number of lineages remaining
            unsigned n = getNumLineages();
                if (i == 1) {
                    n--;
                }
            assert(n > 1);
            
            // Draw n-1 internal node heights and store in vector heights
            vector<double> heights(n - 1, 0.0);
            
            double rho = G::_sim_rho;
            double birth_rate = G::_sim_lambda;
            double death_rate = G::_sim_mu;
            double exp_death_minus_birth = exp(death_rate - birth_rate);
            double phi = 0.0;
            phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
            phi += (death_rate - birth_rate)*exp_death_minus_birth;
            phi /= (exp_death_minus_birth - 1.0);
            for (unsigned i = 0; i < n - 2; i++) {
                double u = rng->uniform();
                double y = u/(1.0 + birth_rate*rho*(1.0 - u));
                if (birth_rate > death_rate) {
                    y = log(phi - u*rho*birth_rate);
                    y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                    y /= (death_rate - birth_rate);
                }
                heights[i] = y;
            }
            heights[n-2] = 1.0;
            sort(heights.begin(), heights.end());
            
            // Waiting time to next speciation event is first height
            // scaled so that max height is 1 - cum_height
            double t = heights[0]*(1.0 - cum_height);
            cum_height += t;
            
            if (i == 1) {
            logf << a;
            logf << "\t" << heights[7] - heights[6] << endl;
            }
        }
        }
#else
        ofstream logf("sim5.log");
        logf << "sample" << "\t" << "increment" << endl;
        // Draw n-1 internal node heights and store in vector heights
        for (unsigned a =0; a < 10; a++) {
        unsigned n = getNumLineages();
        vector<double> heights(n - 1, 0.0);
        
        double rho = G::_sim_rho;
        double birth_rate = G::_sim_lambda;
        double death_rate = G::_sim_mu;
        double exp_death_minus_birth = exp(death_rate - birth_rate);
        double phi = 0.0;
        phi += rho*birth_rate*(exp_death_minus_birth - 1.0);
        phi += (death_rate - birth_rate)*exp_death_minus_birth;
        phi /= (exp_death_minus_birth - 1.0);
        for (unsigned i = 0; i < n - 2; i++) {
            double u = rng->uniform();
            double y = u/(1.0 + birth_rate*rho*(1.0 - u));
            if (birth_rate > death_rate) {
                y = log(phi - u*rho*birth_rate);
                y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                y /= (death_rate - birth_rate);
            }
            heights[i] = y;
        }
        heights[n-2] = 1.0;
        sort(heights.begin(), heights.end());
            
            logf << a << "\t" << heights[4] - heights[3] << endl;
        
        // Now that we have the increments, perform the joins
        //
        // 1   2   3  4  5 n = 5, n - 1 = 4
        // |   |   |  |  |
        // +-+-+   |  |  | i = 0
        //   |     |  |  |
        //   +--+--+  |  | i = 1
        //      |     |  |
        //      +--+--+  | i = 2
        //         |     |
        //         +--+--+ i = 3
            
        }
#endif
    }
#endif

    inline void Forest::buildYuleTree() {
        createTrivialForest();
        unsigned nsteps = G::_ntaxa - 1;
        for (unsigned i = 0; i < nsteps; i++) {
            // Determine number of lineages remaining
            unsigned n = getNumLineages();
            assert(n > 1);
    
            // Waiting time to speciation event is Exponential(rate = n*lambda)
            // u = 1 - exp(-r*t) ==> t = -log(1-u)/r
            double r = G::_sim_lambda*n;
            double u = rng->uniform();
            double t = -log(1.0 - u)/r;
            advanceAllLineagesBy(t);
            joinRandomLineagePair(rng);
        }
        assert(getNumLineages() == 1);
        refreshPreorder();
    }

    inline void Forest::scaleAllEdgeLengthsBy(double scaling_factor) {
        // This function should only be called for complete trees
        assert (getNumLineages() == 1);
        
        // Supplied scaling_factor should be strictly positive
        assert(scaling_factor > 0.0);
        
        for (auto & nd : _nodes) {
            double elen = nd.getEdgeLength();
            nd.setEdgeLength(scaling_factor*elen);
        }
        
#if defined (FOSSILS)
        for (auto & nd : _nodes) {
            nd._accumulated_height *= scaling_factor;
        }
#endif
    }
    
    inline void Forest::advanceAllLineagesBy(double dt) {
        // Add t to the edge length of all lineage root nodes, unless there
        // is just one lineage, in which case do nothing
        unsigned n = getNumLineages();
        if (n > 1) {
            for (auto nd : _lineages) {
                double elen = nd->getEdgeLength() + dt;
                assert(elen >= 0.0 || fabs(elen) < Node::_smallest_edge_length);
                nd->setEdgeLength(elen);
                ++n;
            }
        }
    }

    struct negLogLikeDist {
        negLogLikeDist(unsigned npatterns, unsigned first, const Data::pattern_counts_t & counts, const vector<double> & same, const vector<double> & diff, double v0)
            : _npatterns(npatterns), _first(first), _counts(counts), _same(same), _diff(diff), _v0(v0) {}
        
        double operator()(double const & v) {
            double edgelen = v + _v0;
            double tprob_same = 0.25 + 0.75*exp(-4.0*edgelen/3.0);
            double tprob_diff = 0.25 - 0.25*exp(-4.0*edgelen/3.0);

            double log_like = 0.0;
            for (unsigned p = 0; p < _npatterns; p++) {
                double site_like = 0.25 * (tprob_same * _same[p] + tprob_diff * _diff[p]);
                log_like += log(site_like) * _counts[_first + p];
            }
            
            return -log_like;
        }
        
        private:
            unsigned _npatterns;
            unsigned _first;
            const Data::pattern_counts_t & _counts;
            const vector<double> & _same;
            const vector<double> & _diff;
            double _v0;
    };

    inline double Forest::getLineageHeight(Node* nd) {
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

    inline void Forest::drawBirthDiff(Lot::SharedPtr lot) {
        // birth diff = lambda - mu
        // Gamma(1, n) = Exp(1/n)
        // mean = n
        // for now, n = (G::_lambda - G::_mu) set by user
        double mean = G::_lambda - G::_mu;
        _estimated_birth_difference = lot->gamma(1, mean);
    }

    inline void Forest::drawTurnover(Lot::SharedPtr lot) {
        // turnover = mu / lambda
        // Uniform distribution - must be (0, 1)
        
        _turnover = lot->uniform();
    }

    inline void Forest::calculateLambdaAndMu() {
        _estimated_lambda = _estimated_birth_difference / (1 - _turnover);
        _estimated_mu = (_turnover * _estimated_birth_difference) / (1 - _turnover);
    }

    inline void Forest::drawLambda(Lot::SharedPtr lot) {
        // Yule model; don't estimate mu since it is fixed at 0.0
        // Gamma(1, n) = Exp(1/n)
        // mean = n
        // for now, n = G::_lambda set by user
        _estimated_lambda = lot->gamma(1, G::_lambda);
    }

    inline void Forest::drawRootAge(Lot::SharedPtr lot) {
        // Gamma(1, n) = Exp(1/n)
        // mean = n
        // for now, n = G::_root_age set by user
        _estimated_root_age = lot->gamma(1, G::_root_age);
#if defined (FOSSILS)
        _estimated_root_age = G::_root_age;
#endif
    }

    inline double Forest::getLogLikelihood() {
        double log_likelihood = 0.0;
        for (auto &g:_gene_tree_log_likelihoods) {
            log_likelihood += g;
        }
        return log_likelihood;
    }
    
}
