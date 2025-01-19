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
        //void buildYuleTree();
        void buildBirthDeathTree();
        typedef shared_ptr<Forest> SharedPtr;
        string debugSaveForestInfo() const;
        //POL added above
    
    private:

        void                        clear();
        void                        setForestData(Data::SharedPtr d, bool partials);
        Node *                      findNextPreorder(Node * nd);
        void                        refreshPreorder();
        double                      calcSubsetLogLikelihood(unsigned i);
        void                        addIncrement(Lot::SharedPtr lot);
        double                      joinTaxa(Lot::SharedPtr lot);
        void                        calcPartialArray(Node* new_nd);
        double                      calcTransitionProbability(Node* child, double s, double s_child, unsigned locus);
        
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

        pair<unsigned, unsigned>    chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        string                      makePartialNewick(unsigned precision, bool use_names);
        void                        showForest();
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);

        double                      getTreeHeight();
        double                      getTreeLength();
        double                      getSpeciesTreePrior();
        double                      calcTopologyPrior(unsigned nlineages);
        void                        clearPartials();
    
        Data::SharedPtr             _data;
        vector<Node *>              _lineages;
#if NEWWAY == POLWAY
        vector<Node>                _nodes;
#else
        list<Node>                  _nodes;
#endif
        vector<Node*>               _preorder;
        unsigned                    _first_pattern;
        unsigned                    _npatterns;
        vector<double>              _gene_tree_log_likelihoods;
        unsigned                    _ninternals;
        unsigned                    _nleaves;
        double                      _log_joining_prob;
        vector<pair<double, double>> _increments_and_priors;
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
        _increments_and_priors.clear();
        _nleaves = 0;
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::createTrivialForest() {
        assert(G::_ntaxa > 0);
        assert(G::_ntaxa == G::_taxon_names.size());
        clear();
#if NEWWAY == POLWAY
        unsigned nnodes = 2*G::_ntaxa - 1;
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
#else  //AAM (_nodes is list not fully allocated at start)
        unsigned i = 0;
        unsigned nnodes = 2*G::_ntaxa - 1;
        _lineages.reserve(nnodes);
        for (unsigned i = 0; i < G::_ntaxa; i++) {
            Node * nd = pullNode();
            string taxon_name = G::_taxon_names[i];
            nd->_number = (int)i;
            nd->_name = taxon_name;
            nd->setEdgeLength(0.0);
            nd->_height = 0.0;
            nd->_split.resize(G::_ntaxa);
            nd->_split.setBitAt(i);
            _lineages.push_back(nd);
        }
#endif
        
        refreshPreorder();
    }
    
    inline void Forest::setForestData(Data::SharedPtr d, bool partials) {
        _data = d;
        
        _gene_tree_log_likelihoods.resize(G::_nloci, 0.0);
        
        auto &data_matrix=_data->getDataMatrix();

#if NEWWAY == POLWAY
        unsigned nnodes = 2*G::_ntaxa - 1;
        _nodes.reserve(nnodes);
        _nodes.resize(G::_ntaxa);
#else
        _nodes.resize(G::_ntaxa);
#endif

        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
#if NEWWAY == POLWAY
            _nodes[i]._right_sib=0;
            _nodes[i]._name=" ";
            _nodes[i]._left_child=0;
            _nodes[i]._right_sib=0;
            _nodes[i]._parent=0;
            _nodes[i]._number=i;
            _nodes[i]._edge_length=0.0;
            _nodes[i]._position_in_lineages=i;
            //_nodes[i]._partials.resize(G::_nloci);
            _nodes[i]._partials = nullptr;
            _nodes[i]._name = G::_taxon_names[i];
            _lineages.push_back(&_nodes[i]);
#else
            Node* nd = &*next(_nodes.begin(), i);
            nd->_right_sib=0;
            nd->_name=" ";
            nd->_left_child=0;
            nd->_right_sib=0;
            nd->_parent=0;
            nd->_number=i;
            nd->_edge_length=0.0;
            nd->_position_in_lineages=i;
            _lineages.push_back(nd);
            //nd->_partials.resize(G::_nloci); // _partials contains a vector of partials for each locus
            nd->_partials = nullptr; // _partials contains a vector of partials for each locus
            // replace all spaces with underscores so that other programs do not have
              // trouble parsing tree descriptions
            std::string name = G::_taxon_names[i];
            boost::replace_all(name, " ", "_");
            nd->_name = name;
#endif
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
        _npatterns = _data->getNumPatternsInSubset(i); // TODO: make a vector of patterns by gene?
        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
        _first_pattern = gene_begin_end.first;

        for (auto &nd:_lineages) {

            //temporary!
            if (nd->_partials == nullptr) {
                cerr << "oops" << endl;
            }
            
            assert (nd->_partials != nullptr);
            double log_like = 0.0;
            for (unsigned p=_first_pattern; p<_npatterns + _first_pattern; p++) {
                double site_like = 0.0;
                for (unsigned s=0; s<G::_nstates; s++) {
                    double partial = (*nd->_partials)[p*G::_nstates+s];
                    site_like += 0.25*partial;
                }
                assert(site_like>0);
                log_like += log(site_like)*counts[p];
            }
            _gene_tree_log_likelihoods[i] += log_like;
        }
        
        return _gene_tree_log_likelihoods[i];
    }

    inline void Forest::addIncrement(Lot::SharedPtr lot) {
            unsigned nlineages = (unsigned) _lineages.size();
            double rate = nlineages * G::_lambda;
            
            double increment = -log(1.0 - lot->uniform())/rate;
            
            for (auto &nd:_lineages) {
                nd->_edge_length += increment;
            }
            
            // lorad only works if all topologies the same - then don't include the prior on joins because it is fixed
            double increment_prior = (log(rate)-increment*rate);
                        
            _increments_and_priors.push_back(make_pair(increment, increment_prior));
    }

    inline double Forest::joinTaxa(Lot::SharedPtr lot) {
        double prev_log_likelihood = 0.0;
        
        for (auto &g:_gene_tree_log_likelihoods) {
            prev_log_likelihood += g;
        }
        
        unsigned nlineages = (unsigned) _lineages.size();
        Node *subtree1 = nullptr;
        Node *subtree2 = nullptr;
        
        if (G::_save_memory) {
            double npatterns_total = _data->getNumPatterns();
            for (auto &nd:_lineages) {
                if (nd->_partials == nullptr) {
                    nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                    calcPartialArray(nd);
                }
            }
        }
        
        pair <unsigned, unsigned> t = make_pair (0, 1); // if there is only choice, 0 and 1 will be chosen
        t = chooseTaxaToJoin(nlineages, lot);
        
        subtree1 = _lineages[t.first];
        subtree2 = _lineages[t.second];

        assert (subtree1 != subtree2);
        
        // The commented lines below now done by pullNode()
        //Node nd;
        //_nodes.push_back(nd);
        //Node* new_nd = &_nodes.back();
        //new_nd->_parent=0;
        //new_nd->_number=_nleaves+_ninternals;
        //new_nd->_edge_length=0.0;
        //_ninternals++;
        //new_nd->_right_sib=0;

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
        
        bool calc_likelihood = true;
                
        if (calc_likelihood) {
            calcPartialArray(new_nd);
        }

        for (unsigned index = 0; index < G::_nloci; index++) {
            subtree1->_partials=nullptr; // throw away subtree partials now, no longer needed
            subtree2->_partials=nullptr;
        }
        
        //update node lists
        updateNodeVector(_lineages, subtree1, subtree2, new_nd);

        if (calc_likelihood) {
            for (unsigned index = 0; index<G::_nloci; index++) {
                calcSubsetLogLikelihood(index);
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
        
        calcTopologyPrior((unsigned) _lineages.size()+1);

        return log_weight;
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
        for (int i=0; i < (int) _lineages.size(); i++) {
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
        double sum_height = 0.0;

        // calculate height of lineage
        Node* base_node = _lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
        }
        return sum_height;
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

    inline double Forest::getSpeciesTreePrior() {
        double prior = 0.0;
        for (auto &i:_increments_and_priors) {
            prior += i.second;
        }
        prior += _log_joining_prob;
        return prior;
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
        auto &data_matrix=_data->getDataMatrix();

        for (unsigned i=0; i<G::_nloci; i++) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(i);
            
            if (!new_nd->_left_child) {
                assert (G::_save_memory || G::_start_mode == "sim");
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

        
        double npatterns_total = _data->getNumPatterns();
//        new_nd->_partials = ps.getPartial(G::_nstates*npatterns_total);
        
        auto & parent_partial_array = *(new_nd->_partials);
        for (unsigned i=0; i<G::_nloci; i++) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(i);

            for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

                if (child->_partials == nullptr) {
                    child->_partials = ps.getPartial(npatterns_total*G::_nstates);
                    calcPartialArray(child);
                }
                
                assert (child->_partials != nullptr);
                auto & child_partial_array = *(child->_partials);

                for (unsigned p = 0; p < _npatterns + _first_pattern; p++) {
                    for (unsigned s = 0; s <G::_nstates; s++) {
                        double sum_over_child_states = 0.0;
                        for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                            double child_transition_prob = calcTransitionProbability(child, s, s_child, i);
                            double child_partial = child_partial_array[p*G::_nstates + s_child];
                            sum_over_child_states += child_transition_prob * child_partial;
                        }   // child state loop
                        if (child == new_nd->_left_child)
                            parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                        else
                            parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
                    }   // parent state loop
                }   // pattern loop
            }   // child loop
        }
    }

    inline double Forest::calcTransitionProbability(Node* child, double s, double s_child, unsigned locus) {
        double relative_rate = G::_double_relative_rates[locus];
        assert (relative_rate > 0.0);
        
        double child_transition_prob = 0.0;

        if (G::_model == "JC" ) {
            double expterm = exp(-4.0*(child->_edge_length * relative_rate)/3.0); // TODO: need to include relative rates
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
            double beta_t = 0.5*(child->_edge_length * relative_rate)/phi; // TODO: need to include relative rates


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

#if NEWWAY == POLWAY
    inline void Forest::operator=(const Forest & other) {
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
        _preorder.resize(other._preorder.size());

        _data                      = other._data;
        _first_pattern             = other._first_pattern;
        _npatterns                 = other._npatterns;
        _gene_tree_log_likelihoods = other._gene_tree_log_likelihoods;
        _ninternals                = other._ninternals;
        _nleaves                   = other._nleaves;
        _log_joining_prob          = other._log_joining_prob;
        _increments_and_priors     = other._increments_and_priors;

        // copy tree itself

        for (const Node & othernd : other._nodes) {
            int k = othernd._number;

            if (k > -1) {
                Node * nd = &*next(_nodes.begin(), k);

                // copy parent
                if (othernd._parent) {
                    unsigned parent_number = othernd._parent->_number;
                    Node * parent = &*next(_nodes.begin(), parent_number);
                    nd->_parent = parent;
                }

                // copy left child
                if (othernd._left_child) {
                    unsigned left_child_number = othernd._left_child->_number;
                    Node* left_child = &*next(_nodes.begin(), left_child_number);
                    nd->_left_child = left_child;
                }
                else {
                    nd->_left_child = 0;
                }

                // copy right sibling
                if (othernd._right_sib) {
                    unsigned right_sib_number = othernd._right_sib->_number;
                    Node* right_sib = &*next(_nodes.begin(), right_sib_number);
                    nd->_right_sib = right_sib;
                }
                else
                    nd->_right_sib = 0;

                nd->_number = othernd._number;
                nd->_name = othernd._name;
                nd->_edge_length = othernd._edge_length;
                nd->_position_in_lineages = othernd._position_in_lineages;
                nd->_partials = othernd._partials;
                nd->_split = othernd._split;
                nd->_height = othernd._height;
            } // if k > -1
        }

        unsigned j = 0;
        for (auto othernd : other._lineages) {
            unsigned k = othernd->_number;
            Node * nd = &*next(_nodes.begin(), k);
            _lineages[j] = nd;
            j++;
        }
        
        if (other._preorder.size() > 0) {
            unsigned m = 0;
            for (auto othernd : other._preorder) {
                unsigned n = othernd->_number;
                Node * nd = &*next(_nodes.begin(), n);
                _preorder[m] = nd;
                m++;
            }
        }
    }
#endif

#if NEWWAY == AAMWAY
    inline void Forest::operator=(const Forest & other) {
        _data               = other._data;
        _nodes.clear();
        _nodes.resize(other._nodes.size());
        _lineages.resize(other._lineages.size());
        _preorder.resize(other._preorder.size());
        _first_pattern = other._first_pattern;
        _npatterns = other._npatterns;
        _gene_tree_log_likelihoods = other._gene_tree_log_likelihoods;
        _ninternals = other._ninternals;
        _nleaves = other._nleaves;
        _log_joining_prob = other._log_joining_prob;
        _increments_and_priors = other._increments_and_priors;

        // copy tree itself

        for (auto othernd : other._nodes) {
            //POL _nodes is not in preorder sequence, so comment below is confusing
            // get number of next node in preorder sequence (serves as index of node in _nodes vector)
            int k = othernd._number;

            if (k>-1) {
                Node* nd = &*next(_nodes.begin(), k);

            // copy parent
                if (othernd._parent) {
                    unsigned parent_number = othernd._parent->_number;
                    Node* parent = &*next(_nodes.begin(), parent_number);
                    nd->_parent = parent;
                }

            // copy left child
                if (othernd._left_child) {
                unsigned left_child_number = othernd._left_child->_number;
                    Node* left_child = &*next(_nodes.begin(), left_child_number);
                    nd->_left_child = left_child;
            }
                else {
                    nd->_left_child = 0;
                }

            // copy right sibling
            if (othernd._right_sib) {
                unsigned right_sib_number = othernd._right_sib->_number;
                Node* right_sib = &*next(_nodes.begin(), right_sib_number);
                nd->_right_sib = right_sib;
            }
            else
                nd->_right_sib = 0;

                nd->_number = othernd._number;
                nd->_name = othernd._name;
                nd->_edge_length = othernd._edge_length;
                nd->_position_in_lineages = othernd._position_in_lineages;
                nd->_partials = othernd._partials;
            }
        }

        unsigned j = 0;
        for (auto othernd : other._lineages) {
            unsigned k = othernd->_number;
            Node* nd = &*next(_nodes.begin(), k);
            _lineages[j] = nd;
            j++;
        }
        
        if (other._preorder.size() > 0) {
            unsigned m = 0;
            for (auto othernd : other._preorder) {
                unsigned n = othernd->_number;
                Node* nd = &*next(_nodes.begin(), n);
                _preorder[m] = nd;
                m++;
            }
        }
    }
#endif

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
        
        // Simulate starting sequence at the root node
        Node * nd = *(_lineages.begin());
        int ndnum = nd->_number;
        assert(ndnum < nnodes);
        for (unsigned i = 0; i < nsites; i++) {
            sequences[ndnum][i] = G::multinomialDraw(lot, basefreq);
        }
        
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

        assert(data);
        Data::data_matrix_t & dm = data->getDataMatrixNonConst();
        
#if NEWWAY == POLWAY
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
#else   //AAMWAY
        // Copy sequences to data object
        for (auto & nd : _nodes) {
            if (!nd._left_child) {
                unsigned t = nd._number;
                assert(t < G::_ntaxa);
                
                // Allocate row t of _data's _data_matrix data member
                assert(dm[t].size() == starting_site);
                dm[t].resize(starting_site + nsites);
                            
                // Translate to state codes and copy
                for (unsigned i = 0; i < nsites; i++) {
                    dm[t][starting_site + i] = (Data::state_t)1 << sequences[t][i];
                }
            }
        }
#endif
    }

    inline unsigned Forest::getNumLineages() const {
        return (unsigned)_lineages.size();
    }
    
    inline unsigned Forest::getNumNodes() const {
        return (unsigned)_preorder.size();
    }

    inline Node * Forest::pullNode() {
#if NEWWAY == POLWAY
        unsigned nnodes = (unsigned)_nodes.size();
        _nodes.resize(nnodes + 1);
        auto nditer = next(_nodes.begin(), nnodes);
        Node * new_nd = &(*nditer);
        assert(new_nd->_number == -1);
        new_nd->clear();
        new_nd->_number = _nleaves + _ninternals;
        _ninternals++;
#else  //AAM
        _nodes.push_back(Node());
        Node* new_nd = &_nodes.back();
        //new_nd->_edge_length = 0.0; // should be Node::_smallest_edge_length?
        new_nd->clear();
        if (_nodes.size() <= G::_ntaxa) {
            assert(_nleaves == _nodes.size() - 1);
            new_nd->_number = (int)_nleaves;
            _nleaves++;
        }
        else {
            new_nd->_number = _nleaves + _ninternals;
            _ninternals++;
        }
#endif
        new_nd->_split.resize(G::_ntaxa);

        return new_nd;
    }
    
    inline void Forest::joinRandomLineagePair(Lot::SharedPtr lot) {
        unsigned n = (unsigned)_lineages.size();
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
        // Algorithm from Yang and Rannala. 1997. MBE 14(7):717-724.
        createTrivialForest();
#if 0
        unsigned nsteps = G::_ntaxa - 1;
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
            cum_height += t;
            advanceAllLineagesBy(t);
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

    //inline void Forest::buildYuleTree() {
    //    createTrivialForest();
    //    unsigned nsteps = G::_ntaxa - 1;
    //    for (unsigned i = 0; i < nsteps; i++) {
    //        // Determine number of lineages remaining
    //        unsigned n = getNumLineages();
    //        assert(n > 1);
    //
    //        // Waiting time to speciation event is Exponential(rate = n*lambda)
    //        // u = 1 - exp(-r*t) ==> t = -log(1-u)/r
    //        double r = G::_sim_lambda*n;
    //        double u = rng->uniform();
    //        double t = -log(1.0 - u)/r;
    //        advanceAllLineagesBy(t);
    //        joinRandomLineagePair(rng);
    //    }
    //    assert(getNumLineages() == 1);
    //    refreshPreorder();
    //}

    inline void Forest::scaleAllEdgeLengthsBy(double scaling_factor) {
        // This function should only be called for complete trees
        assert(_lineages.size() == 1);
        
        // Supplied scaling_factor should be strictly positive
        assert(scaling_factor > 0.0);
        
        for (auto & nd : _nodes) {
            double elen = nd.getEdgeLength();
            nd.setEdgeLength(scaling_factor*elen);
        }
    }
    
    inline void Forest::advanceAllLineagesBy(double dt) {
        // Add t to the edge length of all lineage root nodes, unless there
        // is just one lineage, in which case do nothing
        unsigned n = (unsigned)_lineages.size();
        if (n > 1) {
            for (auto nd : _lineages) {
                double elen = nd->getEdgeLength() + dt;
                assert(elen >= 0.0 || fabs(elen) < Node::_smallest_edge_length);
                nd->setEdgeLength(elen);
                ++n;
            }
        }
    }
    
}
