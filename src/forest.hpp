#pragma once

extern proj::PartialStore ps;

namespace proj {

class Forest {
    
    friend class Particle;
    
    public:
        Forest();
        ~Forest();
        Forest(const Forest & other);
    
        void operator=(const Forest & other);
    
    private:
        void                        clear();
        void                        setForestData(Data::SharedPtr d, bool partials);
        Node *                      findNextPreorder(Node * nd);
        void                        refreshPreorder();
        double                      calcSubsetLogLikelihood(unsigned i);
        void                        addIncrement(Lot::SharedPtr lot);
        double                      joinTaxa(Lot::SharedPtr lot);
        void                        calcPartialArray(Node* new_nd);
        double                      calcTransitionProbability(Node* child, double s, double s_child);
        pair<unsigned, unsigned>    chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        string                      makeNewick(unsigned precision, bool use_names);
        string                      makePartialNewick(unsigned precision, bool use_names);
        void                        showForest();
        void                        updateNodeList(list<Node *> & node_list, Node * delnode1, Node * delnode2, Node * addnode);
        void                        updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        double                      getTreeHeight();
        double                      getTreeLength();
        double                      getSpeciesTreePrior();
        double                      calcTopologyPrior(unsigned nlineages);
    
        Data::SharedPtr             _data;
        vector<Node *>              _lineages;
        list<Node>                  _nodes;
        vector<Node*>               _preorder;
        unsigned                    _first_pattern = 0;
        unsigned                    _npatterns;
        vector<double>              _gene_tree_log_likelihoods;
        unsigned                    _ninternals;
        unsigned                    _nleaves;
        double                      _log_joining_prob;
        vector<pair<double, double>> _increments_and_priors;
        
};

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
        _nleaves = G::_ntaxa;
        _log_joining_prob = 0.0;
        _increments_and_priors.clear();
    }

    inline Forest::Forest(const Forest & other) {
        clear();
        *this = other;
    }

    inline void Forest::setForestData(Data::SharedPtr d, bool partials) {
        _data = d;
        
        _gene_tree_log_likelihoods.resize(G::_nloci, 0.0);
        
        auto &data_matrix=_data->getDataMatrix();

        _nodes.resize(G::_ntaxa);
        _lineages.reserve(_nodes.size());
        //create taxa
        for (unsigned i = 0; i < G::_ntaxa; i++) {
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
            nd->_partials.resize(G::_nloci); // _partials contains a vector of partials for each locus
            
            // replace all spaces with underscores so that other programs do not have
              // trouble parsing tree descriptions
            std::string name = G::_taxon_names[i];
            boost::replace_all(name, " ", "_");
            nd->_name = name;
        }

        for (unsigned index = 0; index < G::_nloci; index ++) {
            for (auto &nd:_lineages) {
                Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(index);
                _first_pattern = gene_begin_end.first;
                _npatterns = _data->getNumPatternsInSubset(index);
                
                if (!nd->_left_child) {

                    if (!G::_save_memory || (G::_save_memory && partials)) { // if save memory setting, don't set tip partials yet
//                        for (unsigned i=0; i<G::_nloci; i++) {
                            nd->_partials[index]=ps.getPartial(_npatterns*4);
                            for (unsigned p=0; p<_npatterns; p++) {
                                unsigned pp = _first_pattern+p;
                                for (unsigned s=0; s<G::_nstates; s++) {
                                    Data::state_t state = (Data::state_t)1 << s;
                                    Data::state_t d = data_matrix[nd->_number][pp];
                                    double result = state & d;
                                    (*nd->_partials[index])[p*G::_nstates + s] = (result == 0.0 ? 0.0:1.0);
                                }
                            }
//                        }
                    }
                }
            }
        }
    }

    inline double Forest::calcSubsetLogLikelihood(unsigned i) {
        auto &counts = _data->getPatternCounts();
        _npatterns = _data->getNumPatternsInSubset(i); // TODO: make a vector of patterns by gene?
        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
        _first_pattern = gene_begin_end.first;

        for (auto &nd:_lineages) {
            double log_like = 0.0;
            for (unsigned p=0; p<_npatterns; p++) {
                double site_like = 0.0;
                for (unsigned s=0; s<G::_nstates; s++) {
                    double partial = (*nd->_partials[i])[p*G::_nstates+s];
                    site_like += 0.25*partial;
                }
                assert(site_like>0);
                log_like += log(site_like)*counts[_first_pattern+p];
            }
            _gene_tree_log_likelihoods[i] += log_like;
//            debugLogLikelihood(nd, log_like);
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
            for (unsigned i=0; i<G::_nloci; i++) {
                for (auto &nd:_lineages) {
                    if (nd->_partials[i] == nullptr) {
                        nd->_partials[i] = ps.getPartial(_npatterns*4);
                        calcPartialArray(nd);
                    }
                }
            }
        }
        
        pair <unsigned, unsigned> t = make_pair (0, 1); // if there is only choice, 0 and 1 will be chosen
        t = chooseTaxaToJoin(nlineages, lot);
        
        subtree1 = _lineages[t.first];
        subtree2 = _lineages[t.second];

        assert (subtree1 != subtree2);
        
        //new node is always needed
        Node nd;
        _nodes.push_back(nd);
        Node* new_nd = &_nodes.back();

        new_nd->_parent=0;
        new_nd->_number=_nleaves+_ninternals;
        new_nd->_edge_length=0.0;
        _ninternals++;
        new_nd->_right_sib=0;

        new_nd->_left_child=subtree1;
        subtree1->_right_sib=subtree2;

        subtree1->_parent=new_nd;
        subtree2->_parent=new_nd;

        // calculate new partials
        new_nd->_partials.resize(G::_nloci);
        for (unsigned index = 0; index < G::_nloci; index++) {
            assert (new_nd->_partials[index] == nullptr);
            _npatterns = _data->getNumPatternsInSubset(index);
            new_nd->_partials[index]=ps.getPartial(_npatterns*4);
        }
        assert(new_nd->_left_child->_right_sib);

        if (G::_save_memory) {
            for (auto &nd:_lineages) {
                for (unsigned index = 0; index < G::_nloci; index++) {
                    if (nd->_partials[index] == nullptr) {
                        _npatterns = _data->getNumPatternsInSubset(index);
                        nd->_partials[index] = ps.getPartial(_npatterns*4);
                        calcPartialArray(nd);
                    }
                }
            }
        }
        
        bool calc_likelihood = true;
        
        if (calc_likelihood) {
            calcPartialArray(new_nd);
        }

        for (unsigned index = 0; index < G::_nloci; index++) {
            subtree1->_partials[index]=nullptr; // throw away subtree partials now, no longer needed
            subtree2->_partials[index]=nullptr;
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
               for (unsigned index = 0; index<G::_nloci; index++) {
                   nd._partials[index]=nullptr;
               }
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
                    new_nd->_partials[i]=ps.getPartial(_npatterns*4);
                    for (unsigned p=0; p<_npatterns; p++) {
                        unsigned pp = _first_pattern+p;
                        for (unsigned s=0; s<G::_nstates; s++) {
                            Data::state_t state = (Data::state_t)1 << s;
                            Data::state_t d = data_matrix[new_nd->_number][pp];
                            double result = state & d;
                            (*new_nd->_partials[i])[p*G::_nstates+s]= (result == 0.0 ? 0.0:1.0);
                        }
                    }
                }
            }
        }

        for (unsigned i=0; i<G::_nloci; i++) {
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatternsInSubset(i);
            
            auto & parent_partial_array = *(new_nd->_partials[i]);
            for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {

                if (child->_partials[i] == nullptr) {
                    child->_partials[i] = ps.getPartial(_npatterns*4);
                    calcPartialArray(child);
                }
                assert (child->_partials[i] != nullptr);
                auto & child_partial_array = *(child->_partials[i]);

                for (unsigned p = 0; p < _npatterns; p++) {
                    for (unsigned s = 0; s <G::_nstates; s++) {
                        double sum_over_child_states = 0.0;
                        for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                            double child_transition_prob = calcTransitionProbability(child, s, s_child);
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

    inline double Forest::calcTransitionProbability(Node* child, double s, double s_child) {
        double child_transition_prob = 0.0;

        if (G::_model == "JC" ) {
            double expterm = exp(-4.0*(child->_edge_length)/3.0); // TODO: need to include relative rates
            double prsame = 0.25+0.75*expterm;
            double prdif = 0.25 - 0.25*expterm;

            child_transition_prob = (s == s_child ? prsame : prdif);
            assert (child_transition_prob > 0.0);
            return child_transition_prob;
        }

        if (G::_model == "HKY") { // TODO: add HKY
            double pi_A = G::_base_frequencies[0];
            double pi_C = G::_base_frequencies[1];
            double pi_G = G::_base_frequencies[2];
            double pi_T = G::_base_frequencies[3];

            double pi_j = 0.0;
            double PI_J = 0.0;

            double phi = (pi_A+pi_G)*(pi_C+pi_T)+G::_kappa*(pi_A*pi_G+pi_C*pi_T);
            double beta_t = 0.5*(child->_edge_length)/phi; // TODO: need to include relative rates


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
                const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
                const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
                const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
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
                                        newick += boost::str(boost::format(tip_node_name_format)
                                            % nd->_name
                                            % nd->_edge_length);
                                        } else {
                                        newick += boost::str(boost::format(tip_node_number_format)
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
                                                if (lineage->_edge_length != 0.0) {
                                                    newick += boost::str(boost::format(internal_node_format) % lineage->_edge_length);
                                                }
                                                popped = 0;
                                            }
                                            else {
                                                newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                                                popped = node_stack.top();
                                            }
                                        }
                                        if (popped && popped->_right_sib) {
                                            node_stack.pop();
                                            newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                                            newick += ",";
                                        }
                                }   // leaf node
                            }
                            nd = findNextPreorder(nd);
                        }   // while (subnd)...

                        if (i < _lineages.size() - 1)
                            newick += ",";
                        ++i;
                    }
                    newick += ")";

                    return newick;
                }
            }

    inline string Forest::makePartialNewick(unsigned precision, bool use_names) {
            // this function makes a newick string for a partially constructed tree
            string newick = "(";
            const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
            const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
            const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
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
                            newick += boost::str(boost::format(tip_node_name_format)
                                % nd->_name
                                % nd->_edge_length);
                        } else {
                            newick += boost::str(boost::format(tip_node_number_format)
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
                                    newick += boost::str(boost::format(internal_node_format) % lineage->_edge_length);
                                    popped = 0;
                                }
                                else {
                                    newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                                    popped = node_stack.top();
                                }
                            }
                            if (popped && popped->_right_sib) {
                                node_stack.pop();
                                newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
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
                nd->_partials = othernd._partials; // TODO: will this copy the vector correctly?
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

}
