#pragma once

class Particle;

extern proj::Lot::SharedPtr rng;

namespace proj {

    class ForestPOL {
    
        friend class Particle;
        
        public:
        
            typedef set<Split> treeid_t;
            typedef shared_ptr<ForestPOL> SharedPtr;

            ForestPOL();
            ~ForestPOL();
                
            void        clear();
            string      makeNewick(unsigned precision = 6, bool use_names = true) const;
            void        advanceAllLineagesBy(double dt);
            void        createTrivialForest();
            void        joinRandomLineagePair(Lot::SharedPtr lot);
            unsigned    getNumLineages() const;
            unsigned    getNumNodes() const;
            void        refreshAllPreorders() const;
            void        renumberInternals();
            void        simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites);

        protected:
        
            Node *      pullNode();
            void        stowNode(Node * nd);
            void        joinLineagePair(Node * anc, Node * first, Node * second);
            void        unjoinLineagePair(Node * anc, Node * first, Node * second);
            Node *      findNextPreorder(Node * nd) const;
            void        refreshPreorder(Node::ptr_vect_t & preorder) const;
            void        removeOne(Node::ptr_vect_t & node_vect, Node * del);
            void        addOne(Node::ptr_vect_t & node_vect, Node * add);
            void        removeTwoAddOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        addTwoRemoveOne(Node::ptr_vect_t & node_vect, Node * del1, Node * del2, Node * add);
            void        addTwoRemoveOneAt(Node::ptr_vect_t & node_vect, unsigned pos1, Node * del1, unsigned pos2, Node * del2, Node * add);
            double calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
            
            virtual void operator=(const ForestPOL & other);

            // NOTE: any variables added below here must be copied in operator=
            
            vector<unsigned>        _unused_nodes;
            double                  _forest_height;
            vector<Node>            _nodes;
            Node::ptr_vect_t        _lineages;
                                                
            // Because these can be recalculated at any time, they should not
            // affect the const status of the ForestPOL object
            mutable unsigned                 _next_node_number;
            mutable vector<Node::ptr_vect_t> _preorders;
};
    
    inline ForestPOL::ForestPOL() {
    }

    inline ForestPOL::~ForestPOL() {
        clear();
    }

    inline void ForestPOL::clear() {
        _unused_nodes.clear();
        _forest_height = 0.0;
        _next_node_number = 0;
        _nodes.clear();
        _preorders.clear();
        _lineages.clear();
    }
        
    inline unsigned ForestPOL::getNumLineages() const {
        return (unsigned)_lineages.size();
    }
    
    inline unsigned ForestPOL::getNumNodes() const {
        unsigned n = 0;
        for (auto & preorder : _preorders) {
            n += (unsigned)preorder.size();
        }
        return n;
    }
    
    inline void ForestPOL::joinRandomLineagePair(Lot::SharedPtr lot) {
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
        
        removeTwoAddOne(_lineages, lchild, rchild, ancnd);
    }
    
    inline void ForestPOL::renumberInternals() {
        // First internal node number is the number of leaves
        int next_node_number = (unsigned)G::_taxon_names.size();

        // Renumber internal nodes in postorder sequence for each lineage in turn
        for (auto & preorder : _preorders) {
            for (auto nd : boost::adaptors::reverse(preorder)) {
                if (nd->_left_child) {
                    // nd is an internal node
                    assert(nd->_height != G::_infinity);
                    if (nd->_height > _forest_height)
                        _forest_height = nd->_height;
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
    }
    
    inline string ForestPOL::makeNewick(unsigned precision, bool use_names) const  {
        // Assumes preorders correct
        
        // Place basal polytomy (if there is one) at a height 10% greater than
        // the _forest_height
        double basal_polytomy_height = 0.0; //_forest_height*(1.1);
        bool is_complete_tree = (bool)(_preorders.size() == 1);
        
        //const format basal_subtree_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_name_format( str(format("%%s:%%.%df") % precision) );
        const format tip_node_number_format( str(format("%%d:%%.%df") % precision) );
        const format internal_node_format( str(format("):%%.%df") % precision) );
        
        vector<string> subtree_newicks;
        for (auto preorder : _preorders) {
            string subtree_newick;
            stack<Node *> node_stack;
            for (auto nd : preorder) {
                assert(nd->_number > -1);
                if (nd->_left_child) {
                    subtree_newick += "(";
                    node_stack.push(nd);
                }
                else {
                    double edge_length = nd->_edge_length;
                    if (!nd->_parent) {
                        // Subtree consists of just this one leaf node
                        edge_length += basal_polytomy_height;
                    }
                    if (use_names) {
                        if (precision > 0)
                            subtree_newick += str(format(tip_node_name_format)
                                % nd->_name
                                % edge_length);
                        else
                            subtree_newick += str(format("%s")
                                % nd->_name);
                    } else {
                        if (precision > 0)
                            subtree_newick += str(format(tip_node_number_format)
                                % (nd->_number + 1)
                                % edge_length);
                        else
                            subtree_newick += str(format("%d")
                                % (nd->_number + 1));
                    }
                    if (nd->_right_sib)
                        subtree_newick += ",";
                    else if (nd->_parent) {
                        Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                        double popped_edge_length = popped->_edge_length;
                        while (popped && !popped->_right_sib) {
                            node_stack.pop();
                            if (node_stack.empty()) {
                                if (is_complete_tree) {
                                    subtree_newick += ")";
                                }
                                else {
                                    // This is the root of one of several subtrees, so
                                    // it is important to preserve its edge length
                                    popped_edge_length += basal_polytomy_height;
                                    if (precision > 0)
                                        subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                    else
                                        subtree_newick += ")";
                                }
                                popped = 0;
                            }
                            else {
                                if (precision > 0)
                                    subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                else
                                    subtree_newick += ")";
                                popped = node_stack.top();
                                popped_edge_length = popped->_edge_length;
                            }
                        }
                        if (popped && popped->_right_sib) {
                            node_stack.pop();
                            if (precision > 0) {
                                subtree_newick += str(format(internal_node_format) % popped_edge_length);
                                subtree_newick += ",";
                            }
                            else
                                subtree_newick += "),";
                        }
                    }
                }
            }
            //if (subtree_newick[0] == '(')
            //    subtree_newicks.push_back(str(format(basal_subtree_format) % subtree_newick % basal_polytomy_height));
            //else
                subtree_newicks.push_back(subtree_newick);
        }
        
        string newick;
        if (is_complete_tree)
            newick = subtree_newicks[0];
        else {
            string separator = str(format(":%.5f,") % basal_polytomy_height);
            string insides = join(subtree_newicks, ",");
            newick = str(format("(%s)") % insides);
        }

        return newick;
    }
    
    inline void ForestPOL::advanceAllLineagesBy(double dt) {
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
        
            // Add to to the current forest height
            _forest_height += dt;
        }
    }
    
    inline Node * ForestPOL::pullNode() {
        if (_unused_nodes.empty()) {
            unsigned nleaves = 2*G::_ntaxa - 1;
            throw XProj(str(format("ForestPOL::pullNode tried to return a node beyond the end of the _nodes vector (%d nodes allocated for %d leaves)") % _nodes.size() % nleaves));
        }
        double node_index = _unused_nodes.back();
        _unused_nodes.pop_back();
        
        Node * new_nd = &_nodes[node_index];
        assert(new_nd->_my_index == node_index);
        
        new_nd->clear();
        new_nd->_number = -2;
        new_nd->_split.resize(G::_ntaxa);
        return new_nd;
    }
    
    inline void ForestPOL::stowNode(Node * nd) {
        // Get index of nd in _nodes vector
        int offset = nd->_my_index;
        assert(offset > -1);
        _unused_nodes.push_back(offset);
        nd->clear();
    }
    
    inline void ForestPOL::joinLineagePair(Node * new_nd, Node * subtree1, Node * subtree2) {
        // Note: must call pullNode to obtain new_nd before calling this function
        assert(new_nd);
        assert(subtree1);
        assert(subtree2);
        assert(new_nd->_my_index > -1);
        new_nd->_name        = "anc-" + to_string(new_nd->_my_index);
        new_nd->_left_child  = subtree1;
        new_nd->_edge_length = 0.0;
                
        // Calculate height of the new node
        // (should be the same whether computed via the left or right child)
        double h1 = subtree1->_height + subtree1->_edge_length;
        double h2 = subtree2->_height + subtree2->_edge_length;

        // //temporary! Probably need to reinstate this assert
        //assert(fabs(h1 - h2) < G::_small_enough);

        new_nd->_height = (h1 + h2)/2.0;

        // Finish connecting new trio of nodes
        subtree1->_right_sib = subtree2;
        subtree1->_parent    = new_nd;
        subtree2->_parent    = new_nd;
    }
    
    inline void ForestPOL::unjoinLineagePair(Node * anc, Node * subtree1, Node * subtree2) {
        // Note: be sure to call stowNode for anc after calling this function
        // Reset members set by joinLineagePair function
        anc->_number = -1;
        anc->_name = "";
        anc->_left_child = nullptr;
        anc->_edge_length = 0.0;
        anc->_height = 0.0;
        subtree1->_right_sib = nullptr;
        subtree1->_parent = nullptr;
        subtree2->_parent = nullptr;
    }
    
    inline void ForestPOL::removeOne(Node::ptr_vect_t & v, Node * del) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), del);
        assert(it != v.end());
        v.erase(it);
    }
    
    inline void ForestPOL::addOne(Node::ptr_vect_t & v, Node * add) {
        // Add node
        v.push_back(add);
    }
    
    inline void ForestPOL::removeTwoAddOne(Node::ptr_vect_t & v, Node * del1, Node * del2, Node * add) {
        // Get iterator to first node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), del1);
        assert(it != v.end());
        v.erase(it);
        
        // Get iterator to second node to be deleted and remove from v
        it = find(v.begin(), v.end(), del2);
        assert(it != v.end());
        v.erase(it);
        
        // Add node
        v.push_back(add);
    }
    
    inline void ForestPOL::addTwoRemoveOne(Node::ptr_vect_t & v, Node * add1, Node * add2, Node * rem) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), rem);
        assert(it != v.end());
        v.erase(it);
        
        // Add node2
        v.push_back(add1);
        v.push_back(add2);
    }
    
    inline void ForestPOL::addTwoRemoveOneAt(Node::ptr_vect_t & v, unsigned pos1, Node * add1, unsigned pos2, Node * add2, Node * rem) {
        // Get iterator to node to be deleted and remove from v
        auto it = find(v.begin(), v.end(), rem);
        assert(it != v.end());
        v.erase(it);
        
        // Insert add1 and add2 into v so that they end up in
        // positions pos1 and pos2, respectively
        if (pos1 < pos2) {
            it = v.begin();
            advance(it, pos1);
            v.insert(it, add1);
            
            it = v.begin();
            advance(it, pos2);
            v.insert(it, add2);
        }
        else {
            it = v.begin();
            advance(it, pos2);
            v.insert(it, add2);
            
            it = v.begin();
            advance(it, pos1);
            v.insert(it, add1);
        }
    }
    
    inline void ForestPOL::refreshAllPreorders() const {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _next_node_number = G::_ntaxa;
        _preorders.clear();
        if (_lineages.size() == 0)
            return;
        
        for (auto nd : _lineages) {
            if (nd->_left_child) {
                nd->_number = _next_node_number++;
            }
            
            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to _preorders containing, for now, just the root of the subtree
            _preorders.push_back({nd});
            
            // Now add the nodes above the root in preorder sequence
            Node::ptr_vect_t & preorder_vector = *_preorders.rbegin();
            refreshPreorder(preorder_vector);
        }
    }
    
    inline Node * ForestPOL::findNextPreorder(Node * nd) const {
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

    inline void ForestPOL::refreshPreorder(Node::ptr_vect_t & preorder) const {
        // Assumes preorder just contains the root node when this function is called
        // Also assumes that _next_node_number was initialized prior to calling this function
        assert(preorder.size() == 1);
        
        Node * nd = preorder[0];
        while (true) {
            nd = findNextPreorder(nd);
            if (nd) {
                preorder.push_back(nd);
                if (nd->_left_child)
                    nd->_number = _next_node_number++;
            }
            else
                break;
        }
    }
    
    inline void ForestPOL::operator=(const ForestPOL & other) {
        _forest_height                  = other._forest_height;
        //_next_node_index                = other._next_node_index;
        _unused_nodes                   = other._unused_nodes;
        _next_node_number               = other._next_node_number;
        
        // Create node map: if other._nodes[3]._number = 2, then node_map[2] = 3
        // (i.e. node number 2 is at index 3 in _nodes vector)
        map<unsigned, unsigned> node_map;
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            int other_node_number = other._nodes[i]._number;
            if (other_node_number > -1) {
                node_map[other_node_number] = i;
            }
        }
        
        _nodes.resize(other._nodes.size());
        for (unsigned i = 0; i < other._nodes.size(); ++i) {
            unsigned j = 0;
            int n = 0;
            
            // _left_child
            _nodes[i]._left_child = nullptr;
            if (other._nodes[i]._left_child) {
                n = other._nodes[i]._left_child->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._left_child = &_nodes[j];
            }
            
            // _right_sib
            _nodes[i]._right_sib = nullptr;
            if (other._nodes[i]._right_sib) {
                n = other._nodes[i]._right_sib->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._right_sib = &_nodes[j];
            }

            // _parent
            _nodes[i]._parent = nullptr;
            if (other._nodes[i]._parent) {
                n = other._nodes[i]._parent->_number;
                assert(n > -1);
                j = node_map[n];
                _nodes[i]._parent = &_nodes[j];
            }
            
            _nodes[i]._number      = other._nodes[i]._number;
            _nodes[i]._my_index    = other._nodes[i]._my_index;
            _nodes[i]._name        = other._nodes[i]._name;
            _nodes[i]._edge_length = other._nodes[i]._edge_length;
            _nodes[i]._height      = other._nodes[i]._height;
            _nodes[i]._split       = other._nodes[i]._split;
        }

        // Build _lineages
        _lineages.clear();
        for (auto & nd : _nodes) {
            // Assume that nodes that are being used have _number > -1
            // and nodes serving as the root of a lineage have no parent
            if (nd._number > -1 && !nd._parent)
                _lineages.push_back(&nd);
        }
        assert(_lineages.size() == other._lineages.size());
                
        refreshAllPreorders();
    }

    inline double ForestPOL::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length) {
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
        
    inline void ForestPOL::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites) {
        
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
        unsigned ndnum = nd->_number;
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
}
