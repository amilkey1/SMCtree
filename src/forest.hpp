#pragma once
#include <cmath>

extern proj::PartialStore ps;
std::mutex mtx;

namespace proj {

class Forest {
    
    friend class Particle;
    friend class ForestExtension;
    
    typedef std::shared_ptr<const Forest> ConstSharedPtr;
    
    public:
        Forest();
        ~Forest();
        Forest(const Forest & other);
        string makeNewick(unsigned precision, bool use_names);
        void operator=(const Forest & other);
        
        //POL added below
        void    createTrivialForest();
        void    simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites, double clock_rate);
        double  getLogLikelihood();
        double  getHeightFirstSplit();
        double  getHeightSecondIncr();
        double  getHeightThirdIncr();
    
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
        double joinPriorPrior(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils);
        void calcPartialArray(Node* new_nd);
        double calcPartialArrayLazyJC(Node * new_nd, const Node * lchild, const Node * rchild, double clock_rate) const;
        double calcPartialArrayLazyHKY(Node * new_nd, const Node * lchild, const Node * rchild, double clock_rate) const;
        double calcTransitionProbabilityLazyJC(double s, double s_child, double edge_length, unsigned locus, double clock_rate, unsigned rate_categ) const;
        double calcTransitionProbabilityLazyHKY(double s, double s_child, double edge_length, unsigned locus, double clock_rate, unsigned rate_categ) const;
        
        //POL added below
        Node * pullNode();
        unsigned getNumLineages() const;
        unsigned getNumNodes() const;
        double calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length, double clock_rate);
        //POL added above

        pair<unsigned, unsigned> chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        string makePartialNewick(unsigned precision, bool use_names);
        void showForest();
        void updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);

        double getTreeHeight();
        double getTreeLength();
        double getTreePrior();
        double calcTopologyPrior(unsigned nlineages);
        void   clearPartials();
        double getLineageHeight(Node* nd);
    
        double    drawBirthDeathIncrement(Lot::SharedPtr lot, double age, double estimated_lambda, double estimated_mu, double estimated_root_age);
        pair<bool, vector<bool>>    checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets);
    
        double  getForestHeight() {return _tree_height;}
        PartialStore::partial_t         pullPartial();
    
        void                            setLogLikelihood(double log_likelihood) {_log_likelihood = log_likelihood;}
        void                            addIncrAndJoin(double incr, const Split & lsplit, const Split & rsplit, ForestExtension & gfx);
        void                            advanceAllLineagesBy(double increment);
        void                            refreshAllPreorders();
        
        Data::SharedPtr             _data;
        vector<Node *>              _lineages;
        vector<Node>                _nodes;
        vector<Node*>               _preorder;
        unsigned                    _first_pattern;
        unsigned                    _npatterns;
        vector<double> mutable      _gene_tree_log_likelihoods;
        double                      _log_likelihood;
        unsigned                    _ninternals;
        unsigned                    _nleaves;
        double                      _log_joining_prob;
        vector<double>              _increments;
        vector<pair<Node*, Node*>>  _node_choices;
        double                      _partial_count;
        double                      _weight_correction; // correct for taxon set constraints
        double                      _first_split_prior;
        double                      _tree_height;
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
        _log_likelihood = 0.0;
        _ninternals = 0;
        _log_joining_prob = 0.0;
        _increments.clear();
        _nleaves = 0;
        _node_choices.clear();
        _partial_count = 0;
        _weight_correction = 0.0;
        _tree_height = 0.0;
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
        _nodes.reserve(nnodes);
        _nodes.resize(2*G::_ntaxa - 1);
        
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
            _nodes[i]._position_in_lineages = i;
            _nodes[i]._set_partials = true;
        }
        refreshPreorder();
        _partial_count = 0;
    }
    
    inline void Forest::setForestData(Data::SharedPtr d, bool partials) {
        _data = d;
        
        _gene_tree_log_likelihoods.resize(G::_nloci, 0.0);
        
        auto &data_matrix=_data->getDataMatrix();

        unsigned nnodes = 2*G::_ntaxa - 1;
        _nodes.reserve(nnodes);
        _nodes.resize(2*G::_ntaxa - 1);

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
            _nodes[i]._position_in_lineages=i;
            _nodes[i]._partials = nullptr;
            _nodes[i]._name = G::_taxon_names[i];
            _nodes[i]._set_partials = true;
            _lineages.push_back(&_nodes[i]);
            
            string taxon_name = G::_taxon_names[i];
            _nodes[i]._name = taxon_name;
            _nodes[i]._split.resize(G::_ntaxa);
            _nodes[i]._split.setBitAt(i);
        }
        
        for (unsigned i = G::_ntaxa; i < _nodes.size(); i++) {
            _nodes[i]._use_in_likelihood = false;
        }
        
        _nleaves = G::_ntaxa;

        for (auto &nd:_lineages) {
            mtx.lock();
            nd->_partials=ps.getPartial();
            mtx.unlock();
            
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(0);
            _first_pattern = gene_begin_end.first;
            _npatterns = _data->getNumPatterns();
            
            unsigned npartials_used = 0;
            
            
            if (!nd->_left_child) {

                if (!G::_save_memory || (G::_save_memory && partials)) { // if save memory setting, don't set tip partials yet
                    for (unsigned i=0; i<G::_nloci; i++) {
                        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
                        double first_pattern = gene_begin_end.first;
                        double last_pattern = gene_begin_end.second;
                        double npatterns_in_subset = last_pattern - first_pattern;
                        
                        for (unsigned p=0; p<npatterns_in_subset; p++) {
                            for (unsigned step = 0; step < G::_gamma_rate_cat.size(); step++) {
                                unsigned pp = first_pattern + p;
//                                unsigned start = npartials_used;
//                                unsigned start = step * G::_gamma_rate_cat.size() * npatterns_in_subset + npartials_used;
//                                unsigned pxnstates = p*G::_nstates + start;
                                unsigned pxnstates = npartials_used;
                                
                                for (unsigned s=0; s<G::_nstates; s++) {
                                    Data::state_t state = (Data::state_t)1 << s;
                                    Data::state_t d = data_matrix[nd->_number][pp];
                                    double result = state & d;
                                    (nd->_partials->_v)[pxnstates + s] = (result == 0.0 ? 0.0:1.0);
                                }
                                npartials_used += G::_nstates;
                            }
//                            npartials_used += npatterns_in_subset * G::_nstates;
                        }
                    }
                }
            }
        }
    }

    inline double Forest::calcSubsetLogLikelihood(unsigned i) {

        for (auto &nd:_nodes) {
            nd._use_in_likelihood = false;
        }

        for (auto &nd:_lineages) {
            nd->_use_in_likelihood = true;
        }

        _gene_tree_log_likelihoods[i] = 0.0;

        auto &counts = _data->getPatternCounts();
        _npatterns = _data->getNumPatternsInSubset(i);
        Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
        _first_pattern = gene_begin_end.first;

        for (auto &nd:_lineages) {
            if (nd->_use_in_likelihood) {
                assert (nd->_partials != nullptr); // ignore fossils and fake nodes in likelihood calculations
                double log_like = 0.0;
                unsigned count = 0;
                unsigned skip = G::_gamma_rate_cat.size();
//                for (unsigned p=_first_pattern; p<_npatterns + _first_pattern; p++) { // TODO: need to skip nratecat replicates of each pattern
                for (unsigned p = 0; p < _npatterns * G::_gamma_rate_cat.size(); p++) {
                    double site_like = 0.0;
                    for (unsigned s=0; s<G::_nstates; s++) {
                        double partial = (nd->_partials->_v)[p*G::_nstates+s];
                        site_like += G::_base_frequencies[s]*partial;
                    }
                    if (site_like == 0) {
                        showForest();
                    }
                    assert(site_like>0);
                    log_like += log(site_like)*counts[count];
                    p += skip;
                    count++;
                }
                _gene_tree_log_likelihoods[i] += log_like;
                }
        }
        return _gene_tree_log_likelihoods[i];
    }

    double Forest::drawBirthDeathIncrement(Lot::SharedPtr lot, double age, double estimated_lambda, double estimated_mu, double estimated_root_age) {
        // this function draws an increment but does not add it
        // fossils are age constraints and not part of the tree
        
        // if there is only one lineage left, all extant taxa have been joined and remaining fossils must be added
        
        // birth death
        double cum_height = _tree_height;
        unsigned n = G::_ntaxa;
        
        double birth_rate = estimated_lambda;
        
        double death_rate = estimated_mu;
        
        assert (birth_rate >= death_rate);
        
        double t = 0.0;
        
        assert (n > 1);
        
        double troot = estimated_root_age;
                
        double u = lot->uniform();
            
        double lambda_minus_mu = birth_rate - death_rate;
        
        double b = G::_step;
        double k = b + 1;
        
        double phi = pow((1-u), (1/(n-k-1)));
        
        double log_phi = log(phi);
        double log_lambda = log(birth_rate);
        
        double log_mu = log(death_rate);
        
        double log_x = -1 * lambda_minus_mu * cum_height;
        double log_z = -1 * lambda_minus_mu * troot;

        
        double term1_num = log_lambda + log_phi + log_x;
        double term2_num = log_lambda + log_phi + log_z;
        double term3_num = log_lambda + log_z;
        double term4_num = log_mu + log_x + log_z;

        vector<double> num_values;

        if (G::_mu == 0) {
            num_values.push_back(term1_num);
            num_values.push_back(term2_num);
            num_values.push_back(term3_num);
        }

        else {
            num_values.push_back(term1_num);
            num_values.push_back(term2_num);
            num_values.push_back(term3_num);
            num_values.push_back(term4_num);
        }

        double max_logv = *max_element(num_values.begin(), num_values.end());

        double factored_sum = 0.0;
        unsigned count = 0;
        for (auto & logv : num_values) {
            if (count == 1 || count == 3) {
                factored_sum -= exp(logv - max_logv);
            }
            else {
                factored_sum += exp(logv - max_logv);
            }
            count++;
        }
        double numerator_test = max_logv + log(factored_sum);

        double term1_denom = log_mu + log_phi + log_x;
        double term2_denom = log_mu + log_phi + log_z;
        double term3_denom = log_lambda;
        double term4_denom = log_mu + log_x;
        

        vector<double> denom_values;
        if (G::_mu == 0) {
            denom_values.push_back(term3_denom);
        }
        else {
            denom_values.push_back(term1_denom);
            denom_values.push_back(term2_denom);
            denom_values.push_back(term3_denom);
            denom_values.push_back(term4_denom);
        }

        max_logv = *max_element(denom_values.begin(), denom_values.end());

        factored_sum = 0.0;
        count = 0;
        for (auto & logv : denom_values) {
            if (count == 1 || count == 3) {
                factored_sum -= exp(logv - max_logv);
            }
            else {
                factored_sum += exp(logv - max_logv);
            }
            count++;
        }
        double denominator_test = max_logv + log(factored_sum);

        double new_height = numerator_test - denominator_test;
        
        new_height /= (death_rate - birth_rate);

        assert (new_height > 0);
        assert (new_height <= troot + 0.01);
        assert (new_height <= estimated_root_age + 0.01);
        
        t = new_height - cum_height;
        
        return t;
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
        
        // calculate prior on branch lengths
        for (auto &t:_increments) {
            birth_death_prior += log(t);
        }
        
        double log_joining_prob = 0.0;
        // calculate topology prior
        for (unsigned nlineages=G::_ntaxa; nlineages > 1; nlineages--) {
            log_joining_prob += -log(0.5*nlineages*(nlineages-1));
        }
        
        birth_death_prior += log_joining_prob;
        
        // TODO: if estimating lambda and mu and clock rate and root age, add those priors
        
        assert(birth_death_prior == birth_death_prior);
        
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

    inline double Forest::calcPartialArrayLazyHKY(Node * new_nd, const Node * lchild, const Node * rchild, double clock_rate) const {
        // Get pattern counts
        auto counts = _data->getPatternCounts();
        
        // Determine if there is an edge length extension (this would be the
        // case if new_nd comes from a gene forest extension)
        double lchild_stem_height = lchild->_height + lchild->_edge_length;
        double rchild_stem_height = rchild->_height + rchild->_edge_length;
        assert(fabs(lchild_stem_height - rchild_stem_height) < G::_small_enough);
        
        // Calculate the edge length extension
        double edgelen_extension = new_nd->_height - lchild_stem_height;
        
        // Edge length extension may be slightly negative due to roundoff
        assert(edgelen_extension >= -G::_small_enough);
        if (edgelen_extension < 0.0) {
            edgelen_extension = 0.0;
        }
        
        unsigned n_likelihood_calculations = 1;
            if (G::_plus_G) {
                n_likelihood_calculations = (unsigned) G::_gamma_rate_cat.size();
            }
            double weight = 0.0;
        vector<vector<double>> log_likelihoods(G::_nloci);
        vector<vector<double>> prev_loglikelihoods(G::_nloci);
        vector<double> log_likelihoods_by_pattern(G::_nloci);
        
        unsigned npartials_used = 0.0;
        double log_n_rate_categ = log(G::_gamma_rate_cat.size());

        for (unsigned i=0; i<G::_nloci; i++) {
            _gene_tree_log_likelihoods[i] = 0.0;
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            auto & parent_partial_array = new_nd->_partials->_v;
            unsigned first_pattern = gene_begin_end.first;
            unsigned last_pattern = gene_begin_end.second;
            unsigned npatterns_in_subset = last_pattern - first_pattern;

            for (unsigned p = 0; p < npatterns_in_subset; p++) {

            for (unsigned step = 0; step < n_likelihood_calculations; step++) {
                for (const Node * child : {lchild, rchild})  {
                    assert(child->_partials);
                    auto & child_partial_array = child->_partials->_v;
//                        unsigned start = step * G::_gamma_rate_cat.size() * npatterns_in_subset + npartials_used;
                        unsigned start = npartials_used;
//                        unsigned pxnstates = p*G::_nstates + start;
                        unsigned pxnstates = start;

                        for (unsigned s = 0; s < G::_nstates; s++) {
                            double sum_over_child_states = 0.0;
                            for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                                double child_transition_prob = calcTransitionProbabilityLazyHKY(s, s_child, child->_edge_length + edgelen_extension, i, clock_rate, step);
                                double child_partial = child_partial_array[pxnstates + s_child];
                                                        
                                sum_over_child_states += child_transition_prob * child_partial;
                            }   // child state loop
                            
                            parent_partial_array[pxnstates + s] *= sum_over_child_states;
                        }   // parent state loop
                    }   // pattern loop
                    npartials_used += G::_nstates; // this is the total number of partial steps the step took up
                }
            }
        
                // Compute the ratio of after to before likelihoods
                //TODO: make more efficient
                double prev_loglike = 0.0;
                double curr_loglike = 0.0;
                
                auto & newnd_partial_array = new_nd->_partials->_v;
                auto & lchild_partial_array = lchild->_partials->_v;
                auto & rchild_partial_array = rchild->_partials->_v;
                
                npartials_used = 0;
            
            for (unsigned p = 0; p < npatterns_in_subset; p++) {
                vector<double> prev_log_likelihoods_for_step(G::_gamma_rate_cat.size());
                 vector<double> log_likelihoods_for_step(G::_gamma_rate_cat.size());
                 unsigned pp = first_pattern + p;
                
                for (unsigned step = 0; step < G::_gamma_rate_cat.size(); step++) {
                     unsigned start = npartials_used;
                     unsigned pxnstates = npartials_used;

                    
                //unsigned count = counts[pp];
                double left_sitelike = 0.0;
                double right_sitelike = 0.0;
                double newnd_sitelike = 0.0;
        #if defined (UNROLL_LOOPS)
                // loop 0
                unsigned s = 0;
                left_sitelike += G::_base_frequencies[s]*lchild_partial_array[pxnstates + s];
                right_sitelike += G::_base_frequencies[s]*rchild_partial_array[pxnstates + s];
                newnd_sitelike += G::_base_frequencies[s]*newnd_partial_array[pxnstates + s];

                // loop 1
                s = 1;
                left_sitelike += G::_base_frequencies[s]*lchild_partial_array[pxnstates + s];
                right_sitelike += G::_base_frequencies[s]*rchild_partial_array[pxnstates + s];
                newnd_sitelike += G::_base_frequencies[s]*newnd_partial_array[pxnstates + s];

                // loop 2
                s = 2;
                left_sitelike += G::_base_frequencies[s]*lchild_partial_array[pxnstates + s];
                right_sitelike += G::_base_frequencies[s]*rchild_partial_array[pxnstates + s];
                newnd_sitelike += G::_base_frequencies[s]*newnd_partial_array[pxnstates + s];

                // loop 3
                s = 3;
                left_sitelike += G::_base_frequencies[s]*lchild_partial_array[pxnstates + s];
                right_sitelike += G::_base_frequencies[s]*rchild_partial_array[pxnstates + s];
                newnd_sitelike += G::_base_frequencies[s]*newnd_partial_array[pxnstates + s];

        #else
                for (unsigned s = 0; s < G::_nstates; s++) {
                    left_sitelike += G::_base_frequencies[s]*(lchild_partial_array)[pxnstates + s];
                    right_sitelike += G::_base_frequencies[s]*(rchild_partial_array)[pxnstates + s];
                    newnd_sitelike += G::_base_frequencies[s]*(newnd_partial_array)[pxnstates + s];
                }
        #endif
                    prev_log_likelihoods_for_step[step] += log(left_sitelike);
                    prev_log_likelihoods_for_step[step] += log(right_sitelike);
                    log_likelihoods_for_step[step] += log(newnd_sitelike);
                    
                    npartials_used += G::_nstates;

            }
                // calculate log sum over all rate categories for the site
                 if (G::_plus_G) {
                           assert (log_likelihoods_for_step.size() == G::_gamma_rate_cat.size());
                           assert (prev_log_likelihoods_for_step.size() == G::_gamma_rate_cat.size());
                           double curr_sum = (G::calcLogSum(log_likelihoods_for_step) - log_n_rate_categ) * counts[pp];
                           double prev_sum = (G::calcLogSum(prev_log_likelihoods_for_step) - log_n_rate_categ) * counts[pp];
                         _gene_tree_log_likelihoods[i] += curr_sum;
                     weight += curr_sum - prev_sum;
                     }
                 else {
                         assert (log_likelihoods_for_step.size() == 1);
                         assert (prev_log_likelihoods_for_step.size() == 1);
                         double curr_sum = log_likelihoods_for_step[0];
                         double prev_sum = prev_log_likelihoods_for_step[0];
                       _gene_tree_log_likelihoods[i] += curr_sum;
                     weight += curr_sum - prev_sum;
                 }
             }
             }
         return weight;
     }

    inline double Forest::calcPartialArrayLazyJC(Node * new_nd, const Node * lchild, const Node * rchild, double clock_rate) const {
        // Get pattern counts
        auto counts = _data->getPatternCounts();
        
        // Determine if there is an edge length extension (this would be the
        // case if new_nd comes from a gene forest extension)
        double lchild_stem_height = lchild->_height + lchild->_edge_length;
        double rchild_stem_height = rchild->_height + rchild->_edge_length;
        assert(fabs(lchild_stem_height - rchild_stem_height) < G::_small_enough);
        
        // Calculate the edge length extension
        double edgelen_extension = new_nd->_height - lchild_stem_height;
        
        // Edge length extension may be slightly negative due to roundoff
        assert(edgelen_extension >= -G::_small_enough);
        if (edgelen_extension < 0.0) {
            edgelen_extension = 0.0;
        }

        unsigned n_likelihood_calculations = 1;
        if (G::_plus_G) {
            n_likelihood_calculations = (unsigned) G::_gamma_rate_cat.size();
        }
        double weight = 0.0;
        vector<vector<double>> log_likelihoods(G::_nloci);
        vector<vector<double>> prev_loglikelihoods(G::_nloci);
        vector<double> log_likelihoods_by_pattern(G::_nloci);
        
        unsigned npartials_used = 0.0;
        double log_n_rate_categ = log(G::_gamma_rate_cat.size());

        
        for (unsigned i=0; i<G::_nloci; i++) {
            _gene_tree_log_likelihoods[i] = 0.0;
            Data::begin_end_pair_t gene_begin_end = _data->getSubsetBeginEnd(i);
            auto & parent_partial_array = new_nd->_partials->_v;
            unsigned first_pattern = gene_begin_end.first;
            unsigned last_pattern = gene_begin_end.second;
            unsigned npatterns_in_subset = last_pattern - first_pattern;

            for (unsigned p = 0; p < npatterns_in_subset; p++) {

        for (unsigned step = 0; step < n_likelihood_calculations; step++) {
            for (const Node * child : {lchild, rchild})  {
                assert(child->_partials);
                auto & child_partial_array = child->_partials->_v;

                double pr_same = calcTransitionProbabilityLazyJC(0, 0, child->_edge_length + edgelen_extension, i, clock_rate, step);
                double pr_diff = calcTransitionProbabilityLazyJC(0, 1, child->_edge_length + edgelen_extension, i, clock_rate, step);
//                for (unsigned p = first_pattern; p < last_pattern; p++) {
//                for (unsigned p = 0; p < npatterns_in_subset; p++) {
                    unsigned start = npartials_used;
//                    unsigned start = step * G::_gamma_rate_cat.size() * npatterns_in_subset + npartials_used;
//                    unsigned pxnstates = p*G::_nstates + start;
                unsigned pxnstates = start;

    #if defined (UNROLL_LOOPS)
                    // unroll parent loop
                    assert (G::_nstates == 4);
                    unsigned s = 0;
                    double sum_over_child_states = 0.0;

                        // child state subloop 0
                    unsigned s_child = 0;
                    double child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_same * child_partial;

                        // child state subloop 1
                    s_child = 1;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;


                        // child state subloop 2
                    s_child = 2;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 3
                    s_child = 3;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                    parent_partial_array[pxnstates + s] *= sum_over_child_states;

                    s = 1;
                    sum_over_child_states = 0.0;
                    // child state subloop 0
                    s_child = 0;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 1
                    s_child = 1;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_same * child_partial;

                        // child state subloop 2
                    s_child = 2;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 3
                    s_child = 3;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                    parent_partial_array[pxnstates + s] *= sum_over_child_states;

                    s = 2;
                    sum_over_child_states = 0.0;
                    // child state subloop 0
                    s_child = 0;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 1
                    s_child = 1;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 2
                    s_child = 2;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_same * child_partial;

                        // child state subloop 3
                    s_child = 3;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                    parent_partial_array[pxnstates + s] *= sum_over_child_states;

                    s = 3;
                    sum_over_child_states = 0.0;
                    // child state subloop 0
                    s_child = 0;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 1
                    s_child = 1;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 2
                    s_child = 2;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_diff * child_partial;

                        // child state subloop 3
                    s_child = 3;
                    child_partial = child_partial_array[pxnstates + s_child];

                    sum_over_child_states += pr_same * child_partial;

                    parent_partial_array[pxnstates + s] *= sum_over_child_states;
    #else
                    for (unsigned s = 0; s < G::_nstates; s++) {
                        double sum_over_child_states = 0.0;
                        for (unsigned s_child = 0; s_child < G::_nstates; s_child++) {
                            double child_transition_prob = (s == s_child ? pr_same : pr_diff);
                            double child_partial = (child_partial_array)[pxnstates + s_child];
//                            cout << "child transition prob = " << child_transition_prob << endl;
//                            cout << "child partial = " << child_partial << endl;

                            sum_over_child_states += child_transition_prob * child_partial;
//                            cout << "sum over child states = " << sum_over_child_states << endl;
                        }   // child state loop

                        (parent_partial_array)[pxnstates + s] *= sum_over_child_states;
                    }   // parent state loop
    #endif
                }   // pattern loop
            npartials_used += G::_nstates; // this is the total number of partial steps the step took up
            }
            }
        
        // Compute the ratio of after to before likelihoods
        //TODO: make more efficient
            double prev_loglike = 0.0;
            double curr_loglike = 0.0;

            
        auto & newnd_partial_array = new_nd->_partials->_v;
        auto & lchild_partial_array = lchild->_partials->_v;
        auto & rchild_partial_array = rchild->_partials->_v;
            
                npartials_used = 0;
            for (unsigned p = 0; p < npatterns_in_subset; p++) {
                vector<double> prev_log_likelihoods_for_step(G::_gamma_rate_cat.size());
                vector<double> log_likelihoods_for_step(G::_gamma_rate_cat.size());
                unsigned pp = first_pattern + p; // TODO: or just first_pattern?

                for (unsigned step = 0; step < G::_gamma_rate_cat.size(); step++) {
                unsigned start = npartials_used;
//            unsigned start = step * G::_gamma_rate_cat.size() * npatterns_in_subset + npartials_used;
//            unsigned pxnstates = p*G::_nstates + start;
                    unsigned pxnstates = npartials_used;
            
            //unsigned count = counts[pp];
            double left_sitelike = 0.0;
            double right_sitelike = 0.0;
            double newnd_sitelike = 0.0;
#if defined (UNROLL_LOOPS)
            // loop 0
            unsigned s = 0;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];

            // loop 1
            s = 1;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];

            // loop 2
            s = 2;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];

            // loop 3
            s = 3;
            left_sitelike += 0.25*lchild_partial_array[pxnstates + s];
            right_sitelike += 0.25*rchild_partial_array[pxnstates + s];
            newnd_sitelike += 0.25*newnd_partial_array[pxnstates + s];

#else
            for (unsigned s = 0; s < G::_nstates; s++) {
                left_sitelike += 0.25*(lchild_partial_array)[pxnstates + s];
                right_sitelike += 0.25*(rchild_partial_array)[pxnstates + s];
                newnd_sitelike += 0.25*(newnd_partial_array)[pxnstates + s];
            }
#endif
                    
                    prev_log_likelihoods_for_step[step] += log(left_sitelike);
                    prev_log_likelihoods_for_step[step] += log(right_sitelike);
                    log_likelihoods_for_step[step] += log(newnd_sitelike);
                    
                    npartials_used += G::_nstates;
        }
                // calculate log sum over all rate categories for the site
                if (G::_plus_G) {
                          assert (log_likelihoods_for_step.size() == G::_gamma_rate_cat.size());
                          assert (prev_log_likelihoods_for_step.size() == G::_gamma_rate_cat.size());
                          double curr_sum = (G::calcLogSum(log_likelihoods_for_step) - log_n_rate_categ) * counts[pp];
                          double prev_sum = (G::calcLogSum(prev_log_likelihoods_for_step) - log_n_rate_categ) * counts[pp];
                        _gene_tree_log_likelihoods[i] += curr_sum;
                    weight += curr_sum - prev_sum;
                    }
                else {
                        assert (log_likelihoods_for_step.size() == 1);
                        assert (prev_log_likelihoods_for_step.size() == 1);
                        double curr_sum = log_likelihoods_for_step[0];
                        double prev_sum = prev_log_likelihoods_for_step[0];
                      _gene_tree_log_likelihoods[i] += curr_sum;
                    weight += curr_sum - prev_sum;
                }
            }
            }
        return weight;
    }

    inline double Forest::calcTransitionProbabilityLazyJC(double s, double s_child, double edge_length, unsigned locus, double clock_rate, unsigned rate_categ) const {
        double child_transition_prob = 0.0;
        double relative_rate = G::_double_relative_rates[locus];
        assert (relative_rate > 0.0);
        assert (rate_categ < G::_gamma_rate_cat.size());
        double gamma_rate = G::_gamma_rate_cat[rate_categ];

            if (s == s_child) {
                child_transition_prob = 0.25 + 0.75*exp(-4.0 * clock_rate * edge_length * gamma_rate * relative_rate / 3.0);
            }
            
            else {
                child_transition_prob = 0.25 - 0.25*exp(-4.0 * edge_length * clock_rate * gamma_rate * relative_rate / 3.0);
            }
            return child_transition_prob;
    }

    inline double Forest::calcTransitionProbabilityLazyHKY(double s, double s_child, double edge_length, unsigned locus, double clock_rate, unsigned rate_categ) const {
        double relative_rate = G::_double_relative_rates[locus];
        assert (relative_rate > 0.0);
        assert (rate_categ < G::_gamma_rate_cat.size());
        double gamma_rate = G::_gamma_rate_cat[rate_categ];

        double child_transition_prob = 0.0;

        double pi_A = G::_base_frequencies[0];
        double pi_C = G::_base_frequencies[1];
        double pi_G = G::_base_frequencies[2];
        double pi_T = G::_base_frequencies[3];

        double pi_j = 0.0;
        double PI_J = 0.0;

        double phi = (pi_A+pi_G)*(pi_C+pi_T)+G::_kappa*(pi_A*pi_G+pi_C*pi_T);
        double beta_t = 0.5*(edge_length * relative_rate * gamma_rate * clock_rate )/phi;

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
        assert (child_transition_prob > 0.0);
        return child_transition_prob;
    }


    inline pair<bool, vector<bool>> Forest::checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets) {
        bool valid = false;
        vector<bool> valid_taxsets;
                
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
                valid = true; // TODO: what if real count is 1? in that case, the fossil has no constraint? just make that not an option?
            }
            valid_taxsets.push_back(valid);
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
        
        valid_taxsets.push_back(valid);
        
        bool at_least_one_valid = false;
        for (unsigned v=0; v<valid_taxsets.size(); v++) {
            if (valid_taxsets[v]) {
                at_least_one_valid = true;
                break;
            }
        }
        
        return make_pair(at_least_one_valid, valid_taxsets);
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
        _partial_count = other._partial_count;
        _weight_correction = other._weight_correction;
        _first_split_prior = other._first_split_prior;
        _tree_height = other._tree_height;
        _log_likelihood = other._log_likelihood;

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
            nd._position_in_lineages = othernd._position_in_lineages;
            nd._partials             = othernd._partials;
            nd._split                = othernd._split;
            nd._height               = othernd._height;
            nd._set_partials         = othernd._set_partials;
            nd._use_in_likelihood    = othernd._use_in_likelihood;
            
            // Sanity check
//            assert(nd._number >= 0);
//            assert(nd._number < _nodes.size());
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
        
        _preorder.resize(other._preorder.size());
        if (other._preorder.size() > 0) {
            unsigned m = 0;
            for (auto &othernd : other._preorder) {
                unsigned n = othernd->_number;
                Node * nd = &_nodes[n];
                _preorder[m] = nd;
                m++;
            }
        }
    }

    inline double Forest::calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length, double clock_rate) {
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
//        double kappa = 1.0;
        double kappa = G::_kappa;
        double betat = 0.5*_relrate*edge_length*clock_rate/((pi[0] + pi[2])*(pi[1] + pi[3]) + kappa*(pi[0]*pi[2] + pi[1]*pi[3]));
        
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
        
    inline void Forest::simulateData(Lot::SharedPtr lot, Data::SharedPtr data, unsigned starting_site, unsigned nsites, double clock_rate) {
        
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
        vector<double> basefreq = G::_base_frequencies;
//        vector<double> basefreq = {0.25, 0.25, 0.25, 0.25};
//        if (G::_comphet != G::_infinity) {
//            // Draw 4 Gamma(G::_comphet, 1) variates
//            double A = lot->gamma(G::_comphet, 1.0);
//            double C = lot->gamma(G::_comphet, 1.0);
//            double G = lot->gamma(G::_comphet, 1.0);
//            double T = lot->gamma(G::_comphet, 1.0);
//            double total = A + C + G + T;
//            basefreq[0] = A/total;
//            basefreq[1] = C/total;
//            basefreq[2] = G/total;
//            basefreq[3] = T/total;
//        }
        
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
                    cum_prob += calcSimTransitionProbability(from_state, to_state, basefreq, site_relrate*nd->_edge_length, clock_rate);
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
//#endif

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

    inline double Forest::getHeightFirstSplit() {
        assert (_increments.size() > 0);
        return _increments[0];
    }

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

    inline PartialStore::partial_t Forest::pullPartial() {
        lock_guard<mutex> guard(mutex);
        PartialStore::partial_t ptr;
        
        // Grab one partial from partial storage
        ptr = ps.getPartial();
        return ptr;
    }

    inline void Forest::advanceAllLineagesBy(double increment) {
        // Add dt to the edge length of all lineage root nodes, unless
        // 1. there is just one lineage or
        // 2. this is a gene forest extension and the
        //    lineage root node belongs to the parent,
        // in which case do nothing
        unsigned n = (unsigned)_lineages.size();
        if (n > 1) {
            for (auto &nd : _lineages) {
                double edge_len = nd->_edge_length + increment;
                assert(edge_len >= 0.0 || fabs(edge_len) < Node::_smallest_edge_length);
                nd->_edge_length = edge_len;
            }
        
            // Add to to the current forest height
            _tree_height += increment;
        }
    }

    inline double Forest::getLogLikelihood() {
        double log_likelihood = 0.0;
        for (auto &g:_gene_tree_log_likelihoods) {
            log_likelihood += g;
        }
        return log_likelihood;
    }

    inline void Forest::refreshAllPreorders() {
        // For each subtree stored in _lineages, create a vector of node pointers in preorder sequence
        _preorder.clear();
        if (_lineages.size() == 0) {
            return;
        }

        for (auto & nd : _lineages) {

            // lineage is a Node::ptr_vect_t (i.e. vector<Node *>)
            // lineage[0] is the first node pointer in the preorder sequence for this lineage
            // Add a new vector to _preorders containing, for now, just the root of the subtree
            _preorder.push_back({nd});

            // Now add the nodes above the root in preorder sequence
//            Node::ptr_vect_t & preorder_vector = *_preorder.rbegin();
            refreshPreorder();
        }
    }
    
}
