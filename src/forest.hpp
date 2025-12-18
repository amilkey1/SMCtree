#pragma once
#include <cmath>

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
        double getLogLikelihood();
        double getHeightFirstSplit();
        double getHeightSecondIncr();
        double getHeightThirdIncr();
    
#if defined (INCREMENT_COMPARISON_TEST)
    void incrementComparisonTest();
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
        double joinPriorPrior(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils);
        void calcPartialArray(Node* new_nd);
        double calcTransitionProbability(Node* child, double s, double s_child, unsigned locus);
        
        //POL added below
        Node * pullNode();
        unsigned getNumLineages() const;
        unsigned getNumNodes() const;
        double calcSimTransitionProbability(unsigned from, unsigned to, const vector<double> & pi, double edge_length);
        //POL added above

        pair<unsigned, unsigned> chooseTaxaToJoin(double s, Lot::SharedPtr lot);
        string makePartialNewick(unsigned precision, bool use_names);
        void showForest();
        void updateNodeVector(vector<Node *> & node_vector, Node * delnode1, Node * delnode2, Node * addnode);
        void drawBirthDiff(Lot::SharedPtr lot);
        void drawTurnover(Lot::SharedPtr lot);
        void drawRootAge(Lot::SharedPtr lot, double max_fossil_age);
        void calculateLambdaAndMu();
        void drawLambda(Lot::SharedPtr lot);

        double getTreeHeight();
        double getTreeLength();
        double getTreePrior();
        double calcTopologyPrior(unsigned nlineages);
        void clearPartials();
        double getLineageHeight(Node* nd);
        void addBirthDeathTreeIncrement(Lot::SharedPtr lot);
    
#if defined (FOSSILS)
        void    addBirthDeathIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name);
        bool    checkForValidTaxonSet(vector<TaxSet> taxset, vector<TaxSet> unused_taxsets);
#endif
    
        Data::SharedPtr             _data;
        vector<Node *>              _lineages;
        vector<Node>                _nodes;
        vector<Node*>               _preorder;
        unsigned                    _first_pattern;
        unsigned                    _npatterns;
        vector<double>              _gene_tree_log_likelihoods;
        unsigned                    _ninternals;
        unsigned                    _nleaves;
        double                      _log_joining_prob;
        vector<double>              _increments;
        vector<pair<Node*, Node*>>  _node_choices;
        double                      _estimated_lambda;
        double                      _estimated_mu;
        double                      _estimated_root_age;
        double                      _estimated_birth_difference;
        double                      _turnover;
        double                      _partial_count;
        map<string, double>         _taxset_ages;
        double                      _clock_rate;
        double                      _weight_correction; // correct for taxon set constraints
        double                      _first_split_height;
        double                      _first_split_prior;
        double                      _second_incr;
        double                      _third_incr;
    
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
        _weight_correction = 0.0;
        _first_split_height = 0.0;
        _second_incr = 0.0;
        _third_incr = 0.0;
        _clock_rate = 1.0;
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
        }
        
        // lorad only works if all topologies the same - then don't include the prior on joins because it is fixed
        double rate = 0.0; // TODO: need to modify this for birth-death
        rate = getNumLineages() * birth_rate;
    }

#if defined (INCREMENT_COMPARISON_TEST)
    inline void Forest::addBirthDeathTreeIncrementTest(Lot::SharedPtr lot) {
        // birth death
        
        bool yule = false;
        
        double beast_test1 = 0.0;
        double beast_test2 = 0.0;
        double ntaxa = 5;
        double a = G::_mu / G::_lambda;
        double r = G::_lambda - G::_mu;
        double ca = 1-a;
        double height = 100.0;
        double mrh = -r * height;
        double emrh = exp(mrh);
        double rho = 1.0;
        
        double erh = exp(r * height);
        
        beast_test1 = -(ntaxa - 2) * r * ca * emrh / (emrh - 1.0) / (emrh - 1.0 + ca);
        
        beast_test2 = -(ntaxa - 2) * ca / height / (r * height +ca);
        
        double tmp = a * emrh;
        double zDeriv = tmp == 0.0 ? 0.0 : r * tmp / log(1.0 - tmp);
        double beast_test3 = -2 * zDeriv - r;
        
//        cout << beast_test1 << endl;
//        cout << beast_test2 << endl;
//        cout << beast_test3 << endl;
//
//        double beast_test4 = log(r * ca * (rho + ca / (erh - 1)));
//        cout << beast_test4 << endl;
        
        height = 40.0;
        double n = 5;
        double mu = 0.5;
        double lamda = 5.0;
        double k = 1;
        double u = lot->uniform();
        double phi = (1 - u)*exp(1/(n-k-1));
        double tprev = 10.0;
        a = exp((lamda - mu) * tprev);
        double troot = 40.0;
        double b = exp((lamda - mu) * troot);
        
        double incr_test = log(lamda * phi * (b - a) + lamda * a - 1 / b) - log(mu * phi * (b - a) + lamda * a * b - 1);
        cout << incr_test << endl;
        
        if (yule) {
                ofstream logf("smc3.log");
                logf << "sample" << "\t" << "increment" << endl;
                double cum_height = 0.0;
                unsigned n = getNumLineages();
            
                double t = 0.0;
                for (unsigned a = 0; a < 100000; a++) {
                    bool old = false;
                
                    for (unsigned b = 0; b < 3; b++) {
                        if (!old) {
                            assert (n > 1);

                            if (n > 1) {

                                // Draw n-1 internal node heights and store in vector heights
                                vector<double> heights(n - 1, 0.0);

    //                    for (unsigned i = 0; i < 1; i++) {
    //                        for (unsigned i = 0; i < n - 2; i++) {
                                
    //                            bool birth_death = false;
    //                            if (birth_death) {
    //                            double phi = 0.0;
    //
    //                            double u = lot->uniform();
    //
    //                            phi += G::_lambda - G::_mu * exp((G::_mu - G::_lambda) * (_estimated_root_age - cum_height));
    //                            phi /= (1 - exp((G::_mu - G::_lambda) * (_estimated_root_age - cum_height)));
    //
    //
    //                            double s = 0.0;
    //                            double inner_term = (u * G::_lambda - phi) / (u * G::_mu - phi);
    //                            s = -1 * log(inner_term) / (G::_lambda - G::_mu);
    //                            s += cum_height;
    //
    //                            assert (s > 0);
    //                            heights[i] = s;
    //                            }
                                
    //                            else {
                            double troot = 2.0;
                            double u = lot->uniform();
                            double k = b + 1;
                            n = 5;
                            double a = pow((1-u), (1 / (n-k-1)));
                            double phi = exp(-1 * G::_lambda * (troot - cum_height));
                            t = (-1 / G::_lambda) * log(phi + a * (1 - phi));
                        
    //                        double k = b + 1;
    //
    //                        n = 5;
    //                        double a = pow((1-u), (1 / (n-k)));
    //                        double b = exp(-1 * G::_lambda * cum_height);
    //                        double c = exp(-1 * G::_lambda * troot);
    //                        double d = -1 / G::_lambda;
    //
    //                        double s = d * log(a*b - a*c + c);
    //
    //                        assert (s > cum_height);
    //                        t = s - cum_height;
    //                            }
    //                        }

    //                        heights[n-2] = _estimated_root_age; // TODO: is this right?
    //                        sort(heights.begin(), heights.end());
    //
    //                    t = heights[0] - cum_height;
                    }
                    }
                    
    //        double cum_height = getLineageHeight(_lineages.back());
            // Determine number of lineages remaining
    //        unsigned n = getNumLineages();
    //        assert(n > 1);
    //
    //                n -= b;
            
            // Draw n-1 internal node heights and store in vector heights
                
                    if (old) {
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
                        t = heights[0]*(_estimated_root_age - cum_height);
                        
                        assert (t > 0.0);
                    }
                
                    if (b < 2) {
                        cum_height += t;
    //                    n -= 1;
                    }
                    else {
                        cum_height += t;
                        logf << a << "\t";
                        logf << cum_height << endl;
                        cum_height = 0.0;
//                        logf << a << "\t";
//                        logf << t << endl;
                        n = getNumLineages();
                    }
            }
            }
//            cout << "x";
        }
        else {
            
            ofstream logf("smc4.log");
            logf << "sample" << "\t" << "increment" << endl;
            double cum_height = 0.0;
            unsigned n = getNumLineages();
        
            double t = 0.0;
            for (unsigned a = 0; a < 100000; a++) {
            
                for (unsigned b = 0; b < 4; b++) {
                    G::_lambda = 20.0;
                    G::_mu = 10.0;
                    
                        assert (n > 1);

                        if (n > 1) {

                            // Draw n-1 internal node heights and store in vector heights
                            vector<double> heights(n - 1, 0.0);
                            
//                            double troot = 1.0;
                            double troot = 40.0;
//                            _estimated_root_age = 2.0;
                            
                            double lambda_minus_mu = G::_lambda - G::_mu;
                            
//                            double untransformed_cum_height = cum_height / _estimated_root_age;
                            
                            double x = exp(-1*lambda_minus_mu*cum_height);
//                            double x = exp(-1*lambda_minus_mu*untransformed_cum_height);
                            double z = exp(-1*lambda_minus_mu*troot);
                            
                            double log_x = -1 * lambda_minus_mu * cum_height;
                            double log_z = -1 * lambda_minus_mu * troot;

                            double u = lot->uniform();
                            
                            double k = b + 1;
                            n = getNumLineages();
                            
                            double phi = pow((1-u), (1/(n-k-1)));
                            double C = phi * (x-z) / (G::_lambda - G::_mu * x);
                            
//                            double C = pow((1-u), (1/(n-k-1))) * (x - z) / (G::_lambda - G::_mu * x);
                            
                            double numerator = C * G::_lambda + z;
                            double denominator = C * G::_mu + 1;
                            
                            double numerator1 = phi * (x-z) / (G::_lambda - G::_mu * x) * G::_lambda;
                            double numerator2 = z;
                            
//                            double new_height = log(numerator / denominator);
                            
                            // caculate each term
                            double log_phi = log(1 - u) / (n - k - 1);
                            
                            double log_lambda = log(G::_lambda);
                            double log_mu = log(G::_mu);
                            
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
                            
//                            double numerator_test = G::calcLogSum(num_values);
                            
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
                                
//                            double denominator_test = G::calcLogSum(denom_values);
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
                            
                            double new_height_test = numerator_test - denominator_test;
                            new_height_test /= (G::_mu - G::_lambda);
                            
                            cout << new_height_test << endl;
                            
                            double test1 = log(G::_lambda * phi * x - G::_lambda * phi * z + G::_lambda * z  - G::_mu * x * z);
                            double test2 = log(G::_mu * phi * x - G::_mu * phi * z + G::_lambda - G::_mu * x);
                            
                            double new_height = log(G::_lambda * phi * x - G::_lambda * phi * z + G::_lambda * z - G::_mu * x * z) - log(G::_mu * phi * x - G::_mu * phi * z + G::_lambda - G::_mu * x);
                            
                            new_height /= (G::_mu - G::_lambda);
                            
                            assert (new_height > 0.0);
                            
                            t = new_height - cum_height;
                            
//                            t = new_height - untransformed_cum_height;
//
//                            t *= _estimated_root_age;
                            
                            assert (t > 0.0);
                }
                           
                if (b < 3) {
                    cum_height += t;
                }
                else {
                    cum_height += t;
                    logf << a << "\t";
                    logf << cum_height << endl;
                    cum_height = 0.0;
//                    logf << a << "\t";
//                    logf << t << endl;
                    n = getNumLineages();
                }
        }
        }
        cout << "x";
        }
    
    }
#endif

#if defined (FOSSILS)
    void Forest::addBirthDeathIncrementFossil(Lot::SharedPtr lot, double age, string fossil_name) {
        // fossils are age constraints and not part of the tree
        
        bool fossil_added = false;
        
        // if there is only one lineage left, all extant taxa have been joined and remaining fossils must be added
        
        // birth death
        double cum_height = _tree_height;
        unsigned n = G::_ntaxa;
        
        double birth_rate = _estimated_lambda;
        
        double death_rate = _estimated_mu;
        
        assert (birth_rate >= death_rate);
        
        double t = 0.0;
        
        assert (n > 1);
        
        double troot = _estimated_root_age;
                
        double u = lot->uniform();
            
        double lambda_minus_mu = G::_lambda - G::_mu;
        
        double b = G::_step;
        double k = b + 1;
        
        double x = exp(-1*lambda_minus_mu*cum_height);
        double z = exp(-1*lambda_minus_mu*troot);
        
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
        new_height /= (G::_mu - G::_lambda);

        assert (new_height > 0);
        assert (new_height <= troot + 0.01);
        assert (new_height <= _estimated_root_age + 0.01);
        
        t = new_height - cum_height;
        
# if defined (INCREMENT_COMPARISON_TEST)
        if (b == 1) {
            _second_incr = t;
        }
        else if (b == 2) {
            _third_incr = t;
        }
#endif
        
        for (auto &nd:_lineages) {
            nd->_edge_length += t;
        }
        
        _tree_height += t;
        
        _increments.push_back(t);
                
        assert (!fossil_added);
    }
#endif

    inline void Forest::addIncrement(Lot::SharedPtr lot) {
#if defined (INCREMENT_COMPARISON_TEST)
        addBirthDeathTreeIncrementTest(lot);
#else
        addBirthDeathTreeIncrement(lot);
#endif
    }

    inline double Forest::joinPriorPrior(double prev_log_likelihood, Lot::SharedPtr lot, vector<TaxSet> &taxset, vector<TaxSet> &unused_taxset, vector<TaxSet> &taxset_no_fossils, vector<TaxSet> &unused_taxset_no_fossils, vector<Fossil> &particle_fossils) {
        _weight_correction = 0.0;
        bool filter = false;
        
        if (G::_ruv || G::_start_mode == "sim") {
            if (_lineages.size() == G::_ntaxa) {
                _first_split_height = _lineages[0]->_edge_length;
            }
        }
        // find the new_nd from the previous step and accumulate height if needed
        
        // if lineages.back() is a fossil, go to the previous node
        // TODO: don't need to do any of this if not actually including fossils in the tree
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
        
        bool fossil_constraint = false;
        string match_string = "FOSSIL";
        double fossil_age = -1;
        bool fossil_age_is_violated = false;
        
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
                if (_valid_taxsets[t]) {
                    node_choices.resize(node_choices.size() + 1);
                    unsigned real_count = 0;
                    for (auto &n:taxset_no_fossils[t]._species_included) {
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
            for (auto &t:unused_taxset_no_fossils) {
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
            
            double weight_posterior = 0;
            
            for (auto &s:set_probabilities) { //s * s  because you have ex (3/4) probability of choosing a taxset * 3 taxa within that taxset
                weight_posterior += (s * s / total_choices);
                s /= total_choices;
            }
            
            double nlineages = (double) _lineages.size();
            double weight_prior = 1 / ((nlineages * (nlineages - 1)) / 2);
            weight_posterior = 1 / weight_posterior;
            
            _weight_correction = weight_prior / weight_posterior;
            
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
            
            // node choice index corresponds to the taxset chosen
            unsigned size1 = (unsigned) taxset.size();
            if (node_choice_index < size1) { // only existing taxon sets will have associated fossils; there is no branch length constraint if the node chosen is not in a taxon set
                vector<string> species_included = taxset[node_choice_index]._species_included;
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
                        // fossil name is s
                        //
                        break;
                    }
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
            // TODO: check for taxsets here
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

            if (new_nd->_set_partials) {
                // calculate new partials
                assert (new_nd->_partials == nullptr);
                if (G::_start_mode != "sim") {
                    double npatterns_total = _data->getNumPatterns();
                    new_nd->_partials = ps.getPartial(G::_nstates*npatterns_total);
                }
                assert(new_nd->_left_child->_right_sib);

                if (G::_save_memory && G::_start_mode != "sim") {
                    double npatterns_total = _data->getNumPatterns();
                    new_nd->_partials = ps.getPartial(npatterns_total*G::_nstates);
                    
                    for (auto &nd:_lineages) {
                        if (nd->_partials == nullptr) {
                            nd->_partials = ps.getPartial(npatterns_total * G::_nstates);
                            calcPartialArray(nd);
                        }
                    }
                }
                
                if (G::_start_mode != "sim") {
                    calcPartialArray(new_nd);
                }
                
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
            
            if (new_nd->_set_partials && G::_start_mode != "sim") {
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
            // go through all taxsets in existence and look for chosen taxa, then update all of them
            // because there might be multiple taxsets with the same taxon names
            
            string name1 = subtree1->_name;
            string name2 = subtree2->_name;
            
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
                    taxset_no_fossils[i]._species_included.erase(remove(taxset_no_fossils[i]._species_included.begin(), taxset_no_fossils[i]._species_included.end(), subtree1->_name));
                    taxset_no_fossils[i]._species_included.erase(remove(taxset_no_fossils[i]._species_included.begin(), taxset_no_fossils[i]._species_included.end(), subtree2->_name));
                    taxset_no_fossils[i]._species_included.push_back(new_nd->_name);
                    
                    if (taxset_no_fossils[i]._species_included.size() == 1) {
                        update_unused = true;
                        string name = taxset_no_fossils[i]._name;
                        _taxset_ages[name] = _tree_height;
                        taxset_no_fossils.erase(taxset_no_fossils.begin() + i);
                    }
                }
            }
            
            for (unsigned i=0; i<update_these_unused_taxsets.size(); i++) {
                if (update_these_unused_taxsets[i] == true) {
                    unused_taxset_no_fossils[i]._species_included.erase(remove(unused_taxset_no_fossils[i]._species_included.begin(), unused_taxset_no_fossils[i]._species_included.end(), subtree1->_name));
                    unused_taxset_no_fossils[i]._species_included.erase(remove(unused_taxset_no_fossils[i]._species_included.begin(), unused_taxset_no_fossils[i]._species_included.end(), subtree2->_name));
                    unused_taxset_no_fossils[i]._species_included.push_back(new_nd->_name);
                }
            }
            
            // if a taxset is down to one taxon, find corresponding unused taxset and replace
            if (update_unused) {
                vector<unsigned> updateable_unused;
                vector<unsigned> updateable_unused_sizes;
                string new_name = new_nd->_name;
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
                    taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), subtree1->_name));
                    taxset[i]._species_included.erase(remove(taxset[i]._species_included.begin(), taxset[i]._species_included.end(), subtree2->_name));
                    taxset[i]._species_included.push_back(new_nd->_name);
                    
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
                        _taxset_ages[name] = _tree_height;
                        taxset.erase(taxset.begin() + i);
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
            
            // if a taxset is down to one taxon, find corresponding unused taxset and replace
            if (update_unused) {
                vector<unsigned> updateable_unused;
                vector<unsigned> updateable_unused_sizes;
                string new_name = new_nd->_name;
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

            for (unsigned index = 0; index<G::_nloci; index++) {
                calcSubsetLogLikelihood(index);
            }
            
            double new_log_likelihood = 0.0;
            for (auto &g:_gene_tree_log_likelihoods) {
                new_log_likelihood += g;
            }
            
            double log_weight = new_log_likelihood - prev_log_likelihood + _weight_correction;
                   
           if (G::_save_memory) {
               for (auto &nd:_nodes) {
                   nd._partials = nullptr;
               }
           }
            
            calcTopologyPrior(getNumLineages() +1);
            
            _valid_taxsets.clear();
        
        if (fossil_constraint) {
            assert (fossil_age != -1);
            // figure out if associated branch length has violated the fossil constraint
            // node must be at least as deep as the fossil, which sets the minimum age for the group
            if (getLineageHeight(_lineages.back()->_left_child) < fossil_age) {
                fossil_age_is_violated = true;
            }
        }
        
        if (fossil_age_is_violated) {
            return -1 * G::_infinity; // fossil constraint has been violated and particle weight is 0
        }
        else {
            if (G::_start_mode == "sim") {
                return 1.0;
            }
            else {
                return log_weight;
            }
        }
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
        
        // https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon)/10%3A_Introduction_to_Birth-Death_Models/10.02%3A_The_Birth-Death_Model
        // Raup 1985
        
        // starting from the root, go increment by increment and calculate the following parameters:
        // TODO: do fossil increments count? - no fossils for now
        
        double E = _estimated_mu / _estimated_lambda;
        double r = _estimated_lambda - _estimated_mu;
        double t = 0.0;
        unsigned count = (unsigned) _increments.size() - 1;
        
        unsigned nbranches = (unsigned) _increments.size();
//        unsigned nbranches = ntips - 1 + (unsigned) G::_fossils.size() * 2; // TODO: not sure - but I think each fossil is associated with an extra branch length - or should it be combined with existing branch length?
        
        for (unsigned i = 1; i < nbranches + 1; i++) {
            // i is number of species in existence
            // i = 1, 2, 3, ... nbranches
//        for (unsigned i=1; i<G::_fossils.size() + G::_ntaxa + 1; i++) {
            // i is number of species in existence
            // density increment is the probability of having exactly i species at time t
            
            // TODO: is it okay to calculate the prior root forwards when the tree was drawn tips backwards?
            
            // if starting at lineage at the root, prob of 1 lineage at time 0 = 1, so no need to include this
            
//            if (count == (unsigned) _increments.size()) {
//                t = 0;
//                // start by calculating the prob of 1 lineage at time 0
//            }
//            else {
                t += _increments[count]; // TODO: = or += ? - I think t is height not increment, so +=
//            }
            double a = (E*(exp(r)*t - 1)) / (exp(r)*t - E);
            double B = a / E;
            
            double density_increment = 0.0;
            if (birth_death_prior == 0) {
                density_increment = a; // starting with a single lineage
                birth_death_prior = log(density_increment);
            }
            else {
                density_increment = (1-a)*(1-B)*pow(B, i-1);
                birth_death_prior += log(density_increment);
            }
            
//            double log_prob_density_increment = log((1-a)*(1-B)*pow(B, i-1));
            
//            birth_death_prior += log_prob_density_increment;
            if (count == 0) {
                _first_split_prior = density_increment;
            }
            count--;
        }
        
//        birth_death_prior = log(birth_death_prior);
        
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
//# if defined (FOSSILS)
//                                child_transition_prob = calcTransitionProbabilityFossil(child, s, s_child, i);
//#else
                                child_transition_prob = calcTransitionProbability(child, s, s_child, i);
//#endif
                                double child_partial = child_partial_array[p*G::_nstates + s_child];
                                sum_over_child_states += child_transition_prob * child_partial;
                            }   // child state loop
//                            if (!skip) {
                            if (a == 0) {
                                parent_partial_array[p*G::_nstates+s] = sum_over_child_states;
                            }
                            else {
                                assert (a == 1);
                                parent_partial_array[p*G::_nstates+s] *= sum_over_child_states;
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
            double expterm = exp(-4.0*(child->_edge_length * _clock_rate * relative_rate)/3.0);
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
            double beta_t = 0.5*(child->_edge_length * _clock_rate * relative_rate)/phi;

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
                valid = true; // TODO: what if real count is 1? in that case, the fossil has no constraint? just make that not an option?
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
        _taxset_ages = other._taxset_ages;
        _clock_rate = other._clock_rate;
        _weight_correction = other._weight_correction;
        _first_split_height = other._first_split_height;
        _first_split_prior = other._first_split_prior;
        _second_incr = other._second_incr;
        _third_incr = other._third_incr;
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
            nd._position_in_lineages = othernd._position_in_lineages;
            nd._partials             = othernd._partials;
            nd._split                = othernd._split;
            nd._height               = othernd._height;
            nd._set_partials         = othernd._set_partials;
            nd._use_in_likelihood    = othernd._use_in_likelihood;
            
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

#if defined (INCREMENT_COMPARISON_TEST)
    inline void Forest::incrementComparisonTest() {
        // Algorithm from Yang and Rannala. 1997. MBE 14(7):717-724.
        createTrivialForest();
#if 0
        ofstream logf("sim3.log");
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
            logf << "\t" << heights[2] - heights[1] << endl;
            }
        }
        }
#else
        ofstream logf("sim3.log");
        logf << "sample" << "\t" << "increment" << endl;
        // Draw n-1 internal node heights and store in vector heights
        for (unsigned a =0; a < 1; a++) {
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
                
//                double troot = 1.0;
                double troot = 2.0;
                
                bool yule = false;
                
                if (yule) {
                    double y = (-1 / G::_lambda) * log(1 - u * (1 - exp(-1 * G::_lambda * troot)));
                }
                
                else {
                    double y = u/(1.0 + birth_rate*rho*(1.0 - u));
                    if (birth_rate > death_rate) {
                        y = log(phi - u*rho*birth_rate);
                        y -= log(phi - u*rho*birth_rate + u*(birth_rate - death_rate));
                        y /= (death_rate - birth_rate);
                    }
                    heights[i] = y;
                }
            }
            heights[n-2] = 1.0;
            sort(heights.begin(), heights.end());
            
            logf << a << "\t" << heights[2] - heights[1] << endl;
        
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
        cout << "x";
#endif
    }
#endif

    inline double Forest::getHeightFirstSplit() {
        return _first_split_height;
    }

    inline double Forest::getHeightSecondIncr() {
        return _second_incr;
    }

    inline double Forest::getHeightThirdIncr() {
        return _third_incr;
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

    inline void Forest::drawRootAge(Lot::SharedPtr lot, double max_fossil_age) {
        bool done = false;
        while (!done) {
            // Gamma(1, n) = Exp(1/n)
            // mean = n
            // for now, n = G::_root_age set by user
            _estimated_root_age = lot->gamma(1, G::_root_age);
#if defined (FOSSILS)
            if (_estimated_root_age > max_fossil_age) {
                done = true;
            }
#else
            done = true;
#endif
        }
    }

    inline double Forest::getLogLikelihood() {
        double log_likelihood = 0.0;
        for (auto &g:_gene_tree_log_likelihoods) {
            log_likelihood += g;
        }
        return log_likelihood;
    }
    
}
