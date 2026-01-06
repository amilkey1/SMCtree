#pragma once

namespace proj {

    class ForestExtension {
        public:
            
                            ForestExtension();
                                
            void            dock(const Forest::SharedPtr gf, PartialStore::partial_t partial, Lot::SharedPtr lot);
            void            dockSim(const Forest::SharedPtr gf, Lot::SharedPtr lot);
            void            undock();

            G::uint_pair_t  chooseNodesToJoin() const;
            double          getProposedDelta() const;
            double          getHeight() const;
            double          getLogWeight() const;
            const Node *    getProposedAnc() const;
            const Node *    getProposedLChild() const;
            const Node *    getProposedRChild() const;
                        
            void            addIncrement(double increment);
            void            joinPriorPrior();
                    
            PartialStore::partial_t getExtensionPartial();

        private:

            Forest::ConstSharedPtr        _docked_gene_forest;
            double                        _log_weight;
            double                        _proposed_delta;
            Node                          _proposed_anc;
            const Node *                  _proposed_lchild;
            const Node *                  _proposed_rchild;
            Lot::SharedPtr                _lot;
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
    
    inline G::uint_pair_t ForestExtension::chooseNodesToJoin() const {
//        double nsubtrees = node_indices.size();
        double nsubtrees = _docked_gene_forest->_lineages.size();
        
        assert (nsubtrees>1);
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

    inline void ForestExtension::joinPriorPrior() {
        // Choose the two nodes to join
         G::uint_pair_t chosen = chooseNodesToJoin();
        
        // Get pointers to the two nodes to join
        _proposed_lchild = _docked_gene_forest->_lineages[chosen.first];
        _proposed_rchild = _docked_gene_forest->_lineages[chosen.second];
        
        // Set _proposed_anc's split to union of the two child splits
        _proposed_anc._split.resize(G::_ntaxa);
        _proposed_anc._split += _proposed_lchild->_split;
        _proposed_anc._split += _proposed_rchild->_split;

        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }
        else {
            // Compute partial likelihood array of ancestral node
            _log_weight = _docked_gene_forest->calcPartialArrayLazy(&_proposed_anc, _proposed_lchild, _proposed_rchild);
        }
        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }

        assert(!isnan(_log_weight));
        assert(!isinf(_log_weight));
    }
    
    inline PartialStore::partial_t ForestExtension::getExtensionPartial() {
        return _proposed_anc._partials;
    }

 }


