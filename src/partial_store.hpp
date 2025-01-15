#pragma once

namespace proj {

    struct Partial {
        Partial(unsigned g, unsigned n);
        ~Partial();
        unsigned        _g; // the gene
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
    };

    inline Partial::Partial(unsigned g, unsigned n) {
        _g = g;
        _v.resize(n);
    }
    
    inline Partial::~Partial() {
    }

    class PartialStore {
    
        public:
                                    PartialStore();
                                    ~PartialStore();
                                    
            typedef shared_ptr<Partial>          partial_t;
            
            typedef vector<partial_t>            vect_partial_t;
            typedef vector<vect_partial_t>       leaf_partials_t;
            typedef vector<vect_partial_t>       storage_t;
            
            void            setNGenes(unsigned ngenes);
            partial_t       getPartial(unsigned gene);
            void            putPartial(unsigned gene, partial_t partial);
            void            setNElements(unsigned nelements, unsigned gene);
            unsigned        getNElements(unsigned gene) const {return _nelements[gene];}
            
        private:
        
            vector<unsigned> _nelements;
            storage_t        _storage;
    };

    inline PartialStore::PartialStore() {
    }

    inline PartialStore::~PartialStore() {
        _nelements.clear();
        _storage.clear();
    }
    
    inline void PartialStore::setNGenes(unsigned ngenes) {
        // Should be called before any partials are stored
        assert(_nelements.empty());
        assert(_storage.empty());
        
        // Resize both containers
        _nelements.resize(ngenes);
        _storage.resize(ngenes);
    }

    inline void PartialStore::setNElements(unsigned nelements, unsigned gene) {
        assert(_nelements.size() > gene);
        _nelements[gene] = nelements;
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned gene) {
        // Check to make sure supplied value of gene is valid
        assert(_nelements.size() > gene);
        assert(_nelements[gene] > 0);
        
        partial_t partial;
        if (_storage[gene].empty()) {
            // No stored partials for this gene, so allocate one
            partial = partial_t(new Partial(gene, _nelements[gene]));
        }
        else {
            size_t n = _storage[gene].size();
            partial = _storage[gene].at(n-1);
            _storage[gene].pop_back();
        }
        return partial;
    }
    
    inline void PartialStore::putPartial(unsigned gene, partial_t partial) {
        // Check to make sure supplied value of gene is valid
        assert(_nelements.size() > gene);
        assert(_nelements[gene] > 0);

        // Store the partial for later
        assert(partial->_v.size() == _nelements[gene]);
        partial->_v.assign(_nelements[gene], 0.0);
        _storage[gene].push_back(partial);
    }
        
 }
