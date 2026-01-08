#pragma once

namespace proj {

    struct Partial {
        Partial(unsigned n);
        ~Partial();
        vector<double>  _v; // the partial array: length = _nstates*<no. patterns>
    };

    inline Partial::Partial(unsigned n) {
        _v.resize(n);
        _v.assign(n, 1.0); // TODO: be careful - check this works for HKY too
    }

    inline Partial::~Partial() {
    }

    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
            typedef std::shared_ptr<Partial>          partial_t;
        
        void            setNLoci(unsigned nloci);
//        partial_t       getPartial(unsigned locus);
        void            putPartial(partial_t partial);
        void            setNElements(unsigned nelements);
        unsigned        getNElements() const {return _nelements;}

        
        typedef partial_t           vect_partial_t;
        typedef vect_partial_t      storage_t;
        
        unsigned         _total_partials_created;
        unsigned         _total_elements_created;
        unsigned         _total_partials_reused;
        unsigned         _total_elements_reused;
        
        PartialStore::partial_t getPartial();
        
        void resetPartial();
        private:
            storage_t        _storage;
            unsigned _nelements;
    };

    inline PartialStore::PartialStore() {
        //std::cout << "Constructing a partial store" << std::endl;
        _total_partials_created = 0;
        _total_elements_created = 0;
        _total_partials_reused = 0;
        _total_elements_reused = 0;
    }

    inline PartialStore::~PartialStore() {
        //std::cout << "Destroying a partial store" << std::endl;
        _nelements = 0;
        _storage = nullptr;
    }


    inline void PartialStore::setNLoci(unsigned nloci) {
        // Should be called before any partials are stored
        assert(_nelements == 0);
        assert(_storage == nullptr);
    }

    inline void PartialStore::setNElements(unsigned nelements) {
        _nelements = nelements;
    }

    inline PartialStore::partial_t PartialStore::getPartial() {
        // Check to make sure supplied value of locus is valid
        assert (_nelements > 0);
        
        partial_t partial;
        if (_storage == nullptr) {
            // No stored partials for this locus, so allocate one
            partial = partial_t(new Partial(_nelements));
            _total_partials_created++;
            _total_elements_created += _nelements;
        }
        else {
            partial = _storage;
            _storage = nullptr;
            _total_partials_reused++;
            _total_elements_reused += _nelements;
        }
        
        return partial;
    }
        
    inline void PartialStore::putPartial(partial_t partial) {
        // Check to make sure supplied value of locus is valid
        assert (_nelements > 0);

        // Store the partial for later
        assert(partial->_v.size() == _nelements);
        partial->_v.assign(_nelements, 1.0);
        _storage = partial;
    }
}
