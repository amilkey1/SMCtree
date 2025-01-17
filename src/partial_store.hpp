#pragma once

namespace proj {
    class PartialStore {
        public:
            PartialStore();
            ~PartialStore();
            typedef std::shared_ptr < std::vector<double> >   partial_t;
        
//            typedef vector<partial_t>                         partials_t;
//            typedef std::shared_ptr <std::vector<std::vector<double> > >     partials_t;
            partial_t getPartial(unsigned nelements);
//            void setnelements(unsigned nelements) {_nelements = nelements;}
        void resetPartial();
//        private:
//            unsigned _nelements;
    };

    inline PartialStore::PartialStore() {
        //std::cout << "Constructing a partial store" << std::endl;
//        _nelements=0;
    }

    inline PartialStore::~PartialStore() {
        //std::cout << "Destroying a partial store" << std::endl;
    }

    inline PartialStore::partial_t PartialStore::getPartial(unsigned nelements) {
        assert(nelements>0);
        return partial_t(new std::vector<double> (nelements));
//        return partials_t(new vector < vector<double> (nelements) > (G::_nloci));
    }
}
