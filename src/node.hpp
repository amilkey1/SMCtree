#pragma once    

namespace proj {

    class TreeManip;
    class Likelihood;
    class Updater;
    class Forest;
    class GeneForest;
    class SpeciesForest;
    class Particle;

    class Node {

        friend class TreeManip;
        friend class Likelihood;
        friend class Updater;
        friend class Forest;
        friend class GeneForest;
        friend class SpeciesForest;
        friend class Particle;

        public:
                                        Node();
                                        ~Node();

                    typedef vector<Node *>  ptr_vect_t;
                    //  0. number of lineages
                    //  1. gene index
                    //  2. species within gene
                    //  3. vector<Node *> lineage roots for gene/species combination
                    typedef tuple<unsigned, unsigned, G::species_t, Node::ptr_vect_t>  species_tuple_t;
        
                    Node *              getParent()                 {return _parent;}
                    const Node *        getParent() const           {return _parent;}
                    void                setParent(Node * nd)        {_parent = nd;}

                    Node *              getLeftChild()              {return _left_child;}
                    const Node *        getLeftChild() const        {return _left_child;}
                    void                setLeftChild(Node * nd)      {_left_child = nd;}

                    Node *              getRightSib()               {return _right_sib;}
                    const Node *        getRightSib() const         {return _right_sib;}
                    void                setRightSib(Node * nd)      {_right_sib = nd;}

                    string              getName()                   {return _name;}
                    const string        getName() const             {return _name;}

                    int                 getNumber() const           {return _number;}
                    int                 getMyIndex() const          {return _my_index;}
                    Split               getSplit()                  {return _split;}
        
                    bool                isSelected()                {return _flags & Flag::Selected;}
                    void                select()                    {_flags |= Flag::Selected;}
                    void                deselect()                  {_flags &= ~Flag::Selected;}

                    bool                isSelPartial()              {return _flags & Flag::SelPartial;}
                    void                selectPartial()             {_flags |= Flag::SelPartial;}
                    void                deselectPartial()           {_flags &= ~Flag::SelPartial;}

                    bool                isSelTMatrix()              {return _flags & Flag::SelTMatrix;}
                    void                selectTMatrix()             {_flags |= Flag::SelTMatrix;}
                    void                deselectTMatrix()           {_flags &= ~Flag::SelTMatrix;}

                    bool                isAltPartial()              {return _flags & Flag::AltPartial;}
                    void                setAltPartial()             {_flags |= Flag::AltPartial;}
                    void                clearAltPartial()           {_flags &= ~Flag::AltPartial;}

                    bool                isAltTMatrix()              {return _flags & Flag::AltTMatrix;}
                    void                setAltTMatrix()             {_flags |= Flag::AltTMatrix;}
                    void                clearAltTMatrix()           {_flags &= ~Flag::AltTMatrix;}
                    
                    void                flipTMatrix()               {isAltTMatrix() ? clearAltTMatrix() : setAltTMatrix();}
                    void                flipPartial()               {isAltPartial() ? clearAltPartial() : setAltPartial();}

                    double              getEdgeLength() const       {return _edge_length;}
                    void                setEdgeLength(double v);
                    
                    const G::species_t &    getSpecies() const {return _species;}
                    G::species_t &          getSpecies() {return _species;}
                                        
                    void                            setSpecies(const G::species_t other);
                    void                            setSpeciesToUnion(const G::species_t left, const G::species_t right);
                                        
                    void                            revertSpecies();
                    bool                            canRevertSpecies() const;
        
                    double              getHeight() const       {return _height;}
                    void                setHeight(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers()             {_left_child = _right_sib = _parent = 0;}
                    
                    static string       taxonNameToSpeciesName(string taxon_name);

                    static void         setSpeciesMask(G::species_t & mask, unsigned nspecies);
                    static void         setSpeciesBits(G::species_t & to_species, const G::species_t & from_species, bool init_to_zero_first);
                    static void         unsetSpeciesBits(G::species_t & to_species, const G::species_t & from_species);
                    static void         setSpeciesBit(G::species_t & to_species, unsigned i, bool init_to_zero_first);
                    static void         unsetSpeciesBit(G::species_t & to_species, unsigned i);
                                        
            static const double _smallest_edge_length;

        private:
        
            enum Flag {
                Selected   = (1 << 0),
                SelPartial = (1 << 1),
                SelTMatrix = (1 << 2),
                AltPartial = (1 << 3),
                AltTMatrix = (1 << 4)
            };

            void                clear();
            
            // Adding data members? Be sure to get them copied in Forest::operator=(const Forest & other)

            Node *          _left_child;
            Node *          _right_sib;
            Node *          _parent;
            int             _number;
            int             _my_index;
            string          _name;
            double          _edge_length;
            
            // distance from node to any leaf
            double          _height;
            
            // Bitset of species (indices) compatible with this node
            G::species_t    _species;
                        
            Split           _split;
            int             _flags;
            
            PartialStore::partial_t _partial;
    };
        
    inline Node::Node() {
        clear();
    }

    inline Node::~Node() {
    }

    inline void Node::clear() {
        _flags = 0;
        clearPointers();
        _number = -1;
        //_my_index should not be cleared because clear is called by Forest::stowNode
        _name = "";
        _edge_length = _smallest_edge_length;
        _height = 0.0;
        _species = (G::species_t)0;
        _partial = nullptr;
    }

    inline void Node::setSpeciesMask(G::species_t & mask, unsigned nspecies) {
        mask = (G::species_t)0;
        for (unsigned i = 0; i < nspecies; ++i) {
            mask |= ((G::species_t)1 << i);
        }
    }
    
    inline void Node::setSpeciesBits(G::species_t & to_species, const G::species_t & from_species, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (G::species_t)0;

        // Copy bits in from_species to to_species
        to_species |= from_species;
    }
    
    inline void Node::unsetSpeciesBits(G::species_t & to_species, const G::species_t & from_species) {
        // Zero from_species' bits in to_species
        to_species &= ~from_species;
    }
    
    inline void Node::setSpeciesBit(G::species_t & to_species, unsigned i, bool init_to_zero_first) {
        if (init_to_zero_first)
            to_species = (G::species_t)0;
            
        // Set ith bit in to_species
        to_species |= ((G::species_t)1 << i);
    }
    
    inline void Node::unsetSpeciesBit(G::species_t & to_species, unsigned i) {
        // Unset ith bit in to_species
        to_species &= ~((G::species_t)1 << i);
    }
    
    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    inline void Node::setHeight(double h) {
        _height = h;
    }

    inline unsigned Node::countChildren() const {
        unsigned n_children = 0;
        for (Node * child = _left_child; child; child=child->_right_sib) {
            n_children ++;
        }
        return n_children;
    }
    
    inline void Node::setSpecies(const G::species_t spp) {
        _species = spp;
    }
    
    inline void Node::setSpeciesToUnion(const G::species_t left, const G::species_t right) {
        // Internal nodes can have species that are unions and hence have size > 1, but
        // leaf nodes should be assigned to just a single species. This function should only
        // be called for internal nodes.
        assert(_left_child);
        _species = left;
        Node::setSpeciesBits(_species, right, /*init_to_zero_first*/false);
    }
        
    inline string Node::taxonNameToSpeciesName(string tname) {
        vector<string> before_after;
        split(before_after, tname, boost::is_any_of("^"));
        if (before_after.size() != 2)
            throw XProj(format("Expecting taxon name to conform to taxon^species pattern: %s") % tname);
        return before_after[1];
    }

}
