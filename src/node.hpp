#pragma once

namespace proj {

    class Likelihood;
    class Forest;
    class Particle;

    class Node {
        friend class Likelihood;
        friend class Forest;
        friend class Particle;

        public:
                                        Node();
                                        ~Node();

                    Node *              getParent()                 {return _parent;}
                    Node *              getLeftChild()              {return _left_child;}
                    Node *              getRightSib()               {return _right_sib;}
                    int                 getNumber()                 {return _number;}
                    std::string         getName()                   {return _name;}
                    Split               getSplit()                  {return _split;}

                    double              getEdgeLength()             {return _edge_length;}
                    void                setEdgeLength(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers()             {_left_child = _right_sib = _parent = 0;}
            void                resetNode();
                                        
            static const double _smallest_edge_length;

            typedef std::vector<Node>    Vector;
            typedef std::vector<Node *>  PtrVector;
        
        private:
        
            enum Flag {
                Selected   = (1 << 0),
                SelPartial = (1 << 1),
                SelTMatrix = (1 << 2),
                AltPartial = (1 << 3),
                AltTMatrix = (1 << 4)
            };

            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            std::string         _name;
            double              _edge_length;
            Split               _split;
            int                 _flags;
//            PartialStore::partial_t _partial;
            PartialStore::partials_t _partials;
            unsigned            _position_in_lineages;
    };
    
    
    inline Node::Node() {
        //std::cout << "Creating Node object" << std::endl;
        clear();
    }

    inline Node::~Node() {
        //std::cout << "Destroying Node object" << std::endl;
    }

    inline void Node::clear() {
        _flags = 0;
        clearPointers();
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
//        _partials->clear();
    }

    inline void Node::setEdgeLength(double v) {
        _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);
    }

    inline unsigned Node::countChildren() const {
        unsigned n_children = 0;
        for (Node * child = _left_child; child; child=child->_right_sib) {
            n_children ++;
        }
        return n_children;
    }

    inline void Node::resetNode() {
        _parent=0;
        _right_sib=0;
    }
}

