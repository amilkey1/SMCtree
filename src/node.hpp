#pragma once

namespace proj {

    class Likelihood;
    class Forest;
    class Forest;
    class Particle;

    class Node {
        friend class Likelihood;
        friend class Forest;
        friend class Forest;
        friend class Particle;

        public:
                                        Node();
                                        ~Node();

                    typedef vector<Node *>  ptr_vect_t;

                    Node *              getParent()                 {return _parent;}
                    Node *              getLeftChild()              {return _left_child;}
                    Node *              getRightSib()               {return _right_sib;}
                    int                 getNumber()                 {return _number;}
                    std::string         getName()                   {return _name;}
                    Split               getSplit()                  {return _split;}

                    double              getEdgeLength()             {return _edge_length;}
                    void                setEdgeLength(double v);
                    double              getHeight() const           {return _height;}
                    void                setHeight(double v);
        
                    unsigned            countChildren() const;

                    void                clearPointers();
                    void                resetNode();
                    
                    string              asString() const;
                                                            
            static const double _smallest_edge_length;
        
        private:
        
            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            std::string         _name;
            double              _edge_length;
            Split               _split;
//            PartialStore::partials_t _partials;
            PartialStore::partial_t _partials;
            unsigned            _position_in_lineages;

            int                 _my_index;
            
            // distance from node to any leaf
            double          _height;
    };
    
    
    inline Node::Node() {
        //std::cout << "Creating Node object" << std::endl;
        clear();
    }

    inline Node::~Node() {
        //std::cout << "Destroying Node object" << std::endl;
    }

    inline void Node::clear() {
        clearPointers();
        _number = -1;
        _name = "";
        _edge_length = _smallest_edge_length;
        _height = 0.0;
        _partials.clear();
    }
    
    inline void Node::clearPointers() {
        _left_child = _right_sib = _parent = nullptr;
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
    
    inline string Node::asString() const {
        string s = str(format("Node %d:\n") % _number);
        s += str(format("  _name                 = %s\n") % _name);
        s += str(format("  _edge_length          = %.9f\n") % _edge_length);
        s += str(format("  _height               = %.9f\n") % _height);
        s += str(format("  _my_index             = %d\n") % _my_index);
        s += str(format("  _position_in_lineages = %d\n") % _position_in_lineages);
        s += str(format("  _split                = %s\n") % _split.createPatternRepresentation());
        s += str(format("  _partials.size()      = %d\n") % _partials.size());
        if (_left_child)
            s += str(format("  _left_child           = Node %d\n") % _left_child->_number);
        else
            s += "  _left_child           = NULL\n";
        if (_right_sib)
            s += str(format("  _right_sib            = Node %d\n") % _right_sib->_number);
        else
            s += "  _right_sib            = NULL\n";
        if (_parent)
            s += str(format("  _parent               = Node %d\n") % _parent->_number);
        else
            s += "  _parent               = NULL\n";
        return s;
    }
}

