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
            
            string              saveNodeInfo(string prefix = "") const;
                                                    
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
            PartialStore::partial_t _partials;
            unsigned            _position_in_lineages;
            bool                _set_partials; // true if node is included in likelihood calculation
            bool                _use_in_likelihood;
        
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
        //_partials.clear();
        _partials = nullptr;
        _set_partials = true;
        _use_in_likelihood = true;
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
    
    inline string Node::saveNodeInfo(string prefix) const {
        string s = str(format("%sNode %d:\n") % prefix % _number);
        s += str(format("%s  _number               = %d\n") % prefix % _number);
        s += str(format("%s  _name                 = %s\n") % prefix % _name);
        s += str(format("%s  _edge_length          = %.9f\n") % prefix % _edge_length);
        s += str(format("%s  _height               = %.9f\n") % prefix % _height);
        s += str(format("%s  _position_in_lineages = %d\n") % prefix % _position_in_lineages);
        s += str(format("%s  _split                = %s\n") % prefix %  _split.createPatternRepresentation());
        if (_partials)
            s += str(format("%s  _partials             = %d\n") % prefix % _partials->size());
        else
            s += str(format("%s  _partials             = nullptr\n") % prefix);
        if (_left_child)
            s += str(format("%s  _left_child           = Node %d\n") % prefix % _left_child->_number);
        else
            s += str(format("%s  _left_child           = nullptr\n") % prefix);
        if (_right_sib)
            s += str(format("%s  _right_sib            = Node %d\n") % prefix % _right_sib->_number);
        else
            s += str(format("%s  _right_sib            = nullptr\n") % prefix);
        if (_parent)
            s += str(format("%s  _parent               = Node %d\n") % prefix % _parent->_number);
        else
            s += str(format("%s  _parent               = nullptr\n") % prefix);
        return s;
    }
}
