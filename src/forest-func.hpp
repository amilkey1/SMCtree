#pragma once

namespace proj {

// put some forest functions in here to avoid circular include issues
    
#if defined(LAZY_COPYING)
    inline void Forest::addIncrAndJoin(double incr, const Split & lsplit, const Split & rsplit, ForestExtension & gfx) {
        // Identify the two nodes to join
        _increments.push_back(incr);
        
        Node * first_node = nullptr;
        Node * second_node = nullptr;
        for (auto nd : _lineages) {
            const Split & s = nd->_split;
            if (s == lsplit) {
                first_node = nd;
            }
            else if (s == rsplit) {
                second_node = nd;
            }
            if (first_node && second_node)
                break;
        }
        
        assert(first_node);
        assert(second_node);
        
        // Increment all lineages
        advanceAllLineagesBy(incr);
        
        // get next available node to serve as ancestral node
        Node * anc_node = &_nodes[G::_ntaxa + _ninternals];

        assert (anc_node->_parent==0);
        assert (anc_node->_number == -1);
        assert (anc_node->_right_sib == 0);
        
        anc_node->_number=G::_ntaxa+_ninternals;
        anc_node->_edge_length=0.0;
         _ninternals++;
        anc_node->_name += "node-" + to_string(anc_node->_number);

        anc_node->_left_child=first_node;
        first_node->_right_sib=second_node;
        
        first_node->_parent = anc_node;
        second_node->_parent = anc_node;
        
        anc_node->_height = _tree_height;

        // Set anc_node split to union of the two child splits
        anc_node->_split.resize(G::_ntaxa);
        anc_node->_split += first_node->_split;
        anc_node->_split += second_node->_split;
        
        // Set partial to supplied (already-calculated) partial
        anc_node->_partials = gfx.getExtensionPartial();

        // Fix up _lineages
        updateNodeVector(_lineages, first_node, second_node, anc_node);
        refreshAllPreorders();
    }
#endif

}

