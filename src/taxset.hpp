#pragma once

namespace proj {

    struct TaxSet {
    
        TaxSet(string name, vector<string> species_names)
                : _name(name), _species_included(species_names) {}
        string         _name;
        vector<string> _species_included;
    };

}
