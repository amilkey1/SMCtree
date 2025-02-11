#pragma once

namespace proj {

    struct Fossil {
    
        Fossil(string name, double lower, double upper, double age)
                : _name(name), _lower(lower), _upper(upper), _age(age) {}
        string         _name;
        double         _lower;
        double         _upper;
        double         _age;
    };

}
