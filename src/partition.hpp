#pragma once    

using boost::trim;

namespace proj {

    class Partition {
        public:
            typedef match_results<string::const_iterator>::const_reference         regex_match_t;
            typedef tuple<unsigned, unsigned, unsigned, unsigned>                  subset_range_t;
            typedef vector<subset_range_t>                                         partition_t;
            typedef vector<DataType>                                               datatype_vect_t;
            typedef vector<unsigned>                                               subset_sizes_vect_t;
            typedef vector<string>                                                 subset_names_vect_t;
            typedef std::shared_ptr<Partition>                                     SharedPtr;

                                                        Partition();
                                                        ~Partition();
                                                                
            unsigned                                    getNumSites() const;
            unsigned                                    getNumSubsets() const;
            string                                      getSubsetName(unsigned subset) const;
        
            const partition_t &                         getSubsetRangeVect() const;
        
            unsigned                                    findSubsetByName(const string & subset_name) const;
            unsigned                                    findSubsetForSite(unsigned site_index) const;
            bool                                        siteInSubset(unsigned site_index, unsigned subset_index) const;
            DataType                                    getDataTypeForSubset(unsigned subset_index) const;
            const datatype_vect_t &                     getSubsetDataTypes() const;
        
            unsigned                                    numSitesInSubset(unsigned subset_index) const;
            subset_sizes_vect_t                         calcSubsetSizes() const;

            void                                        defaultPartition(unsigned nsites = numeric_limits<unsigned>::max());
            void                                        parseSubsetDefinition(string & s);
            void                                        finalize(unsigned nsites);

            void                                        clear();

        private:

            int                                         extractIntFromRegexMatch(regex_match_t s, unsigned min_value);
            void                                        addSubsetRange(unsigned subset_index, string range_definition);
            void                                        addSubset(unsigned subset_index, string subset_definition);

            unsigned                                    _num_sites;
            unsigned                                    _num_subsets;
            subset_names_vect_t                         _subset_names;
            partition_t                                 _subset_ranges;
            datatype_vect_t                             _subset_data_types;

            const unsigned                              _infinity;
    };
    
    inline Partition::Partition() : _infinity(numeric_limits<unsigned>::max()) {
        clear();
    }

    inline Partition::~Partition() {
    }

    inline unsigned Partition::getNumSites() const {
        return _num_sites;
    }
    
    inline unsigned Partition::getNumSubsets() const {
        return _num_subsets;
    }
    
    inline string Partition::getSubsetName(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_names[subset];
    }
    
    inline const Partition::partition_t & Partition::getSubsetRangeVect() const {
        return _subset_ranges;
    }
    
    inline DataType Partition::getDataTypeForSubset(unsigned subset_index) const {
        assert(subset_index < _subset_data_types.size());
        return _subset_data_types[subset_index];
    }

    inline const vector<DataType> & Partition::getSubsetDataTypes() const {
        return _subset_data_types;
    }    

    inline unsigned Partition::findSubsetByName(const string & subset_name) const {    
        auto iter = find(_subset_names.begin(), _subset_names.end(), subset_name);
        if (iter == _subset_names.end())
            throw XProj(format("Specified subset name \"%s\" not found in partition") % subset_name);
        return (unsigned)distance(_subset_names.begin(),iter);
    }    

    inline unsigned Partition::findSubsetForSite(unsigned site_index) const {    
        for (auto & t : _subset_ranges) {
            unsigned begin_site = get<0>(t);
            unsigned end_site = get<1>(t);
            unsigned stride = get<2>(t);
            unsigned site_subset = get<3>(t);
            bool inside_range = site_index >= begin_site && site_index <= end_site;
            if (inside_range && (site_index - begin_site) % stride == 0)
                return site_subset;
        }
        throw XProj(format("Site %d not found in any subset of partition") % (site_index + 1));
    }    
    
    inline bool Partition::siteInSubset(unsigned site_index, unsigned subset_index) const {    
        unsigned which_subset = findSubsetForSite(site_index);
        return (which_subset == subset_index ? true : false);
    }    
    
    inline unsigned Partition::numSitesInSubset(unsigned subset_index) const {    
        unsigned nsites = 0;
        for (auto & t : _subset_ranges) {
            unsigned begin_site = get<0>(t);
            unsigned end_site = get<1>(t);
            unsigned stride = get<2>(t);
            unsigned site_subset = get<3>(t);
            if (site_subset == subset_index) {
                unsigned n = end_site - begin_site + 1;
                nsites += (unsigned)(floor(n/stride)) + (n % stride == 0 ? 0 : 1);
            }
        }
        return nsites;
    }    
    
    inline vector<unsigned> Partition::calcSubsetSizes() const {    
        assert(_num_sites > 0); // only makes sense to call this function after subsets are defined
        vector<unsigned> nsites_vect(_num_subsets, 0);
        for (auto & t : _subset_ranges) {
            unsigned begin_site = get<0>(t);
            unsigned end_site = get<1>(t);
            unsigned stride = get<2>(t);
            unsigned site_subset = get<3>(t);
            unsigned hull = end_site - begin_site + 1;
            unsigned n = (unsigned)(floor(hull/stride)) + (hull % stride == 0 ? 0 : 1);
            nsites_vect[site_subset] += n;
        }
        return nsites_vect;
    }    
    
    inline void Partition::clear() {    
        _num_sites = 0;
        _num_subsets = 1;
        _subset_data_types.clear();
        _subset_data_types.push_back(DataType());
        _subset_names.clear();
        _subset_names.push_back("default");
        _subset_ranges.clear();
        _subset_ranges.push_back(make_tuple(1, _infinity, 1, 0));
    }    
    
    inline void Partition::parseSubsetDefinition(string & s) {    
        vector<string> v;
        
        // first separate part before colon (stored in v[0]) from part after colon (stored in v[1])
        split(v, s, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj("Expecting exactly one colon in partition subset definition");

        string before_colon = v[0];
        string subset_definition = v[1];

        // now see if before_colon contains a data type specification in square brackets
        const char * pattern_string = R"((.+?)\s*(\[(\S+?)\])*)";
        regex re(pattern_string);
        smatch match_obj;
        bool matched = regex_match(before_colon, match_obj, re);
        if (!matched) {
            throw XProj(format("Could not interpret \"%s\" as a subset label with optional data type in square brackets") % before_colon);
        }
        
        // match_obj always yields 2 strings that can be indexed using the operator[] function
        // match_obj[0] equals entire subset label/type string (e.g. "rbcL[codon:standard]")
        // match_obj[1] equals the subset label (e.g. "rbcL")
        
        // Two more elements will exist if the user has specified a data type for this partition subset
        // match_obj[2] equals data type inside square brackets (e.g. "[codon:standard]")
        // match_obj[3] equals data type only (e.g. "codon:standard")
        
        string subset_name = match_obj[1].str();
        DataType dt;    // nucleotide by default
        string datatype = "nucleotide";
        if (match_obj.size() == 4 && match_obj[3].length() > 0) {
            datatype = match_obj[3].str();
            boost::to_lower(datatype);

            // check for comma plus genetic code in case of codon
            regex re(R"(codon\s*,\s*(\S+))");
            smatch m;
            if (regex_match(datatype, m, re)) {
                dt.setCodon();
                string genetic_code_name = m[1].str();
                dt.setGeneticCodeFromName(genetic_code_name);
            }
            else if (datatype == "codon") {
                dt.setCodon();  // assumes standard genetic code
            }
            else if (datatype == "protein") {
                dt.setProtein();
                }
            else if (datatype == "nucleotide") {
                dt.setNucleotide();
                }
            else if (datatype == "standard") {
                dt.setStandard();
                }
            else {
                throw XProj(format("Datatype \"%s\" specified for subset(s) \"%s\" is invalid: must be either nucleotide, codon, protein, or standard") % datatype % subset_name);
                }
            }

        // Remove default subset if there is one
        unsigned end_site = get<1>(_subset_ranges[0]);
        if (_num_subsets == 1 && end_site == _infinity) {
            _subset_names.clear();
            _subset_data_types.clear();
            _subset_ranges.clear();
        }
        else if (subset_name == "default") {
            throw XProj("Cannot specify \"default\" partition subset after already defining other subsets");
        }
        _subset_names.push_back(subset_name);
        _subset_data_types.push_back(dt);
        _num_subsets = (unsigned)_subset_names.size();
        trim(subset_definition);
        addSubset(_num_subsets - 1, subset_definition);
    }    
    
    inline void Partition::addSubset(unsigned subset_index, string subset_definition) {    
        vector<string> parts;
        split(parts, subset_definition, boost::is_any_of(","));
        for (auto subset_component : parts) {
            addSubsetRange(subset_index, subset_component);
        }
    }    
    
    inline void Partition::addSubsetRange(unsigned subset_index, string range_definition) {    
        // match patterns like these: "1-.\3" "1-1000" "1001-."
        const char * pattern_string = R"(\s*(\d+)\s*(-\s*([0-9.]+)(\\\s*(\d+))*)*)";
        regex re(pattern_string);
        smatch match_obj;
        bool matched = regex_match(range_definition, match_obj, re);
        if (!matched) {
            throw XProj(format("Could not interpret \"%s\" as a range of site indices") % range_definition);
        }
        
        // match_obj always yields 6 strings that can be indexed using the operator[] function
        // match_obj[0] equals entire site_range (e.g. "1-.\3")
        // match_obj[1] equals beginning site index (e.g. "1")
        // match_obj[2] equals everything after beginning site index (e.g. "-.\3")
        // match_obj[3] equals "" or ending site index (e.g. ".")
        // match_obj[4] equals "" or everything after ending site index (e.g. "\3")
        // match_obj[5] equals "" or step value (e.g. "3")
        int ibegin = extractIntFromRegexMatch(match_obj[1], 1);
        int iend   = extractIntFromRegexMatch(match_obj[3], ibegin);
        int istep  = extractIntFromRegexMatch(match_obj[5], 1);
        
        // record the triplet
        _subset_ranges.push_back(make_tuple(ibegin, iend, istep, subset_index));
        
        // determine last site in subset
        unsigned last_site_in_subset = iend - ((iend - ibegin) % istep);
        if (last_site_in_subset > _num_sites) {
            _num_sites = last_site_in_subset;
        }
    }    
    
    inline int Partition::extractIntFromRegexMatch(regex_match_t s, unsigned min_value) {    
        int int_value = min_value;
        if (s.length() > 0) {
            string str_value = s.str();
            try {
                int_value = stoi(str_value);
            }
            catch(invalid_argument) {
                throw XProj(format("Could not interpret \"%s\" as a number in partition subset definition") % s.str());
            }
            
            // sanity check
            if (int_value < (int)min_value) {
                throw XProj(format("Value specified in partition subset definition (%d) is lower than minimum value (%d)") % int_value % min_value);
            }
        }
        return int_value;
    }    
    
    inline void Partition::finalize(unsigned nsites) {    
        if (_num_sites == 0) {
            defaultPartition(nsites);
            return;
        }

        // First sanity check:
        //   nsites is the number of sites read in from a data file;
        //   _num_sites is the maximum site index specified in any partition subset.
        //   These two numbers should be the same.
        if (_num_sites != nsites) {
            throw XProj(format("Number of sites specified by the partition (%d) does not match actual number of sites (%d)") % _num_sites % nsites);
        }
        
        // Second sanity check: ensure that no sites were left out of all partition subsets
        // Third sanity check: ensure that no sites were included in more than one partition subset
        vector<int> tmp(nsites, -1);   // begin with -1 for all sites
        for (auto & t : _subset_ranges) {
            unsigned begin_site  = get<0>(t);
            unsigned end_site    = get<1>(t);
            unsigned stride  = get<2>(t);
            unsigned site_subset = get<3>(t);
            for (unsigned s = begin_site; s <= end_site; s += stride) {
                if (tmp[s-1] != -1)
                    throw XProj("Some sites were included in more than one partition subset");
                else
                    tmp[s-1] = site_subset;
            }
        }
        if (find(tmp.begin(), tmp.end(), -1) != tmp.end()) {
            throw XProj("Some sites were not included in any partition subset");
        }
        tmp.clear();
    }    
    
    inline void Partition::defaultPartition(unsigned nsites) {    
        clear();
        _num_sites = nsites;
        _num_subsets = 1;
        _subset_ranges[0] = make_tuple(1, nsites, 1, 0);
    }    
    
}   
