#pragma once

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);

namespace proj {

    class Data;
    class Model;
    class QMatrix;

    class GeneticCode {

        friend class Data;
        friend class Model;
        friend class QMatrix;

        public:

            typedef map<int, int>               genetic_code_map_t;
            typedef map<char, unsigned>         amino_acid_map_t;
            typedef vector<unsigned>            amino_acid_vect_t;
            typedef vector<string>              codon_vect_t;
            typedef vector<char>                amino_acid_symbol_vect_t;
            typedef map<string, string>         genetic_code_definitions_t;
            typedef vector<string>              genetic_code_names_t;
            
        
                                                GeneticCode();
                                                GeneticCode(string name);
                                                ~GeneticCode();
        
            string                              getGeneticCodeName() const;
            void                                useGeneticCode(string name);
        
            unsigned                            getNumNonStopCodons() const;
            int                                 getStateCode(int triplet_index) const;
            char                                getAminoAcidAbbrev(unsigned aa_index) const;

            void                                copyCodons(vector<string> & codon_vect) const;
            void                                copyAminoAcids(vector<unsigned> & aa_vect) const;

            static genetic_code_names_t         getRecognizedGeneticCodeNames();
            static bool                         isRecognizedGeneticCodeName(const string & name);
            static void                         ensureGeneticCodeNameIsValid(const string & name);

        private:

            void                                buildGeneticCodeTranslators();

            string                         _genetic_code_name;
        
            genetic_code_map_t                  _genetic_code_map;
            amino_acid_map_t                    _amino_acid_map;
        
            amino_acid_vect_t                   _amino_acids;
            codon_vect_t                        _codons;
        
            const amino_acid_symbol_vect_t      _all_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
            const vector<string>      _all_codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

            static genetic_code_definitions_t   _definitions;

        public:

            typedef std::shared_ptr<GeneticCode> SharedPtr;
    };

    // member function bodies go here
    
    inline GeneticCode::GeneticCode() {
        useGeneticCode("standard");
    }

    inline GeneticCode::GeneticCode(string name) {
        useGeneticCode(name);
    }

    inline GeneticCode::~GeneticCode() {
    }

    inline string GeneticCode::getGeneticCodeName() const {
        return _genetic_code_name;
    }

    inline void GeneticCode::useGeneticCode(string name) {
        _genetic_code_name = name;
        buildGeneticCodeTranslators();
    }

    inline unsigned GeneticCode::getNumNonStopCodons() const {
        return (unsigned)_codons.size();
    }

    inline int GeneticCode::getStateCode(int triplet_index) const {
        int retval = 0;
        if (_genetic_code_map.count(triplet_index) == 0)
            throw XProj(format("GeneticCode::getStateCode failed because key %d does not exist in _genetic_code_map") % triplet_index);
        else {
            retval = _genetic_code_map.at(triplet_index);
        }
        return retval;
    }
    
    inline char GeneticCode::getAminoAcidAbbrev(unsigned aa_index) const {
        return _all_amino_acids[aa_index];
    }
    
    inline void GeneticCode::copyCodons(vector<string> & codon_vect) const {
        codon_vect.resize(_codons.size());
        copy(_codons.begin(), _codons.end(), codon_vect.begin());
    }
    
    inline void GeneticCode::copyAminoAcids(vector<unsigned> & aa_vect) const {
        aa_vect.resize(_amino_acids.size());
        copy(_amino_acids.begin(), _amino_acids.end(), aa_vect.begin());
    }

    inline void GeneticCode::buildGeneticCodeTranslators() {
        _amino_acid_map.clear();
        for (unsigned i = 0; i < 20; ++i) {
            char aa = _all_amino_acids[i];
            _amino_acid_map[aa] = i;
        }

        ensureGeneticCodeNameIsValid(_genetic_code_name);
        string gcode_aa = _definitions[_genetic_code_name];  // e.g. "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
        
        int k = 0;
        int state_code = 0;
        _codons.clear();
        _amino_acids.clear();
        _genetic_code_map.clear();
        for (char ch : gcode_aa) {
            if (ch != '*') {
                _genetic_code_map[k] = state_code++;
                _codons.push_back(_all_codons[k]);
                _amino_acids.push_back(_amino_acid_map[ch]);
            }
            ++k;
        }
    }

    inline GeneticCode::genetic_code_names_t GeneticCode::getRecognizedGeneticCodeNames() {
        genetic_code_names_t names;
        for (auto it = _definitions.begin(); it != _definitions.end(); ++it)
            names.push_back(it->first);
        return names;
    }
    
    inline bool GeneticCode::isRecognizedGeneticCodeName(const string & name) {
        string lcname = name;
        boost::to_lower(lcname);
        return (_definitions.find(lcname) != _definitions.end());
    }
   
    inline void GeneticCode::ensureGeneticCodeNameIsValid(const string & name) {
        if (!isRecognizedGeneticCodeName(name)) {
            auto valid_genetic_code_names = getRecognizedGeneticCodeNames();
            output("Recognized genetic codes:\n", 2);
            for (string name : valid_genetic_code_names) {
                output(format("  %s\n") % name, 2);
            }
            output("\n", 2);
            throw XProj(format("%s is not a recognized genetic code") % name);
        }
    }
    
}
