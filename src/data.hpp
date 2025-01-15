#pragma once

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);

namespace proj {

    class Proj;

    class Data {
        
        friend class Proj;
        friend class Forest;
        friend class GeneForest;
        friend class Particle;
        
        public:
            typedef vector<string>                        taxon_names_t;
            typedef unsigned long long                    state_t;
            typedef vector<state_t>                       pattern_vect_t;
            typedef vector<state_t>                       monomorphic_vect_t;
            typedef vector<int>                           partition_key_t;
            typedef map<pattern_vect_t,unsigned>          pattern_map_t;
            typedef vector<pattern_vect_t>                data_matrix_t;
            typedef vector<pattern_map_t>                 pattern_map_vect_t;
            typedef vector<double>                        pattern_counts_t;
            typedef vector<unsigned>                      subset_end_t;
            typedef vector<unsigned>                      npatterns_vect_t;
            typedef map<pattern_vect_t,vector<unsigned> > pattern_origsite_map_t;
            typedef map<unsigned,vector<unsigned> >       orig_site_lookup_t;
            typedef pair<unsigned, unsigned>              begin_end_pair_t;
            typedef std::shared_ptr<Data>                 SharedPtr;

                                                        Data();
                                                        ~Data();
                                                        
            Partition::SharedPtr                        getPartition();
            void                                        setPartition(Partition::SharedPtr partition);

            void                                        getDataFromFile(const string filename);
            void                                        writeDataToFile(const string filename);

            unsigned                                    getNumSubsets() const;
            string                                      getSubsetName(unsigned subset) const;

            unsigned                                    getNumTaxa() const;

            unsigned                                    setTaxonNames(const vector<string> & names);
            const taxon_names_t &                       getTaxonNames() const;

            void                                        copyTaxonNames(taxon_names_t & dest) const;

            unsigned                                    getNumPatterns() const;
            npatterns_vect_t                            calcNumPatternsVect() const;
            unsigned                                    getNumPatternsInSubset(unsigned subset) const;
            unsigned                                    getNumStatesForSubset(unsigned subset) const;
            unsigned                                    getTotalNumberOfSites() const;
            unsigned                                    calcSeqLen() const;
            unsigned                                    calcSeqLenInSubset(unsigned subset) const;
            const data_matrix_t &                       getDataMatrix() const;
            begin_end_pair_t                            getSubsetBeginEnd(unsigned subset) const;
            const pattern_counts_t &                    getPatternCounts() const;
            const monomorphic_vect_t &                  getMonomorphic() const;
            const partition_key_t &                     getPartitionKey() const;
            
            char                                        stateCodeToDNALetter(state_t s) const;
            string                                      getPatternAsString(unsigned i) const;

            void                                        clear();
            
            void                                        memoryReport(ofstream & memf) const;

        private:

            data_matrix_t &                             getDataMatrixNonConst();
            unsigned                                    storeTaxonNames(NxsTaxaBlock * taxaBlock, unsigned taxa_block_index);
            unsigned                                    storeData(unsigned ntax, unsigned nchar_before, NxsCharactersBlock * charBlock, NxsCharactersBlock::DataTypesEnum datatype);
            unsigned                                    buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets);
            void                                        updatePatternMap(Data::pattern_vect_t & pattern, unsigned subset);
            void                                        compressPatterns();

            Partition::SharedPtr                        _partition;
            pattern_counts_t                            _pattern_counts;
            monomorphic_vect_t                          _monomorphic;
            partition_key_t                             _partition_key;
            pattern_map_vect_t                          _pattern_map_vect;
            taxon_names_t                               _taxon_names;
            data_matrix_t                               _data_matrix;
            subset_end_t                                _subset_end;
            pattern_origsite_map_t                      _orig_site_map;
            orig_site_lookup_t                          _orig_site_lookup;
            unsigned                                    _cum_nchar;
    };

    inline Data::Data() {
        clear();
    }

    inline Data::~Data() {
    }
    
    inline void Data::setPartition(Partition::SharedPtr partition) {
        _partition = partition;
    }    

    inline Partition::SharedPtr Data::getPartition() {    
        return _partition;
    }

    inline unsigned Data::getNumSubsets() const {
        return (_partition ? _partition->getNumSubsets() : 1);
    }
    
    inline string Data::getSubsetName(unsigned subset) const {
        return _partition ? _partition->getSubsetName(subset) : string("default");
    }    

    inline const Data::partition_key_t & Data::getPartitionKey() const {    
        return _partition_key;
    }
    
    inline const Data::pattern_counts_t & Data::getPatternCounts() const {
        return _pattern_counts;
    }
    
    inline const Data::monomorphic_vect_t & Data::getMonomorphic() const {
        return _monomorphic;
    }

    inline const Data::taxon_names_t & Data::getTaxonNames() const {
        return _taxon_names;
    }

    inline void Data::copyTaxonNames(Data::taxon_names_t & dest) const {
        dest.resize(_taxon_names.size());
        for (unsigned i = 0; i < _taxon_names.size(); ++i) {
            dest[i] = _taxon_names[i];
            boost::replace_all(dest[i], " ", "_");
        }
    }

    inline const Data::data_matrix_t & Data::getDataMatrix() const {
        return _data_matrix;
    }

    inline Data::data_matrix_t & Data::getDataMatrixNonConst() {
        return _data_matrix;
    }

    inline Data::begin_end_pair_t Data::getSubsetBeginEnd(unsigned subset) const {
        assert(_subset_end.size() > subset);
        if (subset == 0)
            return make_pair(0, _subset_end[0]);
        else
            return make_pair(_subset_end[subset-1], _subset_end[subset]);
    }    

    inline void Data::clear() {    
        _partition_key.clear();
        _pattern_counts.clear();
        _monomorphic.clear();
        _pattern_map_vect.clear();
        _taxon_names.clear();
        _data_matrix.clear();
        _subset_end.clear();
    }    

    inline unsigned Data::getNumPatterns() const {    
        if (_data_matrix.size() > 0)
            return (unsigned)_data_matrix[0].size();
        else
            return 0;
    }    

    inline Data::npatterns_vect_t Data::calcNumPatternsVect() const {    
        unsigned nsubsets = (unsigned)_subset_end.size();
        vector<unsigned> num_patterns_vect(nsubsets, 0);
        for (unsigned s = 0; s < nsubsets; s++)
            num_patterns_vect[s] = getNumPatternsInSubset(s);
        return num_patterns_vect;
    }    
    
    inline unsigned Data::getNumStatesForSubset(unsigned subset) const {    
        DataType data_type = _partition->getDataTypeForSubset(subset);
        return data_type.getNumStates();
    }    

    inline unsigned Data::getNumPatternsInSubset(unsigned subset) const {    
        assert(_subset_end.size() > subset);
        return (unsigned)_subset_end[subset] - (subset == 0 ? 0 : _subset_end[subset-1]);
    }    
    
    inline unsigned Data::getNumTaxa() const {    
        return (unsigned)_taxon_names.size();
    }    

    inline unsigned Data::getTotalNumberOfSites() const {
        return _cum_nchar;
    }

    inline unsigned Data::calcSeqLen() const {
        unsigned total_sites = accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
        assert(total_sites == _cum_nchar);
        return total_sites;
    }

    inline unsigned Data::calcSeqLenInSubset(unsigned subset) const {
        begin_end_pair_t s = getSubsetBeginEnd(subset);
        return accumulate(_pattern_counts.begin() + s.first, _pattern_counts.begin() + s.second, 0);
    }    
    
    inline unsigned Data::buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets) {    
        pattern_vect_t pattern(ntaxa);

        _orig_site_map.clear();
        _pattern_map_vect.clear();
        _pattern_map_vect.resize(nsubsets);
        
        const Partition::partition_t & tuples = _partition->getSubsetRangeVect();
        for (auto & t : tuples) {
            unsigned site_begin  = get<0>(t);
            unsigned site_end    = get<1>(t);
            unsigned site_skip   = get<2>(t);
            unsigned site_subset = get<3>(t);
            for (unsigned site = site_begin; site <= site_end; site += site_skip) {
                // Copy site into pattern
                for (unsigned taxon = 0; taxon < ntaxa; ++taxon) {
                    pattern[taxon] = _data_matrix[taxon][site-1];
                }
                
                // Add this pattern to _pattern_map_vect element corresponding to subset site_subset
                updatePatternMap(pattern, site_subset);
                
                // Update the map that stores original site numbers for each pattern
                auto it0 = _orig_site_map.find(pattern);
                if (it0 == _orig_site_map.end()) {
                    _orig_site_map[pattern] = {site};
                }
                else {
                    it0->second.push_back(site);
                }
            }
        }
        
        // Tally total number of patterns across all subsets
        unsigned npatterns = 0;
        for (auto & map : _pattern_map_vect) {
            npatterns += (unsigned)map.size();
        }
        
        return npatterns;
    }    

    inline void Data::updatePatternMap(Data::pattern_vect_t & pattern, unsigned subset) {
        // If pattern is not already in pattern_map, insert it and set value to 1.
        // If it does exist, increment its current value.
        // (see item 24, p. 110, in Meyers' Efficient STL for more info on the technique used here)
        pattern_map_t::iterator lowb = _pattern_map_vect[subset].lower_bound(pattern);
        if (lowb != _pattern_map_vect[subset].end() && !(_pattern_map_vect[subset].key_comp()(pattern, lowb->first))) {
            // this pattern has already been seen
            lowb->second += 1;
        }
        else
            {
            // this pattern has not yet been seen
            _pattern_map_vect[subset].insert(lowb, pattern_map_t::value_type(pattern, 1));
        }
    }

    inline void Data::compressPatterns() {
        // Perform sanity checks
        if (_data_matrix.empty())
            throw XProj("Attempted to compress an empty data matrix");
        
        unsigned ntaxa = (unsigned)_data_matrix.size();
        unsigned seqlen = (unsigned)_data_matrix[0].size();
        
        // Finalize partition
        unsigned nsubsets = getNumSubsets();
        _subset_end.resize(nsubsets);
        _partition->finalize(seqlen);

        // Compact the data, storing it in _pattern_map_vect
        unsigned npatterns = buildSubsetSpecificMaps(ntaxa, seqlen, nsubsets);
        _pattern_counts.assign(npatterns, 0);
        _monomorphic.assign(npatterns, 0);
        _partition_key.assign(npatterns, -1);

        // Rebuild _data_matrix to hold compact data, storing counts in _pattern_counts
        _data_matrix.resize(ntaxa);
        for (auto & row : _data_matrix) {
            row.resize(npatterns);
        }
        
        _orig_site_lookup.clear();

        unsigned p = 0; 
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            for (auto & pc : _pattern_map_vect[subset]) {
                _pattern_counts[p] = pc.second; // record how many sites have pattern p
                _partition_key[p] = subset;     // record the subset to which pattern p belongs
                
                // Make it possible to lookup all original site numbers
                // corresponding to a particular pattern index
                _orig_site_lookup[p] = _orig_site_map[pc.first];
                
                state_t constant_state = pc.first[0];
                unsigned t = 0;
                for (auto sc : pc.first) {
                    assert(sc > 0);
                    constant_state &= sc;
                    _data_matrix[t][p] = sc;
                    ++t;
                }
                // constant_state equals 0 if polymorphic or state code of state present if monomorphic
                _monomorphic[p] = constant_state;
                ++p;
            }   
            
            _subset_end[subset] = p;

            // Everything for this subset has been transferred to _data_matrix and _pattern_counts,
            // so we can now free this memory
            _pattern_map_vect[subset].clear();
        }
        
        // No longer need this map
        _orig_site_map.clear();
    }    

    inline unsigned Data::setTaxonNames(const vector<string> & names) {
        _taxon_names.resize(names.size());
        copy(names.begin(), names.end(), _taxon_names.begin());
        unsigned ntax = (unsigned)_taxon_names.size();
        _data_matrix.resize(ntax);
        return ntax;
    }
    
    inline unsigned Data::storeTaxonNames(NxsTaxaBlock * taxaBlock, unsigned taxa_block_index) {    
        unsigned ntax = 0;
        if (taxa_block_index == 0) {
            // First taxa block encountered in the file
            _taxon_names.clear();
            for (auto s : taxaBlock->GetAllLabels())
                _taxon_names.push_back(s);
            ntax = (unsigned)_taxon_names.size();
            _data_matrix.resize(ntax);
        }
        else {
            // Second (or later) taxa block encountered in the file
            // Check to ensure taxa block is identical to the first one
            for (auto s : taxaBlock->GetAllLabels()) {
                if (_taxon_names[ntax++] != s)
                    throw XProj(format("Taxa block %d in data file is not identical to first taxa block read") % (taxa_block_index+1));
            }
        }
        
        return ntax;
    }
    
    inline unsigned Data::storeData(unsigned ntax, unsigned nchar_before, NxsCharactersBlock * charBlock, NxsCharactersBlock::DataTypesEnum datatype) {
        unsigned seqlen = 0;
        
        // Find the data type for the partition subset containing the first site in this NxsCharactersBlock
        // Assumes that all sites in any given NxsCharactersBlock have the same type (i.e. mixed not allowed)
        assert(_partition);
        unsigned subset_index = _partition->findSubsetForSite(nchar_before + 1); // remember that sites begin at 1, not 0, in partition definitions
        DataType dt = _partition->getDataTypeForSubset(subset_index);

        // Determine number of states and bail out if data type not handled
        // 1 = standard, 2 = dna, 3 = rna, 4 = nucleotide, 5 = protein, 6 = continuous, 7 = codon, 8 = mixed
        NxsCharactersBlock * block = charBlock;
        if (datatype == NxsCharactersBlock::dna || datatype == NxsCharactersBlock::rna || datatype == NxsCharactersBlock::nucleotide) {
            if (dt.isCodon()) {
                // Create a NxsCharactersBlock containing codons rather than nucleotides
                block = NxsCharactersBlock::NewCodonsCharactersBlock(
                    charBlock,
                    true,   // map partial ambiguities to completely missing (note: false is not yet implemented in NCL)
                    true,   // gaps to missing
                    true,   // inactive characters treated as missing
                    NULL,   // if non-NULL, specifies the indices of the positions in the gene
                    NULL);  // if non-NULL, specifies a pointer to a NxsCharactersBlock that contains all non-coding positions in gene
            }
            else {
                if (!dt.isNucleotide())
                    throw XProj(format("Partition subset has data type \"%s\" but data read from file has data type \"nucleotide\"") % dt.getDataTypeAsString());
            }
        }
        else if (datatype == NxsCharactersBlock::protein) {
            if (!dt.isProtein())
                throw XProj(format("Partition subset has data type \"%s\" but data read from file has data type \"protein\"") % dt.getDataTypeAsString());
        }
        else if (datatype == NxsCharactersBlock::standard) {
            if (!dt.isStandard())
                throw XProj(format("Partition subset has data type \"%s\" but data read from file has data type \"standard\"") % dt.getDataTypeAsString());
            assert(charBlock->GetSymbols());
            string symbols = string(charBlock->GetSymbols());
            dt.setStandardNumStates((unsigned)symbols.size());
        }
        else {
            // ignore block because data type is not one that is supported
            return nchar_before;
        }
        
        unsigned num_states = dt.getNumStates();
        
        // Make sure all states can be accommodated in a variable of type state_t   
        unsigned bits_in_state_t = 8*sizeof(state_t);
        if (num_states > bits_in_state_t)
            throw XProj(format("This program can only process data types with fewer than %d states") % bits_in_state_t);   
        
        // Copy data matrix from NxsCharactersBlock object to _data_matrix
        // Loop through all taxa, processing one row from block for each taxon
        for (unsigned t = 0; t < ntax; ++t) {

            const NxsDiscreteStateRow & row = block->GetDiscreteMatrixRow(t);
            if (seqlen == 0)
                seqlen = (unsigned)row.size();
            _data_matrix[t].resize(nchar_before + seqlen);
            
            // Loop through all sites/characters in row corresponding to taxon t
            unsigned k = nchar_before;
            for (int raw_state_code : row) {
                // For codon model, raw_state_code ranges from 0-63, but deletion of stop codons means fewer state codes
                state_t state = numeric_limits<state_t>::max(); // complete ambiguity, all bits set
                bool complete_ambiguity = (!dt.isCodon() && raw_state_code == (int)num_states);
                bool all_missing_or_gaps = (raw_state_code < 0);
                if ((!complete_ambiguity) && (!all_missing_or_gaps)) {
                    int state_code = raw_state_code;
                    if (dt.isCodon())
                        state_code = dt.getGeneticCode()->getStateCode(raw_state_code);

                    if (state_code < (int)num_states) {
                        state = (state_t)1 << state_code;
                    }
                    else {
                        // incomplete ambiguity (NCL state code > num_states)
                        const NxsDiscreteDatatypeMapper      * mapper = block->GetDatatypeMapperForChar(k - nchar_before);
                        const set<NxsDiscreteStateCell> & state_set = mapper->GetStateSetForCode(raw_state_code);
                        state = 0;
                        for (auto s : state_set) {
                             state |= (state_t)1 << s;
                        }
                    }
                }
                _data_matrix[t][k++] = state;
            }
        }
        
        return seqlen;
    }    

    inline void Data::writeDataToFile(const string filename) {
        // Creates NEXUS data file with specified filename
        
        // Gather information
        unsigned ntax = (unsigned)_taxon_names.size();
        unsigned nchar = (unsigned)accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
        unsigned longest_taxon_name = 0;
        for (auto nm : _taxon_names) {
            if (nm.size() > longest_taxon_name)
                longest_taxon_name = (unsigned)nm.size();
        }
        const format name_format( str(format("    %%%ds ") % longest_taxon_name) );
        
        ofstream nexf(filename);
        
        nexf << "#NEXUS\n\n";
        
        nexf << "begin data;\n";
        nexf << "  dimensions ntax=" << ntax << " nchar=" << nchar << ";\n";
        nexf << "  format datatype=dna gap=- missing=?;\n";
        nexf << "  matrix\n";
        
        for (unsigned t = 0; t < _taxon_names.size(); t++) {
            nexf << str(format(name_format) % _taxon_names[t]);
            
            unsigned pattern = 0;
            for (unsigned g = 0; g < _partition->getNumSubsets(); g++) {
                // Determine whether this taxon will have all missing data
                // for this subset
                bool taxon_all_missing = false;
                if (G::_occupancy < 1.0) {
                    double u = rng->uniform();
                    if (u > G::_occupancy)
                        taxon_all_missing = true;
                }
            
                unsigned npatterns = getNumPatternsInSubset(g);
                for (unsigned j = 0; j < npatterns; j++) {
                    char dna_letter = '?';

                    if (!taxon_all_missing) {
                        state_t s = _data_matrix[t][pattern];
                        bool is_A = ((state_t)1 << 0 == s);
                        bool is_C = ((state_t)1 << 1 == s);
                        bool is_G = ((state_t)1 << 2 == s);
                        bool is_T = ((state_t)1 << 3 == s);
                        
                        if (is_A && is_C && is_G && is_T)
                            dna_letter = 'N';
                        else if (is_A && is_C && is_T)
                            dna_letter = 'H';
                        else if (is_C && is_G && is_T)
                            dna_letter = 'B';
                        else if (is_A && is_C && is_G)
                            dna_letter = 'V';
                        else if (is_A && is_G && is_T)
                            dna_letter = 'D';
                        else if (is_A && is_G)
                            dna_letter = 'R';
                        else if (is_C && is_T)
                            dna_letter = 'Y';
                        else if (is_A && is_C)
                            dna_letter = 'M';
                        else if (is_G && is_T)
                            dna_letter = 'K';
                        else if (is_C && is_G)
                            dna_letter = 'S';
                        else if (is_A && is_T)
                            dna_letter = 'W';
                        else if (is_A)
                            dna_letter = 'A';
                        else if (is_C)
                            dna_letter = 'C';
                        else if (is_G)
                            dna_letter = 'G';
                        else if (is_T)
                            dna_letter = 'T';
                        assert(dna_letter != '?');
                    }
                    
                    unsigned pattern_count = _pattern_counts[pattern];
                    for (unsigned k = 0; k < pattern_count; k++) {
                        nexf << dna_letter;
                    }
                    pattern++;
                }
            }
            nexf << endl;
        }
        
        nexf << "  ;\n";
        nexf << "end;\n";
        
        nexf.close();
    }
    
    inline void Data::getDataFromFile(const string filename) {
        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for documentation
        //
        // -1 means "process all blocks found" (this is a bit field and -1 fills the bit field with 1s)
        // Here are the bits (and nexus blocks) that are defined:
        //     enum NexusBlocksToRead
        //     {
        //         NEXUS_TAXA_BLOCK_BIT = 0x01,
        //         NEXUS_TREES_BLOCK_BIT = 0x02,
        //         NEXUS_CHARACTERS_BLOCK_BIT = 0x04,
        //         NEXUS_ASSUMPTIONS_BLOCK_BIT = 0x08,
        //         NEXUS_SETS_BLOCK_BIT = 0x10,
        //         NEXUS_UNALIGNED_BLOCK_BIT = 0x20,
        //         NEXUS_DISTANCES_BLOCK_BIT = 0x40,
        //         NEXUS_UNKNOWN_BLOCK_BIT = 0x80
        //     };
        //MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
        MultiFormatReader nexusReader(-1, NxsReader::IGNORE_WARNINGS);
        
        // Both of these needed to suppress "storing read block" messages
        // see NxsReader::statusMessage in nxsreader.cpp
        nexusReader.SetAlwaysReportStatusMessages(false);
        nexusReader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);

        try {
            nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {    
            nexusReader.DeleteBlocksFromFactories();
            throw;
        }   

        // Commit to storing new data
        clear();

        // Ensure that Data::setPartition was called before reading data
        assert(_partition);

        int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
        if (numTaxaBlocks == 0)
            throw XProj("No taxa blocks were found in the data file");
            
        _cum_nchar = 0;
        for (int i = 0; i < numTaxaBlocks; ++i) {
            NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
            unsigned ntax = storeTaxonNames(taxaBlock, i);
            const unsigned numCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
            for (unsigned j = 0; j < numCharBlocks; ++j) {
                NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
                NxsCharactersBlock::DataTypesEnum datatype = charBlock->GetOriginalDataType();
                _cum_nchar += storeData(ntax, _cum_nchar, charBlock, datatype);
            }
        }   

        // No longer any need to store raw data from nexus file
        nexusReader.DeleteBlocksFromFactories();

        // Compress _data_matrix so that it holds only unique patterns (counts stored in _pattern_counts)
        if (_data_matrix.empty()) {
            output(format("No data were stored from the file \"%s\"\n") % filename, 2);
            clear();
        }
        else {
            compressPatterns();
        }
    }
    
    inline char Data::stateCodeToDNALetter(state_t s) const {
        bool is_A = ((state_t)1 << 0 == s);
        bool is_C = ((state_t)1 << 1 == s);
        bool is_G = ((state_t)1 << 2 == s);
        bool is_T = ((state_t)1 << 3 == s);
        
        char dna_letter = '?';
        if (is_A && is_C && is_G && is_T)
            dna_letter = 'N';
        else if (is_A && is_C && is_T)
            dna_letter = 'H';
        else if (is_C && is_G && is_T)
            dna_letter = 'B';
        else if (is_A && is_C && is_G)
            dna_letter = 'V';
        else if (is_A && is_G && is_T)
            dna_letter = 'D';
        else if (is_A && is_G)
            dna_letter = 'R';
        else if (is_C && is_T)
            dna_letter = 'Y';
        else if (is_A && is_C)
            dna_letter = 'M';
        else if (is_G && is_T)
            dna_letter = 'K';
        else if (is_C && is_G)
            dna_letter = 'S';
        else if (is_A && is_T)
            dna_letter = 'W';
        else if (is_A)
            dna_letter = 'A';
        else if (is_C)
            dna_letter = 'C';
        else if (is_G)
            dna_letter = 'G';
        else if (is_T)
            dna_letter = 'T';
        assert(dna_letter != '?');
        return dna_letter;
    }
    
    inline string Data::getPatternAsString(unsigned i) const {
        ostringstream pattern_as_string;
        for (unsigned t = 0; t < _taxon_names.size(); t++) {
            char dna_letter = stateCodeToDNALetter(_data_matrix[t][i]);
            pattern_as_string << dna_letter;
        }
        return pattern_as_string.str();
    }
    
    inline void Data::memoryReport(ofstream & memf) const {
        memf << "\nData memory report:\n\n";
        unsigned long total_nbytes = 0;
        unsigned ntaxa = getNumTaxa();
        unsigned nsubsets = getNumSubsets();
        memf << str(format("  Number of taxa: %d\n") % ntaxa);
        memf << str(format("  Number of subsets: %d\n") % nsubsets);
        memf << str(format("  Number of sites: %d\n\n") % _cum_nchar);
        
        memf << str(format("  %12s %12s %12s\n") % "subset" % "nstates" % "npatterns");
        memf << str(format("  %12s %12s %12s\n") % " -----------" % " -----------" % " -----------");
        unsigned total_npatterns = 0;
        for (unsigned subset = 0; subset < nsubsets; ++subset) {
            unsigned npatterns = getNumPatternsInSubset(subset);
            unsigned nstates = getNumStatesForSubset(subset);
            memf << str(format("  %12d %12d %12d\n") % subset % nstates % npatterns);
            total_npatterns += npatterns;
        }
        memf << str(format("  %12s %12s %12s\n") % " -----------" % " -----------" % " -----------");
        memf << str(format("  %12s %12s %12d\n\n") % " " % " " % total_npatterns);
        
        unsigned element_size = (unsigned)sizeof(unsigned);
        unsigned nelements = (unsigned)_pattern_counts.size();
        unsigned nbytes = nelements*element_size;
        total_nbytes += nbytes;
        memf << str(format("  _pattern_counts:     %d bytes = %d patterns * %d bytes\n") % nbytes % nelements % element_size);

        element_size = (unsigned)sizeof(state_t);
        nelements = (unsigned)_monomorphic.size();
        nbytes = nelements*element_size;
        total_nbytes += nbytes;
        memf << str(format("  _monomorphic:        %d bytes = %d patterns * %d bytes\n") % nbytes % nelements % element_size);
        
        element_size = (unsigned)sizeof(int);
        nelements = (unsigned)_partition_key.size();
        nbytes = nelements*element_size;
        total_nbytes += nbytes;
        memf << str(format("  _partition_key:      %d bytes = %d patterns * %d bytes\n") % nbytes % nelements % element_size);

        element_size = (unsigned)sizeof(state_t);
        nelements = ntaxa*_cum_nchar;
        nbytes = nelements*element_size;
        //total_nbytes += nbytes;
        memf << str(format("  _data_matrix (pre):  %d bytes = %d taxa * %d sites * %d bytes\n") % nbytes % ntaxa % _cum_nchar % element_size);
        
        element_size = (unsigned)sizeof(state_t);
        nelements = ntaxa*total_npatterns;
        nbytes = nelements*element_size;
        total_nbytes += nbytes;
        memf << str(format("  _data_matrix (post): %d bytes = %d taxa * %d patterns * %d bytes\n") % nbytes % ntaxa % total_npatterns % element_size);
        
        //TODO: add other data members
        
        memf << str(format("\n  Total megabytes: %.5f\n") % (1.0*total_nbytes/1048576));
        
    }
    
}
