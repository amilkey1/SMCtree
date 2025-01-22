//#define NDEBUG
#include <cassert>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>   // shared_ptr
#include <new>      // used by Mallocator: bad_alloc, bad_array_new_length
#include <numeric>
#include <queue>
#include <regex>
#include <set>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/program_options.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "ncl/nxsmultiformat.h"

using namespace std;
using boost::format;
using boost::str;
using boost::algorithm::join;

#include "conditionals.hpp"
#include "xproj.hpp"
#include "lot.hpp"
#include "split.hpp"
#include "g.hpp"
#include "stopwatch.hpp"
#include "genetic-code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "partial_store.hpp"
#include "node.hpp"
#include "forest.hpp"
#include "particle.hpp"
//#include "forestpol.hpp"
#include "proj.hpp"

proj::PartialStore ps;

using namespace proj;

void output(format & fmt, unsigned level) {
    if (G::_verbosity > 0 && level <= G::_verbosity)
        cout << str(fmt);
}

void output(string msg, unsigned level) {
    if (G::_verbosity > 0 && level <= G::_verbosity)
        cout << msg;
}

Lot::SharedPtr   rng(new Lot());
StopWatch        stopwatch;

bool     G::_simulating        = false;
bool     G::_debugging         = false;
unsigned G::_nthreads          = 1;
unsigned G::_nstates           = 4;
unsigned G::_ntaxa             = 0;
unsigned G::_nloci             = 0;
double   G::_lambda            = 1.0;
double   G::_mu                = 1.0;
double   G::_kappa             = 1.0;
double   G::_root_age          = 1.0;
string   G::_string_base_frequencies = "";
vector<double> G::_base_frequencies;
string   G::_string_relative_rates = "";
vector<double> G::_double_relative_rates;
string G::_proposal = "prior-prior";
bool G::_est_lambda = "false";
bool G::_est_mu = "false";
bool G::_est_root_age = "false";

double   G::_small_enough      = 0.00001;
double   G::_infinity          = numeric_limits<double>::infinity();
double   G::_negative_infinity = -numeric_limits<double>::infinity();

unsigned  G::_rnseed = 1;
unsigned  G::_verbosity = 1;

string    G::_program_name = "smctree";
unsigned  G::_major_version = 0;
unsigned  G::_minor_version = 0;
string    G::_model = "JC";

unsigned  G::_sim_ntaxa    = 4;
double    G::_sim_lambda   = 1.0;
double    G::_sim_mu       = 0.0;
double    G::_sim_rho      = 1.0;
double    G::_sim_root_age = 1.0;
string    G::_sim_filename_prefix = "sim";

double    G::_occupancy  = 1.0;
double    G::_comphet    = G::_infinity;
double    G::_asrv_shape = G::_infinity;

vector<string>         G::_taxon_names;
vector<string>         G::_locus_names;
vector<unsigned>       G::_nsites_per_locus;
map<unsigned,unsigned> G::_nexus_taxon_map;
map<unsigned, double>  G::_relrate_for_locus;
bool                   G::_save_memory;

string                  G::_start_mode = "smc";
string                  G::_filename = "data.nex";
unsigned                G::_nparticles = 500;

const double Node::_smallest_edge_length=1.0e-12;

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required in order to use infinity()");

GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = { 
                             // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

int main(int argc, const char * argv[]) {

    Proj proj;
    bool normal_termination = true;
    try {
        proj.processCommandLineOptions(argc, argv);
        
        StopWatch sw;
        sw.start();
        proj.run();
        double total_seconds = sw.stop();
        output(format("\nTotal time: %.5f seconds\n") % total_seconds, 1);
    }
    catch(std::exception & x) {
        cerr << str(format("Exception: %s\n") % x.what());
        cerr << "Aborted.\n";
        normal_termination = false;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        normal_termination = false;
    }
    
    return 0;
}
