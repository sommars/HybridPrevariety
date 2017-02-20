#include <ppl.hh>
#include <boost/dynamic_bitset.hpp>
#include "Thread_Pool_defs.hh"
#include "/home/jeff/Desktop/Software/soplex-2.2.1/src/soplex.h"

using namespace std;
using namespace soplex;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > RandomSimplices(int n);
