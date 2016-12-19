#include <ppl.hh>
#include <boost/dynamic_bitset.hpp>
#include "Thread_Pool_defs.hh"

using namespace std;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > RandomSimplices(int n);
