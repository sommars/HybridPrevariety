#include <ppl.hh>
#include <iostream>
#include <list>
#include <vector>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > RandomSimplices(int n);
