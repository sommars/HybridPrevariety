// Defines some basic polynomial systems to be used in testing.
#include "prevariety_types.h"
#include "prevariety_inlines.h"

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced);
// Returns the support set for the cyclic n roots, or the reduced cyclic n roots

//------------------------------------------------------------------------------
vector<vector<vector<int> > > RandomSimplices(int n);
// Returns the support set for a set of n random simplices of dimension n.

//------------------------------------------------------------------------------
vector<vector<vector<int> > > Viviani(int n);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > HiddenRay(int n);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > FourByFourMinors(int n);


