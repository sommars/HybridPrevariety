// Main function in this unit is NewHull, which takes a set of points and
// returns a vector of cones. These cones are what are used to initialize
// the prevariety algorithm.

#include "prevariety_util.h"

//------------------------------------------------------------------------------
vector<Cone> NewHull(
   vector<vector<int> > &Points,
   vector<double> &VectorForOrientation,
   bool Verbose);

//------------------------------------------------------------------------------
void FindFacets(Hull &H);

//------------------------------------------------------------------------------
void FindEdges(Hull &H);

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull &H);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > &Points);
