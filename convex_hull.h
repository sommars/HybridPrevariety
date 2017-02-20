#include "prevariety_util.h"

//------------------------------------------------------------------------------
vector<Cone> NewHull(vector<vector<int> > Points, vector<double> VectorForOrientation, bool Verbose);

//------------------------------------------------------------------------------
void FindFacets(Hull &H);

//------------------------------------------------------------------------------
void FindEdges(Hull &H);

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull H);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points);
