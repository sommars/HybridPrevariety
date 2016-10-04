#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

class Cone {
	public:
		C_Polyhedron HOPolyhedron;
		C_Polyhedron ClosedPolyhedron;
		vector<set<int> > ClosedIntersectionIndices;
		set<int> PolytopesVisited;
};

class Edge {
	public:
		Cone EdgeCone;
		bool HasEdgeCone;
		set<int> PointIndices;
		set<int> NeighborIndices;
};

class Facet {
	public:
		set<int> PointIndices;
};

class Hull {
	public:
		vector<vector<int> > Points;
		map<vector<int>,int> PointToIndexMap;
		map<int,vector<int> > IndexToPointMap;
		vector<Edge> Edges;
		vector<Facet> Facets;
		C_Polyhedron CPolyhedron;
		int AffineDimension;
		int SpaceDimension;
};

//------------------------------------------------------------------------------
double GetPolyhedralIntersectionTime();

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron &ph1, C_Polyhedron &ph2);

//------------------------------------------------------------------------------
Cone IntersectCones(Cone C1, Cone C2);

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
vector<vector<int> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void FindFacets(Hull &H);

//------------------------------------------------------------------------------
void FindEdges(Hull &H);

//------------------------------------------------------------------------------
Constraint StrictInequalityToNonStrictInequality(Constraint c);

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull H);

//------------------------------------------------------------------------------
int InnerProduct(vector<int> V1, vector<int> V2);

//------------------------------------------------------------------------------
double DoubleInnerProduct(vector<int> V1, vector<double> V2);

//------------------------------------------------------------------------------
vector<vector<int> > FindInitialForm(vector<vector<int> > &Points, vector<int> &Vector);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point);

//------------------------------------------------------------------------------
void PrintCPolyhedron(C_Polyhedron ph, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<C_Polyhedron> phs, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
set<int> IntersectSets(set<int> S1, set<int> S2);

//------------------------------------------------------------------------------
bool SetsDoIntersect(set<int> &S1, set<int> &S2);
