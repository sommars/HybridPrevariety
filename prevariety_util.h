#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

class Cone {
	public:
		NNC_Polyhedron HOPolyhedron;
		NNC_Polyhedron ClosedPolyhedron;
		vector<set<int> > ClosedIntersectionIndices;
		vector<int> PolytopesVisited;
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
		NNC_Polyhedron CPolyhedron;
		int AffineDimension;
		int SpaceDimension;
		vector<Constraint> Eqs;
};

//------------------------------------------------------------------------------
inline NNC_Polyhedron IntersectCones(NNC_Polyhedron &ph1, NNC_Polyhedron &ph2) {
	NNC_Polyhedron ph = ph1;
	ph.add_constraints(ph2.minimized_constraints());
	ph.affine_dimension();
	return ph;
};

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
vector<vector<int> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<int> > Points, vector<double> VectorForOrientation);

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
inline vector<vector<int> > FindInitialForm(vector<vector<int> > &Points, vector<int> &Vector) {
	/*
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};
	vector<vector<int> > InitialForm;

	vector<int> *Point;
	Point = &(Points[0]);
	InitialForm.push_back(*Point);
	int MinimalIP = 0;
	for (size_t i = 0; i != (*Point).size(); i++)
		MinimalIP += Vector[i] * (*Point)[i];

	for (size_t i = 1; i != Points.size(); i++) {
		Point = &(Points[i]);
		int IP = 0;
		for (size_t j = 0; j != (*Point).size(); j++)
			IP += Vector[j] * (*Point)[j];
		if (MinimalIP > IP) {
			MinimalIP = IP;
			InitialForm.clear();
			InitialForm.push_back(*Point);
		} else if (IP == MinimalIP)
			InitialForm.push_back(*Point);
	};

	return InitialForm;
};

//------------------------------------------------------------------------------
NNC_Polyhedron FindCPolyhedron(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point);

//------------------------------------------------------------------------------
void PrintCPolyhedron(NNC_Polyhedron ph, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<NNC_Polyhedron> phs, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
inline set<int> IntersectSets(set<int> S1, set<int> S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	set<int> Result;
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			Result.insert(*S1Itr);
			S1Itr++;
			S2Itr++;
		};
	};
	return Result;
};

//------------------------------------------------------------------------------
inline bool SetsDoIntersect(set<int> &S1, set<int> &S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			return true;
		};
	};
	return false;
};

//------------------------------------------------------------------------------
set<int> PreintersectWalkPolytope(int HullIndex, Cone NewCone, vector<Hull> &Hulls);
