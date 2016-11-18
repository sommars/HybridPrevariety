#include <ppl.hh>
#include <iostream>
#include <list>
#include <mutex>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

struct Cone {
	NNC_Polyhedron HOPolyhedron;
	vector<set<int> > IntersectionIndices; // 1 means that there is an intersection, 0 that there isn't
	vector<int> PolytopesVisited; // 1 means visited, 0 means unvisited
	int PolytopesVisitedCount;
};

struct Edge {
	Cone EdgeCone;
	bool HasEdgeCone;
	set<int> PointIndices;
	set<int> NeighborIndices;
};

struct Facet {
	set<int> PointIndices;
};

struct Hull {
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

struct ConeWithIndicator {
	vector<int> RayIndices;
	bool IsMaximal;
};

inline bool operator== (const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
    return  s1.RayIndices == s2.RayIndices;
}

inline bool operator< (const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
    return  s1.RayIndices < s2.RayIndices;
}

struct TropicalPrevariety {
	map<vector<int>, int> RayToIndexMap;
	vector<set<ConeWithIndicator > > ConeTree;
	vector<vector<int> > Pretropisms; // Consider ripping out
	map<int, Generator> IndexToGenMap;
	vector<int> FVector; // Consider ripping out
};

struct ThreadJob {
	mutable mutex M;
	vector<vector<Cone> > SharedCones;
	ThreadJob(vector<vector<Cone> > ConeVector): SharedCones(ConeVector) {};
	ThreadJob(const ThreadJob& TJ) {
		lock_guard<mutex> lock(TJ.M);
		SharedCones = TJ.SharedCones;
	};
};

//------------------------------------------------------------------------------
inline NNC_Polyhedron IntersectCones(NNC_Polyhedron &ph1, NNC_Polyhedron &ph2) {
	NNC_Polyhedron ph = ph1;
	ph.add_constraints(ph2.constraints());
	ph.affine_dimension();
	return ph;
};

//------------------------------------------------------------------------------
inline C_Polyhedron IntersectCones(C_Polyhedron &ph1, C_Polyhedron &ph2) {
	C_Polyhedron ph = ph1;
	ph.add_constraints(ph2.constraints());
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
inline set<int> IntersectSets(set<int> &S1, set<int> &S2) {
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

