#include "polynomial_systems.h"

//------------------------------------------------------------------------------
struct BitsetWithCount {
	boost::dynamic_bitset<> Indices; // 1 means that there is an intersection, 0 that there isn't
	int Count;
};

//------------------------------------------------------------------------------
struct Cone {
	vector<BitsetWithCount> RelationTables;
	BitsetWithCount PolytopesVisited;
	C_Polyhedron HOPolyhedron;
};

//------------------------------------------------------------------------------
struct Edge {
	set<int> PointIndices;
	set<int> NeighborIndices;
};

//------------------------------------------------------------------------------
struct Facet {
	set<int> PointIndices;
};

//------------------------------------------------------------------------------
struct Hull {
	vector<vector<int> > Points;
	map<vector<int>,int> PointToIndexMap;
	map<int,vector<int> > IndexToPointMap;
	vector<Edge> Edges;
	vector<Cone> Cones;
	vector<Facet> Facets;
	C_Polyhedron CPolyhedron;
	int AffineDimension;
	int SpaceDimension;
};

//------------------------------------------------------------------------------
struct ConeWithIndicator {
	vector<int> RayIndices;
	bool IsMaximal;
};

//------------------------------------------------------------------------------
inline bool operator== (const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
    return  s1.RayIndices == s2.RayIndices;
}

//------------------------------------------------------------------------------
inline bool operator< (const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
    return  s1.RayIndices < s2.RayIndices;
}

//------------------------------------------------------------------------------
struct TropicalPrevariety {
	map<vector<int>, int> RayToIndexMap;
	vector<set<ConeWithIndicator > > ConeTree;
	vector<vector<int> > Pretropisms; // Consider ripping out
	map<int, Generator> IndexToGenMap;
	vector<int> FVector; // Consider ripping out
};

//------------------------------------------------------------------------------
struct ThreadJob {
	mutable mutex M;
	vector<list<Cone> > SharedCones;
	ThreadJob(vector<list<Cone> > ConeVector): SharedCones(ConeVector) {};
	ThreadJob(const ThreadJob& TJ) {
		lock_guard<mutex> lock(TJ.M);
		SharedCones = TJ.SharedCones;
	};
};

//------------------------------------------------------------------------------
inline NNC_Polyhedron IntersectCones(NNC_Polyhedron ph1, NNC_Polyhedron &ph2) {
	ph1.add_constraints(ph2.constraints());
	return ph1;
};

//------------------------------------------------------------------------------
inline C_Polyhedron IntersectCones(C_Polyhedron ph1, C_Polyhedron &ph2) {
	ph1.add_constraints(ph2.constraints());
	ph1.affine_dimension();
	return ph1;
};

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g, bool KnockOffLastTerm);

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
vector<Cone> NewHull(vector<vector<int> > Points, vector<double> VectorForOrientation, bool Verbose);

//------------------------------------------------------------------------------
void FindFacets(Hull &H);

//------------------------------------------------------------------------------
void FindEdges(Hull &H);

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
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point);

//------------------------------------------------------------------------------
void PrintPointForPython(vector<int> Point);

//------------------------------------------------------------------------------
inline void PrintPoint(vector<bool> Point) {
	vector<bool>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
};


//------------------------------------------------------------------------------
void PrintPointsForPython(vector<vector<int> > Points);

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

//------------------------------------------------------------------------------
inline BitsetWithCount IntersectRTs(BitsetWithCount R1, BitsetWithCount &R2) {
	R1.Indices = R1.Indices&=R2.Indices;
	R1.Count = R1.Indices.count();
	return R1;
};

//------------------------------------------------------------------------------
vector<vector<vector<int> > > ParseToSupport(string Input);
