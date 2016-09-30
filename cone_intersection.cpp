#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "polynomial_systems.h"
#include "prevariety_util.h"
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

double ContainmentTime, IntersectionTime, ParallelTime, TestTime, GetConesTime, PolytopePickingTime;
int ConeIntersectionCount;

//------------------------------------------------------------------------------
vector<Cone> WalkPolytope(int HullIndex, Cone &NewCone, vector<Hull> &Hulls) {
	clock_t BeginTestTime = clock();
	Recycle_Input dummy;
	Hull *HIndex;
	HIndex = &Hulls[HullIndex];
	set<int> *ClosedIntersectionIndices;
	ClosedIntersectionIndices = &NewCone.ClosedIntersectionIndices[HullIndex];
	TestTime += double(clock() - BeginTestTime);
	//take a random vector from half open cone
	vector<int> RandomVector(Hulls[0].SpaceDimension, 0);
	Generator_System gs = NewCone.HOPolyhedron.minimized_generators();
	for (Generator_System::const_iterator i = gs.begin(),gs_end = gs.end(); i != gs_end; ++i) {
		if (!(*i).is_ray() && !(*i).is_line()) {
			continue;
		};
		for (size_t j = 0; j != Hulls[0].SpaceDimension; j++) {
			stringstream s;
			s << (*i).coefficient(Variable(j));
			int ToAppend;
			istringstream(s.str()) >> ToAppend;
			RandomVector[j] += ToAppend;
		};
	};
	//take initial form using that vector.
	vector<vector<int> > InitialForm = FindInitialForm((*HIndex).Points, RandomVector);

	set<int> InitialIndices;
	for (vector<vector<int> >::iterator InitialFormItr=InitialForm.begin();
	InitialFormItr != InitialForm.end(); InitialFormItr++) {
		InitialIndices.insert((*HIndex).PointToIndexMap[*InitialFormItr]);
	};

	// List of indices of edges to visit.
	vector<int> EdgesToTest;
	for(size_t EdgeIndex = 0; EdgeIndex != (*HIndex).Edges.size(); EdgeIndex++) {
		if (SetsDoIntersect(InitialIndices, (*HIndex).Edges[EdgeIndex].PointIndices)
		&& ( find((*ClosedIntersectionIndices).begin(), (*ClosedIntersectionIndices).end(), EdgeIndex) != (*ClosedIntersectionIndices).end())
		) {
			EdgesToTest.push_back(EdgeIndex);
		};
	};
	
	vector<Cone> NewCones;
	// Explore edge skeleton
	set<int> PretropGraphEdges;
	set<int> NotPretropGraphEdges;
	int EdgeToTestIndex;
	Edge *EdgeToTest;
	Constraint_System cs2 = NewCone.HOPolyhedron.minimized_constraints();
	while(!EdgesToTest.empty()) {
		EdgeToTestIndex = EdgesToTest.back();
		EdgesToTest.pop_back();
		clock_t IntBegin = clock();

		EdgeToTest = &(*HIndex).Edges[EdgeToTestIndex];
		Constraint_System cs1 = (*EdgeToTest).EdgeCone.ClosedPolyhedron.minimized_constraints();
		
		for (Constraint_System::const_iterator i = cs2.begin(),
		cs1_end = cs2.end(); i != cs1_end; ++i) {
			cs1.insert(*i);
		};
		
		Cone TempCone;
		TempCone.ClosedPolyhedron = NNC_Polyhedron(cs1);
		TempCone.ClosedPolyhedron.affine_dimension();
		TempCone.ClosedPolyhedron.minimized_constraints();
		ConeIntersectionCount++;
		IntersectionTime += double(clock() - IntBegin);
		if (TempCone.ClosedPolyhedron.affine_dimension() > 0) {
			PretropGraphEdges.insert(EdgeToTestIndex);
			
			// HOIntersectionIndices are ONLY useful so that we know whether or not to
			// bother at this point.
			clock_t IntBegin = clock();
			Constraint_System cs3 = (*EdgeToTest).EdgeCone.HOPolyhedron.minimized_constraints();
		
			for (Constraint_System::const_iterator i = cs2.begin(),
			cs1_end = cs2.end(); i != cs1_end; ++i) {
				cs3.insert(*i);
			};
			TempCone.HOPolyhedron = NNC_Polyhedron(cs3);
		
			TempCone.HOPolyhedron.affine_dimension();
			TempCone.HOPolyhedron.minimized_constraints();
			ConeIntersectionCount++;
			IntersectionTime += double(clock() - IntBegin);
	

			if (TempCone.HOPolyhedron.affine_dimension() > 0) {
				vector<set<int> > InitialSet1 (NewCone.ClosedIntersectionIndices.size());
				TempCone.ClosedIntersectionIndices = InitialSet1;

				for (size_t i = 0; i != NewCone.ClosedIntersectionIndices.size(); i++) {
					if (find(NewCone.PolytopesVisited.begin(), NewCone.PolytopesVisited.end(), i) == NewCone.PolytopesVisited.end()) {
						TempCone.ClosedIntersectionIndices[i] = IntersectSets((*EdgeToTest).EdgeCone.ClosedIntersectionIndices[i], NewCone.ClosedIntersectionIndices[i]);
					};
				};
				NewCones.push_back(TempCone);
			};
			set<int>::iterator NeighborItr;
			for(NeighborItr=(*EdgeToTest).NeighborIndices.begin(); NeighborItr!=(*EdgeToTest).NeighborIndices.end(); NeighborItr++) {
				int Neighbor = *NeighborItr;
				if (( find((*ClosedIntersectionIndices).begin(), (*ClosedIntersectionIndices).end(), Neighbor) != (*ClosedIntersectionIndices).end())
				&& ( find(PretropGraphEdges.begin(), PretropGraphEdges.end(), Neighbor) == PretropGraphEdges.end() )
				&& ( find(NotPretropGraphEdges.begin(), NotPretropGraphEdges.end(), Neighbor) == NotPretropGraphEdges.end() )
				&& ( find(EdgesToTest.begin(), EdgesToTest.end(), Neighbor) == EdgesToTest.end() )) {
					EdgesToTest.push_back(Neighbor);
				};
			};
		} else {
			NotPretropGraphEdges.insert(EdgeToTestIndex);
		};

	};
	return NewCones;
}

//------------------------------------------------------------------------------
void DynamicEnumerate(Cone &C, vector<Hull> &Hulls, vector<vector<int> > &Pretropisms) {
	clock_t BeginTime = clock();
	// Figure out which polytope we want to pick 
	int SmallestInt = 10000000; // Lazy.
	int SmallestIndex = -1;
	for (size_t i = 0; i != C.ClosedIntersectionIndices.size(); i++) {
		if ((C.ClosedIntersectionIndices[i].size() < SmallestInt)
		&&  (find(C.PolytopesVisited.begin(), C.PolytopesVisited.end(), i) == C.PolytopesVisited.end())) {
			SmallestInt = C.ClosedIntersectionIndices[i].size();
			SmallestIndex = i;
		};
	};
	if (SmallestIndex == -1) {
		cout << "Internal error: DynamicEnumerate had a value of -1 for SmallestIndex" << endl;
		cin.get();
	};
	C.PolytopesVisited.insert(SmallestIndex);
	
	PolytopePickingTime += double(clock() - BeginTime);
	vector<Cone> ResultCones = WalkPolytope(SmallestIndex, C, Hulls);
	for (size_t i = 0; i != ResultCones.size(); i++) {
		ResultCones[i].PolytopesVisited = C.PolytopesVisited;
	};
	if (ResultCones.size() == 0) {
		return;
	};

	if (ResultCones[0].PolytopesVisited.size() == Hulls.size()) {
		for (size_t i = 0; i != ResultCones.size(); i++) {
			Generator_System gs = ResultCones[i].HOPolyhedron.minimized_generators();
			for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
				if ((*gsi).is_point() or (*gsi).is_closure_point()) {
					continue;
				};
				vector<int> Pt = GeneratorToPoint(*gsi);
				if ( find(Pretropisms.begin(), Pretropisms.end(), Pt) == Pretropisms.end() ) {
				//	gv.push_back(gen);
				Pretropisms.push_back(GeneratorToPoint(*gsi));
				};
			};
		};
	} else {
		for (size_t i = 0; i != ResultCones.size(); i++) {
			DynamicEnumerate(ResultCones[i], Hulls, Pretropisms);
		};
	};
};

//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	clock_t StartTime = clock();
	if (argc != 3) {
		cout << "Internal error: expected two arguments." << endl;
		return 1;
	}
	int n = atoi(argv[1]);

	string SystemName = string(argv[2]);
	vector<vector<vector<int> > > PolynomialSystemSupport;
	if (SystemName == "reducedcyclicn") {
		PolynomialSystemSupport = CyclicN(n, true);
	} else if (SystemName == "cyclicn") {
		PolynomialSystemSupport = CyclicN(n, false);
	} else if (SystemName == "random") {
		PolynomialSystemSupport = RandomSimplices(n);
	} else {
		cout << "Internal error: only supported systems are: "
			<< "reducedcyclicn, cyclicn, or random." << endl;
		return 1;
	};
	
	vector<Hull> Hulls;
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++) {
		Hulls.push_back(NewHull(PolynomialSystemSupport[i]));
	}
	
	double HullTime = double(clock() - StartTime);
	clock_t AlgorithmStartTime = clock();
	double PreintersectTime = 0;

	for(int i = 0; i != Hulls.size(); i++){
		for(int j = 0; j != Hulls[i].Edges.size(); j++){
			vector<set<int> > InitialSet1 (Hulls.size());
			vector<set<int> > InitialSet2 (Hulls.size());
			Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices = InitialSet1;
		};
	};
	clock_t PreintTimeStart = clock();
	int TotalInt = 0;
	int NonInt = 0;
	// TODO: switch this to explore the edge skeleton method.
	for(int i = 0; i != Hulls.size(); i++){
		int ExpectedDim = Hulls[0].Edges[0].EdgeCone.ClosedPolyhedron.affine_dimension() - 1;
		vector<Edge> Edges1 = Hulls[i].Edges;
		for (int k = 0; k != Hulls[i].Edges.size(); k++) {
			Hulls[i].Edges[k].EdgeCone.PolytopesVisited.insert(i);
		};

		for(int j = i+1; j != Hulls.size(); j++){
			vector<Edge> Edges2 = Hulls[j].Edges;
		
			for(int k = 0; k != Edges1.size(); k++){
				for(int l = 0; l != Edges2.size(); l++){
					// Intersect pairs of closed cones. Osserman/Payne applies
					if (IntersectCones(Edges1[k].EdgeCone.ClosedPolyhedron, Edges2[l].EdgeCone.ClosedPolyhedron).affine_dimension() >= ExpectedDim) {
						Hulls[i].Edges[k].EdgeCone.ClosedIntersectionIndices[j].insert(l);
						Hulls[j].Edges[l].EdgeCone.ClosedIntersectionIndices[i].insert(k);
					} else {
						NonInt++;
					};
					// Intersect pairs of half open cones. Osserman/Payne does not apply.
					TotalInt++;
				};
			};
		};
		printf("Finished level %d of pre-intersections.\n", i);
	};
	cout << "Total Intersections: " << TotalInt << ", Non Intersections: " << NonInt << endl;
	PreintersectTime = double(clock() - PreintTimeStart);
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;

	// This is one way to do it. It seems reasonable to pick sum, median, mean, min...
	int SmallestInt = 1000000;
	int SmallestIndex = -1;
	for (size_t i = 0; i != Hulls.size(); i++) {
		int TestValue = 0;
		for (size_t j = 0; j != Hulls[i].Edges.size(); j++) {
			TestValue += Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices.size();
		};
		if (TestValue < SmallestInt) {
			TestValue = SmallestInt;
			SmallestIndex = i;
		};
	};

	if (SmallestIndex == -1) {
		cout << "Internal error: DynamicEnumerate had a value of -1 for SmallestIndex" << endl;
		cin.get();
	};
	// TODO: Kill line below when not testing.
	//SmallestIndex = 0;
	vector<vector<int> > Pretropisms;
	for (size_t i = 0; i != Hulls[SmallestIndex].Edges.size(); i++) {
		DynamicEnumerate(Hulls[SmallestIndex].Edges[i].EdgeCone, Hulls, Pretropisms);
	};

	clock_t CleanupStart = clock();
	vector<Cone>::iterator PolyItr;
	double AlgTime = double(clock() - AlgorithmStartTime);
	// This is parsing and displaying all of the pretropisms.
	sort(Pretropisms.begin(), Pretropisms.end());
	PrintPoints(Pretropisms);
	cout << "Number of pretropisms found: " << Pretropisms.size() << endl;
	cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
	cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Parallel time: " << ParallelTime / CLOCKS_PER_SEC << endl;
	cout << "Cleanup time: " << double(clock() - CleanupStart) / CLOCKS_PER_SEC << endl;
	cout << "GetCones time: " << GetConesTime / CLOCKS_PER_SEC << endl;
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	cout << "Alg time: " << AlgTime / CLOCKS_PER_SEC << endl;
	cout << "Test time: " << TestTime / CLOCKS_PER_SEC << endl;
	cout << "Polytope Picking time: " << PolytopePickingTime / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;
	
}
