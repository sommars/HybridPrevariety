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

double IntersectionTime, TestTime, PolytopePickingTime;
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
	for (Generator_System::const_iterator i = 
	NewCone.HOPolyhedron.minimized_generators().begin(), gs_end = 
	NewCone.HOPolyhedron.minimized_generators().end(); i != gs_end; ++i) {
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
	while(!EdgesToTest.empty()) {
		EdgeToTestIndex = EdgesToTest.back();
		EdgesToTest.pop_back();
		clock_t IntBegin = clock();

		EdgeToTest = &(*HIndex).Edges[EdgeToTestIndex];
		Constraint_System cs1 = (*EdgeToTest).EdgeCone.ClosedPolyhedron.minimized_constraints();
		
		for (Constraint_System::const_iterator i = NewCone.HOPolyhedron.minimized_constraints().begin(),
		cs1_end = NewCone.HOPolyhedron.minimized_constraints().end(); i != cs1_end; ++i) {
			cs1.insert(*i);
		};
		
		Cone TempCone;
		TempCone.ClosedPolyhedron = C_Polyhedron(cs1);
		TempCone.ClosedPolyhedron.affine_dimension();
		ConeIntersectionCount++;
		IntersectionTime += double(clock() - IntBegin);
		if (TempCone.ClosedPolyhedron.affine_dimension() > 0) {
			TempCone.ClosedPolyhedron.minimized_constraints();
			PretropGraphEdges.insert(EdgeToTestIndex);
			
			clock_t IntBegin = clock();
			Constraint_System cs3 = (*EdgeToTest).EdgeCone.HOPolyhedron.minimized_constraints();
		
			for (Constraint_System::const_iterator i = NewCone.HOPolyhedron.minimized_constraints().begin(),
			cs1_end = NewCone.HOPolyhedron.minimized_constraints().end(); i != cs1_end; ++i) {
				cs3.insert(*i);
			};
			TempCone.HOPolyhedron = C_Polyhedron(cs3);
		
			TempCone.HOPolyhedron.affine_dimension();
			ConeIntersectionCount++;
			IntersectionTime += double(clock() - IntBegin);


			if (TempCone.HOPolyhedron.affine_dimension() > 0) {
				TempCone.HOPolyhedron.minimized_constraints();
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
//			cout << ResultCones[i].HOPolyhedron.minimized_generators() << endl;
			for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
				if ((*gsi).is_point() or (*gsi).is_closure_point()) {
					continue;
				};
				if ((*gsi).coefficient(Variable((*gsi).space_dimension() - 1)) != 0) {
					continue;
				};
				
				vector<int> Pt = GeneratorToPoint(*gsi);
				if (find(Pretropisms.begin(), Pretropisms.end(), Pt) == Pretropisms.end()) {
					Pretropisms.push_back(Pt);
				};
				
/*
				// There has to be a better way to do this. This works for now.
				// Eventually, when switched to closed cones of dimension n+1, this
				// won't be necessary.
				Generator_System TempGS;
				TempGS.insert((*gsi));
				//
				Linear_Expression LE;
				for (size_t j = 0; j < (*gsi).space_dimension(); j++) {
					LE += Variable(j) * (*gsi).coefficient(Variable(j));
				};
				TempGS.insert(point(LE));
				C_Polyhedron RayTestPoly(TempGS);
				if (not ResultCones[i].HOPolyhedron.contains(RayTestPoly)) {
					continue;
				}
				
				if ((*gsi).is_ray()) {
					Pretropisms.push_back(GeneratorToPoint(*gsi));
				} else if ((*gsi).is_line()) {
					vector<int> PositivePretropism = GeneratorToPoint(*gsi);
					Pretropisms.push_back(PositivePretropism);
					vector<int> NegativePretropism;
					for (size_t j = 0; j != PositivePretropism.size(); j++) {
						NegativePretropism.push_back((-1) * PositivePretropism[j]);
					};
					Pretropisms.push_back(NegativePretropism);
				};
				*/
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

	srand ( time(NULL) );
	vector<Hull> Hulls;
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++) {
		Hulls.push_back(NewHull(PolynomialSystemSupport[i]));
	}
	
	
	
	
	
	
	
	/*
	vector<Cone> ConeVector;
	
	if (true == true) {
		vector<vector<Cone> > Cones;
		int HullIndex = 0;
		vector<Hull>::iterator it;
		for (it=Hulls.begin(); it != Hulls.end(); it++) {
			vector<Edge> Edges = (*it).Edges;
			vector<Cone> HullCones;
			vector<Edge>::iterator itr;

			for (itr=Edges.begin(); itr != Edges.end(); itr++) {
				HullCones.push_back((*itr).EdgeCone);
			};
		
			if (HullIndex == 0) {
				ConeVector = HullCones;
			} else {
				Cones.push_back(HullCones);
			};			
			HullIndex++;
		};

		cout << "ConeVector count: " << ConeVector.size() << endl;
		cout << "Cones count: " << Cones.size() << endl;
		
		//Iterate through Cones
		vector<vector<Cone> >::iterator ConesItr;
		int TreeLevel = 1;
		for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
			vector<Cone> TestCones = *ConesItr;
			vector<Cone> NewCones;
			//Iterate through Cones
			vector<Cone>::iterator ConeItr;
			for (ConeItr=ConeVector.begin(); ConeItr != ConeVector.end(); ConeItr++) {
				//Iterate through TestCones
				Cone Cone1 = *ConeItr;
				vector<Cone>::iterator TestConesItr;
				for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
					C_Polyhedron Temp = IntersectCones((*TestConesItr).HOPolyhedron, Cone1.HOPolyhedron);
					ConeIntersectionCount++;
					if (Temp.affine_dimension() == 0) {
						continue;
					};
					Cone NewCone;
					NewCone.HOPolyhedron = Temp;
					NewCones.push_back(NewCone);
				};
			};
			ConeVector = NewCones;
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", TreeLevel, Cones.size(), ConeVector.size(),ConeIntersectionCount);
			TreeLevel++;
		};
		
		for (size_t ii = 0; ii != ConeVector.size(); ii++) {
			cout << ConeVector[ii].HOPolyhedron.minimized_generators() << endl;
		};
		cin.get();
	}
*/











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
					if (
					(IntersectCones(Edges1[k].EdgeCone.ClosedPolyhedron, Edges2[l].EdgeCone.ClosedPolyhedron).affine_dimension() >= ExpectedDim)
					||
					  (IntersectCones(Edges1[k].EdgeCone.HOPolyhedron, Edges2[l].EdgeCone.HOPolyhedron).affine_dimension() >= 1)
					  ) {
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
	cout << "TOTAL CONES: " << Hulls[SmallestIndex].Edges.size() << endl;
	clock_t DynamicEnumerateTime = clock();
	for (size_t i = 0; i != Hulls[SmallestIndex].Edges.size(); i++) {
		cout << i << endl;
		int StartCount = ConeIntersectionCount;
		clock_t OneConeStart = clock();
		DynamicEnumerate(Hulls[SmallestIndex].Edges[i].EdgeCone, Hulls, Pretropisms);
		cout << "Ints: " << ConeIntersectionCount - StartCount << endl;
		cout << double(clock() - OneConeStart) / CLOCKS_PER_SEC << endl << endl;
	};

	vector<Cone>::iterator PolyItr;
	double AlgTime = double(clock() - AlgorithmStartTime);
	sort(Pretropisms.begin(), Pretropisms.end());
	PrintPoints(Pretropisms);
	cout << "Number of pretropisms found: " << Pretropisms.size() << endl;
	cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	cout << "Alg time: " << AlgTime / CLOCKS_PER_SEC << endl;
	cout << "Test time: " << TestTime / CLOCKS_PER_SEC << endl;
	cout << "Polytope Picking time: " << PolytopePickingTime / CLOCKS_PER_SEC << endl;
	cout << "Dynamic Enumerate time: " << double(clock() - DynamicEnumerateTime) / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;
}
