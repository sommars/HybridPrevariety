#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "polynomial_systems.h"
#include "prevariety_util.h"
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <mutex>
#include "Thread_Pool_defs.hh"
#include <chrono>
#include <thread>


using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

double IntersectionTime, TestTime, CleanupTime;
int ConeIntersectionCount;

bool c_poly_compare (const Cone& lhs, const Cone& rhs) {
  return rhs.HOPolyhedron.affine_dimension() < lhs.HOPolyhedron.affine_dimension();
};

//------------------------------------------------------------------------------
inline vector<Cone> WalkPolytope(int HullIndex, Cone &NewCone, vector<Hull> &Hulls) {
	Hull *HIndex;
	HIndex = &Hulls[HullIndex];
	set<int> *ClosedIntersectionIndices = &NewCone.ClosedIntersectionIndices[HullIndex];
	//take a random vector from half open cone
	// List of indices of edges to visit.
	set<int> EdgesToTest = *ClosedIntersectionIndices;
	vector<Cone> NewCones;
	// Explore edge skeleton
	set<int> PretropGraphEdges;
	set<int> NotPretropGraphEdges;
	int EdgeToTestIndex;
	Edge *EdgeToTest;
	while(!EdgesToTest.empty()) {
		EdgeToTestIndex = *(EdgesToTest.begin());
		EdgesToTest.erase(EdgesToTest.begin());
		clock_t IntBegin = clock();

		EdgeToTest = &(*HIndex).Edges[EdgeToTestIndex];

		Cone TempCone;
		TempCone.HOPolyhedron = (*EdgeToTest).EdgeCone.HOPolyhedron;

		TempCone.HOPolyhedron.add_constraints(NewCone.HOPolyhedron.constraints());
		TempCone.HOPolyhedron.affine_dimension();
		ConeIntersectionCount++;
		IntersectionTime += double(clock() - IntBegin);
		if (TempCone.HOPolyhedron.affine_dimension() > 0) {
			PretropGraphEdges.insert(EdgeToTestIndex);
			
			vector<set<int> > InitialSet1 (NewCone.ClosedIntersectionIndices.size());
			TempCone.ClosedIntersectionIndices = InitialSet1;
			for (size_t i = 0; i != NewCone.ClosedIntersectionIndices.size(); i++) {
			//	if (NewCone.PolytopesVisited[i] == 0)
					TempCone.ClosedIntersectionIndices[i] = IntersectSets((*EdgeToTest).EdgeCone.ClosedIntersectionIndices[i], NewCone.ClosedIntersectionIndices[i]);
			};
			NewCones.push_back(TempCone);
		};
	};
	
	return NewCones;
}

//------------------------------------------------------------------------------
inline vector<Cone> DynamicEnumerate(Cone &C, vector<Hull> &Hulls) {
	// Figure out which polytope we want to pick
	int SmallestInt = 10000000; // Lazy.
	int SmallestIndex = -1;
	for (size_t i = 0; i != C.ClosedIntersectionIndices.size(); i++) {
		if ((C.PolytopesVisited[i] == 0)
		&& (C.ClosedIntersectionIndices[i].size() < SmallestInt)) {
			SmallestInt = C.ClosedIntersectionIndices[i].size();
			SmallestIndex = i;
		};
	};
	for (size_t i = 0; i != C.PolytopesVisited.size(); i++) {
		if (C.PolytopesVisited[i] == 0) {
			SmallestIndex = i;
			break;
		}
	};

	if (SmallestIndex == -1) {
		cout << "Internal error: DynamicEnumerate had a value of -1 for SmallestIndex" << endl;
		cin.get();
	};
	C.PolytopesVisited[SmallestIndex] = 1;
	C.PolytopesVisitedCount++;
	vector<Cone> ResultCones = WalkPolytope(SmallestIndex, C, Hulls);

	for (size_t i = 0; i != ResultCones.size(); i++) {
		ResultCones[i].PolytopesVisited = C.PolytopesVisited;
		ResultCones[i].PolytopesVisitedCount = C.PolytopesVisitedCount;
	};
	return ResultCones;
};

//------------------------------------------------------------------------------
void ParseToPrevariety(vector<int> &ConeRayIndices, int ConeDim, TropicalPrevariety &Output, bool CouldBeMaximal) {
	while (Output.ConeTree.size() < ConeDim) {
		set<ConeWithIndicator> Temp;
		Output.ConeTree.push_back(Temp);
		Output.FVector.push_back(0);
	};

	ConeWithIndicator CWI;
	CWI.IsMaximal = CouldBeMaximal;
	CWI.RayIndices = ConeRayIndices;

	set<ConeWithIndicator>::iterator CWIIt = find(Output.ConeTree[ConeDim - 1].begin(), Output.ConeTree[ConeDim - 1].end(), CWI);
	if (CWIIt == Output.ConeTree[ConeDim - 1].end()) {
		Output.ConeTree[ConeDim - 1].insert(CWI);
	} else {
		if ((*CWIIt).IsMaximal) {
			ConeWithIndicator CWITemp;
			CWITemp.RayIndices = ConeRayIndices;
			CWITemp.IsMaximal = false;
			Output.ConeTree[ConeDim - 1].erase(*CWIIt);
			Output.ConeTree[ConeDim - 1].insert(CWITemp);
		};
		if (ConeDim == 3) {
		
	//	cout << "AAA" << endl;
	//	PrintPoint(ConeRayIndices);
};
		return;
	};
	if (ConeDim == ConeRayIndices.size()) {
		// This is the easy case. We can recurse exactly as we would want to.
		if (CWIIt == Output.ConeTree[ConeDim - 1].end()) {
			if (ConeDim == 1)
				return;
			for (size_t i = 0; i != ConeRayIndices.size(); i++) {
				vector<int> ToRecurse;
				for (size_t j = 0; j != ConeRayIndices.size(); j++) {
					if (i != j)
						ToRecurse.push_back(ConeRayIndices[j]);
				};
				ParseToPrevariety(ToRecurse, ConeDim - 1, Output, false);
			};
		}
	} else {
/*		vector<Generator> Gens;
		Generator_System gs;
		for (size_t i = 0; i != ConeRayIndices.size(); i++) {
			Gens.push_back(Output.IndexToGenMap[ConeRayIndices[i]]);
			gs.insert(Output.IndexToGenMap[ConeRayIndices[i]]);
		};
		Linear_Expression LE;
		for (size_t i = 0; i != Gens[0].space_dimension(); i++) {
			LE += Variable(i) * 0;
		};
		Generator ZeroPoint = point(LE); // Is this necessary? Could this just be point(0)?
		gs.insert(ZeroPoint);
		NNC_Polyhedron InitialPoly(gs);
		
		vector<Constraint> Inequalities;
		Constraint_System Equalities;
		for (Constraint_System::const_iterator i = InitialPoly.constraints().begin(),
		cs_end = InitialPoly.constraints().end(); i != cs_end; ++i) {
			if ((*i).is_equality()) {
				Equalities.insert(*i);
			} else {
				Inequalities.push_back(*i);
			};
		};
		for (size_t i = 0; i != Inequalities.size(); i++) {
			Constraint_System cs = Equalities;
			for (size_t j = 0; j != Inequalities.size(); j++) {
				if (i != j)
					cs.insert(Inequalities[j]);
			};
			NNC_Polyhedron TempPoly(cs);
			vector<int> TempConeRayIndices;
			for (Generator_System::const_iterator gsi = TempPoly.generators().begin(), gs_end = TempPoly.generators().end(); gsi != gs_end; ++gsi) {
				if ((*gsi).is_point() or (*gsi).is_closure_point())
					continue;
				// TODO!!!! If a value isn't defined, the map outputs 0 instead of throwing an error!
				TempConeRayIndices.push_back(Output.RayToIndexMap[GeneratorToPoint(*gsi)]);
			};
			sort(TempConeRayIndices.begin(), TempConeRayIndices.end());
			ParseToPrevariety(TempConeRayIndices, TempPoly.affine_dimension(), Output, false);
		};*/
	};
};

//------------------------------------------------------------------------------
void ENUMERATETEST(vector<Hull> Hulls, int ProcessID, int ProcessCount, vector<ThreadJob> &ThreadJobs, vector<int> &BoredProcesses, TropicalPrevariety &Output, mutex &BPmtx, mutex &Outputmtx) {
	vector<vector<int> > Pretropisms;
	Cone C;
	bool HasCone = false;
	ThreadJob *TJ;
	while (true) {
		int j = ProcessID;
		while (not HasCone) {
			TJ = &ThreadJobs[j % ProcessCount];
			(*TJ).M.lock();
			// If work stealing happens, we really want to steal the less traveled cones,
			// not the more traveled cones. Need to rework that...somehow.
			for (size_t i = (*TJ).SharedCones.size() - 1; i !=-1; i--) {
				if ((*TJ).SharedCones[i].size() > 0) {
					C = (*TJ).SharedCones[i].back();
					(*TJ).SharedCones[i].pop_back();
					HasCone = true;
					// This case happens when the process tried to steal from every process
					// possible but did not succeed, which made it bored. It added itself
					// to BoredProcesses, but then continued looking and found a cone it
					// could have. This made it no longer a bored process.
					if (j > (ProcessID + ProcessCount)) {
						BPmtx.lock();
						BoredProcesses[ProcessID] = 0;
						BoredProcesses[ProcessCount] -= 1;
						BPmtx.unlock();
					};
					break;
				};
			};
			(*TJ).M.unlock();
			if (HasCone)
				break;
			// This case means that we spun through all of the other threads' queues
			// and they were all empty. Now we need to make the process bored.
			if (j == ProcessCount + ProcessID) {
				BPmtx.lock();
				BoredProcesses[ProcessID] = 0;
				BoredProcesses[ProcessCount] += 1;
				BPmtx.unlock();
			};
			if ((j % ProcessCount) == ProcessID) {
				if (BoredProcesses[ProcessCount] == ProcessCount)
					return;
			};
			j++;
		};
		vector<Cone> ResultCones = DynamicEnumerate(C, Hulls);
		
		// If there are remaining new cones, give them to the job queue.
		if (ResultCones.size() > 0) {
			int Index = -1;
			for (size_t i = 0; i != ResultCones[0].PolytopesVisited.size(); i++) {
				if (ResultCones[0].PolytopesVisited[i] == 1)
					Index++;
			};
			if (ResultCones[0].PolytopesVisitedCount == Hulls.size() - 1) {
				Outputmtx.lock();
				clock_t TimeStart = clock();
				for (size_t i = 0; i != ResultCones.size(); i++) {
					vector<int> ConeRayIndices;
					int ConeDim = ResultCones[i].HOPolyhedron.affine_dimension();
					Generator_System pregs = ResultCones[i].HOPolyhedron.generators();
					Generator_System gs;
					gs.insert(point(Variable(0)*0));
					for (Generator_System::const_iterator gsi = pregs.begin(), gs_end = pregs.end(); gsi != gs_end; ++gsi) {
						if ((*gsi).is_point() or (*gsi).is_closure_point())
							continue;
						gs.insert(*gsi);
					};
					gs = NNC_Polyhedron(gs).minimized_generators();
					for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
						if ((*gsi).is_point() or (*gsi).is_closure_point() or (*gsi).is_line())
							continue;							
						vector<int> Ray = GeneratorToPoint(*gsi);
						map<vector<int>, int>::iterator GSIt;
						GSIt = Output.RayToIndexMap.find(Ray);
						if (GSIt == Output.RayToIndexMap.end()) {
							int NewIndex = Output.Pretropisms.size();
							ConeRayIndices.push_back(NewIndex);
							Output.RayToIndexMap[Ray] = NewIndex;
							Output.IndexToGenMap[NewIndex] = (*gsi);
							Output.Pretropisms.push_back(Ray);
							cout << (*gsi) << endl;
						} else {
							ConeRayIndices.push_back((*GSIt).second);
						};
					};
					sort(ConeRayIndices.begin(), ConeRayIndices.end());
					
/*					if (ResultCones[i].ClosedPolyhedron.affine_dimension() == 3) {
						cout << ResultCones[i].ClosedPolyhedron.affine_dimension() << endl;
						cout << ResultCones[i].ClosedPolyhedron.generators() << endl;
						cout << ResultCones[i].ClosedPolyhedron.constraints() << endl << endl;
						PrintPoint(ConeRayIndices);
					};
*/
					ParseToPrevariety(ConeRayIndices, ConeDim, Output, true);
				};
				CleanupTime += double(clock() - TimeStart);
				Outputmtx.unlock();
				HasCone = false;
			} else {
				// If there is at least one cone, hold onto it for the next round.
				C = ResultCones.back();
				ResultCones.pop_back();
				if (ResultCones.size() > 0) {
					TJ = &ThreadJobs[ProcessID];
					(*TJ).M.lock();
					for (size_t i = 0; i != ResultCones.size(); i++) {
						(*TJ).SharedCones[Index].push_back(ResultCones[i]);
					};
					(*TJ).M.unlock();
				};
			};
		} else
			HasCone = false;
	};
};

//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	clock_t StartTime = clock();
	if (argc != 4) {
		cout << "Internal error: expected three arguments." << endl;
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
	
	double RandomSeed = time(NULL);
	cout <<fixed<< "Random seed value: " << RandomSeed << endl;
//	RandomSeed = 1479237790.000000; // FAVORITE 1479237711.000000; // 1478813622.000000 ALSO 1479237711.000000 ALSO 1479237790.000000 Reducedcyclicn 8 weird f vector? // 1477423487.000000 Reducedcyclicn 9 wrong # pretrop (seems to work now!)
	srand(RandomSeed);
	
	vector<Hull> Hulls;
	vector<double> VectorForOrientation;
	for (size_t i = 0; i != PolynomialSystemSupport[0][0].size(); i++) {
		double dd = rand();
		VectorForOrientation.push_back(dd);
	};
	
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++) {
		Hulls.push_back(NewHull(PolynomialSystemSupport[i], VectorForOrientation));
	}
	double HullTime = double(clock() - StartTime);

	double PreintersectTime = 0;

	for(size_t i = 0; i != Hulls.size(); i++){
		for(size_t j = 0; j != Hulls[i].Edges.size(); j++){
			vector<set<int> > InitialSet1 (Hulls.size());
			vector<set<int> > InitialSet2 (Hulls.size());
			Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices = InitialSet1;
		};
	};
	
	clock_t PreintTimeStart = clock();
	int TotalInt = 0;
	int NonInt = 0;
	// TODO: parallelize this.
	for(int i = 0; i != Hulls.size(); i++){
		vector<Edge> Edges1 = Hulls[i].Edges;
		for (size_t k = 0; k != Hulls[i].Edges.size(); k++) {
			for (size_t l = 0; l != Hulls.size(); l++)
				Hulls[i].Edges[k].EdgeCone.PolytopesVisited.push_back(0);
			Hulls[i].Edges[k].EdgeCone.PolytopesVisited[i] = 1;
			Hulls[i].Edges[k].EdgeCone.PolytopesVisitedCount = 0;
		};

		for(size_t j = i+1; j != Hulls.size(); j++){
			vector<Edge> Edges2 = Hulls[j].Edges;
			for(size_t k = 0; k != Edges1.size(); k++){
				for(size_t l = 0; l != Edges2.size(); l++){
					if (IntersectCones(Edges1[k].EdgeCone.HOPolyhedron, Edges2[l].EdgeCone.HOPolyhedron).affine_dimension() >= 1) {
//if (true) {
						// Intersect pairs of half open cones. Osserman/Payne does not apply.
						Hulls[i].Edges[k].EdgeCone.ClosedIntersectionIndices[j].insert(l);
						Hulls[j].Edges[l].EdgeCone.ClosedIntersectionIndices[i].insert(k);
						TotalInt++;
					} else {
						NonInt++;
						TotalInt++;
					};
				};
			};
		};
		printf("Finished level %d of pre-intersections.\n", i);
	};
	cout << "Total Intersections: " << TotalInt << ", Non Intersections: " << NonInt << endl;
	PreintersectTime = double(clock() - PreintTimeStart);
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	
	
	
	
/*
	Constraint_System cz;
	cz.insert(Variable(2) - Variable(6) == 0);
	cz.insert(Variable(1) - Variable(5) == 0);
	cz.insert(Variable(0) - Variable(4) == 0);
	cz.insert(Variable(3) == 0);
	cz.insert(-1*Variable(0) > 0);
	cz.insert(Variable(0) - Variable(2) > 0);
	cz.insert(-1*Variable(1) + Variable(2) > 0);
	
	
	NNC_Polyhedron JEFF(cz);
	
//THIS ONE!C - G = 0, B - F = 0, A - E = 0, D = 0, -A > 0, A - C > 0, -B + C > 0
//C - G = 0, B - F = 0, A - E = 0, D = 0, B - C > 0, A >= 0, -B > 0

	
	
	
	vector<Cone> ConeVector;
	
	
		// Start by initializing the objects.
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
	for (size_t i = 0; i != Cones.size(); i++) {
		cout << "Level " << i << endl;
		vector<Cone> ResultCones;
		for (size_t j = 0; j != ConeVector.size(); j++) {
			vector<Cone> Temp = WalkPolytope(i+1, ConeVector[j], Hulls);
			for (size_t k = 0; k != Temp.size(); k++) {
				ResultCones.push_back(Temp[k]);
			};
		};
		ConeVector = ResultCones;
	};

	int ThreeCount = 0;
	vector<Cone> TestCones;
	vector<int> OutVec;
	for (size_t i = 0; i != ConeVector.size(); i++) {
		int Dim = ConeVector[i].HOPolyhedron.affine_dimension();
		while (OutVec.size() <= Dim) {
			OutVec.push_back(0);
		};
		OutVec[Dim] += 1;
		
	};
	cout << "THREECOUNT " << ThreeCount << endl;
	PrintPoint(OutVec);
	cin.get();
	


	int SmallestInta = 1000000;
	int SmallestIndexa = -1;
	for (size_t i = 0; i != Hulls.size(); i++) {
		int TestValue = 0;
		for (size_t j = 0; j != Hulls[i].Edges.size(); j++) {
			for (size_t k = 0; k != Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices.size(); k++) {
				TestValue += Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices[k].size();
			};
		};
		if (TestValue < SmallestInta) {
			SmallestInta = TestValue;
			SmallestIndexa = i;
		};
	};
	if (SmallestIndexa == -1) {
		cout << "Internal error: value of -1 for SmallestIndex" << endl;
		return 1;
	};
*/	



/*
	vector<Cone> StartCones;
	for (size_t i = 0; i != Hulls[SmallestIndexa].Edges.size(); i++) {
		StartCones.push_back(Hulls[SmallestIndexa].Edges[i].EdgeCone);
	};

	for (size_t i = 0; i != Hulls.size() - 1; i++) {
		vector<Cone> TempCones;
		for(size_t j = 0; j != StartCones.size(); j++) {
			bool THISONE = (IntersectCones(JEFF, StartCones[j].HOPolyhedron).affine_dimension() > 0);
			if (THISONE) {
				cout << StartCones[j].ClosedPolyhedron.minimized_constraints() << endl;
				cout << StartCones[j].HOPolyhedron.minimized_constraints() << endl;
				
				
				PrintPoint(StartCones[j].PolytopesVisited);
				for(size_t as = 0; as != StartCones[j].ClosedIntersectionIndices.size(); as++)
					PrintPoint(StartCones[j].ClosedIntersectionIndices[as]);
				cout << "AAAA " << j << endl;
				cout << StartCones[j].HOPolyhedron.minimized_generators() << endl;
			};
			vector<Cone> OutputCones = DynamicEnumerate(StartCones[j], Hulls);
			if (THISONE) {
			cout << OutputCones.size() << endl;
			for (size_t qq = 0; qq != OutputCones.size(); qq++) {
				cout << OutputCones[qq].ClosedPolyhedron.affine_dimension() << endl;
				cout << OutputCones[qq].ClosedPolyhedron.minimized_constraints() << endl;
			};
			
			cout << endl;
			cin.get();
			};
			for (size_t k = 0; k != OutputCones.size(); k++)
				TempCones.push_back(OutputCones[k]);
		}
		StartCones = TempCones;
	};

	int ThreeCount2 = 0;
	vector<Cone> TestCones2;
	vector<int> OutVec2;
	for (size_t i = 0; i != StartCones.size(); i++) {
		int Dim = StartCones[i].ClosedPolyhedron.affine_dimension();
		while (OutVec2.size() <= Dim) {
			OutVec2.push_back(0);
		};
		OutVec2[Dim] += 1;
		if (StartCones[i].ClosedPolyhedron.affine_dimension() == 3) {
		cout << StartCones[i].ClosedPolyhedron.minimized_constraints() << endl;
		ThreeCount2++;
		TestCones2.push_back(StartCones[i]);
		};
	};
	for (size_t i = 0; i != TestCones2.size(); i++) {
		for (size_t j = i+1; j != TestCones2.size(); j++) {
			int Dimm = IntersectCones(TestCones2[i].ClosedPolyhedron, TestCones2[j].ClosedPolyhedron).affine_dimension();
			if (Dimm == 3) 
				cout << "AAAA" << endl;
			
		};
	};
	cout << "THREECOUNT " << ThreeCount2 << endl;
	PrintPoint(OutVec2);
	cin.get();
	/*
	for (size_t i = 0; i != TestCones.size(); i++) {
		bool FoundMatch = false;
		for (size_t j = 0; j != StartCones.size(); j++) {
			if (IntersectCones(TestCones[i].ClosedPolyhedron, StartCones[j].ClosedPolyhedron).affine_dimension() == 3) {
				FoundMatch = true;
				break;
			};
		};	
		if (FoundMatch == false) {
			cout << i << endl;
			cout << TestCones[i].HOPolyhedron.minimized_generators() << endl;
			cout << TestCones[i].HOPolyhedron.minimized_constraints() << endl;
		};
	};
	cin.get();
	*/
	/*
	vector<Cone> ConeVector;
		// Start by initializing the objects.
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
			for (int i = 0; i != TestCones.size(); i++) {
				if (find(Cone1.ClosedIntersectionIndices[TreeLevel].begin(), Cone1.ClosedIntersectionIndices[TreeLevel].end(), i) == Cone1.ClosedIntersectionIndices[TreeLevel].end()) continue;
				Cone NewCone;
				NewCone.ClosedPolyhedron = IntersectCones(TestCones[i].ClosedPolyhedron, Cone1.ClosedPolyhedron);
				ConeIntersectionCount++;
				int ExpectedDim = Cone1.ClosedPolyhedron.affine_dimension() - 1;
				if (NewCone.ClosedPolyhedron.affine_dimension() == 0) {
					continue;
				};
				NewCone.HOPolyhedron = IntersectCones(TestCones[i].HOPolyhedron, Cone1.HOPolyhedron);
				if (NewCone.HOPolyhedron.affine_dimension() == 0) {
					continue;
				};
				
				for (size_t j = 0; j != Cone1.ClosedIntersectionIndices.size(); j++) {
					NewCone.ClosedIntersectionIndices.push_back(IntersectSets(Cone1.ClosedIntersectionIndices[j], TestCones[i].ClosedIntersectionIndices[j]));
				};
				NewCones.push_back(NewCone);
			};
		};
		ConeVector = NewCones;
		printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", TreeLevel, Cones.size(), ConeVector.size(),ConeIntersectionCount);
		TreeLevel++;
	};
	vector<vector<int> > PretropismsTemp;
	int ThreeDimCount = 0;
	for (size_t i = 0; i != ConeVector.size(); i++) {
		Generator_System gs = ConeVector[i].ClosedPolyhedron.minimized_generators();
		cout << ConeVector[i].ClosedPolyhedron.affine_dimension() << endl;
		cout << ConeVector[i].ClosedPolyhedron.minimized_constraints() << endl;
		cout << gs << endl;
		cout << "TEST" << endl;
		if (ConeVector[i].ClosedPolyhedron.affine_dimension() == 3) ThreeDimCount++;
		for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
			if ((*gsi).is_point() or (*gsi).is_closure_point()) {
				continue;
			};
			vector<int> Pretrop = GeneratorToPoint(*gsi);
			if (find(PretropismsTemp.begin(), PretropismsTemp.end(), Pretrop) == PretropismsTemp.end() ) {
				PretropismsTemp.push_back(Pretrop);
			};
		};
	};
	//sort(PretropismsTemp.begin(), PretropismsTemp.end());

	//PrintPoints(PretropismsTemp);
	cout << "Number of pretropisms found: " << PretropismsTemp.size() << endl;

	
	cout << "Three dim count: " << ThreeDimCount << endl;
	cin.get();
	
	
	
	*/
	
	
	
	
	
	
	
	
	
	
	// This is one way to do it. It seems reasonable to pick sum, median, mean, min...
	int SmallestInt = 1000000;
	int SmallestIndex = -1;
	for (size_t i = 0; i != Hulls.size(); i++) {
		int TestValue = 0;
		for (size_t j = 0; j != Hulls[i].Edges.size(); j++) {
			for (size_t k = 0; k != Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices.size(); k++) {
				TestValue += Hulls[i].Edges[j].EdgeCone.ClosedIntersectionIndices[k].size();
			};
		};
		if (TestValue < SmallestInt) {
			SmallestInt = TestValue;
			SmallestIndex = i;
		};
	};
	if (SmallestIndex == -1) {
		cout << "Internal error: value of -1 for SmallestIndex" << endl;
		return 1;
	};
	SmallestIndex = 0;
	int TotalProcessCount = atoi(argv[3]);
	if (TotalProcessCount > thread::hardware_concurrency()) {
		cout << "Internal error: hardware_concurrency = " << thread::hardware_concurrency() << " but TotalProcessCount = " << TotalProcessCount << endl;
		return 1;
	};
	
	vector<ThreadJob> ThreadJobs;
	for (size_t i = 0; i != TotalProcessCount; i++) {
		vector<vector<Cone> > SharedCones;
		for (size_t i = 0; i != Hulls.size() - 1; i++) {
			vector<Cone> Temp;
			SharedCones.push_back(Temp);
		};
		ThreadJob TJ(SharedCones);
		ThreadJobs.push_back(TJ);
	};

	for (size_t i = 0; i != Hulls[SmallestIndex].Edges.size(); i++) {
		ThreadJobs[i % TotalProcessCount].SharedCones[0].push_back(Hulls[SmallestIndex].Edges[i].EdgeCone);
	};
	
	mutex BPmtx;
	mutex Outputmtx;
	TropicalPrevariety Output;
	vector<int> BoredProcesses (TotalProcessCount + 1, 0);
	clock_t AlgorithmStartTime = clock();
	typedef function<void()> work_type;
	Thread_Pool<work_type> thread_pool(TotalProcessCount);
	for (size_t i = 0; i != TotalProcessCount; i++) {
		work_type work = bind(ENUMERATETEST, Hulls, i, TotalProcessCount, ref(ThreadJobs), ref(BoredProcesses), ref(Output), ref(BPmtx), ref(Outputmtx));
		thread_pool.submit(make_threadable(work));
	}
	// Wait for all workers to complete.
	thread_pool.finalize();
	
	while (true) {
		if (BoredProcesses[TotalProcessCount] == TotalProcessCount)
			break;
		else
			this_thread::sleep_for(chrono::milliseconds(10));
	};
	vector<int> TotalVec;
	for (size_t i = 0; i != Output.ConeTree.size(); i++) {
		TotalVec.push_back(Output.ConeTree[i].size());
		set<ConeWithIndicator>::iterator CWIIt;
		for (CWIIt=Output.ConeTree[i].begin(); CWIIt != Output.ConeTree[i].end(); CWIIt++) {
							
			if ((*CWIIt).IsMaximal) {
				Output.FVector[i]++;

			} else {
//			PrintPoint((*CWIIt).RayIndices);
			};
		};
	};
	cout << "MAXIMAL CONES------------------" << endl;
		for (size_t i = 0; i != Output.ConeTree.size(); i++) {
		cout << i << " " << Output.ConeTree[i].size() <<endl;
		set<ConeWithIndicator>::iterator CWIIt;
		for (CWIIt=Output.ConeTree[i].begin(); CWIIt != Output.ConeTree[i].end(); CWIIt++) {
			PrintPoint((*CWIIt).RayIndices);
							
			if ((*CWIIt).IsMaximal) {

			} else {
			};
		};
		cout << endl;
	};
//	cin.get();
	for(map<vector<int>,int >::iterator itt = Output.RayToIndexMap.begin(); itt != Output.RayToIndexMap.end(); ++itt) {
		//PrintPoint(itt->first);
	//	cout << itt->second << "\n";
}
	
	
	sort(Output.Pretropisms.begin(), Output.Pretropisms.end());
	PrintPoints(Output.Pretropisms);
	
	cout << "Maximal cone count--------" << endl;
	PrintPoint(Output.FVector);
	PrintPoint(TotalVec);
	cout << "Number of pretropisms found: " << Output.Pretropisms.size() << endl;
	cout << "Cleanup time: " << CleanupTime / CLOCKS_PER_SEC << endl;
	cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	cout << "Test time: " << TestTime / CLOCKS_PER_SEC << endl;
//	cout << "Time after init: " << TimeAfterInit / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;

}
