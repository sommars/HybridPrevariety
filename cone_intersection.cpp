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

double IntersectionTime, TestTime;
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
	vector<int> RandomVector(Hulls[0].SpaceDimension, 0);
	for (Generator_System::const_iterator i = 
	NewCone.HOPolyhedron.generators().begin(), gs_end = 
	NewCone.HOPolyhedron.generators().end(); i != gs_end; ++i) {
		if (!(*i).is_ray()) {
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
	set<int> EdgesToTest;
	for(size_t EdgeIndex = 0; EdgeIndex != (*HIndex).Edges.size(); EdgeIndex++) {
		if (SetsDoIntersect(InitialIndices, (*HIndex).Edges[EdgeIndex].PointIndices)
		&& (find((*ClosedIntersectionIndices).begin(), (*ClosedIntersectionIndices).end(), EdgeIndex) != (*ClosedIntersectionIndices).end()))
			EdgesToTest.insert(EdgeIndex);
	};
	
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
		TempCone.ClosedPolyhedron = (*EdgeToTest).EdgeCone.ClosedPolyhedron;

		TempCone.ClosedPolyhedron.add_constraints(NewCone.HOPolyhedron.constraints());
		TempCone.ClosedPolyhedron.affine_dimension();
		ConeIntersectionCount++;
		IntersectionTime += double(clock() - IntBegin);
		if (TempCone.ClosedPolyhedron.affine_dimension() > 0) {
			PretropGraphEdges.insert(EdgeToTestIndex);

			clock_t IntBegin = clock();
			TempCone.HOPolyhedron = (*EdgeToTest).EdgeCone.HOPolyhedron;
			TempCone.HOPolyhedron.add_constraints(NewCone.HOPolyhedron.constraints());

			TempCone.HOPolyhedron.affine_dimension();
			ConeIntersectionCount++;
			IntersectionTime += double(clock() - IntBegin);


			if (TempCone.HOPolyhedron.affine_dimension() > 0) {
				clock_t IntBegin = clock();
				IntersectionTime += double(clock() - IntBegin);
				
				vector<set<int> > InitialSet1 (NewCone.ClosedIntersectionIndices.size());
				TempCone.ClosedIntersectionIndices = InitialSet1;
				for (size_t i = 0; i != NewCone.ClosedIntersectionIndices.size(); i++) {
					if (NewCone.PolytopesVisited[i] == 0)
						TempCone.ClosedIntersectionIndices[i] = IntersectSets((*EdgeToTest).EdgeCone.ClosedIntersectionIndices[i], NewCone.ClosedIntersectionIndices[i]);
				};
				TempCone.ConeHistory = NewCone.ConeHistory;
				TempCone.ConeHistory.push_back((*EdgeToTest).EdgeCone.HOPolyhedron);
				NewCones.push_back(TempCone);
			};
			set<int>::iterator NeighborItr;
			for(NeighborItr=(*EdgeToTest).NeighborIndices.begin(); NeighborItr!=(*EdgeToTest).NeighborIndices.end(); NeighborItr++) {
				int Neighbor = *NeighborItr;
				if (( find((*ClosedIntersectionIndices).begin(), (*ClosedIntersectionIndices).end(), Neighbor) != (*ClosedIntersectionIndices).end())
				&& ( find(PretropGraphEdges.begin(), PretropGraphEdges.end(), Neighbor) == PretropGraphEdges.end() )
				&& ( find(NotPretropGraphEdges.begin(), NotPretropGraphEdges.end(), Neighbor) == NotPretropGraphEdges.end() )
				&& ( find(EdgesToTest.begin(), EdgesToTest.end(), Neighbor) == EdgesToTest.end() )) {
					EdgesToTest.insert(Neighbor);
				};
			};
		} else {
			NotPretropGraphEdges.insert(EdgeToTestIndex);
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
	if (SmallestIndex == -1) {
		cout << "Internal error: DynamicEnumerate had a value of -1 for SmallestIndex" << endl;
		cin.get();
	};
	C.PolytopesVisited[SmallestIndex] = 1;
	vector<Cone> ResultCones = WalkPolytope(SmallestIndex, C, Hulls);

	for (size_t i = 0; i != ResultCones.size(); i++)
		ResultCones[i].PolytopesVisited = C.PolytopesVisited;
	return ResultCones;
};

//------------------------------------------------------------------------------
void ParseToPrevariety(vector<int> &ConeRayIndices, int ConeDim, TropicalPrevariety &Output, bool CouldBeMaximal) {
	if (ConeDim == ConeRayIndices.size()) {
		// This is the easy case. We can recurse exactly as we would want to.
		while (Output.ConeTree.size() < ConeDim) {
			set<ConeWithIndicator > Temp;
			Output.ConeTree.push_back(Temp);
			Output.FVector.push_back(0);
		};
		ConeWithIndicator CWI;
		CWI.IsMaximal = CouldBeMaximal;
		CWI.RayIndices = ConeRayIndices;
		set<ConeWithIndicator>::iterator CWIIt = find(Output.ConeTree[ConeDim - 1].begin(), Output.ConeTree[ConeDim - 1].end(), CWI);
		if (CWIIt == Output.ConeTree[ConeDim - 1].end()) {
			Output.ConeTree[ConeDim - 1].insert(CWI);
			
			if (ConeDim == 1)
				return;	
			for (size_t i = 0; i != ConeRayIndices.size(); i++) {
				// There has to be a better way to do this...
				vector<int> ToRecurse;
				for (size_t j = 0; j != ConeRayIndices.size(); j++) {
					if (i != j)
						ToRecurse.push_back(ConeRayIndices[j]);
				};
				ParseToPrevariety(ToRecurse, ConeDim - 1, Output, false);
			};
		} else {
			if ((*CWIIt).IsMaximal) {
				ConeWithIndicator CWITemp;
				CWITemp.RayIndices = ConeRayIndices;
				CWITemp.IsMaximal = false;
				Output.ConeTree[ConeDim - 1].erase(*CWIIt);
				Output.ConeTree[ConeDim - 1].insert(CWITemp);
			};
			return;
		}
	} else {
		// TODO!!
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
			if (Index == Hulls.size() - 1) {
				Outputmtx.lock();
				for (size_t i = 0; i != ResultCones.size(); i++) {
					vector<int> ConeRayIndices;
					int ConeDim = ResultCones[i].HOPolyhedron.affine_dimension();
					Generator_System gs = ResultCones[i].HOPolyhedron.generators();
					for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
						if ((*gsi).is_point() or (*gsi).is_closure_point())
							continue;
						vector<int> Ray = GeneratorToPoint(*gsi);
						map<vector<int>, int>::iterator GSIt;
						GSIt = Output.RayToIndexMap.find(Ray);
						if (GSIt == Output.RayToIndexMap.end()) {
							int NewIndex = Output.RayToIndexMap.size();
							ConeRayIndices.push_back(NewIndex);
							Output.RayToIndexMap[Ray] = NewIndex;
							Output.Pretropisms.push_back(Ray);
						} else {
							ConeRayIndices.push_back((*GSIt).second);
						};
					};
					sort(ConeRayIndices.begin(), ConeRayIndices.end());
					ParseToPrevariety(ConeRayIndices, ConeDim, Output, true);
				};
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
	srand(RandomSeed);
	
	vector<Hull> Hulls;
	vector<double> VectorForOrientation;
	for (size_t i = 0; i != PolynomialSystemSupport[0][0].size(); i++) {
		double dd = rand();
		VectorForOrientation.push_back(dd);
	};
	
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++)
		Hulls.push_back(NewHull(PolynomialSystemSupport[i], VectorForOrientation));
	
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
		int ExpectedDim = Hulls[0].Edges[0].EdgeCone.ClosedPolyhedron.affine_dimension() - 1;
		vector<Edge> Edges1 = Hulls[i].Edges;
		for (size_t k = 0; k != Hulls[i].Edges.size(); k++) {
			for (size_t l = 0; l != Hulls.size(); l++)
				Hulls[i].Edges[k].EdgeCone.PolytopesVisited.push_back(0);
			Hulls[i].Edges[k].EdgeCone.PolytopesVisited[i] = 1;
		};

		for(size_t j = i+1; j != Hulls.size(); j++){
			vector<Edge> Edges2 = Hulls[j].Edges;
			for(size_t k = 0; k != Edges1.size(); k++){
				for(size_t l = 0; l != Edges2.size(); l++){
					if (IntersectCones(Edges1[k].EdgeCone.ClosedPolyhedron, Edges2[l].EdgeCone.ClosedPolyhedron).affine_dimension() >= ExpectedDim) {
						// Intersect pairs of closed cones. Osserman/Payne applies
						Hulls[i].Edges[k].EdgeCone.ClosedIntersectionIndices[j].insert(l);
						Hulls[j].Edges[l].EdgeCone.ClosedIntersectionIndices[i].insert(k);
						TotalInt++;
					} else if (IntersectCones(Edges1[k].EdgeCone.HOPolyhedron, Edges2[l].EdgeCone.HOPolyhedron).affine_dimension() >= 1) {
						// Intersect pairs of half open cones. Osserman/Payne does not apply.
						Hulls[i].Edges[k].EdgeCone.ClosedIntersectionIndices[j].insert(l);
						Hulls[j].Edges[l].EdgeCone.ClosedIntersectionIndices[i].insert(k);
						TotalInt++;
						TotalInt++;
					} else {
						NonInt++;
						TotalInt++;
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
	//Hulls[SmallestIndex].Edges[0].EdgeCone.HOPolyhedron = Hulls[SmallestIndex].Edges[0].EdgeCone.ClosedPolyhedron;
	//SharedCones[0].push_back(Hulls[SmallestIndex].Edges[0].EdgeCone);
	
	
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
			this_thread::sleep_for(chrono::milliseconds(100));
	};
	
	for (size_t i = 0; i != Output.ConeTree.size(); i++) {
		
		set<ConeWithIndicator>::iterator CWIIt;
		for (CWIIt=Output.ConeTree[i].begin(); CWIIt != Output.ConeTree[i].end(); CWIIt++) {
			if ((*CWIIt).IsMaximal)
				Output.FVector[i]++;
		};
	};
	sort(Output.Pretropisms.begin(), Output.Pretropisms.end());
	PrintPoints(Output.Pretropisms);
	cout << "Maximal cone count--------" << endl;
	PrintPoint(Output.FVector);
	cout << "Number of pretropisms found: " << Output.Pretropisms.size() << endl;
//	cout << "Cleanup time: " << double(clock() - CleanupStart) / CLOCKS_PER_SEC << endl;
	cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	cout << "Test time: " << TestTime / CLOCKS_PER_SEC << endl;
//	cout << "Time after init: " << TimeAfterInit / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;

}
