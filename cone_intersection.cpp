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

double IntersectionTime, SetIntersectTime, CleanupTime, AffineDimTime, ConeCopyTime;
int ConeIntersectionCount;

//------------------------------------------------------------------------------
inline list<Cone> DoCommonRefinement(int HullIndex, Cone &NewCone, vector<vector<Cone> > &HullCones) {
	vector<Cone> *HIndex;
	HIndex = &HullCones[HullIndex];
	BitsetWithCount *RT = &NewCone.RelationTables[HullIndex];
	list<Cone> Result;
	Cone *ConeToTest;
	for (boost::dynamic_bitset<>::size_type j = 0; j != (*RT).Indices.size(); j++) {
		if (!(*RT).Indices[j])
			continue;
		ConeToTest = &(*HIndex)[j];

		
		clock_t IntBegin = clock();
		C_Polyhedron HOPolyhedron = IntersectCones((*ConeToTest).HOPolyhedron, NewCone.HOPolyhedron);
		IntersectionTime += double(clock() - IntBegin);
		
		clock_t IntBegin2 = clock();
		HOPolyhedron.affine_dimension();
		AffineDimTime += double(clock() - IntBegin2);
		ConeIntersectionCount++;
		if (HOPolyhedron.affine_dimension() > 0) {
			Cone TestCone;
			TestCone.HOPolyhedron = HOPolyhedron;
			TestCone.PolytopesVisited.Count = NewCone.PolytopesVisited.Count;
			TestCone.RelationTables.resize(HullCones.size());
			TestCone.PolytopesVisited = NewCone.PolytopesVisited;
			
			
			for (size_t i = 0; i != NewCone.RelationTables.size(); i++) {
				if (!NewCone.PolytopesVisited.Indices[i]) {
					clock_t SetIntersectBegin = clock();
					TestCone.RelationTables[i] = IntersectRTs((*ConeToTest).RelationTables[i], NewCone.RelationTables[i]);
					SetIntersectTime += double(clock() - SetIntersectBegin);
				}
			};
			
			Result.push_back(TestCone);
		};
	};
	return Result;
}

//------------------------------------------------------------------------------
inline list<Cone> DynamicEnumerate(Cone &C, vector<vector<Cone> > &HullCones) {
	// Figure out which polytope we want to pick
	int SmallestInt = 10000000;
	int SmallestIndex = -1;
	for (size_t i = 0; i != C.RelationTables.size(); i++) {
		if ((!C.PolytopesVisited.Indices[i])
		&& (C.RelationTables[i].Count < SmallestInt)) {
			SmallestInt = C.RelationTables[i].Count;
			SmallestIndex = i;
		};
	};

	if (SmallestIndex == -1) {
		cout << "Internal error: DynamicEnumerate had a value of -1 for SmallestIndex" << endl;
		cin.get();
	};
	C.PolytopesVisited.Indices[SmallestIndex] = 1;
	C.PolytopesVisited.Count++;

	return DoCommonRefinement(SmallestIndex, C, HullCones);
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
void ThreadEnum(vector<vector<Cone> > HullCones, int ProcessID, int ProcessCount, vector<ThreadJob> &ThreadJobs, vector<int> &BoredProcesses, TropicalPrevariety &Output, mutex &BPmtx, mutex &Outputmtx) {
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
					clock_t BEGIN = clock();
					C = (*TJ).SharedCones[i].back();
					ConeCopyTime += double(clock() - BEGIN);
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
		list<Cone> ResultCones = DynamicEnumerate(C, HullCones);

		// If there are remaining new cones, give them to the job queue.
		if (ResultCones.size() > 0) {
			int Index = ResultCones.front().PolytopesVisited.Count;
			
			// The cones have visited all of the polytopes.
			if (Index == HullCones.size()) {
				Outputmtx.lock();
				clock_t TimeStart = clock();
				list<Cone>::iterator i;
				for (i = ResultCones.begin(); i != ResultCones.end(); i++) {
					
					vector<int> ConeRayIndices;
					int ConeDim = (*i).HOPolyhedron.affine_dimension();
					Generator_System pregs = (*i).HOPolyhedron.generators();
					Generator_System gs;
					gs.insert(point(Variable(0)*0));
					for (Generator_System::const_iterator gsi = pregs.begin(), gs_end = pregs.end(); gsi != gs_end; ++gsi) {
						if ((*gsi).is_point() or (*gsi).is_closure_point())
							continue;
						gs.insert(*gsi);
					};
					gs = C_Polyhedron(gs).minimized_generators(); // Should be a better way to do this...
					for (Generator_System::const_iterator gsi = gs.begin(), gs_end = gs.end(); gsi != gs_end; ++gsi) {
						if ((*gsi).is_point() or (*gsi).is_closure_point() or (*gsi).is_line())
							continue;
						vector<int> Ray = GeneratorToPoint(*gsi);
						if (Ray.back() != 0)
							continue;
						Ray.pop_back();
						map<vector<int>, int>::iterator GSIt;
						GSIt = Output.RayToIndexMap.find(Ray);
						if (GSIt == Output.RayToIndexMap.end()) {
							int NewIndex = Output.Pretropisms.size();
							//ConeRayIndices.push_back(NewIndex);
							Output.RayToIndexMap[Ray] = NewIndex;
							//Output.IndexToGenMap[NewIndex] = (*gsi);
							Output.Pretropisms.push_back(Ray);
							//cout << (*gsi) << endl;
						} else {
						//	ConeRayIndices.push_back((*GSIt).second);
						};
					};
					//sort(ConeRayIndices.begin(), ConeRayIndices.end());
					//ParseToPrevariety(ConeRayIndices, ConeDim, Output, true);
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
					list<Cone>::iterator i;
					for (i = ResultCones.begin(); i != ResultCones.end(); i++) {
						(*TJ).SharedCones[Index - 1].push_back(*i);
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
	
	bool Verbose;
	
	int TotalProcessCount;
	vector<vector<vector<int> > > PolynomialSystemSupport;
	if (argc != 4) {
		if (Verbose) {
			cout << "Expected three arguments. Waiting for user input system..." << endl;
		}
		TotalProcessCount = 1;
		string input;
		cin >> input;
		PolynomialSystemSupport = ParseToSupport(input);
		Verbose = false;
	} else {
		int n = atoi(argv[1]);
		string SystemName = string(argv[2]);
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
		TotalProcessCount = atoi(argv[3]);
		if (TotalProcessCount > thread::hardware_concurrency()) {
			cout << "Internal error: hardware_concurrency = " << thread::hardware_concurrency() << " but TotalProcessCount = " << TotalProcessCount << endl;
			return 1;
		};
		Verbose = true;
	};

	double RandomSeed = time(NULL);
	if (Verbose) {
		cout << fixed << "Random seed value: " << RandomSeed << endl;
	};
	srand(RandomSeed);
	
	vector<vector<Cone> > HullCones;
	vector<double> VectorForOrientation;
	for (size_t i = 0; i != PolynomialSystemSupport[0][0].size(); i++) {
		double dd = rand();
		VectorForOrientation.push_back(dd);
	};
	
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++) {
		HullCones.push_back(NewHull(PolynomialSystemSupport[i], VectorForOrientation, Verbose));
	}
	double HullTime = double(clock() - StartTime);

	double PreintersectTime = 0;

	for(size_t i = 0; i != HullCones.size(); i++){
		for(size_t j = 0; j != HullCones[i].size(); j++){
			for(size_t k = 0; k != HullCones.size(); k++) {
				BitsetWithCount RT;
				RT.Indices.resize(HullCones[k].size());
				RT.Count = 0;
				HullCones[i][j].RelationTables.push_back(RT);
			};
		};
	};
	clock_t PreintTimeStart = clock();
	int TotalInt = 0;
	int NonInt = 0;
	// TODO: parallelize this.
	for(int i = 0; i != HullCones.size(); i++){
		vector<Cone> Cones1 = HullCones[i];
		for (size_t k = 0; k != HullCones[i].size(); k++) {
			HullCones[i][k].PolytopesVisited.Indices.resize(HullCones.size());
			HullCones[i][k].PolytopesVisited.Indices[i] = 1;
			HullCones[i][k].PolytopesVisited.Count = 1;
		};
		for(size_t j = i+1; j != HullCones.size(); j++){
			vector<Cone> Cones2 = HullCones[j];
			for(size_t k = 0; k != Cones1.size(); k++){
				for(size_t l = 0; l != Cones2.size(); l++){
					if (IntersectCones(Cones1[k].HOPolyhedron, Cones2[l].HOPolyhedron).affine_dimension() >= 1) {
						HullCones[i][k].RelationTables[j].Indices[l] = 1;
						HullCones[j][l].RelationTables[i].Indices[k] = 1;
						HullCones[i][k].RelationTables[j].Count++;
						HullCones[j][l].RelationTables[i].Count++;
					} else {
						NonInt++;
					};
					TotalInt++;
				};
			};
		};
		if (Verbose) {
			printf("Finished level %d of pre-intersections.\n", i);
		};
	};
	PreintersectTime = double(clock() - PreintTimeStart);
	if (Verbose) {
		cout << "Total Intersections: " << TotalInt << ", Non Intersections: " << NonInt << endl;
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	};
	
	// This is one way to do it. It seems reasonable to pick sum, median, mean, min...
	int SmallestInt = 1000000;
	int SmallestIndex = -1;
	for (size_t i = 0; i != HullCones.size(); i++) {
		int TestValue = 0;
		for (size_t j = 0; j != HullCones[i].size(); j++) {
			for (size_t k = 0; k != HullCones[i][j].RelationTables.size(); k++) {
				TestValue += HullCones[i][j].RelationTables[k].Count;
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
	
	vector<ThreadJob> ThreadJobs;
	for (size_t i = 0; i != TotalProcessCount; i++) {
		vector<list<Cone> > SharedCones;
		for (size_t i = 0; i != HullCones.size() - 1; i++) {
			list<Cone> Temp;
			SharedCones.push_back(Temp);
		};
		ThreadJob TJ(SharedCones);
		ThreadJobs.push_back(TJ);
	};

	for (size_t i = 0; i != HullCones[SmallestIndex].size(); i++) {
		ThreadJobs[i % TotalProcessCount].SharedCones[0].push_back(HullCones[SmallestIndex][i]);
	};
	
	mutex BPmtx;
	mutex Outputmtx;
	TropicalPrevariety Output;
	vector<int> BoredProcesses (TotalProcessCount + 1, 0);
	clock_t AlgorithmStartTime = clock();
	typedef function<void()> work_type;
	Thread_Pool<work_type> thread_pool(TotalProcessCount);
	for (size_t i = 0; i != TotalProcessCount; i++) {
		work_type work = bind(ThreadEnum, HullCones, i, TotalProcessCount, ref(ThreadJobs), ref(BoredProcesses), ref(Output), ref(BPmtx), ref(Outputmtx));
		thread_pool.submit(make_threadable(work));
	}
	// Wait for all workers to complete.
	thread_pool.finalize();
	
	//There has to be a better way to do this...
	while (true) {
		if (BoredProcesses[TotalProcessCount] == TotalProcessCount)
			break;
		else
			this_thread::sleep_for(chrono::milliseconds(10));
	};
	
	cout << "Cone intersection count: " << ConeIntersectionCount << endl;
	
	double AlgorithmTotalTime = double(clock() - AlgorithmStartTime);
	clock_t TimeAfterEndOfAlg = clock();
	vector<int> TotalVec;
	for (size_t i = 0; i != Output.ConeTree.size(); i++) {
		TotalVec.push_back(Output.ConeTree[i].size());
		set<ConeWithIndicator>::iterator CWIIt;
		for (CWIIt=Output.ConeTree[i].begin(); CWIIt != Output.ConeTree[i].end(); CWIIt++) {
			if ((*CWIIt).IsMaximal)
				Output.FVector[i]++;
		};
	};
	if (Verbose) {
		cout << "MAXIMAL CONES------------------" << endl;
	};
	for (size_t i = 0; i != Output.ConeTree.size(); i++) {
		if (Verbose) {
			cout << i << " " << Output.ConeTree[i].size() <<endl;
		};
		set<ConeWithIndicator>::iterator CWIIt;
		for (CWIIt=Output.ConeTree[i].begin(); CWIIt != Output.ConeTree[i].end(); CWIIt++) {
			if (Verbose) {
				PrintPoint((*CWIIt).RayIndices);
			};
			if ((*CWIIt).IsMaximal) {

			} else {
			};
		};
	};	
	
	sort(Output.Pretropisms.begin(), Output.Pretropisms.end());
	if (Verbose) {
		PrintPoints(Output.Pretropisms);
	} else {
		PrintPointsForPython(Output.Pretropisms);
	};
	
	if (Verbose) {
		cout << "Maximal cone count--------" << endl;
		PrintPoint(Output.FVector);
		PrintPoint(TotalVec);
		cout << "Number of pretropisms found: " << Output.Pretropisms.size() << endl;
		cout << "AffineDim time: " << AffineDimTime / CLOCKS_PER_SEC << endl;
		cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
		cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
		cout << "Cleanup time: " << CleanupTime / CLOCKS_PER_SEC << endl;
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
		cout << "Set intersect time: " << SetIntersectTime / CLOCKS_PER_SEC << endl;
		cout << "Number of intersections: " << ConeIntersectionCount << endl;
		cout << "Time after end of alg: " << double(clock() - TimeAfterEndOfAlg) / CLOCKS_PER_SEC << endl;
		cout << "Algorithm total time: " << AlgorithmTotalTime / CLOCKS_PER_SEC << endl;
		cout << "Cone copy time: " << ConeCopyTime / CLOCKS_PER_SEC << endl;
		cout << "Alg time - cone int time - set int time: " << (AlgorithmTotalTime - IntersectionTime - SetIntersectTime) / CLOCKS_PER_SEC << endl;
	};
}
