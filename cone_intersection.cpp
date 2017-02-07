#include "soplex_test.h"

//double IntersectionTime, SetIntersectTime, CleanupTime, AffineDimTime, ConeCopyTime;
int ConeIntersectionCount;
double SoplexTime;
bool UsePPL;

//------------------------------------------------------------------------------
list<Cone> DoCommonRefinement(int HullIndex, Cone &NewCone, vector<vector<Cone> > &HullCones, MySoPlex &MySop) {
	vector<Cone> *HIndex;
	HIndex = &HullCones[HullIndex];
	BitsetWithCount *RT = &NewCone.RelationTables[HullIndex];
	list<Cone> Result;
	Cone *ConeToTest;
	for (boost::dynamic_bitset<>::size_type j = 0; j != RT->Indices.size(); j++) {
		if (!RT->Indices[j])
			continue;
		ConeToTest = &(*HIndex)[j];

		bool SkipIntersection = false;
		vector<BitsetWithCount> RelationTables(HullCones.size());
		for (size_t i = 0; i != NewCone.RelationTables.size(); i++) {
			if (!NewCone.PolytopesVisited.Indices[i]) {
				RelationTables[i] = IntersectRTs(ConeToTest->RelationTables[i], NewCone.RelationTables[i]);
				if (RelationTables[i].Count == 0) {
					SkipIntersection = true;
					break;
				};
			};
		};
		
		if (SkipIntersection)
			continue;
		ConeIntersectionCount++;
		Cone TestCone;
		bool ShouldContinue;
		if (UsePPL) {
			TestCone.HOPolyhedron = NewCone.HOPolyhedron;
			TestCone.HOPolyhedron.add_constraints(ConeToTest->Constraints);
			ShouldContinue = TestCone.HOPolyhedron.is_empty();
		} else {
			MySop.ClearLP();
		//	for (size_t i = 0; i != NewCone.PolytopesVisitedIndices.size(); i++) {
		//		MySop.AddRowsToLP(HullCones[NewCone.PolytopesVisitedIndices[i]][NewCone.ConesVisitedIndices[i]].Rows);
	//		};
			MySop.AddRowsToLP(NewCone.Rows);
			MySop.AddRowsToLP(HullCones[HullIndex][j].Rows);
			ShouldContinue = !MySop.SoplexSaysIntersect();
			if (!ShouldContinue) {
				TestCone.MIP = NewCone.MIP;
				TestCone.MIP.add_constraints(ConeToTest->Constraints);
//				TestCone.HOPolyhedron = NewCone.HOPolyhedron;
//				TestCone.HOPolyhedron.add_constraints(ConeToTest->Constraints);
				// Reduce the MIP.
				TestCone.MIP.is_satisfiable();
//				TestCone.Rows = ConstraintSystemToSoplexRows(TestCone.HOPolyhedron.constraints());
	
				for (MIP_Problem::const_iterator c = TestCone.MIP.constraints_begin(); c != TestCone.MIP.constraints_end();c++) {
					DSVector NewRow(c->space_dimension());
					for (size_t i = 0; i < c->space_dimension(); i++) {
						double val = raw_value(c->coefficient(Variable(i))).get_ui();
						if (c->coefficient(Variable(i)) < 0)
							val *= -1;
						NewRow.add(i, val);
					};
					if (c->is_equality()) {
						TestCone.Rows.add(LPRow(0, NewRow, 0));
					} else {
						stringstream s;
						s << c->inhomogeneous_term();
						int ToAppend;
						istringstream(s.str()) >> ToAppend;
						TestCone.Rows.add(LPRow(-1 * ToAppend, NewRow, infinity));
					};
				};
				
			};
		};
		if (ShouldContinue)
			continue;
		TestCone.PolytopesVisitedIndices = NewCone.PolytopesVisitedIndices;
		TestCone.ConesVisitedIndices = NewCone.ConesVisitedIndices;
		TestCone.PolytopesVisitedIndices.push_back(HullIndex);
		TestCone.ConesVisitedIndices.push_back(j);
		
		TestCone.RelationTables = RelationTables;
		TestCone.PolytopesVisited.Count = NewCone.PolytopesVisited.Count;
		TestCone.PolytopesVisited = NewCone.PolytopesVisited;

		Result.push_back(TestCone);
	};
	return Result;
}

//------------------------------------------------------------------------------
inline list<Cone> DynamicEnumerate(Cone &C, vector<vector<Cone> > &HullCones, MySoPlex &MySop) {
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

	return DoCommonRefinement(SmallestIndex, C, HullCones, MySop);
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
		if (CWIIt->IsMaximal) {
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
void ThreadEnum(vector<vector<Cone> > HullCones, int ProcessID, int ProcessCount, vector<ThreadJob> &ThreadJobs, vector<int> &BoredProcesses, TropicalPrevariety &Output, mutex &BPmtx, mutex &Outputmtx, MySoPlex MySop) {
	vector<vector<int> > Pretropisms;
	Cone C;
	bool HasCone = false;
	dimension_type dimtype = HullCones[0][0].HOPolyhedron.space_dimension();
	ThreadJob *TJ;
	while (true) {
		int j = ProcessID;
		while (not HasCone) {
			int StartIndex;
			int EndIndex;
			int Incrementer;
			TJ = &ThreadJobs[j % ProcessCount];
			if (j == ProcessID) {
				StartIndex = TJ->SharedCones.size() - 1;
				EndIndex = -1;
				Incrementer = -1;
			} else {
				StartIndex = 0;
				EndIndex = TJ->SharedCones.size();
				Incrementer = 1;
			};
			TJ->M.lock();
			// If work stealing happens, we really want to steal the less traveled cones,
			// not the more traveled cones. Need to rework that...somehow.
			for (size_t i = StartIndex; i !=EndIndex; ) {
				if (TJ->SharedCones[i].size() > 0) {
			//		clock_t BEGIN = clock();
					C = TJ->SharedCones[i].back();
			//		ConeCopyTime += double(clock() - BEGIN);
					TJ->SharedCones[i].pop_back();
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
				i+=Incrementer;
			};
			TJ->M.unlock();
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
				if (BoredProcesses[ProcessCount] == ProcessCount) {
					SoplexTime += MySop.SoplexTime;
					return;
				};
			};
			j++;
		};
		list<Cone> ResultCones = DynamicEnumerate(C, HullCones, MySop);

		// If there are remaining new cones, give them to the job queue.
		if (ResultCones.size() > 0) {
			int Index = ResultCones.front().PolytopesVisited.Count;
			
			// The cones have visited all of the polytopes.
			if (Index == HullCones.size()) {
				Outputmtx.lock();
			//	clock_t TimeStart = clock();
				list<Cone>::iterator i;
				for (i = ResultCones.begin(); i != ResultCones.end(); i++) {
					vector<int> ConeRayIndices;
					Constraint_System cs;
					C_Polyhedron HOPolyhedron(dimtype);
					if (UsePPL) {
						HOPolyhedron = i->HOPolyhedron;
					} else {
						for (size_t qq = 0; qq != i->PolytopesVisitedIndices.size(); qq++) {
							HOPolyhedron.add_constraints(HullCones[i->PolytopesVisitedIndices[qq]][i->ConesVisitedIndices[qq]].HOPolyhedron.constraints());
						};
					};
					int ConeDim = HOPolyhedron.affine_dimension();
					//cout << (*i).HOPolyhedron.generators() << endl;
					for (Generator_System::const_iterator gsi = HOPolyhedron.generators().begin(), gs_end = HOPolyhedron.generators().end(); gsi != gs_end; ++gsi) {
						// cout << (*gsi) << endl;
						if (gsi->is_point() or gsi->is_closure_point() or gsi->coefficient(Variable(gsi->space_dimension() -1)) != 0)
							continue;
						vector<int> Ray = GeneratorToPoint(*gsi, true);
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
		//		CleanupTime += double(clock() - TimeStart);
				Outputmtx.unlock();
				HasCone = false;
			} else {
				// If there is at least one cone, hold onto it for the next round.
				C = ResultCones.back();
				ResultCones.pop_back();
				if (ResultCones.size() > 0) {
					TJ = &ThreadJobs[ProcessID];
					TJ->M.lock();
					list<Cone>::iterator i;
					for (i = ResultCones.begin(); i != ResultCones.end(); i++) {
						TJ->SharedCones[Index - 1].push_back(*i);
					};
					TJ->M.unlock();
				};
			};
		} else
			HasCone = false;
	};
};

//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	bool Verbose;

	double RandomSeed = time(NULL);
	if (Verbose) {
		cout << fixed << "Random seed value: " << RandomSeed << endl;
	};
	srand(RandomSeed);
	
	int TotalProcessCount;
	vector<vector<vector<int> > > PolynomialSystemSupport;
	if (argc != 5) {
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
			cout << "Internal error: hardware_concurrency = "
				<< thread::hardware_concurrency()
				<<" but TotalProcessCount = " << TotalProcessCount << endl;
			return 1;
		};
		Verbose = true;
		UsePPL = atoi(argv[4]);
	};
	
	vector<vector<Cone> > HullCones;
	vector<double> VectorForOrientation;
	for (size_t i = 0; i != PolynomialSystemSupport[0][0].size(); i++) {
		VectorForOrientation.push_back(rand());
	};
	
	for (size_t i = 0; i != PolynomialSystemSupport.size(); i++) {
		HullCones.push_back(NewHull(PolynomialSystemSupport[i], VectorForOrientation, Verbose));
	}

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
	int Dim = HullCones[0][0].HOPolyhedron.space_dimension();
	DSVector dummycol(0);
	LPColSetReal Cols;
	for (size_t i = 0; i != Dim; i++)
		Cols.add(LPCol(0, dummycol, infinity, -infinity));

	MySoPlex MySop;
	MySop.Columns = Cols;
	DVector farkasx(MySop.numRowsReal());
	// TODO: parallelize this.
	bool ConesIntersect;
	for(int i = 0; i != HullCones.size(); i++){
		for (size_t k = 0; k != HullCones[i].size(); k++) {
			HullCones[i][k].PolytopesVisited.Indices.resize(HullCones.size());
			HullCones[i][k].PolytopesVisited.Indices[i] = 1;
			HullCones[i][k].PolytopesVisited.Count = 1;
		};
		for(size_t j = i+1; j != HullCones.size(); j++){
			for(size_t k = 0; k != HullCones[i].size(); k++){
				for(size_t l = 0; l != HullCones[j].size(); l++){
					if (UsePPL) {
						ConesIntersect = !HullCones[i][k].HOPolyhedron.is_disjoint_from(HullCones[j][l].HOPolyhedron);
					} else {
						MySop.ClearLP();
						MySop.AddRowsToLP(HullCones[i][k].Rows);
						MySop.AddRowsToLP(HullCones[j][l].Rows);
						ConesIntersect = MySop.SoplexSaysIntersect();
						// Could do dual farkas test here
					};
					if (ConesIntersect) {
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
	double PreintersectTime = double(clock() - PreintTimeStart);
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
		HullCones[SmallestIndex][i].PolytopesVisitedIndices.push_back(SmallestIndex);
		HullCones[SmallestIndex][i].ConesVisitedIndices.push_back(i);
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
		work_type work = bind(ThreadEnum, HullCones, i, TotalProcessCount, ref(ThreadJobs), ref(BoredProcesses), ref(Output), ref(BPmtx), ref(Outputmtx), MySop);
		thread_pool.submit(Parma_Polyhedra_Library::make_threadable(work));
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
	
	//double AlgorithmTotalTime = double(clock() - AlgorithmStartTime);
	//clock_t TimeAfterEndOfAlg = clock();
	/*
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
	*/
	sort(Output.Pretropisms.begin(), Output.Pretropisms.end());
	if (Verbose) {
		PrintPoints(Output.Pretropisms);
	} else {
		PrintPointsForPython(Output.Pretropisms);
	};
	cout << "Soplex time: " << SoplexTime / CLOCKS_PER_SEC << endl;
	if (Verbose) {
		//cout << "Maximal cone count--------" << endl;
		//PrintPoint(Output.FVector);
		//PrintPoint(TotalVec);
		cout << "Pre intersections: " << TotalInt << endl;
		cout << "Alg intersections: " << ConeIntersectionCount << endl;
		cout << "Total intersections: " << TotalInt + ConeIntersectionCount << endl;
		cout << "Number of pretropisms found: " << Output.Pretropisms.size() << endl;
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;/*
		cout << "AffineDim time: " << AffineDimTime / CLOCKS_PER_SEC << endl;
		cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
		cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
		cout << "Cleanup time: " << CleanupTime / CLOCKS_PER_SEC << endl;
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
		cout << "Set intersect time: " << SetIntersectTime / CLOCKS_PER_SEC << endl;
		cout << "Time after end of alg: " << double(clock() - TimeAfterEndOfAlg) / CLOCKS_PER_SEC << endl;
		cout << "Algorithm total time: " << AlgorithmTotalTime / CLOCKS_PER_SEC << endl;
		cout << "Cone copy time: " << ConeCopyTime / CLOCKS_PER_SEC << endl;
		cout << "Alg time - cone int time - set int time: " << (AlgorithmTotalTime - IntersectionTime - SetIntersectTime) / CLOCKS_PER_SEC << endl;*/
	};
}
