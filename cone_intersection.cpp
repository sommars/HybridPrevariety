#include "soplex_test.h"

//double IntersectionTime, SetIntersectTime, CleanupTime, AffineDimTime, ConeCopyTime;
int ConeIntersectionCount;

//------------------------------------------------------------------------------
list<Cone> DoCommonRefinement(
   int HullIndex,
   Cone &NewCone,
   vector<vector<Cone> > &HullCones,
   MySoPlex &MySop)
{
   // Perform common refinement for the specified cone and set of cones.
   vector<Cone> *HIndex;
   HIndex = &HullCones[HullIndex];
   BitsetWithCount *RT = &NewCone.RelationTables[HullIndex];
   list<Cone> Result;
   Cone *ConeToTest;
   for (boost::dynamic_bitset<>::size_type j = 0; j != RT->Indices.size(); j++)
   {
      if (!RT->Indices[j])
         continue;
      ConeToTest = &(*HIndex)[j];

      bool SkipIntersection = false;
      vector<BitsetWithCount> RelationTables(HullCones.size());
      for (size_t i = 0; i != NewCone.RelationTables.size(); i++)
      {
         if (!NewCone.PolytopesVisited.Indices[i])
         {
            RelationTables[i] = IntersectRTs(
               ConeToTest->RelationTables[i], NewCone.RelationTables[i]);
            if (RelationTables[i].Count == 0)
            {
               SkipIntersection = true;
               break;
            };
         };
      };
      
      if (SkipIntersection)
         continue;
      ConeIntersectionCount++;
      Cone TestCone = NewCone;
      TestCone.HOPolyhedron.add_constraints(ConeToTest->HOPolyhedron.constraints());
      if (TestCone.HOPolyhedron.is_empty())
         continue;
      TestCone.RelationTables = RelationTables;

      Result.push_back(TestCone);
   };
   return Result;
}

//------------------------------------------------------------------------------
inline list<Cone> DynamicEnumerate(
   Cone &C, vector<vector<Cone> > &HullCones, MySoPlex &MySop)
{
   // Figure out which polytope we want to visit next
   bool HaveCandidatePolytope = false;
   int SmallestInt;
   int SmallestIndex;
   for (size_t i = 0; i != C.RelationTables.size(); i++)
   {
      if ((!C.PolytopesVisited.Indices[i])
      && ((!HaveCandidatePolytope) || (C.RelationTables[i].Count<SmallestInt)))
      {
         SmallestInt = C.RelationTables[i].Count;
         SmallestIndex = i;
         HaveCandidatePolytope = true;
      };
   };

   if (!HaveCandidatePolytope)
   {
      cout << "Internal error: ";
      cout << "DynamicEnumerate did not find a next polytope to visit" << endl;
      cin.get();
   };
   C.PolytopesVisited.Indices[SmallestIndex] = 1;
   C.PolytopesVisited.Count++;

   return DoCommonRefinement(SmallestIndex, C, HullCones, MySop);
};

//------------------------------------------------------------------------------
void ParseToPrevariety(
   vector<int> &ConeRayIndices,
   int ConeDim,
   TropicalPrevariety &Output,
   bool CouldBeMaximal)
{
   // Dynamically parses output of algorithm to create full prevariety.
   while (Output.ConeTree.size() < ConeDim)
   {
      set<ConeWithIndicator> Temp;
      Output.ConeTree.push_back(Temp);
      Output.FVector.push_back(0);
   };

   ConeWithIndicator CWI;
   CWI.IsMaximal = CouldBeMaximal;
   CWI.RayIndices = ConeRayIndices;

   set<ConeWithIndicator>::iterator CWIIt = find(
      Output.ConeTree[ConeDim - 1].begin(),
      Output.ConeTree[ConeDim - 1].end(), CWI);
      
   if (CWIIt == Output.ConeTree[ConeDim - 1].end())
      Output.ConeTree[ConeDim - 1].insert(CWI);
   else
   {
      if (CWIIt->IsMaximal)
      {
         ConeWithIndicator CWITemp;
         CWITemp.RayIndices = ConeRayIndices;
         CWITemp.IsMaximal = false;
         Output.ConeTree[ConeDim - 1].erase(*CWIIt);
         Output.ConeTree[ConeDim - 1].insert(CWITemp);
      };
      return;
   };
   if (ConeDim == ConeRayIndices.size())
   {
      // This is the easy case. We can recurse exactly as we would want to.
      if (CWIIt == Output.ConeTree[ConeDim - 1].end())
      {
         if (ConeDim == 1)
            return;
         for (size_t i = 0; i != ConeRayIndices.size(); i++)
         {
            vector<int> ToRecurse;
            for (size_t j = 0; j != ConeRayIndices.size(); j++)
            {
               if (i != j)
                  ToRecurse.push_back(ConeRayIndices[j]);
            };
            ParseToPrevariety(ToRecurse, ConeDim - 1, Output, false);
         };
      }
   } else {
/*      vector<Generator> Gens;
      Generator_System gs;
      for (size_t i = 0; i != ConeRayIndices.size(); i++)
      {
         Gens.push_back(Output.IndexToGenMap[ConeRayIndices[i]]);
         gs.insert(Output.IndexToGenMap[ConeRayIndices[i]]);
      };
      Linear_Expression LE;
      for (size_t i = 0; i != Gens[0].space_dimension(); i++)
         LE += Variable(i) * 0;
      Generator ZeroPoint = point(LE); // Is this necessary? Could this just be point(0)?
      gs.insert(ZeroPoint);
      NNC_Polyhedron InitialPoly(gs);
      
      vector<Constraint> Inequalities;
      Constraint_System Equalities;
      for (Constraint_System::const_iterator i = InitialPoly.constraints().begin(),
      cs_end = InitialPoly.constraints().end(); i != cs_end; ++i)
      {
         if ((*i).is_equality())
            Equalities.insert(*i);
         else
            Inequalities.push_back(*i);
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
void ThreadEnum(
   vector<vector<Cone> > HullCones,
   int ProcessID,
   int ProcessCount,
   vector<ThreadQueue> &ThreadQueues,
   vector<int> &BoredProcesses,
   TropicalPrevariety &Output,
   mutex &BPmtx,
   mutex &Outputmtx,
   MySoPlex MySop)
{
   // Manages work stealing and dishing out jobs from thread queues.
   vector<vector<int> > Pretropisms;
   Cone C;
   bool HasCone = false;
   ThreadQueue *TQ;
   stringstream s;
   while (true)
   {
      int j = ProcessID;
      while (not HasCone)
      {
         int StartIndex;
         int EndIndex;
         int Incrementer;
         TQ = &ThreadQueues[j % ProcessCount];
         if (j == ProcessID)
         {
            StartIndex = TQ->SharedCones.size() - 1;
            EndIndex = -1;
            Incrementer = -1;
         } else {
            StartIndex = 0;
            EndIndex = TQ->SharedCones.size();
            Incrementer = 1;
         };
         TQ->M.lock();
         // If work stealing happens, we want to steal the less traveled cones,
         // not the more traveled cones.
         for (size_t i = StartIndex; i !=EndIndex; )
         {
            if (TQ->SharedCones[i].size() > 0)
            {
               C = TQ->SharedCones[i].back();
               TQ->SharedCones[i].pop_back();
               HasCone = true;
               // This case happens when the process tried to steal from every
               // process possible but did not succeed, making it bored.
               // It added itself to BoredProcesses, but then continued looking
               // and found a cone it could have. This made it no longer a
               // bored process.
               if (j > (ProcessID + ProcessCount))
               {
                  BPmtx.lock();
                  BoredProcesses[ProcessID] = 0;
                  BoredProcesses[ProcessCount] -= 1;
                  BPmtx.unlock();
               };
               break;
            };
            i+=Incrementer;
         };
         TQ->M.unlock();
         if (HasCone)
            break;
         // This case means that we spun through each of the other ThreadQueues
         // and they were all empty. Now we need to make the process bored.
         if (j == ProcessCount + ProcessID)
         {
            BPmtx.lock();
            BoredProcesses[ProcessID] = 0;
            BoredProcesses[ProcessCount] += 1;
            BPmtx.unlock();
         };
         if ((j % ProcessCount) == ProcessID)
         {
            if (BoredProcesses[ProcessCount] == ProcessCount)
            {
               ofstream myfile;
               myfile.open("output" + to_string(ProcessID) + ".txt");
               myfile << s.rdbuf();
               myfile.close();
               return;
            };
         };
         j++;
      };
      list<Cone> ResultCones = DynamicEnumerate(C, HullCones, MySop);

      // If there are remaining new cones, give them to the job queue.
      if (ResultCones.size() > 0)
      {
         int Index = ResultCones.front().PolytopesVisited.Count;
         
         // The cones have visited all of the polytopes.
         if (Index == HullCones.size())
         {
            list<Cone>::iterator i;
            for (i = ResultCones.begin(); i != ResultCones.end(); i++)
            {
               vector<int> ConeRayIndices;
               int ConeDim = i->HOPolyhedron.affine_dimension() - 1;
               if (ConeDim == 0)
                  continue;
               s << endl;
               s << ConeDim << endl;
               for (Generator_System::const_iterator gsi = 
                       i->HOPolyhedron.generators().begin(),
                       gs_end = i->HOPolyhedron.generators().end();
                   gsi != gs_end;
                   ++gsi)
               {
                  if (gsi->is_point() 
                  || gsi->is_closure_point() 
                  || gsi->coefficient(Variable(gsi->space_dimension() -1)) != 0)
                     continue;
                  Generator G = *gsi;
                  vector<int> Ray = GeneratorToPoint(G, true);
                  s << "{ ";
                     for (vector<int>::iterator it=Ray.begin();
                          it != Ray.end();
                          it++)
                     {
                        s << (*it) << " ";
                     }
                     s << "}" << endl;
                  if (gsi->is_line())
                  {
                     // If it's a line, print the other ray
                     s << "{ ";
                     for (vector<int>::iterator it=Ray.begin();
                          it != Ray.end();
                          it++)
                     {
                        s << (*it) * (-1) << " ";
                     }
                     s << "}" << endl;
                  }
               };
            };
            HasCone = false;
         } else {
            // If there is at least one cone, hold onto it for the next round.
            C = ResultCones.back();
            ResultCones.pop_back();
            if (ResultCones.size() > 0)
            {
               TQ = &ThreadQueues[ProcessID];
               TQ->M.lock();
               list<Cone>::iterator i;
               for (i = ResultCones.begin(); i != ResultCones.end(); i++)
                  TQ->SharedCones[Index - 1].push_back(*i);
               TQ->M.unlock();
            };
         };
      } else
         HasCone = false;
   };
};

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
   // Main program for computing tropical prevarieties.
   bool Verbose = true;

   double RandomSeed = time(NULL);
   if (Verbose)
      cout << fixed << "Random seed value: " << RandomSeed << endl;
   srand(RandomSeed);
   
   int TotalProcessCount;
   vector<vector<vector<int> > > PolynomialSystemSupport;
   if (argc != 4)
   {
      if (Verbose)
      {
         cout << "Expected three arguments. ";
         cout << "Waiting for user input system..." << endl;
      };
      TotalProcessCount = 1;
      string input;
      cin >> input;
      PolynomialSystemSupport = ParseToSupport(input);
      Verbose = false;
   } else {
      int n = atoi(argv[1]);
      string SystemName = string(argv[2]);
      if (SystemName == "reducedcyclicn")
         PolynomialSystemSupport = CyclicN(n, true);
      else if (SystemName == "cyclicn")
         PolynomialSystemSupport = CyclicN(n, false);
      else if (SystemName == "random")
         PolynomialSystemSupport = RandomSimplices(n);
      else {
         cout << "Internal error: only supported systems are: ";
         cout << "reducedcyclicn, cyclicn, or random." << endl;
         return 1;
      };
      TotalProcessCount = atoi(argv[3]);
      if (TotalProcessCount > thread::hardware_concurrency())
      {
         cout << "Internal error: hardware_concurrency = "
            << thread::hardware_concurrency()
            <<" but TotalProcessCount = " << TotalProcessCount << endl;
         return 1;
      };
      Verbose = true;
   };
   
   vector<vector<Cone> > HullCones;
   vector<double> VectorForOrientation;
   for (size_t i = 0; i != PolynomialSystemSupport[0][0].size(); i++)
      VectorForOrientation.push_back(rand());
   
   for (size_t i = 0; i != PolynomialSystemSupport.size(); i++)
      HullCones.push_back(NewHull(PolynomialSystemSupport[i], VectorForOrientation, Verbose));

   for(size_t i = 0; i != HullCones.size(); i++)
   {
      for(size_t j = 0; j != HullCones[i].size(); j++)
      {
         for(size_t k = 0; k != HullCones.size(); k++)
         {
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
   MySoPlex MySop;
   // TODO: parallelize this.
   for(int i = 0; i != HullCones.size(); i++)
   {
      for (size_t k = 0; k != HullCones[i].size(); k++)
      {
         HullCones[i][k].PolytopesVisited.Indices.resize(HullCones.size());
         HullCones[i][k].PolytopesVisited.Indices[i] = 1;
         HullCones[i][k].PolytopesVisited.Count = 1;
      };
      for(size_t j = i+1; j != HullCones.size(); j++)
      {
         for(size_t k = 0; k != HullCones[i].size(); k++)
         {
            for(size_t l = 0; l != HullCones[j].size(); l++)
            {
               if (!HullCones[i][k].HOPolyhedron.is_disjoint_from(
                   HullCones[j][l].HOPolyhedron))
               {
                  HullCones[i][k].RelationTables[j].Indices[l] = 1;
                  HullCones[j][l].RelationTables[i].Indices[k] = 1;
                  HullCones[i][k].RelationTables[j].Count++;
                  HullCones[j][l].RelationTables[i].Count++;
               } else
                  NonInt++;
               TotalInt++;
            };
         };
      };
      if (Verbose)
         printf("Finished level %d of pre-intersections.\n", i);
   };
   double PreintersectTime = double(clock() - PreintTimeStart);
   if (Verbose)
   {
      cout << "Total Intersections: " << TotalInt << endl;
      cout << ", Non Intersections: " << NonInt << endl;
      PreintersectTime = PreintersectTime / CLOCKS_PER_SEC;
      cout << "Preintersection time: " << PreintersectTime << endl;
   };
   
   
   // This is one way to do it. It seems reasonable to pick sum, median, mean...
   int SmallestInt = 1000000;
   int SmallestIndex = -1;
   for (size_t i = 0; i != HullCones.size(); i++)
   {
      int TestValue = 0;
      for (size_t j = 0; j != HullCones[i].size(); j++)
      {
         for (size_t k = 0; k != HullCones[i][j].RelationTables.size(); k++)
            TestValue += HullCones[i][j].RelationTables[k].Count;
      };
      if (TestValue < SmallestInt)
      {
         SmallestInt = TestValue;
         SmallestIndex = i;
      };
   };
   if (SmallestIndex == -1)
   {
      cout << "Internal error: value of -1 for SmallestIndex" << endl;
      return 1;
   };
   
   vector<ThreadQueue> ThreadQueues;
   for (size_t i = 0; i != TotalProcessCount; i++)
   {
      vector<list<Cone> > SharedCones;
      for (size_t i = 0; i != HullCones.size() - 1; i++)
      {
         list<Cone> Temp;
         SharedCones.push_back(Temp);
      };
      ThreadQueue TQ(SharedCones);
      ThreadQueues.push_back(TQ);
   };

   for (size_t i = 0; i != HullCones[SmallestIndex].size(); i++)
      ThreadQueues[i % TotalProcessCount].SharedCones[0].push_back(
         HullCones[SmallestIndex][i]);
   
   mutex BPmtx;
   mutex Outputmtx;
   TropicalPrevariety Output;
   vector<int> BoredProcesses (TotalProcessCount + 1, 0);
   clock_t AlgorithmStartTime = clock();
   typedef function<void()> work_type;
   Thread_Pool<work_type> thread_pool(TotalProcessCount);
   for (size_t i = 0; i != TotalProcessCount; i++)
   {
      work_type work = bind(
         ThreadEnum,
         HullCones,
         i,
         TotalProcessCount,
         ref(ThreadQueues),
         ref(BoredProcesses),
         ref(Output),
         ref(BPmtx),
         ref(Outputmtx),
         MySop);
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
   
   sort(Output.Pretropisms.begin(), Output.Pretropisms.end());
   if (Verbose)
      PrintPoints(Output.Pretropisms);
   else
      PrintPointsForPython(Output.Pretropisms);
   if (Verbose)
   {
      //cout << "Maximal cone count--------" << endl;
      //PrintPoint(Output.FVector);
      //PrintPoint(TotalVec);
      cout << "Pre intersections: " << TotalInt << endl;
      cout << "Alg intersections: " << ConeIntersectionCount << endl;
      cout << "Total intersections: " << TotalInt + ConeIntersectionCount << endl;
      cout << "Number of pretropisms found: " << Output.Pretropisms.size() << endl;
      cout << "Preintersection time: " << PreintersectTime << endl;/*
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
