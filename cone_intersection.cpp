#include "soplex_test.h"

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
int ThreadEnum(
   vector<vector<Cone> > HullCones,
   int ProcessID,
   int ProcessCount,
   vector<ThreadQueue> &ThreadQueues,
   vector<int> &BoredProcesses,
   TropicalPrevariety &Output,
   mutex &BPMtx,
   mutex &OutputMtx,
   MySoPlex MySop)
{
   // Manages work stealing and dishing out jobs from thread queues.
   vector<vector<int> > Pretropisms;
   Cone C;
   bool HasCone = false;
   ThreadQueue *TQ;
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
                  BPMtx.lock();
                  BoredProcesses[ProcessID] = 0;
                  BoredProcesses[ProcessCount] -= 1;
                  BPMtx.unlock();
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
            BPMtx.lock();
            BoredProcesses[ProcessID] = 0;
            BoredProcesses[ProcessCount] += 1;
            BPMtx.unlock();
         };
         if ((j % ProcessCount) == ProcessID)
         {
            if (BoredProcesses[ProcessCount] == ProcessCount)
            {
               return 1;
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
               int ConeDim = i->HOPolyhedron.affine_dimension() - 1;
               if (ConeDim == 0)
                  continue;
               ConeWithIndicator CWI;
               CWI.IsMaximal = true;
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
                  OutputMtx.lock();
                  clock_t FindStart = clock();
                  map<vector<int>, int>::iterator itr = Output.RayToIndexMap.find(Ray);
                  if (itr == Output.RayToIndexMap.end())
                  {
                     Output.RayToIndexMap[Ray] = Output.RayToIndexMap.size();
                     CWI.RayIndices.insert(Output.RayToIndexMap.size());
                  } else
                     CWI.RayIndices.insert(itr->second);
                  OutputMtx.unlock();
                  
                  if (gsi->is_line()) {
                     // If it's a line, we have to add both rays that make it up
                     for(size_t j = 0; j != Ray.size(); j++)
                        Ray[j] = -1 * Ray[j];
                     
                     OutputMtx.lock();
                     map<vector<int>, int>::iterator itr2 = Output.RayToIndexMap.find(Ray);
                     if (itr2 == Output.RayToIndexMap.end())
                     {
                        Output.RayToIndexMap[Ray] = Output.RayToIndexMap.size();
                        CWI.RayIndices.insert(Output.RayToIndexMap.size());
                     } else
                        CWI.RayIndices.insert(itr2->second);
                     OutputMtx.unlock();
                  }
               };
               OutputMtx.lock();
               while (Output.ConeTree.size() < (ConeDim))
               {
                  vector<ConeWithIndicator> Temp;
                  Output.ConeTree.push_back(Temp);
               };
               Output.ConeTree[ConeDim - 1].push_back(CWI);
               OutputMtx.unlock();
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
      cout << "Non Intersections: " << NonInt << endl;
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
   
   mutex BPMtx;
   mutex OutputMtx;
   TropicalPrevariety Output;
   vector<int> BoredProcesses (TotalProcessCount + 1, 0);
   clock_t AlgorithmStartTime = clock();
   Thread_Pool<function<int()>> thread_pool(TotalProcessCount);
   for (size_t i = 0; i != TotalProcessCount; i++)
   {
      thread_pool.submit(make_threadable(bind(
         ThreadEnum,
         HullCones,
         i,
         TotalProcessCount,
         ref(ThreadQueues),
         ref(BoredProcesses),
         ref(Output),
         ref(BPMtx),
         ref(OutputMtx),
         MySop)));
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
   clock_t MarkingTimeStart = clock();
   Output.MarkMaximalCones();
   double MarkingTime = double(clock() - MarkingTimeStart) / CLOCKS_PER_SEC;
   Output.PrintRayToIndexMap();
   Output.PrintMaximalCones();
   if (Verbose)
   {
      cout << "Intersections for building RT: " << TotalInt << endl;
      cout << "Alg intersections: " << ConeIntersectionCount << endl;
      cout << "Total intersections: " << TotalInt + ConeIntersectionCount << endl;
      cout << "Preintersection time: " << PreintersectTime << endl;
      cout << "Marking time: " << MarkingTime << endl;
   };
}
