#include "process_output.h"

//------------------------------------------------------------------------------
void StreamRayToIndexMap(TropicalPrevariety &TP, stringstream &s)
{
   // Streams ray to index map. Necessary for interpreting output.
   if (TP.RayToIndexMap.size() != 0)
      s << "------ Rays  ------" << endl;
   for(map<vector<int>, int>::iterator itr = TP.RayToIndexMap.begin();
       itr != TP.RayToIndexMap.end();
       ++itr)
   {
      s << itr->second << ": { ";
      for (size_t i = 0; i != itr->first.size(); i++)
         s << itr->first[i] << " ";
      s << "}" << endl;
   }
   return;
}

//------------------------------------------------------------------------------
void ParallelMarking(TropicalPrevariety &TP, mutex &ConeMtx)
{
   for(size_t i = TP.ConeTree.size() - 2; i != -1; i--)
   {
      for(size_t j = 0; j != TP.ConeTree[i].size(); j++)
      {
         ConeMtx.lock();
         if (TP.ConeTree[i][j].Status != 2)
         {
            ConeMtx.unlock();
            continue;
         };
         ConeMtx.unlock();
         TP.ConeTree[i][j].Status = 2;
         
         bool WeKnowConeIsNotMaximal = false;
         for(size_t k = TP.ConeTree.size() - 1; k != i; k--)
         {
            for(size_t l = 0; l != TP.ConeTree[k].size(); l++)
            {
               if (TP.ConeTree[k][l].Status == 0)
                  continue;
                  
               if (Set1IsSubsetOfSet2(TP.ConeTree[i][j].RayIndices, TP.ConeTree[k][l].RayIndices))
               {
                  WeKnowConeIsNotMaximal = true;
                  TP.ConeTree[i][j].Status = 0;
               }
               if (WeKnowConeIsNotMaximal)
                  break;
            }
            if (WeKnowConeIsNotMaximal)
               break;
         }
         // TODO!!!! FIX THIS! Currently not ignoring non-maximal cones. easy fix
         if (!WeKnowConeIsNotMaximal)
            TP.ConeTree[i][j].Status = 1;
      };
   };
   return;
}

//------------------------------------------------------------------------------
void MarkMaximalCones(TropicalPrevariety &TP, int ProcessCount)
{
   if (TP.ConeTree.size() < 2)
      return;
   {
      Thread_Pool<function<void()>> thread_pool(ProcessCount);
      mutex ConeMtx;
      for (size_t i = 0; i != ProcessCount; i++)
      {
         thread_pool.submit(bind(
            ParallelMarking,
            ref(TP),
            ref(ConeMtx)));
      }
      
      // Wait for all workers to complete.
      thread_pool.finalize();
   }
   return;
}

