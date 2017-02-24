#include "prevariety_types.h"

//------------------------------------------------------------------------------
inline bool Set1IsSubsetOfSet2(set<int> &S1, set<int> &S2)
{
   // Tests if S1 is a subset of S2;
   set<int>::iterator S1Itr = S1.begin();
   set<int>::iterator S2Itr = S2.begin();
   while ((S1Itr != S1.end()) && (S2Itr != S2.end()))
   {
      if (*S1Itr < *S2Itr)
         return false;
      else if (*S2Itr<*S1Itr)
         ++S2Itr;
      else 
      {
         S1Itr++;
         S2Itr++;
      };
   };
   
   // If S2 is at the end and there still are elements in S1, then false
   if ((S1Itr != S1.end()) && (S2Itr == S2.end()))
      return false;
   else
      return true;
}

//------------------------------------------------------------------------------
void TropicalPrevariety::MarkMaximalCones(void)
{
   if (ConeTree.size() < 2)
      return;
      
   // All of the highest dimensional cones must be maximal,
   // so should start at second highest dimension
   for(size_t i = ConeTree.size() - 2; i != -1; i--)
   {
      for(size_t j = 0; j != ConeTree[i].size(); j++)
      {
         bool ShouldMoveToNextCone = false;
         for(size_t k = ConeTree.size() - 1; k != i; k--)
         {
            for(size_t l = 0; l != ConeTree[k].size(); l++)
            {
               if (Set1IsSubsetOfSet2(ConeTree[i][j].RayIndices, ConeTree[k][l].RayIndices))
               {
                  ShouldMoveToNextCone = true;
                  ConeTree[i][j].IsMaximal = false;
               }
               if (ShouldMoveToNextCone)
                  break;
            }
            if (ShouldMoveToNextCone)
               break;
         }
      };
   };
}

//------------------------------------------------------------------------------
void TropicalPrevariety::PrintRayToIndexMap(void)
{
   // Prints ray to index map. Necessary for interpreting output.
   if (RayToIndexMap.size() != 0)
      cout << "------ Rays  ------"<< endl;
   for(map<vector<int>, int>::iterator itr = RayToIndexMap.begin();
       itr != RayToIndexMap.end();
       ++itr)
   {
      cout << itr->second << ": { ";
      for (size_t i = 0; i != itr->first.size(); i++)
         cout << itr->first[i] << " ";
      cout << "}" << endl;
   }  
}

//------------------------------------------------------------------------------
void TropicalPrevariety::PrintMaximalCones(void)
{
   // Prints maximal cones from prevariety object.
   for (size_t i = 0; i != ConeTree.size(); i++)
   {
      if (ConeTree[i].size() > 0)
         cout << "------ Cones of dimension " << i + 1 << " ------"<< endl;
      for (size_t j = 0; j != ConeTree[i].size(); j++)
      {
         if (!ConeTree[i][j].IsMaximal)
            continue;
         set<int>::iterator it;
         cout << "{ ";
         for (it=ConeTree[i][j].RayIndices.begin();
              it != ConeTree[i][j].RayIndices.end();
              it++)
            cout << (*it) << " ";
         cout << "}" << endl;
      };
   };
}


