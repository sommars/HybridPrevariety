#include "convex_hull.h"

//------------------------------------------------------------------------------
vector<Cone> NewHull(
   vector<vector<int> > &Points,
   vector<double> &VectorForOrientation, 
   bool Verbose)
{
   // For a given set of points and orientation vector, this function computes
   // the set of half open edge cones for the polytope defined by these points.
   Hull H;
   H.CPolyhedron = FindCPolyhedron(Points);
   map<double,vector<int> > DoubleToPt;
   vector<double> IPs;
   for (size_t i = 0; i != Points.size(); i++)
   {
      double IP = DoubleInnerProduct(Points[i], VectorForOrientation);
      IPs.push_back(IP);
      DoubleToPt[IP] = Points[i];
   };
   sort(IPs.begin(), IPs.end());
   for (size_t i = 0; i != IPs.size(); i++)
      H.Points.push_back(DoubleToPt[IPs[i]]);

   H.AffineDimension = H.CPolyhedron.affine_dimension();
   H.SpaceDimension = H.CPolyhedron.space_dimension();
      
   // Create PointToIndexMap
   for (size_t i = 0; i != H.Points.size(); i++)
   {
      vector<int> Point = H.Points[i];
      H.PointToIndexMap[Point]=i;
      H.IndexToPointMap[i]=Point;
   };

   FindFacets(H);
   FindEdges(H);
   for (size_t i = 0; i != H.Points.size(); i++)
   {
      vector<Constraint> Constraints;
      vector<int> Pt = H.Points[i];
      int PtIndex = H.PointToIndexMap[Pt];
      // Go through all of the edges. If the edge is not on the facet.
      for (size_t j = 0; j != H.Edges.size(); j++)
      {
         if (H.Edges[j].PointIndices.find(PtIndex) 
             != H.Edges[j].PointIndices.end())
         {
            set<int>::iterator PtIter;
            vector<int> OtherPt;
            int OtherPtIndex;
            for (PtIter=H.Edges[j].PointIndices.begin();
                 PtIter != H.Edges[j].PointIndices.end();
                 PtIter++)
            {
               OtherPtIndex = (*PtIter);
               if (OtherPtIndex != PtIndex)
                  OtherPt = H.IndexToPointMap[(*PtIter)];
            };
            
            Linear_Expression LE;
            for (size_t k = 0; k != Pt.size(); k++)
               LE += (OtherPt[k] - Pt[k]) * Variable(k);
            
            // Manufacture the constraint from the edge.
            Constraint c;
            if (OtherPtIndex > PtIndex)
               c = (LE >= 0);
            else
               c = (LE > 0);
            Constraints.push_back(c);
         };
      };
      // If no non-strict inequalities exist, then you may simply throw the cone
      // away. This process produces just one half open cone for each edge.
      for (size_t j = 0; j != Constraints.size(); j++)
      {
         // Take one of your non-strict inequalities for the cone.
         Constraint TempConstraint = Constraints[j];
         if (!TempConstraint.is_strict_inequality())
         {
            // Introduce it as an equation to get a lower dimensional cone,
            // which you output.
            Constraints[j] = InequalityToEquation(TempConstraint);

            // Now we raise to a higher dimension to represent as a standard
            // linear programming problem.
            Constraint_System csHigherDim;
            csHigherDim.insert(Constraint(Variable(H.SpaceDimension) >= 1));
            for (size_t k = 0; k != Constraints.size(); k++)
            {
               if (Constraints[k].is_strict_inequality())
               {
                  Linear_Expression LE;
                  for (size_t l = 0; l < Constraints[k].space_dimension(); l++)
                     LE += Constraints[k].coefficient(Variable(l))*Variable(l);
                  LE += -1 * Variable(Constraints[k].space_dimension());
                  csHigherDim.insert(Constraint(LE >= 0));
               } else
                  csHigherDim.insert(Constraints[k]);
            };
            
            Cone NewCone;
            NewCone.HOPolyhedron = C_Polyhedron(csHigherDim);
            NewCone.HOPolyhedron.minimized_constraints();
            NewCone.HOPolyhedron.minimized_generators();
            NewCone.HOPolyhedron.affine_dimension();

            H.Cones.push_back(NewCone);
            //Introduce it as a strict inequality to describe the rest.
            // Call recursively.
            Constraints[j] = InequalityToStrictInequality(TempConstraint);
         };
      };
   };

   if (Verbose)
   {
      cout << "Convex hull------------------------" << endl;
      PrintPoints(H.Points);
      cout << "Affine dimension: " << H.AffineDimension << endl;
      cout << "Space dimension: " << H.SpaceDimension << endl;
      cout << "Number of edges: " << H.Edges.size() << endl;
      cout << "Number of facets: " << H.Facets.size() << endl << endl;
   };
   return H.Cones;
}

//------------------------------------------------------------------------------
void FindFacets(Hull &H)
{
   // Find the Facets by shooting rays at the polytopes
   // spin through all of the constraints. shoot each constraint at the polytope
   Constraint_System cs = H.CPolyhedron.minimized_constraints();
   for (Constraint_System::const_iterator i = cs.begin(), cs_end = cs.end();
        i != cs_end; 
        ++i)
   {
      if (!i->is_inequality())
         continue;
      Constraint C = *i;
      vector<int> Pt = ConstraintToPoint(C);
      vector<vector<int> > FacetPts = FindInitialForm(H.Points, Pt);
      
      Facet F;
      for (vector<vector<int> >::iterator itr = FacetPts.begin();
           itr != FacetPts.end();
           itr++)
      {
         F.PointIndices.insert(H.PointToIndexMap[*itr]);
      };
      H.Facets.push_back(F);
   }
}

//------------------------------------------------------------------------------
void FindEdges(Hull &H)
{
   // Find the set of edges. This could be done in a better way, but this works.
   vector<vector<int> > CandidateEdges = FindCandidateEdges(H);

   int Dim = H.CPolyhedron.affine_dimension() - 1;

   vector<vector<int> >::iterator itr;
   for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++)
   {
      vector<int> CandidateEdge = (*itr);
      int Point1 = CandidateEdge[0];
      int Point2 = CandidateEdge[1];
   
      int FacetCount = 0;
      for (vector<Facet>::iterator FacetIt = H.Facets.begin();
           FacetIt != H.Facets.end(); 
           FacetIt++)
      {
         set<int> PtIndices = FacetIt->PointIndices;
         bool Point1IsInFacet = PtIndices.find(Point1) != PtIndices.end();
         bool Point2IsInFacet = PtIndices.find(Point2) != PtIndices.end();
         if (Point1IsInFacet and Point2IsInFacet)
            FacetCount++;
      };

      if (FacetCount >= Dim)
      {
         Edge NewEdge;
         NewEdge.PointIndices.insert(Point1);
         NewEdge.PointIndices.insert(Point2);
         H.Edges.push_back(NewEdge);
      }
   };
   
   // After all of the edges have been generated, 
   // fill out all of the neighbors on all of the edges.
   for (size_t Edge1Index = 0; Edge1Index != H.Edges.size(); Edge1Index++)
   {
      Edge Edge1 = H.Edges[Edge1Index];

      for (size_t Edge2Index = 0; Edge2Index != H.Edges.size(); Edge2Index++)
      {
         if (Edge1Index == Edge2Index)
            continue;

         int IntersectionCount = IntersectSets(
            Edge1.PointIndices, H.Edges[Edge2Index].PointIndices).size();
         if (IntersectionCount == 1)
         {
            H.Edges[Edge1Index].NeighborIndices.insert(Edge2Index);
            H.Edges[Edge2Index].NeighborIndices.insert(Edge1Index);
         } else if (IntersectionCount > 1)
         {
            cout << "Internal Error: ";
            cout << "Two edges intersect at more than one point" << endl;
            cin.get();
         }
      };
   };
}

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull &H)
{
   // Helper function for FindEdges. Again, there _has_ to be a better way to
   // do this.
   vector<vector<int> > CandidateEdges;
   int n = H.Points.size();
   vector<int> d(n);
   for (size_t i = 0; i != d.size(); ++i)
      d[i] = i;
   do
   {
      if (d[0] < d[1])
      {
         vector<int> CandidateEdgeIndices;
         CandidateEdgeIndices.push_back(d[0]);
         CandidateEdgeIndices.push_back(d[1]);
         CandidateEdges.push_back(CandidateEdgeIndices);
      };
      reverse(d.begin()+2, d.end());
   } while (next_permutation(d.begin(), d.end()));
   return CandidateEdges;
}

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > &Points)
{
   // Converts an input set of pts to a PPL C_Polyhedron
   Generator_System gs;
   for (vector<vector<int> >::iterator itr=Points.begin();
        itr != Points.end();
        itr++)
   {
      Linear_Expression LE;
      for (size_t i = 0; i != itr->size(); i++)
         LE += Variable(i) * ((*itr)[i]);
      gs.insert(point(LE));
   };
   return C_Polyhedron(gs);
}

