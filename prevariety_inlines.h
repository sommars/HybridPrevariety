// Various functions that ought to be inlined

//------------------------------------------------------------------------------
inline bool operator== (
   const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
   // Test equality of two ConeWithIndicator
   return  s1.RayIndices == s2.RayIndices;
}

//------------------------------------------------------------------------------
inline bool operator< (
   const ConeWithIndicator & s1, const ConeWithIndicator & s2)
{
   // Test < for two ConeWithIndicator
   return  s1.RayIndices < s2.RayIndices;
}

//------------------------------------------------------------------------------
inline NNC_Polyhedron IntersectCones(NNC_Polyhedron ph1, NNC_Polyhedron &ph2)
{
   // Intersect NNC_Polyhedron
   ph1.add_constraints(ph2.constraints());
   return ph1;
};

//------------------------------------------------------------------------------
inline C_Polyhedron IntersectCones(C_Polyhedron ph1, C_Polyhedron &ph2)
{
   // Intersect C_Polyhedron
   ph1.add_constraints(ph2.constraints());
   ph1.affine_dimension();
   return ph1;
};

//------------------------------------------------------------------------------
inline LPRowSetReal ConstraintSystemToSoplexRows(const Constraint_System &cs)
{
   // Converts a PPL Constraint_System to an LPRowSetReal for consumption
   // by SoPlex.
   LPRowSetReal Rows;
   for (Constraint_System::const_iterator c = cs.begin(), cs_end = cs.end();
        c != cs_end;
        c++)
   {
      DSVector NewRow(c->space_dimension());
      for (size_t i = 0; i < c->space_dimension(); i++)
      {
         double val = raw_value(c->coefficient(Variable(i))).get_ui();
         if (c->coefficient(Variable(i)) < 0)
            val *= -1;
         NewRow.add(i, val);
      };
      if (c->is_equality())
         Rows.add(LPRow(0, NewRow, 0));
      else
      {
         stringstream s;
         s << c->inhomogeneous_term();
         int ToAppend;
         istringstream(s.str()) >> ToAppend;
         Rows.add(LPRow(-1 * ToAppend, NewRow, infinity));
      };
   };
   return Rows;
};

//------------------------------------------------------------------------------
inline vector<vector<int> > FindInitialForm(
   vector<vector<int> > &Points, vector<int> &Vector)
{
   // Computes the initial form of a vector and a set of points.
   if (Points.size() == 0)
      return Points;
   vector<vector<int> > InitialForm;

   vector<int> *Point;
   Point = &(Points[0]);
   InitialForm.push_back(*Point);
   int MinimalIP = 0;
   for (size_t i = 0; i != (*Point).size(); i++)
      MinimalIP += Vector[i] * (*Point)[i];

   for (size_t i = 1; i != Points.size(); i++)
   {
      Point = &(Points[i]);
      int IP = 0;
      for (size_t j = 0; j != (*Point).size(); j++)
         IP += Vector[j] * (*Point)[j];
      if (MinimalIP > IP)
      {
         MinimalIP = IP;
         InitialForm.clear();
         InitialForm.push_back(*Point);
      } else if (IP == MinimalIP)
         InitialForm.push_back(*Point);
   };
   
   return InitialForm;
};

//------------------------------------------------------------------------------
inline set<int> IntersectSets(set<int> &S1, set<int> &S2)
{
   // Computes and returns the intersection of two sets.
   set<int>::iterator S1Itr = S1.begin();
   set<int>::iterator S2Itr = S2.begin();
   set<int> Result;
   while ((S1Itr != S1.end()) && (S2Itr != S2.end()))
   {
      if (*S1Itr < *S2Itr)
         ++S1Itr;
      else if (*S2Itr<*S1Itr)
         ++S2Itr;
      else 
      {
         Result.insert(*S1Itr);
         S1Itr++;
         S2Itr++;
      };
   };
   return Result;
};

//------------------------------------------------------------------------------
inline bool SetsDoIntersect(set<int> &S1, set<int> &S2)
{
   // Tests if two sets have a non-trivial intersection
   set<int>::iterator S1Itr = S1.begin();
   set<int>::iterator S2Itr = S2.begin();
   while ((S1Itr != S1.end()) && (S2Itr != S2.end()))
   {
      if (*S1Itr < *S2Itr)
         ++S1Itr;
      else if (*S2Itr<*S1Itr)
         ++S2Itr;
      else
         return true;
   };
   return false;
};

//------------------------------------------------------------------------------
inline BitsetWithCount IntersectRTs(BitsetWithCount R1, BitsetWithCount &R2)
{
   // Intersects two relation tables and correctly establishes the 
   // count property
   R1.Indices = R1.Indices&=R2.Indices;
   R1.Count = R1.Indices.count();
   return R1;
};