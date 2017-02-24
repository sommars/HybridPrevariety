#include "prevariety_util.h"

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator &g, bool KnockOffLastTerm)
{
   // Converts a PPL Generator object to a vector of int
   vector<int> Result;
   int Dim = g.space_dimension();
   if (KnockOffLastTerm)
      Dim = Dim - 1;
   for (size_t i = 0; i < Dim; i++)
   {
      stringstream s;
      s << (g).coefficient(Variable(i));
      int ToAppend;
      istringstream(s.str()) >> ToAppend;
      Result.push_back(ToAppend);
   }
   return Result;
}

//------------------------------------------------------------------------------
Constraint InequalityToStrictInequality(Constraint &c)
{
   // Converts a PPL inequality to a strict inequality
   Linear_Expression LE;
   for (size_t i = 0; i < c.space_dimension(); i++)
      LE += c.coefficient(Variable(i)) * Variable(i);
   Constraint c2(LE > c.inhomogeneous_term());
   return c2;
}

//------------------------------------------------------------------------------
Constraint InequalityToEquation(Constraint &c)
{
   // Converts a PPL inequality to an equation
   Linear_Expression LE;
   for (size_t i = 0; i < c.space_dimension(); i++)
      LE += c.coefficient(Variable(i)) * Variable(i);
   Constraint c2(LE == c.inhomogeneous_term());
   return c2;
}

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint &c)
{
   // Converts a PPL Constraint to a point
   vector<int> Result;
   for (size_t i = 0; i < c.space_dimension(); i++)
   {
      stringstream s;
      s << (c).coefficient(Variable(i));
      int ToAppend;
      istringstream(s.str()) >> ToAppend;
      Result.push_back(ToAppend);
   }
   return Result;
}

//------------------------------------------------------------------------------
double DoubleInnerProduct(vector<int> &V1, vector<double> &V2)
{
   // Computes the inner product of two vectors.
   if (V1.size() != V2.size())
      throw out_of_range("Internal error: DoubleInnerProduct for vectors"
                             "with different sizes");
   double Result = 0;
   for (size_t i = 0; i != V1.size(); i++)
      Result += V1[i] * V2[i];
   return Result;
}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > ParseToSupport(string &Input)
{
   // Takes in strings that look like this:
   //'[[[1,0,0,0][0,1,0,0][0,0,1,0][0,0,0,1]]' + 
   //'[[1,1,0,0][0,1,1,0][1,0,0,1][0,0,1,1]]' + 
   //'[[1,1,1,0][1,1,0,1][1,0,1,1][0,1,1,1]]]'
   // and converts them to support sets. This exists for the Python interface.
   vector<vector<vector<int> > > Result;
   int ParenCount = 0;
   string NewInt;
   vector<vector<int> > Polynomial;
   vector<int> Monomial;
   for (string::iterator it=Input.begin(); it!=Input.end(); ++it)
   {
      if ((*it) == '[')
         ParenCount++;
      else if ((*it) == ']')
      {
         ParenCount--;
         if (NewInt.length() > 0)
         {
            Monomial.push_back(stoi(NewInt));
            NewInt.clear();
         };
         if (ParenCount == 2)
         {
            Polynomial.push_back(Monomial);
            Monomial.clear();
         } else if (ParenCount == 1)
         {
            Result.push_back(Polynomial);
            Polynomial.clear();
         };
      } else if ((*it) == ',')
      {
         if (NewInt.length() > 0)
         {
            Monomial.push_back(stoi(NewInt));
            NewInt.clear();
         };
      } else
         NewInt += (*it);
   };

   return Result;
};
