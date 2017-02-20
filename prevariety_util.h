// Assorted utility functions for computing prevarieties

#include "printer.h"

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator &g, bool KnockOffLastTerm);

//------------------------------------------------------------------------------
Constraint InequalityToStrictInequality(Constraint &c);

//------------------------------------------------------------------------------
Constraint InequalityToEquation(Constraint &c);

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint &c);

//------------------------------------------------------------------------------
double DoubleInnerProduct(vector<int> &V1, vector<double> &V2);

//------------------------------------------------------------------------------
vector<vector<vector<int> > > ParseToSupport(string &Input);
