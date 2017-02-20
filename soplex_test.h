// Experimental unit, needs rewrite if ever moved into production.
// This is a first crack at interfacing SoPlex with the rest of the software

#include "convex_hull.h"

class MySoPlex : public SoPlex {
private:
   SLUFactor m_slu;
public:
   MySoPlex(): SoPlex()
   {
      setIntParam(ALGORITHM, SPxSolver::LEAVE);
      setIntParam(REPRESENTATION, SPxSolver::COLUMN);
      setIntParam(FACTOR_UPDATE_TYPE, SLUFactor::FOREST_TOMLIN);
      setIntParam(VERBOSITY, 0);
      setIntParam(RATIOTESTER, 0);
      setIntParam(PRICER, 0);
      setIntParam(SIMPLIFIER, 0);
      setIntParam(STARTER, 0);
      setRealParam(EPSILON_UPDATE, DEFAULT_EPS_ZERO);
      setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   }
   double SoplexTime;
   LPColSetReal Columns;
   bool SoplexSaysIntersect()
   {
      bool Result;
      SPxSolver::Status stat;
      clock_t StartTime = clock();
      stat = solve();
      SoplexTime += double(clock() - StartTime);
   
      if ( stat == SPxSolver::INFEASIBLE )
         Result = false;
      else
         Result = true;
      return Result;
   };
   void ClearLP()
   {
      clearLPReal();
      addColsReal(Columns);
   };
   void AddRowsToLP(LPRowSetReal RowVec)
   {
      addRowsReal(RowVec);
   };
};
