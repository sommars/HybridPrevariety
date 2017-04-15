// Various print methods that are used or were used at one point.

#include "polynomial_systems.h"

//------------------------------------------------------------------------------
void PrintPoint(vector<int> &Point);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > &Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> &Point);

//------------------------------------------------------------------------------
void PrintPoint(vector<bool> &Point);

//------------------------------------------------------------------------------
void PrintPointForPython(vector<int> &Point);

//------------------------------------------------------------------------------
void PrintPointsForPython(vector<vector<int> > &Points);

//------------------------------------------------------------------------------
void printLP(SPxLP &w);

//------------------------------------------------------------------------------
void PrintMaximalCones(TropicalPrevariety &TP, stringstream &s);

//------------------------------------------------------------------------------
void StreamPoint(set<int> &Point, stringstream &s);

//------------------------------------------------------------------------------
void PrintRT(vector<BitsetWithCount> &RT);
