#include "printer.h"

//------------------------------------------------------------------------------
void PrintPoint(vector<int> &Point)
{
   cout << "{ ";
   for (vector<int>::iterator it=Point.begin(); it != Point.end(); it++) {
      cout << (*it) << " ";
   }
   cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > &Points)
{
   for (vector<vector<int> >::iterator itr=Points.begin(); itr != Points.end(); itr++)
      PrintPoint(*itr);
}

//------------------------------------------------------------------------------
void PrintPoint(set<int> &Point)
{
   set<int>::iterator it;
   cout << "{ ";
   for (it=Point.begin(); it != Point.end(); it++) {
      cout << (*it) << " ";
   }
   cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoint(vector<bool> &Point)
{
   vector<bool>::iterator it;
   cout << "{ ";
   for (it=Point.begin(); it != Point.end(); it++) {
      cout << (*it) << " ";
   }
   cout << "}" << endl;
};

//------------------------------------------------------------------------------
void PrintPointForPython(vector<int> &Point)
{
   cout << "[ ";
   for (vector<int>::iterator it=Point.begin(); it != Point.end(); it++) {
      cout << (*it) << ",";
   }
   cout << "]";
}

//------------------------------------------------------------------------------
void PrintPointsForPython(vector<vector<int> > &Points)
{
   cout << "[";
   for (vector<vector<int> >::iterator itr=Points.begin(); itr != Points.end(); itr++) {
      PrintPointForPython(*itr);
      cout << ",";
   };
   cout << "]";
}
