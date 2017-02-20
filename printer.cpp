#include "printer.h"

//------------------------------------------------------------------------------
void PrintPoint(vector<int> &Point)
{
   cout << "{ ";
   for (vector<int>::iterator it=Point.begin(); it != Point.end(); it++)
      cout << (*it) << " ";
   cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > &Points)
{
   for (vector<vector<int> >::iterator itr=Points.begin();
        itr != Points.end();
        itr++)
      PrintPoint(*itr);
}

//------------------------------------------------------------------------------
void PrintPoint(set<int> &Point)
{
   set<int>::iterator it;
   cout << "{ ";
   for (it=Point.begin(); it != Point.end(); it++)
      cout << (*it) << " ";
   cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoint(vector<bool> &Point)
{
   vector<bool>::iterator it;
   cout << "{ ";
   for (it=Point.begin(); it != Point.end(); it++)
      cout << (*it) << " ";
   cout << "}" << endl;
};

//------------------------------------------------------------------------------
void PrintPointForPython(vector<int> &Point)
{
   // Print a point such that it can be immediately read in by Python
   cout << "[ ";
   for (vector<int>::iterator it=Point.begin(); it != Point.end(); it++)
      cout << (*it) << ",";
   cout << "]";
}

//------------------------------------------------------------------------------
void PrintPointsForPython(vector<vector<int> > &Points)
{
   cout << "[";
   for (vector<vector<int> >::iterator itr=Points.begin(); 
        itr != Points.end(); 
        itr++)
   {
      PrintPointForPython(*itr);
      cout << ",";
   };
   cout << "]";
}

//------------------------------------------------------------------------------
int toint(float r)
{
   // Helper function for printLP
   return *((int*)&r);
};

//------------------------------------------------------------------------------
void printLP(SPxLP &w)
{
   // Prints the SoPlex LP nicely
   cout << "LP has " << w.nRows() << "\trows and\n       ";
   cout << w.nCols() << "\tcolumns" << endl;
   
   int nr=w.nRows();
   int nc=w.nCols();

   for(int i=0;i<nr;i++)
   {
      for(int j=0;j<nc;j++)
      {
         LPRow R;
         w.getRow(i,R);
         cout<<R.rowVector()[j]<<" ";
      }
      cout << endl;
   }
   cout << endl;

   for(int i=0;i<nr;i++)
   {
      for(int j=0;j<nc;j++)
      {
         LPCol C;
         w.getCol(j,C);
         cout << C.colVector()[i]<<" ";
      }
      cout << endl;
   }

   cout << "cols:" << endl;

   for(int j=0;j<nc;j++)
   {
      LPCol C;
      w.getCol(j,C);
      cout << C.lower() << " " << C.upper() << " " << C.obj() << endl;
      cout << toint(C.lower()) << " " << toint(C.upper());
      cout << " " << toint(C.obj()) << endl;
   }

   cout<<"rows:"<<endl;

   for(int i=0;i<nr;i++)
   {
      LPRow R;
      w.getRow(i,R);
      cout << toint(R.lhs()) << " " << toint(R.rhs());
      cout << " " << R.type() << endl;
   }
}
