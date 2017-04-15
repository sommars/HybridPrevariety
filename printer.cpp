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

//------------------------------------------------------------------------------
void PrintMaximalCones(TropicalPrevariety &TP, stringstream &s)
{
   if (TP.ConeTree.size() == 0)
      return;
   
   // Prints maximal cones from prevariety object.
   for (size_t i = 0; i != TP.ConeTree.size() - 1; i++)
   {
      if (TP.ConeTree[i].size() > 0)
         s << "------ Cones of dimension " << i + 1 << " ------"<< endl;
      for (size_t j = 0; j != TP.ConeTree[i].size(); j++)
      {
         if (TP.ConeTree[i][j].Status != 1)
            continue;
         set<int>::iterator it;
         s << "{ ";
         for (it=TP.ConeTree[i][j].RayIndices.begin();
              it != TP.ConeTree[i][j].RayIndices.end();
              it++)
            s << (*it) << " ";
         s << "}" << endl;
      };
   };
   int Dim = TP.ConeTree.size() - 1;
   
   if (Dim < 0)
      return;
      
   s << "------ Cones of dimension " << Dim + 1 << " ------"<< endl;
   
   for (size_t j = 0; j != TP.ConeTree[Dim].size(); j++)
   {
      set<int>::iterator it;
      s << "{ ";
      for (it=TP.ConeTree[Dim][j].RayIndices.begin();
           it != TP.ConeTree[Dim][j].RayIndices.end();
           it++)
         s << (*it) << " ";
      s << "}" << endl;
   };
   
}

//------------------------------------------------------------------------------
void StreamPoint(set<int> &Point, stringstream &s)
{
   set<int>::iterator it;
   s << "{ ";
   for (it = Point.begin(); it != Point.end(); it++)
      s << (*it) << " ";
   s << "}" << endl;
};

//------------------------------------------------------------------------------
void PrintRT(vector<BitsetWithCount> &RT)
{

cout << endl;
   for (size_t i = 0; i != RT.size(); i++)
   {
      cout << RT[i].Indices << endl;
   };
cout << endl;
};
