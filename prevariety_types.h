#include <ppl.hh>
#include <boost/dynamic_bitset.hpp>
#include "Thread_Pool_defs.hh"
#include "/home/jeff/Desktop/Software/soplex-2.2.1/src/soplex.h"

using namespace std;
using namespace soplex;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
struct BitsetWithCount {
   // Used for relation tables and for keeping track which
   // polytopes have been visited by a polytope.
   // As expected, 1 represents true, 0 represents false.
   boost::dynamic_bitset<> Indices;
   int Count;
};

//------------------------------------------------------------------------------
struct Cone {
   // This is the base object used in computation. They hold relation tables
   // and PPL C_Polyhedron.
   vector<BitsetWithCount> RelationTables;
   BitsetWithCount PolytopesVisited;
   C_Polyhedron HOPolyhedron;
   //LPRowSetReal Rows;
};

//------------------------------------------------------------------------------
struct Edge {
   set<int> PointIndices;
   set<int> NeighborIndices;
};

//------------------------------------------------------------------------------
struct Facet {
   set<int> PointIndices;
};

//------------------------------------------------------------------------------
struct Hull {
   vector<vector<int> > Points;
   map<vector<int>,int> PointToIndexMap;
   map<int,vector<int> > IndexToPointMap;
   vector<Edge> Edges;
   vector<Cone> Cones;
   vector<Facet> Facets;
   C_Polyhedron CPolyhedron;
   int AffineDimension;
   int SpaceDimension;
};

//------------------------------------------------------------------------------
struct ConeWithIndicator {
   vector<int> RayIndices;
   bool IsMaximal;
};

//------------------------------------------------------------------------------
struct TropicalPrevariety {
   map<vector<int>, int> RayToIndexMap;
   vector<set<ConeWithIndicator > > ConeTree;
   vector<vector<int> > Pretropisms; // Consider ripping out
   map<int, Generator> IndexToGenMap;
   vector<int> FVector; // Consider ripping out
};

//------------------------------------------------------------------------------
struct ThreadQueue {
   mutable mutex M;
   vector<list<Cone> > SharedCones;
   ThreadQueue(vector<list<Cone> > ConeVector): SharedCones(ConeVector) {};
   ThreadQueue(const ThreadQueue& TQ) {
      lock_guard<mutex> lock(TQ.M);
      SharedCones = TQ.SharedCones;
   };
};

