#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>
#include <string>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}
int ThreeDimCount = 0;
//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < g.space_dimension(); i++) {
		stringstream s;
		s << (g).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
vector<int> GeneratorToIntPoint(Generator g) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < g.space_dimension(); i++) {
		stringstream s;
		s << (g).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
Constraint InequalityToStrictInequality(Constraint c) {
	Linear_Expression LE;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		LE += c.coefficient(Variable(i)) * Variable(i);
	};
	Constraint c2(LE > c.inhomogeneous_term());
	vector<int> CP1 = ConstraintToPoint(c);
	vector<int> CP2 = ConstraintToPoint(c2);
	return c2;
}
//------------------------------------------------------------------------------
Constraint StrictInequalityToNonStrictInequality(Constraint c) {
	Linear_Expression LE;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		LE += c.coefficient(Variable(i)) * Variable(i);
	};
	Constraint c2(LE >= c.inhomogeneous_term());
	vector<int> CP1 = ConstraintToPoint(c);
	vector<int> CP2 = ConstraintToPoint(c2);
	return c2;
}

//------------------------------------------------------------------------------
Constraint InequalityToEquation(Constraint c) {
	Linear_Expression LE;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		LE += c.coefficient(Variable(i)) * Variable(i);
	};
	Constraint c2(LE == c.inhomogeneous_term());
	return c2;
}

//------------------------------------------------------------------------------
vector<vector<int> > GeneratorSystemToPoints(Generator_System gs) {
	vector<vector<int> > Result;
	for (Generator_System::const_iterator i = gs.begin(),
	gs_end = gs.end(); i != gs_end; ++i) {
		Result.push_back(GeneratorToPoint(*i));
	}
	return Result;
}

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		stringstream s;
		s << (c).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<int> > Points, vector<double> VectorForOrientation, bool Verbose) {
	Hull H;
	H.CPolyhedron = FindCPolyhedron(Points);
	
	Constraint_System cstemp = H.CPolyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cstemp.begin(),
	cs1_end = cstemp.end(); i != cs1_end; ++i) {
		if ((*i).is_equality()) {
			H.Eqs.push_back(*i);
		};
	};
	map<double,vector<int> > DoubleToPt;
	vector<double> IPs;
	for (size_t i = 0; i != Points.size(); i++) {
		double IP = DoubleInnerProduct(Points[i], VectorForOrientation);
		IPs.push_back(IP);
		DoubleToPt[IP] = Points[i];
	};
	sort(IPs.begin(), IPs.end());
	for (size_t i = 0; i != IPs.size(); i++) {
		H.Points.push_back(DoubleToPt[IPs[i]]);
	};

	H.AffineDimension = H.CPolyhedron.affine_dimension();
	H.SpaceDimension = H.CPolyhedron.space_dimension();
		
	// Create PointToIndexMap
	for (size_t i = 0; i != H.Points.size(); i++) {
		vector<int> Point = H.Points[i];
		H.PointToIndexMap[Point]=i;
		H.IndexToPointMap[i]=Point;
	};

	FindFacets(H);
	FindEdges(H);
	vector<Cone> HalfOpenCones;
	for (size_t i = 0; i != H.Points.size(); i++) {
		vector<Constraint> Constraints;
		vector<int> Pt = H.Points[i];
		int PtIndex = H.PointToIndexMap[Pt];
		vector<int> Edges;
		vector<int> Facets;
		// Go through all of the edges. If the edge is not on the facet.
		for (size_t j = 0; j != H.Edges.size(); j++) {
			Edge E = H.Edges[j];
			if (E.PointIndices.find(PtIndex) != E.PointIndices.end()) {
				Edges.push_back(j);
				set<int>::iterator PtIter;
				vector<int> OtherPt;
				int OtherPtIndex;
				for (PtIter=E.PointIndices.begin(); PtIter != E.PointIndices.end(); PtIter++) {
					OtherPtIndex = (*PtIter);
					if (OtherPtIndex != PtIndex) {
						OtherPt = H.IndexToPointMap[(*PtIter)];
					};
				};
				
				Linear_Expression LE;
				for (size_t k = 0; k != Pt.size(); k++) {
					LE += (OtherPt[k] - Pt[k]) * Variable(k);
				};
				
				//Here we manufacture the constraint from the edge.
				Constraint c;
				if (OtherPtIndex > PtIndex) {
					c = (LE >= 0);
				} else {
					c = (LE > 0);
				};
				Constraints.push_back(c);
			};
		};
		//If no non-strict inequalities exist, then you may simply throw the cone away.
		//This process produces just one half open cone for each edge.
		for (size_t j = 0; j != Constraints.size(); j++) {
			//Take one of your non-strict inequalities for the cone.
			Constraint TempConstraint = Constraints[j];
			if (!TempConstraint.is_strict_inequality()) {
				//Introduce it as an equation to get a lower dimensional cone, which you output.
				Constraints[j] = InequalityToEquation(TempConstraint);
				Constraint_System cs1;
				for (size_t k = 0; k != Constraints.size(); k++) {
					cs1.insert(Constraints[k]);
				};
				//Introduce it as a strict inequality to describe the rest. Call recursively on this cone.
				Cone NewCone;
				NewCone.HOPolyhedron = NNC_Polyhedron(cs1);
				HalfOpenCones.push_back(NewCone);
				Constraints[j] = InequalityToStrictInequality(TempConstraint);
			};
		};

	};
	for (size_t i = 0; i != HalfOpenCones.size(); i++) {
		for (size_t j = i+1; j != HalfOpenCones.size(); j++) {
//			if (!HalfOpenCones[i].HOPolyhedron.is_disjoint_from(HalfOpenCones[j].HOPolyhedron)) {
			if (IntersectCones(HalfOpenCones[i].HOPolyhedron,HalfOpenCones[j].HOPolyhedron).affine_dimension() > 0) {
				cout << "Internal Error: Two half open cones from the same polytope are not disjoint." << endl;
				cin.get();
			};
		};
	};
	for (size_t i = 0; i != HalfOpenCones.size(); i++) {
		//take a random vector from cone
		vector<int> RandomVector(H.SpaceDimension, 0);
		for (Generator_System::const_iterator ii = HalfOpenCones[i].HOPolyhedron.minimized_generators().begin(), gs_end = HalfOpenCones[i].HOPolyhedron.minimized_generators().end(); ii != gs_end; ++ii) {
			if (!(*ii).is_ray() && !(*ii).is_line()) {
				continue;
			};
			for (size_t j = 0; j != H.SpaceDimension; j++) {
				stringstream s;
				s << (*ii).coefficient(Variable(j));
				int ToAppend;
				istringstream(s.str()) >> ToAppend;
				RandomVector[j] += ToAppend;
			};
		};
		//take initial form using that vector.
		vector<vector<int> > InitialForm = FindInitialForm(H.Points, RandomVector);
		if (InitialForm.size() != 2) {
			cout << "Internal Error: initial form is not of length two. It is of size " << InitialForm.size() << endl;
			PrintPoint(RandomVector);
			cout << endl;
			PrintPoints(InitialForm);
			cin.get();
		};
		set<int> InitialIndices;
		for (vector<vector<int> >::iterator InitialFormItr=InitialForm.begin();
		InitialFormItr != InitialForm.end(); InitialFormItr++) {
			InitialIndices.insert(H.PointToIndexMap[*InitialFormItr]);
		};
		for (size_t k = 0; k != H.Edges.size(); k++) {
			if (H.Edges[k].PointIndices == InitialIndices) {
				H.Edges[k].EdgeCone = HalfOpenCones[i];
				Constraint_System csClosed;
				for (Constraint_System::const_iterator csc = HalfOpenCones[i].HOPolyhedron.minimized_constraints().begin(), cs_end = HalfOpenCones[i].HOPolyhedron.minimized_constraints().end(); csc != cs_end; ++csc) {
					Constraint cc = (*csc);
					if (cc.is_strict_inequality()) {
						csClosed.insert(StrictInequalityToNonStrictInequality(cc));
					} else {
						csClosed.insert(cc);
					};
				};
				
				
				Constraint_System csTest;
				int PtIndex = *(H.Edges[k].PointIndices.begin());
				int PtIndex2 = *(H.Edges[k].PointIndices.rbegin());
				vector<int> Pt11 = H.IndexToPointMap[PtIndex];
				vector<int> Pt22 = H.IndexToPointMap[PtIndex2];
				Linear_Expression LEe;
				for (size_t kk = 0; kk != Pt11.size(); kk++) {
					LEe += (Pt22[kk] - Pt11[kk]) * Variable(kk);
				};
				csTest.insert(LEe == 0);
				for (size_t q = 0; q != H.Edges.size(); q++) {
					if (k == q)
						continue;
					
					int OtherPtIndex = -1;
					if (*(H.Edges[q].PointIndices.begin()) == PtIndex)
						OtherPtIndex = *(H.Edges[q].PointIndices.rbegin());
					
					if (*(H.Edges[q].PointIndices.rbegin()) == PtIndex)
						OtherPtIndex = *(H.Edges[q].PointIndices.begin());
						
					if (OtherPtIndex == -1)
						continue;
					vector<int> Pt1 = H.IndexToPointMap[PtIndex];
					vector<int> Pt2 = H.IndexToPointMap[OtherPtIndex];
//					cout << "HERE!" << endl;
//					PrintPoint(Pt1);
//					PrintPoint(Pt2);
					Linear_Expression LE;
					for (size_t kk = 0; kk != Pt1.size(); kk++) {
						LE += (Pt2[kk] - Pt1[kk]) * Variable(kk);
					};
//					cout << "CONSTRAINT! " << (Constraint(LE >= 0)) << endl;
					csTest.insert(Constraint(LE >= 0));
				};
//				cout << "NEW" << endl;
//				if (NNC_Polyhedron(csTest).minimized_constraints() != NNC_Polyhedron(csClosed).minimized_constraints()) {
//					cout << NNC_Polyhedron(csTest).minimized_constraints() << endl;
//					cout << NNC_Polyhedron(csClosed).minimized_constraints() << endl;
//					cin.get();
//				};
	/*			
				cout << "NEW" << endl;
				cout << NNC_Polyhedron(csTest).minimized_constraints() << endl;
				cout << NNC_Polyhedron(csTest).minimized_generators() << endl;
				cout << NNC_Polyhedron(csTest).affine_dimension() << endl;
		*/		
				
				
				
//				H.Edges[k].EdgeCone.ClosedPolyhedron = NNC_Polyhedron(csClosed);
/*				H.Edges[k].EdgeCone.ClosedPolyhedron = NNC_Polyhedron(csTest);
				H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_constraints();
				H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_generators();
				H.Edges[k].EdgeCone.ClosedPolyhedron.affine_dimension();*/
			/*	cout << "OLD" << endl;
				cout << H.Edges[k].EdgeCone.ClosedPolyhedron.affine_dimension() << endl;
				cout << H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_constraints() << endl;
				cout << H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_generators() << endl;
				cout << "above" << endl;
				cin.get();*/
				
				//H.Edges[k].EdgeCone.HOPolyhedron = H.Edges[k].EdgeCone.ClosedPolyhedron;
				
				//H.Edges[k].EdgeCone.HOPolyhedron = NNC_Polyhedron(csClosed);
				H.Edges[k].EdgeCone.HOPolyhedron.minimized_constraints();
				H.Edges[k].EdgeCone.HOPolyhedron.minimized_generators();
				H.Edges[k].EdgeCone.HOPolyhedron.affine_dimension();
				H.Edges[k].HasEdgeCone = true;
				
				break;
			} else if (k == H.Edges.size() - 1) {
				cout << "Internal Error: edge not matched with cone." << endl;
				cin.get();
			};
		};
	};

	for (size_t i = 0; i != H.Edges.size(); i++) {
		if (!H.Edges[i].HasEdgeCone) {
			cout << "Internal Error: edge never found a cone" << endl;
			cin.get();
		};
	};
	if (Verbose) {
		cout << "Convex hull------------------------" << endl;
		PrintPoints(H.Points);
		cout << "Affine dimension: " << H.AffineDimension << endl;
		cout << "Space dimension: " << H.SpaceDimension << endl;
		cout << "Number of edges: " << H.Edges.size() << endl;
		cout << "Number of facets: " << H.Facets.size() << endl << endl;
	};
	return H;
}

//------------------------------------------------------------------------------
void FindFacets(Hull &H) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		
		if (!(*i).is_inequality()) {
			continue;
		};
		vector<int> Pt = ConstraintToPoint(*i);
		vector<vector<int> > FacetPts = FindInitialForm(H.Points, Pt);
		Facet F;
		vector<vector<int> >::iterator itr;
		
		for (itr=FacetPts.begin(); itr != FacetPts.end(); itr++) {
			F.PointIndices.insert(H.PointToIndexMap[*itr]);
		};
		H.Facets.push_back(F);
	}
}

//------------------------------------------------------------------------------
void FindEdges(Hull &H) {
	vector<vector<int> > CandidateEdges = FindCandidateEdges(H);

	//The number of facets we want is equal to the dimension of the ambient space minus the number of equations -1
	int Dim = H.CPolyhedron.affine_dimension() - 1;

	vector<vector<int> >::iterator itr;
	for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++) {
		vector<int> CandidateEdge = (*itr);
		int Point1 = CandidateEdge[0];
		int Point2 = CandidateEdge[1];
	
		int FacetCount = 0;
		for (vector<Facet>::iterator FacetIt = H.Facets.begin(); FacetIt != H.Facets.end(); FacetIt++) {
			set<int> PtIndices = (*FacetIt).PointIndices;
			bool Point1IsInFacet = PtIndices.find(Point1) != PtIndices.end();
			bool Point2IsInFacet = PtIndices.find(Point2) != PtIndices.end();
			if (Point1IsInFacet and Point2IsInFacet) {
				FacetCount++;
			};
		};

		if (FacetCount >= Dim) {
			Edge NewEdge;
			NewEdge.PointIndices.insert(Point1);
			NewEdge.PointIndices.insert(Point2);
			H.Edges.push_back(NewEdge);
		}
	};
	// After all of the edges have been generated, fill out all of the neighbors on all of the edges.
	for (size_t Edge1Index = 0; Edge1Index != H.Edges.size(); Edge1Index++) {
		Edge Edge1 = H.Edges[Edge1Index];

		for (size_t Edge2Index = 0; Edge2Index != H.Edges.size(); Edge2Index++) {
			if (Edge1Index == Edge2Index) {
				continue;
			};

			int IntersectionCount = IntersectSets(Edge1.PointIndices, H.Edges[Edge2Index].PointIndices).size();
			if (IntersectionCount == 1) {
				H.Edges[Edge1Index].NeighborIndices.insert(Edge2Index);
				H.Edges[Edge2Index].NeighborIndices.insert(Edge1Index);
			} else if (IntersectionCount > 1) {
				cout << "Internal Error: Two edges intersect at more than one point" << endl;
				cin.get();
			}
		};
	};
}

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull H) {
	vector<vector<int> > CandidateEdges;
	int n = H.Points.size();
	vector<int> d(n);
	for (size_t i = 0; i != d.size(); ++i) {
		d[i] = i;
	}
	do {
		if (d[0] < d[1]) {
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
double DoubleInnerProduct(vector<int> V1, vector<double> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
		cin.get();
	};
	double Result = 0;
	for (size_t i = 0; i != V1.size(); i++) {
		Result += V1[i] * V2[i];
	}
	return Result;
}

//------------------------------------------------------------------------------
NNC_Polyhedron FindCPolyhedron(vector<vector<int> > Points) {
	Generator_System gs;
	vector<vector<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		vector<int> Point = *itr;
		Linear_Expression LE;
		for (size_t i = 0; i != Point.size(); i++) {
			LE = LE + Variable(i) * (Point[i]);
		};
		gs.insert(point(LE));
	};
	NNC_Polyhedron ph = NNC_Polyhedron(gs);
	
	return ph;
}

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point) {
	vector<int>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points) {
	vector<vector<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		PrintPoint(*itr);
	};
}

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point) {
	set<int>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPointForPython(vector<int> Point) {
	vector<int>::iterator it;
	cout << "[ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << ",";
	}
	cout << "]";
}

//------------------------------------------------------------------------------
void PrintPointsForPython(vector<vector<int> > Points) {
	vector<vector<int> >::iterator itr;
	cout << "[";
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		PrintPointForPython(*itr);
		cout << ",";
	};
	cout << "]";
}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > ParseToSupport(string Input) {
	vector<vector<vector<int> > > Result;
	int ParenCount = 0;
	string NewInt;
	// Takes in strings that look like this:
	//'[[[1,0,0,0][0,1,0,0][0,0,1,0][0,0,0,1]][[1,1,0,0][0,1,1,0][1,0,0,1][0,0,1,1]][[1,1,1,0][1,1,0,1][1,0,1,1][0,1,1,1]]]'
	vector<vector<int> > Polynomial;
	vector<int> Monomial;
	for (string::iterator it=Input.begin(); it!=Input.end(); ++it) {
		if ((*it) == '[') {
			ParenCount++;
		} else if ((*it) == ']') {
			ParenCount--;
			if (NewInt.length() > 0) {
				Monomial.push_back(stoi(NewInt));
				NewInt.clear();
			};
			if (ParenCount == 2) {
				Polynomial.push_back(Monomial);
				Monomial.clear();
			} else if (ParenCount == 1) {
				Result.push_back(Polynomial);
				Polynomial.clear();
			};
		} else if ((*it) == ',') {
			if (NewInt.length() > 0) {
				Monomial.push_back(stoi(NewInt));
				NewInt.clear();
			};
		} else {
			NewInt += (*it);
		};
	};

	return Result;
};
