#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>
#include <string>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
NNC_Polyhedron IntersectCones(NNC_Polyhedron &ph1, NNC_Polyhedron &ph2) {
	Constraint_System cs1 = ph1.minimized_constraints();
	Constraint_System cs2 = ph2.minimized_constraints();
	for (Constraint_System::const_iterator i = cs1.begin(),
	cs1_end = cs1.end(); i != cs1_end; ++i) {
		cs2.insert(*i);
	}
	Recycle_Input dummy;
	NNC_Polyhedron ph(cs2,dummy);
	ph.affine_dimension();
	return ph;
}

//------------------------------------------------------------------------------
Cone IntersectCones(Cone C1, Cone C2) {
	cout << "Need to reimplement this." << endl;
	cin.get();
/*	
	Constraint_System cs1 = C1.Polyhedron.minimized_constraints();
	Constraint_System cs2 = C2.Polyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs1.begin(),
	cs1_end = cs1.end(); i != cs1_end; ++i) {
		cs2.insert(*i);
	}
	Recycle_Input dummy;
	NNC_Polyhedron ph(cs2,dummy);
	ph.affine_dimension();
	Cone C3;
	C3.Polyhedron = ph;
	*/
//	return C3;
}

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
Hull NewHull(vector<vector<int> > Points) {
	Hull H;
	H.CPolyhedron = FindCPolyhedron(Points);
	
	// We put the points in an order here such that the order can later
	// be used to make a unique sink orientation of the polytope.
	vector<double> VectorForOrientation;
	for (size_t i = 0; i != Points[0].size(); i++) {
		VectorForOrientation.push_back(rand());
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
			if (!HalfOpenCones[i].HOPolyhedron.is_disjoint_from(HalfOpenCones[j].HOPolyhedron)) {
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
				H.Edges[k].EdgeCone.ClosedPolyhedron = NNC_Polyhedron(csClosed);
				H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_constraints();
				H.Edges[k].EdgeCone.ClosedPolyhedron.minimized_generators();
				H.Edges[k].EdgeCone.ClosedPolyhedron.affine_dimension();
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

	Generator_System gs;
	Linear_Expression LE;
	LE = 2*(Variable(0)) + 1*Variable(2) + 1*Variable(3) + 2*Variable(4) + 1*Variable(5) + 1*Variable(6);
	gs.insert(ray(LE));
	gs.insert(point(0*Variable(0)));
	NNC_Polyhedron TestPoly(gs);
	for (size_t i = 0; i != H.Edges.size(); i++) {
		if (!H.Edges[i].HasEdgeCone) {
			cout << "Internal Error: edge never found a cone" << endl;
			cin.get();
		};
	};
	cout << "Convex hull------------------------" << endl;
	PrintPoints(H.Points);
	cout << "Affine dimension: " << H.AffineDimension << endl;
	cout << "Space dimension: " << H.SpaceDimension << endl;
	cout << "Number of edges: " << H.Edges.size() << endl;
	cout << "Number of facets: " << H.Facets.size() << endl << endl;
	
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
int InnerProduct(vector<int> V1, vector<int> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
		cin.get();
	};
	int Result = 0;
	for (size_t i = 0; i != V1.size(); i++) {
		Result += V1[i] * V2[i];
	}
	return Result;
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
vector<vector<int> > FindInitialForm(vector<vector<int> > &Points, vector<int> &Vector) {
	/*
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};
	vector<vector<int> > InitialForm;

	InitialForm.push_back(Points[0]);
	int MinimalIP = InnerProduct(Vector, Points[0]);

	for (size_t i = 1; i != Points.size(); i++) {
		vector<int> Point = Points[i];
		int IP = InnerProduct(Vector, Point);
		if (MinimalIP > IP) {
			MinimalIP = IP;
			InitialForm.clear();
			InitialForm.push_back(Point);
		} else if (IP == MinimalIP) {
			InitialForm.push_back(Point);
		};
	};

	return InitialForm;
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
void PrintCPolyhedron(NNC_Polyhedron ph, bool PrintIf0Dim) {
	cout << "NNC_Polyhedron dimension:" << ph.affine_dimension() << endl;
	cout << "Generators: " << ph.minimized_generators() << endl;
}

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<NNC_Polyhedron> phs, bool PrintIf0Dim) {
	vector<NNC_Polyhedron>::iterator it;
	for (it=phs.begin(); it != phs.end(); it++) {
		PrintCPolyhedron(*it, PrintIf0Dim);
	};
}

//------------------------------------------------------------------------------
set<int> IntersectSets(set<int> S1, set<int> S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	set<int> Result;
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			Result.insert(*S1Itr);
			S1Itr++;
			S2Itr++;
		};
	};
	return Result;
}

//------------------------------------------------------------------------------
bool SetsDoIntersect(set<int> &S1, set<int> &S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			return true;
		};
	};
	return false;
}
