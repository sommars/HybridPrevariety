#include "prevariety_util.h"

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g, bool KnockOffLastTerm) {
	vector<int> Result;
	int Dim = 0;
	if (KnockOffLastTerm)
		Dim = 1;
	for (size_t i = 0; i < g.space_dimension() - Dim; i++) {
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
Constraint InequalityToEquation(Constraint c) {
	Linear_Expression LE;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		LE += c.coefficient(Variable(i)) * Variable(i);
	};
	Constraint c2(LE == c.inhomogeneous_term());
	return c2;
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
vector<Cone> NewHull(vector<vector<int> > Points, vector<double> VectorForOrientation, bool Verbose) {
	Hull H;
	H.CPolyhedron = FindCPolyhedron(Points);
	vector<Cone> OutputCones;
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
		// Go through all of the edges. If the edge is not on the facet.
		for (size_t j = 0; j != H.Edges.size(); j++) {
			if (H.Edges[j].PointIndices.find(PtIndex) != H.Edges[j].PointIndices.end()) {
				set<int>::iterator PtIter;
				vector<int> OtherPt;
				int OtherPtIndex;
				for (PtIter=H.Edges[j].PointIndices.begin(); PtIter != H.Edges[j].PointIndices.end(); PtIter++) {
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

				// Now we raise to a higher dimension to represent as a standard linear programming problem.
				Constraint_System csHigherDim;
				csHigherDim.insert(Constraint(Variable(H.SpaceDimension) >= 1));
				for (size_t k = 0; k != Constraints.size(); k++) {
					if (Constraints[k].is_strict_inequality()) {
						Linear_Expression LE;
						for (size_t l = 0; l < Constraints[k].space_dimension(); l++) {
							LE += Constraints[k].coefficient(Variable(l)) * Variable(l);
						};
						LE += -1 * Variable(Constraints[k].space_dimension());
						csHigherDim.insert(Constraint(LE >= 0));
					} else {
						csHigherDim.insert(Constraints[k]);
					};

				};
				Cone NewCone;
				NewCone.HOPolyhedron = C_Polyhedron(csHigherDim);
				NewCone.HOPolyhedron.minimized_constraints();
				NewCone.HOPolyhedron.minimized_generators();
				NewCone.HOPolyhedron.affine_dimension();
				
				for (Constraint_System::const_iterator cc = NewCone.HOPolyhedron.constraints().begin(), cs_end = NewCone.HOPolyhedron.constraints().end(); cc != cs_end; cc++) {
					DSVector NewRow(cc->space_dimension());
					for (size_t jj = 0; jj < cc->space_dimension(); jj++) {
						stringstream s;
						s << cc->coefficient(Variable(jj));
						int ToAppend;
						istringstream(s.str()) >> ToAppend;
						NewRow.add(jj,(double)ToAppend);
					};
					if (cc->is_equality()) {
						NewCone.Rows.add(LPRow(0, NewRow, 0));
					} else {
						stringstream s;
						s << cc->inhomogeneous_term();
						int ToAppend;
						istringstream(s.str()) >> ToAppend;
						NewCone.Rows.add(LPRow(-1 * ToAppend, NewRow, infinity));
					};
				};
				
				
				/*
				cout << NewCone.Constraints << endl;
				for (MIP_Problem::const_iterator i = NewCone.MP.constraints_begin(); i != NewCone.MP.constraints_end(); i++)
					cout << (*i) << endl;
				cout << endl << endl;
				cin.get();
				*/
				OutputCones.push_back(NewCone);
				//Introduce it as a strict inequality to describe the rest. Call recursively.
				Constraints[j] = InequalityToStrictInequality(TempConstraint);
			};
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
	return OutputCones;
}

//------------------------------------------------------------------------------
void FindFacets(Hull &H) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		
		if (!i->is_inequality()) {
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
			set<int> PtIndices = FacetIt->PointIndices;
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
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points) {
	Generator_System gs;
	vector<vector<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		Linear_Expression LE;
		for (size_t i = 0; i != itr->size(); i++) {
			LE += Variable(i) * ((*itr)[i]);
		};
		gs.insert(point(LE));
	};
	return C_Polyhedron(gs);
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
