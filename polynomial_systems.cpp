#include "polynomial_systems.h"

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced) {
	vector<vector<vector<int> > > System;
	int Length = n;
	if (Reduced == true) {
		Length = n - 1;
	};
	
	for (size_t i = 0; i != n-1; i++) {
		vector<vector<int> > Equation;
		for (size_t j = 0; j != n; j++) {
			vector<int> Monomial;
			set<int> OneSet;
			for (size_t k = 0; k != i+1; k++) {
				OneSet.insert((j+k)%n);
			};
			for (size_t l = 0; l != Length; l++) {
				if (OneSet.find(l) != OneSet.end()) {
					Monomial.push_back(1);
				} else {
					Monomial.push_back(0);
				};
			};
			Equation.push_back(Monomial);
		};
		System.push_back(Equation);
	};
	if (Reduced == false) {
		vector<vector<int> > Equation;
		vector<int> Monomial1;
		vector<int> Monomial2;
		for (size_t i = 0; i != n; i++) {
			Monomial1.push_back(1);
			Monomial2.push_back(0);
		};
		Equation.push_back(Monomial1);
		Equation.push_back(Monomial2);
		System.push_back(Equation);
	};
	return System;
}

//------------------------------------------------------------------------------
vector<vector<vector<int> > > RandomSimplices(int n) {
	vector<vector<vector<int> > > System;
	while (true) {
		vector<vector<int> > Equation;
		for (size_t i = 0; i != n + 1; i++) {
			vector<int> Monomial;
			for (size_t j = 0; j != n; j++) {
				Monomial.push_back(rand() % 10);
			};
			Equation.push_back(Monomial);
		};
		
		Generator_System gs;
		vector<vector<int> >::iterator itr;
		for (itr=Equation.begin(); itr != Equation.end(); itr++) {
			vector<int> Point = *itr;
			Linear_Expression LE;
			for (size_t i = 0; i != Point.size(); i++) {
				LE = LE + Variable(i) * (Point[i]);
			};
			gs.insert(point(LE));
		};
		C_Polyhedron ph = C_Polyhedron(gs);

		if (ph.affine_dimension() == ph.space_dimension()) {
			System.push_back(Equation);
			if (System.size() == n - 1) {
				break;
			};
		};
	};
	return System;
}
