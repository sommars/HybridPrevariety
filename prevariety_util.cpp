#include "prevariety_util.h"

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g, bool KnockOffLastTerm) {
	// Converts a PPL Generator object to a vector of int
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
