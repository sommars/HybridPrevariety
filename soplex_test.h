#include "convex_hull.h"

class MySoPlex : public SoPlex {
private:
	SLUFactor m_slu;
public:
	/// default constructor
	MySoPlex(): SoPlex() {
		setIntParam(ALGORITHM, SPxSolver::LEAVE);
		setIntParam(REPRESENTATION, SPxSolver::COLUMN);
		setIntParam(FACTOR_UPDATE_TYPE, SLUFactor::FOREST_TOMLIN);
		setIntParam(VERBOSITY, 0);
		setIntParam(RATIOTESTER, 0);
		setIntParam(PRICER, 0);
		setIntParam(SIMPLIFIER, 0);
		setIntParam(STARTER, 0);
		setRealParam(EPSILON_UPDATE, DEFAULT_EPS_ZERO);
		setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
	}
	double SoplexTime;
	LPColSetReal Columns;
	bool SoplexSaysIntersect() {
		bool Result;
		SPxSolver::Status stat;
		clock_t StartTime = clock();
		stat = solve();
		SoplexTime += double(clock() - StartTime);
	
		if ( stat == SPxSolver::INFEASIBLE ) {
			Result = false;
		} else {
			Result = true;
		};
		return Result;
	};
	void ClearLP() {
		clearLPReal();
		addColsReal(Columns);
	};
	void AddRowsToLP(LPRowSetReal RowVec) {
		addRowsReal(RowVec);
	//	writeFileReal("dump.lp", NULL, NULL, NULL);
	};
};



//------------------------------------------------------------------------------
//bool SoplexSaysIntersect();

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
static int toint(float r) {
	return *((int*)&r);
}

//------------------------------------------------------------------------------
static void printLP(SPxLP &w) {
	cout << "LP has " << w.nRows() << "\trows and\n       ";
	cout << w.nCols() << "\tcolumns" << endl;
	
	int nr=w.nRows();
	int nc=w.nCols();

	for(int i=0;i<nr;i++) {
		for(int j=0;j<nc;j++) {
			LPRow R;
			w.getRow(i,R);
			cout<<R.rowVector()[j]<<" ";
		}
		cout << endl;
	}
	cout << endl;

	for(int i=0;i<nr;i++) {
		for(int j=0;j<nc;j++) {
			LPCol C;
			w.getCol(j,C);
			cout << C.colVector()[i]<<" ";
		}
		cout << endl;
	}

	cout << "cols:" << endl;

	for(int j=0;j<nc;j++) {
		LPCol C;
		w.getCol(j,C);
		cout << C.lower() << " " << C.upper() << " " << C.obj() << endl;
		cout << toint(C.lower()) << " " << toint(C.upper()) << " " << toint(C.obj()) << endl;
	}

	cout<<"rows:"<<endl;

	for(int i=0;i<nr;i++) {
		LPRow R;
		w.getRow(i,R);
		cout << toint(R.lhs()) << " " << toint(R.rhs()) << " " << R.type() << endl;
	}
}
