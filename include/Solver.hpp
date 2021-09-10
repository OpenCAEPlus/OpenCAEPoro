#pragma once
#include "MAT.hpp"

class Solver
{
	friend class Connection_BB;
	friend class Well;
public:
	void assembleSolverA();
	void clearData();

private:
	int									Dim;
	// for Mat assemble
	std::vector<std::vector<int>>		ColId;
	std::vector<std::vector<double>>	Val;
	std::vector<int>					DiagPtr;
	std::vector<double>					DiagVal;  // just for assembling: intermediate variable
	// for solver
	MAT									A;


	std::vector<double>					b;
	std::vector<double>					u;

};
