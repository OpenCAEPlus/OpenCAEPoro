#include "Solver.hpp"

void Solver::allocate(int dimMax)
{
	MaxDim = dimMax;
	RowCapacity.resize(MaxDim, 0);
	ColId.resize(MaxDim);
	Val.resize(MaxDim);
	DiagPtr.resize(MaxDim, 0);
	DiagVal.resize(MaxDim, 0);
	b.resize(MaxDim, 0);
}

void Solver::allocateColVal()
{
	for (int n = 0; n < MaxDim; n++) {
		ColId[n].reserve(RowCapacity[n]);
		Val[n].reserve(RowCapacity[n]);
	}
}

void Solver::assembleSolverA()
{

}


void Solver::clearData()
{
	for (int i = 0; i < MaxDim; i++) {
		ColId[i].clear();
		Val[i].clear();
	}
	DiagPtr.resize(MaxDim, -1);
	DiagVal.resize(MaxDim, 0);
	b.resize(MaxDim, 0);
}
