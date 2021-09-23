#pragma once
#include "MAT.hpp"

template <typename T>
class Solver
{
	friend class OpenCAEPoro;
	friend class Connection_BB;
	friend class Well;
public:
	void allocate(int dimMax);
	void allocateColVal();
	void assembleSolverA();

	void clearData();

private:
	int									MaxDim;
	int									Dim;
	// for Mat assemble
	std::vector<int>					RowCapacity;
	std::vector<std::vector<int>>		ColId;
	std::vector<std::vector<T>>			Val;
	std::vector<int>					DiagPtr;
	std::vector<T>						DiagVal;  // just for assembling: intermediate variable
	// for solver
	MAT<T>								A;


	std::vector<T>						b;
	std::vector<T>						u;

};


template <typename T>
void Solver<T>::allocate(int dimMax)
{
	MaxDim = dimMax;
	RowCapacity.resize(MaxDim, 0);
	ColId.resize(MaxDim);
	Val.resize(MaxDim);
	DiagPtr.resize(MaxDim, 0);
	DiagVal.resize(MaxDim, 0);
	b.resize(MaxDim, 0);
	u.resize(MaxDim, 0);
}

template <typename T>
void Solver<T>::allocateColVal()
{
	for (int n = 0; n < MaxDim; n++) {
		ColId[n].reserve(RowCapacity[n]);
		Val[n].reserve(RowCapacity[n]);
	}
}

template <typename T>
void Solver<T>::assembleSolverA()
{

}

template <typename T>
void Solver<T>::clearData()
{
	for (int i = 0; i < MaxDim; i++) {
		ColId[i].clear();
		Val[i].clear();
	}
	DiagPtr.resize(MaxDim, -1);
	DiagVal.resize(MaxDim, 0);
	b.resize(MaxDim, 0);
}
