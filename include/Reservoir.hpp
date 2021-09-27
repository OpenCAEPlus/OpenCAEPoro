#pragma once
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Connection_BB.hpp"
#include "WellGroup.hpp"
#include "ParamRead.hpp"


class Reservoir
{
	friend class OpenCAEPoro;
	friend class Summary;
	friend class CAEControl;
public:

	void inputParam(ParamRead& param);
	void setup();

	void init();

	double calCFL(double dt);
	// assemble mat
	template<typename T>
	void allocateMat(Solver<T>& mySolver);

	template<typename T>
	void initAssembleMat(Solver<T>& mySolver);

	void assembleMat(Solver<double>& mysolver, double dt);
	void getSol_IMPES(vector<double>& u);

private:
	Grid					grid;
	Bulk					bulk;
	WellGroup				wellgroup;
	Connection_BB			conn;
};


// allocate memory
template<typename T>
void Reservoir::allocateMat(Solver<T>& mySolver)
{
	mySolver.allocate(conn.getActiveBulkNum() + wellgroup.getWellNum());
	conn.allocateMat(mySolver);
	wellgroup.allocateMat(mySolver);
	mySolver.allocateColVal();
}

template<typename T>
void Reservoir::initAssembleMat(Solver<T>& mySolver)
{
	// initialize ColId and DiagPtr
	conn.initAssembleMat(mySolver);
}