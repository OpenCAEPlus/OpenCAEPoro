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

	OCP_DBL calCFL(OCP_DBL dt);
	// assemble mat
	template<typename T>
	void allocateMat(Solver<T>& mySolver);

	template<typename T>
	void initAssembleMat(Solver<T>& mySolver);

	void assembleMat(Solver<OCP_DBL>& mysolver, OCP_DBL dt);
	void getSol_IMPES(vector<OCP_DBL>& u);

	int checkP();
	bool checkNi() { return bulk.checkNi(); }
	void resetVal01();
	bool checkVe(const OCP_DBL Vlim) { return bulk.checkVe(Vlim); }
	void resetVal02();

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