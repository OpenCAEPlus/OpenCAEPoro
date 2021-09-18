#include "Reservoir.hpp"

void Reservoir::inputParam(ParamRead& param)
{
	grid.inputParam(param.Rs_param);
	bulk.inputParam(param.Rs_param);
	wellgroup.inputParam(param.Well_param);
}

void Reservoir::setup()
{
	grid.setup();
	bulk.setup(grid);
	conn.setup(grid, bulk);
	wellgroup.setup(grid, bulk);
}

// allocate memory
void Reservoir::allocateMat(Solver& mySolver)
{
	mySolver.allocate(conn.getActiveBulkNum() + wellgroup.getWellNum());
	conn.allocateMat(mySolver);
	wellgroup.allocateMat(mySolver);
	mySolver.allocateColVal();
}


void Reservoir::initAssembleMat(Solver& mySolver) 
{
	// initialize ColId and DiagPtr
	conn.initAssembleMat(mySolver);
}

// assemble mat
void Reservoir::assembleMat(Solver& mysolver) 
{
	conn.assembleMat(mysolver, bulk);
	wellgroup.assemblaMat_WB(mysolver, bulk);
}
