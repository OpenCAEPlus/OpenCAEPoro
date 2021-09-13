#include "Reservoir.hpp"

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
