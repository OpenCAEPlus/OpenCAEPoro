#include "Reservoir.hpp"

inline void Reservoir::allocateMat(Solver& mySolver)
{
	mySolver.allocate(conn.getActiveBulkNum() + wellgroup.getWellNum());
	conn.allocateMat(mySolver);
	wellgroup.allocateMat(mySolver);
	mySolver.allocateColVal();
}

// allocate mat memory
inline void Reservoir::initAssembleMat(Solver& mySolver) 
{
	conn.initAssembleMat(mySolver);
}

// assemble mat
inline void Reservoir::assembleMat(Solver& mysolver) 
{
	conn.assembleMat(mysolver, bulk);
	wellgroup.assemblaMat_WB(mysolver, bulk);
}
