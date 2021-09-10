#include "MyReservoir.h"

// allocate mat memory
inline void MyReservoir::initAssembleMat(Solver& mySolver) 
{
	conn.initAssembleMat(mySolver, wellgroup.getWellNum());
}

// assemble mat
inline void MyReservoir::AssembleMat(Solver& mysolver) 
{
	conn.assembleMat(mysolver, bulk);
	wellgroup.assemblaMat_WB(mysolver, bulk);
}
