#include "Reservoir.hpp"

// allocate mat memory
inline void Reservoir::initAssembleMat(Solver& mySolver) 
{
	conn.initAssembleMat(mySolver, wellgroup.getWellNum());
}

// assemble mat
inline void Reservoir::AssembleMat(Solver& mysolver) 
{
	conn.assembleMat(mysolver, bulk);
	wellgroup.assemblaMat_WB(mysolver, bulk);
}
