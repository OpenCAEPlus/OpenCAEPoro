#pragma once
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Connection_BB.hpp"
#include "WellGroup.hpp"




class Reservoir
{
	
	// assemble mat
	void allocateMat(Solver& mySolver);
	void initAssembleMat(Solver& mySolver);
	void assembleMat(Solver& mysolver);


private:
	Grid					grid;
	Bulk					bulk;
	WellGroup				wellgroup;
	Connection_BB			conn;
};
