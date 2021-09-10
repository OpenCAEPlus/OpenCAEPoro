#pragma once
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Connection_BB.hpp"
#include "WellGroup.hpp"




class Reservoir
{
	
	// assemble mat
	void initAssembleMat(Solver& mySolver);
	void AssembleMat(Solver& mysolver);

private:
	Grid					grid;
	Bulk					bulk;
	WellGroup				wellgroup;
	Connection_BB			conn;
};
