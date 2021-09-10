#pragma once
#include "Grid.h"
#include "Bulk.h"
#include "Connection_BB.h"
#include "WellGroup.h"




class MyReservoir
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
