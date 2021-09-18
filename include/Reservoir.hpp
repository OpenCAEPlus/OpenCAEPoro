#pragma once
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Connection_BB.hpp"
#include "WellGroup.hpp"
#include "ParamRead.hpp"


class Reservoir
{
public:

	void inputParam(ParamRead& param);
	void setup();

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
