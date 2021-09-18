#pragma once
#include "Reservoir.hpp"
#include "ParamRead.hpp"

class OpenCAEPoro
{
public:
	void inputParam(ParamRead& param);
	void setup();

private:

	// main component
	Reservoir			reservoir;

	// linear solver
	Solver				solver;

	// I/O
	string				inputfile;

};
