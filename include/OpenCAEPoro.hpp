#pragma once
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "CAEOutput.hpp"
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

	// control
	CAEControl			control;

	// output file
	CAEOutput			output;
};
