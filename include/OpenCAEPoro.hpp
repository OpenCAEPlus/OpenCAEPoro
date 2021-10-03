#pragma once
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "CAEOutput.hpp"
#include "ParamRead.hpp"
#include "Timing.hxx"

class OpenCAEPoro
{
public:
	void inputParam(ParamRead& param);
	void setup();
	void allocateMat();

	void init();

	void run();
	void runIMPES(OCP_DBL& dt);
	void SolveP(OCP_DBL dt);

private:

	// main component
	Reservoir			reservoir;

	// linear solver
	Solver<OCP_DBL>		solver;

	// control
	CAEControl			control;

	// output file
	CAEOutput			output;
};
