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
	void allocateMat();

	void init();

	void run();
	void runIMPES(double& dt);
	void SolveP(double dt);

private:

	// main component
	Reservoir			reservoir;

	// linear solver
	Solver<double>		solver;

	// control
	CAEControl			control;

	// output file
	CAEOutput			output;
};
