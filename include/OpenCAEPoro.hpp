#pragma once
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "OpenCAEOutput.hpp"
#include "ParamRead.hpp"
#include "Method.hpp"
#include "Timing.hxx"



class OpenCAEPoro
{
public:
	void inputParam(ParamRead& param);
	void setup(ParamRead& param);
	void setupSolver();

	void init();
	void run();
	void out();

private:

	// main component
	Reservoir			reservoir;

	// Method
	OCP_IMPES			Impes;

	OCP_FIM				Fim;

	// control
	OCP_Control			control;

	// output file
	OCP_Output			output;
};
