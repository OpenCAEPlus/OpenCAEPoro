#pragma once
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "OpenCAEOutput.hpp"
#include "ParamRead.hpp"
#include "Method.hpp"
#include "Timing.hxx"


/// OpenCAEPoro is the topest structure in our simulator, in which there are Reservoir class,
/// which stroes all the properties of reservoir, and Method class such as IMPES, and OCP_Control class
/// which manages the params of method and time step, and OCP_output class, in which you can output the
/// results you are interested in.
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

	/// the core component.
	Reservoir			reservoir;

	// Method
	OCP_IMPES			Impes;

	OCP_FIM				Fim;

	// control
	OCP_Control			control;

	// output file
	OCP_Output			output;
};
