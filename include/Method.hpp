#pragma once

// OpenCAEPoro header files
#include "Solver.hxx"
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "OpenCAEOutput.hpp"
#include "Timing.hxx"

/// OCP_IMPES is IMPES(implict pressure explict saturation) method in our simulator(OpenCAEPoro),
/// which consists of the functions in Reservoir.
class OCP_IMPES
{
public:
	void setupParam(const string& dir, const string& file);
	void allocateMat(const Reservoir& rs);
	void run(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output);
	void goOneStep(Reservoir& rs, OCP_Control& ctrl);
	void SolveP(Reservoir& rs,OCP_Control& ctrl, const OCP_DBL& dt);

private:
	Solver<OCP_DBL>		        solver;
};


/// (to do) OCP_FIM is FIM(fully implict method) method in our simulator(OpenCAEPoro),
/// which consists of the functions in Reservoir.
class OCP_FIM
{
public:

private:
	Solver<vector<OCP_DBL>>		solver;
};