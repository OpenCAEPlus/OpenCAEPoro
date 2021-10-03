#pragma once
#include "Solver.hxx"
#include "Reservoir.hpp"
#include "OpenCAEControl.hpp"
#include "OpenCAEOutput.hpp"
#include "Timing.hxx"


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


class OCP_FIM
{
public:

private:
	Solver<vector<OCP_DBL>>		solver;
};