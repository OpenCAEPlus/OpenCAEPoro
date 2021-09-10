#pragma once
// #include "ParamReservoir.h"
#include "MyReservoir.h"

//struct CAEPoroParam {
//	string title, dir, filename;
//	ParamReservoir rs_param;
//	ControlParam control_param;
//};


class OpenCAEPoro
{
public:


private:

	// main component
	MyReservoir			reservoir;

	// linear solver
	Solver				solver;

	// I/O
	string				inputfile;


};
