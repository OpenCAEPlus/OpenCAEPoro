#pragma once
#include "Reservoir.hpp"



class OpenCAEPoro
{
public:


private:

	// main component
	Reservoir			reservoir;

	// linear solver
	Solver				solver;

	// I/O
	string				inputfile;


};
