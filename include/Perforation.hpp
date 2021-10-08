#pragma once

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "OpenCAEPoro_consts.hpp"

class Perforation
{
	friend class Well;
public:
	 Perforation() = default;
	 void setState(bool flag) { State = flag; };

private:
	bool					State;
	int						Location;
	OCP_DBL					Depth;
	OCP_DBL					Trans;
	OCP_DBL					P;
	OCP_DBL					WI;
	OCP_DBL					Multiplier;

	mutable OCP_DBL					Xi;			// inj  single phase
	std::vector<OCP_DBL>		qi_lbmol;
	std::vector<OCP_DBL>		transj;
	OCP_DBL					qt_ft3;
};


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/