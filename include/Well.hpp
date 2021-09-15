#pragma once
#include "Perforation.hpp"
#include "Solver.hpp"
#include "Bulk.hpp"
#include "Grid.hpp"
#include "OpenCAEPoro_consts.hpp"
#include <cassert>

using namespace std;

class Well
{
	friend class WellGroup;
public:
	Well() = default;
	void setState(bool flag) { State = flag; };


	// cal Well Index
	void calWI_Peaceman_Vertical(const Bulk& myBulk);

	// Assemble Mat
	void allocateMat(Solver& mySolver);
	void assembleMat_INJ(const Bulk& myBulk, Solver& mySolver);
	void assembleMat_PROD_BLK(const Bulk& myBulk, Solver& mySolver);

private:

	double						Radius;			// well radius
	double						SkinFactor;		// skin factor

	int							Type;			// inj or prod
	int							State;			// open or close
	int							FluidType;		// inj/prod type
	int							Direction;		// direction of well
	int							OptModeInit;	// initial opt mode
	int							OptMode;		// the control mode
	double						OptValue;		// corresponding values
	double						MaxRate;
	double						MaxBHP;
	double						MinBHP;
	std::vector<double>			Zi;				// inj for inj Well, prod for prod well

	double						BHP;			// pressure in reference depth
	double						Depth;			// reference depth
	int							PerfNum;
	std::vector<Perforation>	Perf;
	std::vector<double>			dG;

};


