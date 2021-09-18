#pragma once
#include "Perforation.hpp"
#include "Solver.hpp"
#include "Bulk.hpp"
#include "Grid.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamWell.hpp"
#include <cassert>

using namespace std;

class WellOpt
{
public:
	WellOpt() = default;
	WellOpt(WellOptParam& Optparam);

	int							Type{ -1 };
	int							FluidType{ -1 };
	int							State{ -1 };
	int							OptMode{ -1 };
	double						OptValue;
	double						MaxRate;
	double						MaxBHP;
	double						MinBHP;
	std::vector<double>			Zi;			// inj or prod
};

class Well
{
	friend class WellGroup;
public:
	Well() = default;
	void setState(bool flag) { Opt.State = flag; };

	void setupPerf();
	// cal Well Index
	void calWI_Peaceman_Vertical(const Bulk& myBulk);

	// Assemble Mat
	void allocateMat(Solver& mySolver);
	void assembleMat_INJ(const Bulk& myBulk, Solver& mySolver);
	void assembleMat_PROD_BLK(const Bulk& myBulk, Solver& mySolver);

private:

	double						Radius;			// well radius
	double						Kh;
	double						SkinFactor;		// skin factor
	double						Trans;
	
	string						Name;
	int							I, J;
	int							K1, K2;
	string						Direction;		// direction of well: x, y, z
	int							OptModeInit;	// initial opt mode
	WellOpt						Opt;
	std::vector<WellOpt>		OptSet;

	
	double						BHP;			// pressure in reference depth
	double						Depth;			// reference depth
	int							PerfNum;
	std::vector<Perforation>	Perf;
	std::vector<double>			dG;

};


