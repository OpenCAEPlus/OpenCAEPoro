#pragma once
#include "Perforation.hpp"
#include "Solver.hxx"
#include "Bulk.hpp"
#include "Grid.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamWell.hpp"
#include <cassert>

using namespace std;

class WellOpt
{
	friend class WellGroup;    // temp
	friend class Well;
public:
	WellOpt() = default;
	WellOpt(WellOptParam& Optparam);

private:

	int							Type{ -1 };
	int							FluidType{ -1 };
	bool						State{ false };
	int							OptMode{ -1 };
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

	bool WellState() { return Opt.State; }
	int  WellType() { return Opt.Type; }

	void setup(Grid& myGrid, Bulk& myBulk);
	// cal Well Index
	void calWI_Peaceman_Vertical(const Bulk& myBulk);

	// init
	void init(const Bulk& myBulk);

	// calculate dG
	void calInjdG(const Bulk& myBulk);
	void calProddG(const Bulk& myBulk);

	// calculate rate -- zi, uesd to check well opt mode
	double calInjRate_blk(const Bulk& myBulk);
	double calProdRate_blk(const Bulk& myBulk);

	// check optmode
	void checkOptMode(const Bulk& myBulk);

	// Assemble Mat
	template <typename T>
	void allocateMat(Solver<T>& mySolver);
	void assembleMat_INJ(const Bulk& myBulk, Solver<double>& mySolver, double dt);
	void assembleMat_PROD_BLK(const Bulk& myBulk, Solver<double>& mySolver, double dt);

private:

	double						Radius;			// well radius
	double						Kh;
	double						SkinFactor;		// skin factor
	double						WI;				// connection factor
	
	string						Name;
	int							I, J;
	int							K1, K2;
	string						Direction;		// direction of well: x, y, z
	WellOpt						Opt;
	std::vector<WellOpt>		OptSet;

	
	double						BHP;			// pressure in reference depth
	double						Depth;			// reference depth
	int							PerfNum;
	std::vector<Perforation>	Perf;
	std::vector<double>			dG;

	// production rate and injection rate
	double						WOPR;
	double						WGPR;
	double						WWPR;
	double						WGIR;
	double						WWIR;

};

template <typename T>
void Well::allocateMat(Solver<T>& mySolver)
{
	for (int p = 0; p < PerfNum; p++) {
		mySolver.RowCapacity[Perf[p].Location]++;
	}
}
