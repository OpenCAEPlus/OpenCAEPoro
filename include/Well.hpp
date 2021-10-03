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
	OCP_DBL						MaxRate;
	OCP_DBL						MaxBHP;
	OCP_DBL						MinBHP;
	std::vector<OCP_DBL>			Zi;			// inj or prod
};

class Well
{
	friend class WellGroup;
public:
	Well() = default;

	bool WellState() const { return Opt.State; }
	int  WellType() const { return Opt.Type; }

	void setup(Grid& myGrid, Bulk& myBulk);
	// cal Well Index
	void calWI_Peaceman_Vertical(const Bulk& myBulk);

	// init
	void init(const Bulk& myBulk);

	OCP_DBL calCFL(const Bulk& myBulk, OCP_DBL dt);

	// calculate dG
	void caldG(const Bulk& myBulk);
	void calInjdG(const Bulk& myBulk);
	void calProddG(const Bulk& myBulk);
	// test
	void smoothdG();

	void calTrans(const Bulk& myBulk);
	// cal flux ---- perf: qt_ft3 & qi_lbmol
	void calFlux(const Bulk& myBulk, bool flag = false);
	void massConserve(Bulk& myBulk, OCP_DBL dt);
	

	// calculate rate -- zi, uesd to check well opt mode and calculate well rate
	OCP_DBL calInjRate_blk(const Bulk& myBulk);
	OCP_DBL calProdRate_blk(const Bulk& myBulk);
	void calInjqi_blk(const Bulk& myBulk, OCP_DBL dt);
	void calProdqi_blk(const Bulk& myBulk, OCP_DBL dt);



	// check optmode
	void checkOptMode(const Bulk& myBulk);

	// Assemble Mat
	template <typename T>
	void allocateMat(Solver<T>& mySolver) const;
	void assembleMat_INJ(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, OCP_DBL dt) const;
	void assembleMat_PROD_BLK(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, OCP_DBL dt) const;


	void updatePerfP(){ for (int p = 0; p < PerfNum; p++) Perf[p].P = BHP + dG[p]; }
    int checkP(const Bulk& myBulk);
	int checkCrossFlow(const Bulk& myBulk);

	// show info
	void showPerfStatus();

private:

	OCP_DBL						Radius;			// well radius
	OCP_DBL						Kh;
	OCP_DBL						SkinFactor;		// skin factor
	OCP_DBL						WI;				// connection factor
	
	string						Name;
	int							I, J;
	int							K1, K2;
	string						Direction;		// direction of well: x, y, z
	WellOpt						Opt;
	std::vector<WellOpt>		OptSet;

	
	OCP_DBL						BHP;			// pressure in reference depth
	OCP_DBL						Depth;			// reference depth
	int							PerfNum;
	std::vector<Perforation>	Perf;
	std::vector<OCP_DBL>			dG;
	std::vector<OCP_DBL>			ldG;

	// production rate and injection rate
	std::vector<OCP_DBL>			Qi_lbmol;
	OCP_DBL						WOPR{ 0 }, WOPT{ 0 };
	OCP_DBL						WGPR{ 0 }, WGPT{ 0 };
	OCP_DBL						WWPR{ 0 }, WWPT{ 0 };
	OCP_DBL						WGIR{ 0 }, WGIT{ 0 };
	OCP_DBL						WWIR{ 0 }, WWIT{ 0 };

};

template <typename T>
void Well::allocateMat(Solver<T>& mySolver) const
{
	for (int p = 0; p < PerfNum; p++) {
		mySolver.RowCapacity[Perf[p].Location]++;
	}
}
