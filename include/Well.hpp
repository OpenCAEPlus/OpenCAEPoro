/*! \file    Well.hpp
 *  \brief   Well class declaration
 *  \author  Shizhe Li
 *  \date    Oct/05/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELL_HEADER__
#define __WELL_HEADER__


// Standard header files
#include <cassert>

// OpenCAEPoro header files
#include "Perforation.hpp"
#include "Solver.hxx"
#include "Bulk.hpp"
#include "Grid.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamWell.hpp"

using namespace std;

/// WellOpt describes the operation mode of a well.
/// usually it changes over time, specifically, each attributes could be changed including the well type.
class WellOpt
{
	friend class Well;
public:
	WellOpt() = default;
	WellOpt(const WellOptParam& Optparam);

private:

	USI							Type{ 0 };		///< type of well, Inj or Prod.
	/// indicate which type of fluids will be injected, water, gas, or other solvent.
	/// it's decided by users and only useful for injection well.
	USI							FluidType{ 0 };	
	bool						State{ false }; ///< state of well, close or open.
	USI							OptMode;		///< control mode of well: constant pressure, or constant flow rate of specified fluids.
	/// it gives the upper limit of flow rate of specified fluids if the well is under the control of constant pressure.
	/// it gives the flow rate of specified fluids if the well is under the control of constant flow rate.
	OCP_DBL						MaxRate;	
	/// used for injection well.
	/// it gives the upper limit of well pressure if the well is under the control of constant flow rate.
	/// it gives the pressure of well if the well is under the control of constant pressure. 
	OCP_DBL						MaxBHP;			
	/// used for production well.
	/// it gives the lower limit of well pressure if the well is under the control of constant flow rate.
	/// it gives the pressure of well if the well is under the control of constant pressure. 
	OCP_DBL						MinBHP;
	/// it's decided by users input.
	/// for injection well, it describes the components of injected fluids.
	/// for production well, it gives the the components of fluids which we are interested in.
	vector<OCP_DBL>			Zi;
};

/// Well class defines well, and any operations referred to wells are in it.
/// Well connects to the bulks by perforations, which serve as source and sink.
/// Due to practical difficulties in production, a good treatment for well is important,
/// excellent treatment will make the flow rate in well more stable. 
class Well
{
	friend class WellGroup;
public:
	Well() = default;
	/// return the state of the well, Open or Close.
	bool WellState() const { return Opt.State; }
	/// return the type of well, Inj or Prod.
	USI  WellType() const { return Opt.Type; }

	/// setup the well, it will be called when Grid and Bulk finish setupping.
	void setup(const Grid& myGrid, const Bulk& myBulk);
	/// cal Well Index with Peaceman model for vertical well.
	void calWI_Peaceman_Vertical(const Bulk& myBulk);

	/// guess the well pressure at the beginning of simulation.
	/// usually the pressure equals the ones in topest bulk which connects to the well.
	void init(const Bulk& myBulk);
	/// calculate the CFL number, only parts related to wells are considered.
	OCP_DBL calCFL(const Bulk& myBulk, const OCP_DBL& dt) const;

	/// calculate pressure difference between well and perforations.
	/// it calculates pressure difference between perforations iteratively.
	/// this function can be used in both black oil model and compositional model. stability of this method shoule be tested.
	void caldG(const Bulk& myBulk);
	/// calculate pressure difference between well and perforations for injection.
	void calInjdG(const Bulk& myBulk);
	/// calculate pressure difference between well and perforations for prodcution.
	void calProddG(const Bulk& myBulk);
	// test
	/// try to smooth the dG by average it with dG at last time step.
	/// it's just a test now to make dG more stable.
	void smoothdG();
	/// calculate transmissibility for each phase in perforations.
	void calTrans(const Bulk& myBulk);
	/// calculate the flow rate of moles of components and total flow rate of volume in each perforations.
	void calFlux(const Bulk& myBulk, const bool flag = false);
	/// update moles of components in those bulks who connects to the well.
	void massConserve(Bulk& myBulk, const OCP_DBL& dt) const;
	

	/// calculate flow rate of moles of components for injection well in black oil model,
	/// where pressure in injection well equals minial ones in injection well, which is input by users.
	/// this function is used to check if operation mode of well shoubld be swtched.
	OCP_DBL calInjRate_blk(const Bulk& myBulk);
	/// calculate flow rate of moles of components for production well in black oil model,
	/// where pressure in production well equals minial ones in production well, which is input by users.
	/// this function is used to check if operation mode of well shoubld be swtched.
	OCP_DBL calProdRate_blk(const Bulk& myBulk);
	/// calculate flow rate of moles of components for injection well in black oil model.
	void calInjqi_blk(const Bulk& myBulk, const OCP_DBL& dt);
	/// calculate flow rate of moles of components for production well in black oil model.
	void calProdqi_blk(const Bulk& myBulk, const OCP_DBL& dt);


	/// check if well operation mode would be changed.
	/// constant well pressure would be applied if flow rate is too large.
	/// constant flow rate would be applied if well pressure is outranged.
	void checkOptMode(const Bulk& myBulk);

	// Assemble Mat
	/// allocate memory for matrix.
	template <typename T>
	void allocateMat(Solver<T>& mySolver) const;
	/// assemble matrix, parts related to injection well are included.
	void assembleMat_INJ_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, const OCP_DBL& dt) const;
	/// assemble matrix, parts related to production well are included.
	/// this function can only applied in Black Oil model now.
	void assembleMat_PROD_BLK_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, const OCP_DBL& dt) const;

	/// update pressure in Perforation after well pressure updates.
	void updatePerfP(){ for (USI p = 0; p < PerfNum; p++) Perf[p].P = BHP + dG[p]; }
    /// check if abnormal Pressure occurs.
	OCP_INT checkP(const Bulk& myBulk);
	/// check if crossflow happens.
	OCP_INT checkCrossFlow(const Bulk& myBulk);

	/// display operation mode of well and state of perforations.
	void showPerfStatus() const;

private:

	OCP_DBL						Radius;			///< well radius.
	OCP_DBL						Kh;				///< effective permeability times net thickness of the connection.
	OCP_DBL						SkinFactor;		///< skin factor.
	OCP_DBL						WI;				///< connection transmissibility factor, it can be provided directly from the users.
	
	string						Name;			///< well name
	USI							I;				///< I-index of the well header.
	USI							J;				///< J-index of the well header.
	USI							K1;				///< K-location of upper connecting block in this set of data.
	USI							K2;				///< K-location of lower connecting block in this set of data.
	string						Direction;		///< direction of well: x, y, z.
	WellOpt						Opt;			///< well control parameters, contains current control parameters.
	vector<WellOpt>		OptSet;			///< well control parameters set, contains control parameters in all critical time.

	
	OCP_DBL						BHP;			///< well pressure in reference depth.
	OCP_DBL						Depth;			///< reference depth of well.
	USI							PerfNum;		///< num of perforations belonging to this well.
	vector<Perforation>	Perf;			///< information of perforation belonging to this well.
	vector<OCP_DBL>			dG;			///< difference of pressure between well and perforation: PerfNum.
	vector<OCP_DBL>			ldG;		///< difference of pressure between well and perforation at last time step: PerfNum.

	// production rate and injection rate
	vector<OCP_DBL>			Qi_lbmol;	///< flow rate of moles of component inflowing/outflowing well: num of components.
	OCP_DBL						WOPR{ 0 };		///< well oil production rate.
	OCP_DBL						WOPT{ 0 };		///< well total oil production.
	OCP_DBL						WGPR{ 0 };		///< well gas production rate.
	OCP_DBL						WGPT{ 0 };		///< well total gas production.
	OCP_DBL						WWPR{ 0 };		///< well water production rate.
	OCP_DBL						WWPT{ 0 };		///< well total water production.
	OCP_DBL						WGIR{ 0 };		///< well gas injection rate.
	OCP_DBL						WGIT{ 0 };		///< well total gas injection.
	OCP_DBL						WWIR{ 0 };		///< well water injection rate.
	OCP_DBL						WWIT{ 0 };		///< well total water injection.

};

template <typename T>
void Well::allocateMat(Solver<T>& mySolver) const
{
	for (USI p = 0; p < PerfNum; p++) {
		mySolver.RowCapacity[Perf[p].Location]++;
	}
}

#endif