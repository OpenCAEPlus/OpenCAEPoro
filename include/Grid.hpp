/*! \file    Grid.hpp
 *  \brief   Grid class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __GRID_HEADER__
#define __GRID_HEADER__

// Standard header files
#include <vector>
#include <iostream>

// OpenCAEPoro header files
#include "ParamReservoir.hpp"
#include "OpenCAEPoro_consts.hpp"

using namespace std;

/// GB_Pair contains two variables, which indicates if a grid is active and what its active index is if active.
class GB_Pair
{
public:
	GB_Pair() = default;
	GB_Pair(bool act, OCP_USI i) : activity(act), index(i) {};
	/// return activity of some grid.
	bool getAct() const  { return activity; }
	/// return active index of some grid if active.
	OCP_USI getId() const { return index; }
private:
	bool		activity;		///< activity of grid
	OCP_USI		index;			///< active index of grid if active
};

/// Grid contains basic information of grids of reservoir, the rock properties in grids are
/// also included. all grids are stored here, you can regard it as a database of reservoir.
/// considering the existence of inactive grids(whose volume of pores is too small or effective resource is too little)
/// or switch of activity of grid, Grid class is necessary. Grid class is static while simulating, active grids
/// will be stored in bulks, which is "area" for calculating.
class Grid
{	
	friend class Bulk;
	friend class Connection_BB;
	friend class Well;
public:
	Grid() = default;

	void setup();

	OCP_USI getBulkNum() { return Num; }
	OCP_USI getConnNum() { return ConnNum; }
	OCP_USI getActiveBulkNum() { return ActiveBulkNum; }

	void calDepthV();
	void calActiveBulk(OCP_DBL e1, OCP_DBL e2);		// fill ActiveMap_B2G and ActiveMap_G2B

	void inputParam(const ParamReservoir& rs_param);

	OCP_USI getIndex(USI i, USI j, USI k);

private:
	USI					Nx;					///< num of bulks along x-direction
	USI					Ny;					///< num of bulks along y-direction
	USI					Nz;					///< num of bulks along z-direction
	OCP_USI				Num;				///< num of grids, Nx * Ny * Nz
	OCP_USI				ConnNum;			///< num of connection


	vector<OCP_DBL>		Tops;				///< depth of top face of topest gird: Nx*Ny
	vector<OCP_DBL>		Depth;				///< depth of center of grid: Num.
	vector<OCP_DBL>		Dx;					///< size of gird along the x direction: Num.
	vector<OCP_DBL>		Dy;					///< size of grid along the y direction: Num.
	vector<OCP_DBL>		Dz;					///< size of grid along the z direction: Num.
	vector<OCP_DBL>		V;					///< volume of grid: Num.
	vector<OCP_DBL>		Ntg;				///< net to gross of grid: Num
	vector<OCP_DBL>		Poro;				///< initial porosity of rock: Num
	vector<OCP_DBL>		Kx;					///< Absolute permeability of rock along x direction: Num
	vector<OCP_DBL>		Ky;					///< Absolute permeability of rock along y direction: Num
	vector<OCP_DBL>		Kz;					///< Absolute permeability of rock along z direction: Num

	// Region
	vector<USI>			SATNUM;				///< used to identify SAT region: Num.
	vector<USI>			PVTNUM;             ///< used to identify PVT region in blackoil model: Num.

	OCP_USI				ActiveBulkNum;		///< num of active grid.
	vector<OCP_USI>		ActiveMap_B2G;		///< a index map form active grid to grid: ActiveBulkNum.
	vector<GB_Pair>		ActiveMap_G2B;		///< a index map form grid to active grid: Num.
	
};


#endif


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/