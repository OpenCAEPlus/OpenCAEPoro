/*! \file    Reservoir.hpp
 *  \brief   Reservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/07/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__


#include "Grid.hpp"
#include "Bulk.hpp"
#include "Connection_BB.hpp"
#include "WellGroup.hpp"
#include "ParamRead.hpp"

/// Reservoir is the core component in our simulator, it contains the all reservoir information, 
/// and all operations on it.
/// 
/// Reservoir has four Core components.
/// Grids contains the basic informations of all grids as a database of reservoir.
/// Bulk only stores active grids, which defines the area used for calculation.
/// WellGroup contains the well information, it's used to manage operations related to wells.
/// Connection_BB contains connections between bulks(active grids).
class Reservoir
{
	friend class OpenCAEPoro;
	friend class Summary;
	friend class OCP_Control;
	friend class OCP_IMPES;
public:
	/// input param from internal param data structure, which stores the params from input files.
	void inputParam(ParamRead& param);
	/// setup static information for reservoir with input params.
	void setup();
	/// initialize the reservoir, actually it gives the first step in iterations.
	void init();
	/// calcluate the CFL number, including bulks and wells. 
	OCP_DBL calCFL(const OCP_DBL& dt) const;
	/// allocate memory for linear system, it should be called at the beginning of simulation only once.
	/// it's accessible for both IMPES and FIM.
	template<typename T>
	void allocateMat(Solver<T>& mySolver) const;
	/// assemble the matrix
	/// setup most of sparsity pattern first, and then setup the value only related to the bulks.
	/// finally, assemble the parts related to wells, which will complete the rest sparsity pattern simultaneously
	void assembleMat(Solver<OCP_DBL>& mysolver, const OCP_DBL& dt) const;
	/// get the solution from Solver after the linear system is solved.
	void getSol_IMPES(const vector<OCP_DBL>& u);
	/// check if abnormal pressure occurs including pressure in bulks, wells, perforations.
	/// if so, take corresponding measures and then resolve the linear equations.
	OCP_INT checkP();
	/// check if mole of components occurs
	/// if so, cut the timestep, reset with function resetVal01 and resolve the linear equtions.
	bool checkNi() const { return bulk.checkNi(); }
	/// reset pressure, capillary pressure, flux.
	void resetVal01();
	/// check if relative error between fluids volume and pore volume is too large.
	/// if so, cut the timestep, reset with function resetval02 and resolve the linear equtions.
	bool checkVe(const OCP_DBL& Vlim) { return bulk.checkVe(Vlim); }
	/// reset pressure, capillary pressure, flux, moles of components, volume of pores.
	void resetVal02();

private:
	Grid					grid;			///< Grid class.
	Bulk					bulk;			///< Bulk class.
	WellGroup				wellgroup;		///< WellGroup class.
	Connection_BB			conn;			///< Connection_BB class.
};


// allocate memory
template<typename T>
void Reservoir::allocateMat(Solver<T>& mySolver) const
{
	mySolver.allocate(conn.getActiveBulkNum() + wellgroup.getWellNum());
	conn.allocateMat(mySolver);
	wellgroup.allocateMat(mySolver);
	mySolver.allocateColVal();
}


#endif