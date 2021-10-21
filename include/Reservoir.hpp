/*! \file    Reservoir.hpp
 *  \brief   Reservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__

// OpenCAEPoro header files
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "Grid.hpp"
#include "ParamRead.hpp"
#include "WellGroup.hpp"

/// Reservoir is the core component in our simulator, it contains the all reservoir
/// information, and all operations on it.
///
/// Reservoir has four Core components.
/// Grids contains the basic informations of all grids as a database of reservoir.
/// Bulk only stores active grids, which defines the area used for calculation.
/// WellGroup contains the well information, it's used to manage operations related to
/// wells. Connection_BB contains connections between bulks(active grids).
class Reservoir
{
    friend class OpenCAEPoro;
    friend class OCP_Control;
    friend class OCP_IMPES;
    friend class Summary;
    friend class CriticalInfo;
    friend class DetailInfo;

public:
    /// Input param from internal param data structure, which stores the params from
    /// input files.
    void InputParam(ParamRead& param);
    /// Setup static information for reservoir with input params.
    void Setup();
    /// Initialize the reservoir, actually it gives the first step in iterations.
    void Init();
    /// Calcluate the CFL number, including bulks and wells.
    OCP_DBL CalCFL(const OCP_DBL& dt);
    /// Calcluate the CFL number, including bulks and wells.
    OCP_DBL CalCFL01(const OCP_DBL& dt);
    /// Allocate memory for linear system, it should be called at the beginning of
    /// simulation only once. It's accessible for both IMPES and FIM.
    void AllocateMat(LinearSolver& mySolver) const;
    /// assemble the matrix
    /// Setup most of sparsity pattern first, and then Setup the value only related to
    /// the bulks. finally, assemble the parts related to wells, which will complete the
    /// rest sparsity pattern simultaneously
    void AssembleMat(LinearSolver& mysolver, const OCP_DBL& dt) const;
    /// get the solution from LinearSolver after the linear system is solved.
    void GetSolution_IMPES(const vector<OCP_DBL>& u);
    /// check if abnormal pressure occurs including pressure in bulks, wells,
    /// perforations. if so, take corresponding measures and then resolve the linear
    /// equations.
    OCP_INT CheckP();
    /// check if mole of components occurs
    /// if so, cut the timestep, reset with function ResetVal01 and resolve the linear
    /// equtions.
    bool CheckNi() const { return bulk.CheckNi(); }
    /// reset pressure, capillary pressure, flux.
    void ResetVal();
    /// reset pressure, capillary pressure, moles of componnets, flux.
    void ResetVal01();
    /// check if relative error between fluids volume and pore volume is too large.
    /// if so, cut the timestep, reset with function resetval02 and resolve the linear
    /// equtions.
    bool CheckVe(const OCP_DBL& Vlim) const { return bulk.CheckVe(Vlim); }
    /// reset pressure, capillary pressure, flux, moles of components, volume of pores.
    void ResetVal02();

private:
    Grid          grid;      ///< Grid class.
    Bulk          bulk;      ///< Bulk class.
    WellGroup     wellgroup; ///< WellGroup class.
    Connection_BB conn;      ///< Connection_BB class.

    OCP_DBL cfl; ///< CFL number.
};


#endif /* end if __RESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/