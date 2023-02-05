/*! \file    Solver.hpp
 *  \brief   Solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoro header files
#include "IsothermalSolver.hpp"
#include "OCPOutput.hpp"
#include "ThermalSolver.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

/// Solver class for overall solution methods.
class Solver
{
public:
    /// Setup Solver
    void Setup(Reservoir& rs, const OCPControl& ctrl);
    /// Initialize the reservoir.
    void InitReservoir(Reservoir& rs) const;
    /// Start simulation.
    void RunSimulation(Reservoir& rs, OCPControl& ctrl, OCPOutput& output);

private:
    /// General API
    void GoOneStep(Reservoir& rs, OCPControl& ctrl);

    /// Setup solver for isothermal model
    void SetupIsoT(Reservoir& rs, const OCPControl& ctrl);
    /// Run one time step for isothermal model
    void GoOneStepIsoT(Reservoir& rs, OCPControl& ctrl);

    /// Setup solver for ifThermal model
    void SetupT(Reservoir& rs, const OCPControl& ctrl);
    /// Run one time step for ifThermal model
    void GoOneStepT(Reservoir& rs, OCPControl& ctrl);

private:
    /// Model
    USI OCPModel{ISOTHERMALMODEL};
    /// Solver for isothermal models with fixed T
    IsothermalSolver IsoTSolver;
    /// Solver for ifThermal models with varied T
    ThermalSolver TSolver;
};

#endif /* end if __SOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Shizhe Li           Oct/21/2021      Change from OCPMethod to Solver      */
/*  Chensong Zhang      Oct/27/2021      Rearrange and add comments           */
/*----------------------------------------------------------------------------*/