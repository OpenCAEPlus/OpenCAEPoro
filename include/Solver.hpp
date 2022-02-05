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
#include "FluidSolver.hpp"
#include "OCPOutput.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

/// Solver class for overall solution method.
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
    /// Run one time step.
    void GoOneStep(Reservoir& rs, OCPControl& ctrl);
    /// Setup the solution method.
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);
    /// Before solve: prepare for assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble and Solve: assemble linear system parts together then solve.
    void AssembleSolve(Reservoir& rs, OCPControl& ctrl);
    /// Update reservoir properties after solving for primary variables.
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish the Newton-Raphson iteration.
    bool FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish the current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

private:
    FluidSolver  FSolver;
    LinearSystem LSolver;
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
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/