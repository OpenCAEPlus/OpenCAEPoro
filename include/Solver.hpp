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
#include "LinearSolver.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

class Solver
{
public:
    /// Initialize the reservoir.
    void InitReservoir(Reservoir& rs) const;
    /// Start simulation.
    void RunSimulation(Reservoir& rs, OCPControl& ctrl, OCPOutput& output);
    /// Run one time step.
    void GoOneStep(Reservoir& rs, OCPControl& ctrl);
    /// Setup Mehod
    void SetupMethod(const Reservoir& rs, const OCPControl& ctrl);
    /// Allocate memory for linear system.
    void AllocateMat(const Reservoir& rs);
    /// Setup linear solver params.
    void SetupParamLS(const string& dir, const string& file);

private:
    /// Before solve: prepare for assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble and Solve: assemble linear system parts together then solve.
    void AssembleSolve(Reservoir& rs, OCPControl& ctrl, const OCP_DBL& dt);
    /// Update reservoir properties after solving for primary variables.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
    /// Determine if Newton iteration is finished.
    bool FinishNR();
    /// Finish current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

private:
    FluidSolver  FSolver;
    LinearSolver LSolver;
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