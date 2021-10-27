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

#include "FluidSolver.hpp"
#include "LinearSolver.hpp"

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

class Solver
{
public:
    /// Allocate mat memory.
    void AllocateMat(const Reservoir& rs);
    /// Initialize the reservoir.
    void InitReservoir(Reservoir& rs) const;
    /// Start simulation.
    void RunSimulation(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output);
    /// Run one time step.
    void GoOneStep(Reservoir& rs, OCP_Control& ctrl);

    /// Before solve: prepare for assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Setup linear solver params.
    void SetupParamLS(const string& dir, const string& file);
    /// Assemble and Solve: assemble linear system parts together then solve.
    void AssembleSolve(Reservoir& rs, OCP_Control& ctrl, const OCP_DBL& dt);
    /// Update properties after solving.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
    /// Determine if Newton iteration is finished.
    bool FinishNR();
    /// Finish current time step.
    void FinishStep(Reservoir& rs, OCP_Control& ctrl);

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
/*----------------------------------------------------------------------------*/