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
#include "Solver.hpp"

/// Simulation will go through all time steps and call GoOneStep at each step.
void Solver::RunSimulation(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output)
{
    GetWallTime timer;
    timer.Start();

    USI numTSteps = ctrl.GetNumTSteps();
    output.PrintInfoSched(rs, ctrl, timer.Stop());

    for (USI d = 0; d < numTSteps - 1; d++) {
        rs.ApplyControl(d);
        ctrl.ApplyControl(d);
        ctrl.InitTime(d);
        while (!ctrl.IsCriticalTime(d + 1)) {
            GoOneStep(rs, ctrl);
            output.SetVal(rs, ctrl);
        }
        output.PrintInfoSched(rs, ctrl, timer.Stop());
    }

    ctrl.RecordTotalTime(timer.Stop() / 1000);
}

/// This is one time step of dynamic simulation in an abstract setting.
void Solver::GoOneStep(Reservoir& rs, OCP_Control& ctrl)
{
    OCP_DBL& dt = ctrl.GetCurDt();

    // Prepare for time marching
    Prepare(rs, dt);

    // Time marching with adaptive time stepsize
    while (true) {
        if (dt < MIN_TIME_CURSTEP) OCP_ABORT("Time stepsize is too small!");
        AssembleSolve(rs, ctrl, dt);
        if (!UpdateProperty(rs, dt)) continue;
        if (FinishNR()) break;
    }

    // Finish current time step
    FinishStep(rs, ctrl);
}

void Solver::Prepare(Reservoir& rs, OCP_DBL& dt) { FSolver.Prepare(rs, dt); }

void Solver::AssembleSolve(Reservoir& rs, OCP_Control& ctrl, const OCP_DBL& dt)
{
    // Assemble linear system
    FSolver.AssembleMat(rs, dt);
    // Solve linear system
    FSolver.SolveLinearSystem(rs, ctrl);
}

bool Solver::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{
    return FSolver.UpdateProperty(rs, dt);
}

bool Solver::FinishNR() { return FSolver.FinishNR(); }

void Solver::FinishStep(Reservoir& rs, OCP_Control& ctrl)
{
    FSolver.FinishStep(rs, ctrl);
}

void Solver::SetupParamLS(const string& dir, const string& file)
{
    FSolver.SetupParamLS(dir, file);
}

void Solver::AllocateMat(const Reservoir& rs) { FSolver.AllocateMat(rs); }

void Solver::InitReservoir(Reservoir& rs) const { FSolver.InitReservoir(rs); }

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/