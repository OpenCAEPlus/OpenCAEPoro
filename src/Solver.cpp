/*! \file    Solver.cpp
 *  \brief   Solver class definition
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

void Solver::Setup(Reservoir& rs, const OCPControl& ctrl) { SetupMethod(rs, ctrl); }

/// Initialize the reservoir setting for different solution methods.
void Solver::InitReservoir(Reservoir& rs) const
{
    // Initialize the fluid part
    FSolver.InitReservoir(rs);
}

/// Simulation will go through all time steps and call GoOneStep at each step.
void Solver::RunSimulation(Reservoir& rs, OCPControl& ctrl, OCPOutput& output)
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

    if (rs.bulk.GetMixMode() == EOS_PVTW) {
        cout << "SSMSTA:     " << rs.bulk.GetSSMSTAiters() << endl;
        cout << "NRSTA:      " << rs.bulk.GetNRSTAiters() << endl;
        cout << "SSMSP:      " << rs.bulk.GetSSMSPiters() << endl;
        cout << "NRSP:       " << rs.bulk.GetNRSPiters() << endl;
    }  
    ctrl.RecordTotalTime(timer.Stop() / 1000);
}

/// This is one time step of dynamic simulation in an abstract setting.
void Solver::GoOneStep(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.GetCurDt();

//#ifdef _DEBUG
    //cout << "### DEBUG: " << fixed << ctrl.GetCurTime() << " Days";
    //cout << "  NR: " << ctrl.GetNRiterT() << "  LS: " << ctrl.GetLSiterT() << endl;
//#endif // DEBUG

    // Prepare for time marching
    Prepare(rs, dt);

    // Time marching with adaptive time stepsize
    while (true) {
        if (dt < MIN_TIME_CURSTEP) OCP_ABORT("Time stepsize is too small!");
        AssembleSolve(rs, ctrl);
        if (!UpdateProperty(rs, ctrl)) {
            ctrl.ResetIterNRLS();
            continue;
        }
        if (FinishNR(rs, ctrl)) break;
    }

    // Finish current time step
    FinishStep(rs, ctrl);
}

/// Get ready for assembling the linear system of this time step.
void Solver::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    // Prepare for the fluid part
    FSolver.Prepare(rs, dt);
}

void Solver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    FSolver.SetupMethod(rs, ctrl);
}

/// Assemble linear system and then solve it.
void Solver::AssembleSolve(Reservoir& rs, OCPControl& ctrl)
{
    // Assemble linear system
    FSolver.AssembleMat(rs, ctrl.current_dt);
    // Solve linear system
    FSolver.SolveLinearSystem(rs, ctrl);
}

/// Update properties after solving.
bool Solver::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    // Update for the fluid part
    return FSolver.UpdateProperty(rs, ctrl);
}

/// Clean up Newton-Raphson iteration if there is any.
bool Solver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    // Clean up the fluid part
    return FSolver.FinishNR(rs, ctrl);
}

/// Clean up time step.
void Solver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    // Clean up the fluid part
    FSolver.FinishStep(rs, ctrl);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/