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

void Solver::Setup(Reservoir &rs, const OCPControl &ctrl) { SetupMethod(rs, ctrl); }

/// Initialize the reservoir setting for different solution methods.
void Solver::InitReservoir(Reservoir &rs) const
{
    // Initialize the fluid part
    IsoTSolver.InitReservoir(rs);
}

/// Simulation will go through all time steps and call GoOneStep at each step.
void Solver::RunSimulation(Reservoir &rs, OCPControl &ctrl, OCPOutput &output)
{
    GetWallTime timer;
    timer.Start();
    output.PrintInfoSched(rs, ctrl, timer.Stop());
    USI numTSteps = ctrl.GetNumTSteps();
    for (USI d = 0; d < numTSteps - 1; d++)
    {
        rs.ApplyControl(d);
        ctrl.ApplyControl(d, rs);
        while (!ctrl.IsCriticalTime(d + 1))
        {
            GoOneStep(rs, ctrl);
            output.SetVal(rs, ctrl);
        }
        output.PrintInfoSched(rs, ctrl, timer.Stop());
        if (ctrl.printLevel > 2) {
            // Print Summary and critical information at every TSTEP
            output.PrintInfo();
        }       
        // rs.allWells.ShowWellStatus(rs.bulk);
    }

    if (rs.bulk.GetMixMode() == EOS_PVTW)
    {
        cout << "SSMSTA:     " << setw(12) << rs.bulk.GetSSMSTAiters()
            << setw(15) << rs.bulk.GetSSMSTAiters() * 1.0 / rs.bulk.GetSSMSTAcounts() << endl;
        cout << "NRSTA:      " << setw(12) << rs.bulk.GetNRSTAiters() 
            << setw(15) << rs.bulk.GetNRSTAiters() * 1.0 / rs.bulk.GetNRSTAcounts() << endl;
        cout << "SSMSP:      " << setw(12) << rs.bulk.GetSSMSPiters() 
            << setw(15) << rs.bulk.GetSSMSPiters() * 1.0 / rs.bulk.GetSSMSPcounts() << endl;
        cout << "NRSP:       " << setw(12) << rs.bulk.GetNRSPiters() 
            << setw(15) << rs.bulk.GetNRSPiters() * 1.0 / rs.bulk.GetNRSPcounts() << endl;
        cout << "NRRR:       " << setw(12) << rs.bulk.GetRRiters() 
            << setw(15) << rs.bulk.GetRRiters() * 1.0 / rs.bulk.GetRRcounts() << endl;
    }
    ctrl.RecordTotalTime(timer.Stop() / 1000);
}

/// This is one time step of dynamic simulation in an abstract setting.
void Solver::GoOneStep(Reservoir &rs, OCPControl &ctrl)
{
    OCP_DBL &dt = ctrl.GetCurDt();

    if (ctrl.printLevel > 0) {
        cout << "### DEBUG: " << setprecision(3) << fixed << ctrl.GetCurTime() << " Days";
        cout << "  NR: " << ctrl.GetNRiterT() << "  LS: " << ctrl.GetLSiterT() << "     ";
        cout << "Last dt  " << ctrl.last_dt << " Days" << endl;
    }
    
    // Prepare for time marching
    Prepare(rs, dt);

    // Time marching with adaptive time stepsize
    while (true)
    {
        if (dt < MIN_TIME_CURSTEP)
            OCP_ABORT("Time stepsize is too small!");
        AssembleSolve(rs, ctrl);
        if (!UpdateProperty(rs, ctrl))
        {
            ctrl.ResetIterNRLS();
            continue;
        }
        if (FinishNR(rs, ctrl))
            break;
    }

    // Finish current time step
    FinishStep(rs, ctrl);
}

/// Get ready for assembling the linear system of this time step.
void Solver::Prepare(Reservoir &rs, OCP_DBL &dt)
{
    // Prepare for the fluid part
    IsoTSolver.Prepare(rs, dt);
}

void Solver::SetupMethod(Reservoir &rs, const OCPControl &ctrl)
{
    IsoTSolver.SetupMethod(rs, ctrl);
}

/// Assemble linear system and then solve it.
void Solver::AssembleSolve(Reservoir &rs, OCPControl &ctrl)
{
    // Assemble linear system
    IsoTSolver.AssembleMat(rs, ctrl.current_dt);
    // Solve linear system
    IsoTSolver.SolveLinearSystem(rs, ctrl);
}

/// Update properties after solving.
bool Solver::UpdateProperty(Reservoir &rs, OCPControl &ctrl)
{
    // Update for the fluid part
    return IsoTSolver.UpdateProperty(rs, ctrl);
}

/// Clean up Newton-Raphson iteration if there is any.
bool Solver::FinishNR(Reservoir &rs, OCPControl &ctrl)
{
    // Clean up the fluid part
    return IsoTSolver.FinishNR(rs, ctrl);
}

/// Clean up time step.
void Solver::FinishStep(Reservoir &rs, OCPControl &ctrl)
{
    // Clean up the fluid part
    IsoTSolver.FinishStep(rs, ctrl);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/