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

void Solver::Setup(Reservoir& rs, const OCPControl& ctrl) 
{ 
    OCPModel = ctrl.GetModel();
    switch (OCPModel)
    {
    case ISOTHERMALMODEL:
        SetupIsoT(rs, ctrl);
        break;
    case THERMALMODEL:
        SetupT(rs, ctrl);
        break;
    default:
        OCP_ABORT("WRONG MODEL!");
        break;
    }
}


void Solver::SetupIsoT(Reservoir& rs, const OCPControl& ctrl)
{
    // Setup static infomation for reservoir
    rs.Setup(false);
    IsoTSolver.SetupMethod(rs, ctrl);
}


void Solver::SetupT(Reservoir& rs, const OCPControl& ctrl)
{
    OCP_ABORT("Not Completed Now!");
}

/// Initialize the reservoir setting for different solution methods.
void Solver::InitReservoir(Reservoir& rs) const
{
    switch (OCPModel)
    {
    case ISOTHERMALMODEL:
        InitReservoirIsoT(rs);
        break;
    case THERMALMODEL:
        InitReservoirT(rs);
        break;
    default:
        OCP_ABORT("WRONG MODEL!");
        break;
    }
}


void Solver::InitReservoirIsoT(Reservoir& rs) const
{
    IsoTSolver.InitReservoir(rs);
}


void Solver::InitReservoirT(Reservoir& rs) const
{
    OCP_ABORT("Not Completed Now!");
}


/// Simulation will go through all time steps and call GoOneStep at each step.
void Solver::RunSimulation(Reservoir& rs, OCPControl& ctrl, OCPOutput& output)
{
    GetWallTime timer;
    timer.Start();
    output.PrintInfoSched(rs, ctrl, timer.Stop());
    USI numTSteps = ctrl.GetNumTSteps();
    for (USI d = 0; d < numTSteps - 1; d++) {
        rs.ApplyControl(d);
        ctrl.ApplyControl(d, rs);
        while (!ctrl.IsCriticalTime(d + 1)) {
            GoOneStep(rs, ctrl);
            output.SetVal(rs, ctrl);
            if (ctrl.printLevel >= PRINT_ALL) {
                // Print Summary and critical information at every time step
                output.PrintInfo();
            }
        }
        output.PrintInfoSched(rs, ctrl, timer.Stop());     
        // rs.allWells.ShowWellStatus(rs.bulk);
    }

    if (rs.bulk.GetMixMode() == EOS_PVTW) {
        cout << "SSMSTA:     " << setw(12) << rs.bulk.GetSSMSTAiters() << setw(15)
             << rs.bulk.GetSSMSTAiters() * 1.0 / rs.bulk.GetSSMSTAcounts() << endl;
        cout << "NRSTA:      " << setw(12) << rs.bulk.GetNRSTAiters() << setw(15)
             << rs.bulk.GetNRSTAiters() * 1.0 / rs.bulk.GetNRSTAcounts() << endl;
        cout << "SSMSP:      " << setw(12) << rs.bulk.GetSSMSPiters() << setw(15)
             << rs.bulk.GetSSMSPiters() * 1.0 / rs.bulk.GetSSMSPcounts() << endl;
        cout << "NRSP:       " << setw(12) << rs.bulk.GetNRSPiters() << setw(15)
             << rs.bulk.GetNRSPiters() * 1.0 / rs.bulk.GetNRSPcounts() << endl;
        cout << "NRRR:       " << setw(12) << rs.bulk.GetRRiters() << setw(15)
             << rs.bulk.GetRRiters() * 1.0 / rs.bulk.GetRRcounts() << endl;
    }
    ctrl.RecordTotalTime(timer.Stop() / 1000);
}

/// This is one time step of dynamic simulation in an abstract setting.
void Solver::GoOneStep(Reservoir& rs, OCPControl& ctrl)
{
    
    if (ctrl.printLevel >= PRINT_SOME) {
        cout << "### DEBUG: " << setprecision(3) << fixed << ctrl.GetCurTime()
             << " Days";
        cout << ",  NR: " << ctrl.GetNRiterT() << ",  LS: " << ctrl.GetLSiterT()
             << ",  Last dt: " << ctrl.last_dt << " Days" << endl;
    }

    switch (OCPModel)
    {
    case ISOTHERMALMODEL:
        GoOneStepIsoT(rs, ctrl);
        break;
    case THERMALMODEL:
        GoOneStepT(rs, ctrl);
        break;
    default:
        OCP_ABORT("WRONG MODEL!");
        break;
    }
}


void Solver::GoOneStepIsoT(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.GetCurDt();

    // Prepare for time marching
    IsoTSolver.Prepare(rs, dt);

    // Time marching with adaptive time stepsize
    while (OCP_TRUE) {
        if (dt < MIN_TIME_CURSTEP) OCP_ABORT("Time stepsize is too small!");
        // Assemble linear system
        IsoTSolver.AssembleMat(rs, ctrl);
        // Solve linear system
        IsoTSolver.SolveLinearSystem(rs, ctrl);
        if (!IsoTSolver.UpdateProperty(rs, ctrl)) {
            ctrl.ResetIterNRLS();
            continue;
        }
        if (IsoTSolver.FinishNR(rs, ctrl)) break;
    }

    // Finish current time step
    IsoTSolver.FinishStep(rs, ctrl);
}


void Solver::GoOneStepT(Reservoir& rs, OCPControl& ctrl)
{
    OCP_ABORT("Not Completed Now!");
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/