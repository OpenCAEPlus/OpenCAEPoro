/*! \file    OCP.cpp
 *  \brief   OpenCAEPoro class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCP.hpp"

/// Read from input file and set control and output params.
void OpenCAEPoro::InputParam(ParamRead& param)
{
    reservoir.InputParam(param);
    control.InputParam(param.paramControl);
    output.InputParam(param.paramOutput);
}

/// Call setup processdures for reservoir, output, and linear solver.
void OpenCAEPoro::SetupReservoir(ParamRead& param)
{
    InputParam(param);
    reservoir.Setup();
    output.Setup(reservoir, control);
    SetupSolver();
}

/// Call SetupParm and AllocateMat to prepare the linear solver
void OpenCAEPoro::SetupSolver()
{
    switch (control.method) {
        case FIM: // Fully Implicite Method
            cout << "Applying the FIM Method ..." << endl;
            break;
        case IMPES: // IMplicit Pressure Explicit Composition
            cout << "Applying the IMPES Method ..." << endl;
            impes.SetupParam(control.workDir, control.solveFile);
            impes.AllocateMat(reservoir);
            break;
        default:
            OCP_MESSAGE("Trying to call " << control.method);
            OCP_ABORT("Solution method not supported!");
    }
}

/// Initialize the reservoir class.
void OpenCAEPoro::InitReservoir() { reservoir.Init(); }

/// Call IMPES, FIM, etc for dynamic simulation.
void OpenCAEPoro::RunSimulation()
{
    GetWallTime Timer;
    Timer.Start();

    switch (control.method) {
        case IMPES:
            impes.Run(reservoir, control, output);
            break;
        default:
            OCP_MESSAGE("Trying to call " << control.method);
            OCP_ABORT("Solution method not supported!");
    }

    control.totalTime = Timer.Stop() / 1000;
}

/// Print summary information to cout and SUMMARY.out file.
void OpenCAEPoro::OutputResults()
{
    cout << "=========================================" << endl;
    cout << "Final time:          " << control.current_time << " Days" << endl;
    cout << "Total time steps:    " << control.iterNR_total << endl;
    cout << "Simulation time:     " << control.totalTime << "s" << endl;
    cout << "Total linear steps:  " << control.iterLS_total << endl;
    cout << "Linear solve time:   " << control.timeLS << "s" << endl;

    output.PrintInfo();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/