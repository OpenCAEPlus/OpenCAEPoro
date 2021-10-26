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
void OpenCAEPoro::SetupSimulator(ParamRead& param)
{
    InputParam(param);
    reservoir.Setup();
    output.Setup(reservoir, control);
    SetupLinearSolver();
}

/// Call SetupParm and AllocateMat to prepare the linear solver
void OpenCAEPoro::SetupLinearSolver()
{
    solver.SetupParamLS(control.workDir, control.solveFile);
    solver.AllocateMat(reservoir);
}

/// Initialize the reservoir class.
void OpenCAEPoro::InitReservoir()
{ 
    solver.InitReservoir(reservoir);
}

/// Call IMPEC, FIM, etc for dynamic simulation.
void OpenCAEPoro::RunSimulation()
{
    solver.RunSimulation(reservoir, control, output);
}

/// Print summary information to cout and SUMMARY.out file.
void OpenCAEPoro::OutputResults() const
{
    cout << "=========================================" << endl;
    cout << "Final time:          " << control.current_time << " Days" << endl;
    cout << "Total time steps:    " << control.iterNR_total << endl;
    cout << "Simulation time:     " << control.totalTime << "s" << endl;
    cout << "Total linear steps:  " << control.iterLS_total << endl;
    cout << "Linear solve time:   " << control.timeLS << "s"
         << " (" << 100.0 * control.timeLS / control.totalTime << "%)" << endl;

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