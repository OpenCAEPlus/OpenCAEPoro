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

void OpenCAEPoro::InputParam(ParamRead& param)
{
    reservoir.InputParam(param);
    control.InputParam(param.param_Control);
    output.InputParam(param.param_Output);
}

void OpenCAEPoro::SetupReservoir(ParamRead& param)
{
    InputParam(param);
    reservoir.Setup();
    output.Setup(reservoir, control);
    SetupSolver();
}

void OpenCAEPoro::SetupSolver()
{
    if (control.method == IMPES) {
        impes.SetupParam(control.workDir, control.solveFile);
        impes.AllocateMat(reservoir);
        cout << "IMPES Method Applys !" << endl;
    } else if (control.method == FIM) {
        cout << "FIM Method Applys !" << endl;
    }
}

void OpenCAEPoro::InitReservoir() { reservoir.Init(); }

void OpenCAEPoro::RunSimulation()
{
    GetWallTime Timer;
    Timer.Start();

    switch (control.method) {
        case IMPES:
            impes.Run(reservoir, control, output);
        default:
            break;
    }

    control.totalTime = Timer.Stop() / 1000;
}

void OpenCAEPoro::OutputResults()
{

    cout << endl;
    cout << "Final time:          " << control.current_time << " Days" << endl;
    cout << "Total linear steps:  " << control.iterLS_total << endl;
    cout << "Linear solve time:   " << control.timeLS << "s" << endl;
    cout << "Total time steps:    " << control.iterNR_total << endl;
    cout << "Simulation time:     " << control.totalTime << "s" << endl;
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