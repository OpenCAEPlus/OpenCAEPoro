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

/// Call setup procedures for reservoir, output, and linear solver.
void OpenCAEPoro::SetupSimulator(ParamRead&  param,
                                 const USI&  argc,
                                 const char* options[])
{
    GetWallTime timer;
    timer.Start();

    // Read parameters from input file
    InputParam(param);
    // Read Fast control
    control.SetupFastControl(argc, options);
    // Setup static infomation for reservoir
    reservoir.Setup(output.IfOutputVTK());
    // Setup output for dynamic simulation
    output.Setup(reservoir, control);
    // Setup static information for solver
    solver.Setup(reservoir, control);

    double finalTime = timer.Stop() / 1000;
    if (control.printLevel >= PRINT_MIN) {
        cout << endl
             << "Setup simulation done. Wall time : " << fixed << setprecision(3)
             << finalTime << " Sec" << endl
             << endl;
    }
    control.RecordTotalTime(finalTime);
}

/// Initialize the reservoir class.
void OpenCAEPoro::InitReservoir()
{
    GetWallTime timer;
    timer.Start();

    solver.InitReservoir(reservoir);

    double finalTime = timer.Stop() / 1000;
    if (control.printLevel >= PRINT_MIN) {
        cout << endl
             << "Initialization done. Wall time : " << fixed << setprecision(3)
             << finalTime << " Sesc" << endl;
    }
    control.RecordTotalTime(finalTime);
}

/// Call IMPEC, FIM, AIM, etc for dynamic simulation.
void OpenCAEPoro::RunSimulation()
{
    switch (control.GetMethod()) {
        case IMPEC:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with IMPEC\n" << endl;
            }
            break;
        case FIM:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with FIM\n" << endl;
            }
            break;
        case FIMn:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with FIMn\n" << endl;
            }
            break;
        case AIMc:
            if (control.printLevel >= PRINT_MIN) {
                cout << "\nDynamic simulation with AIMc\n" << endl;
            }
            break;
        default:
            OCP_ABORT("Wrong method type is used!");
    }

    solver.RunSimulation(reservoir, control, output);
}

/// Print summary information on screen and SUMMARY.out file.
void OpenCAEPoro::OutputResults() const
{
    // find an appropriate size for printing times
    int fixWidth =
        MAX(log10(control.current_time), log10(MAX(control.totalSimTime, 1.0))) + 6;

    cout << "==================================================" << endl;

    // print numbers of steps
    cout << "Final time:                  " << right << fixed << setprecision(3)
         << setw(fixWidth) << control.current_time << " (Days)" << endl;
    cout << " - Avg time step size ......." << setw(fixWidth)
         << control.current_time / control.numTstep << " (" << control.numTstep
         << " steps)" << endl;
    cout << " - Avg Newton steps ........." << setw(fixWidth)
         << static_cast<double>(control.iterNR_total) / control.numTstep << " ("
         << control.iterNR_total << " succeeded + " << control.wastedIterNR
         << " wasted)" << endl;
    cout << " - Avg linear steps ........." << setw(fixWidth)
         << static_cast<double>(control.iterLS_total) / control.iterNR_total << " ("
         << control.iterLS_total << " succeeded + " << control.wastedIterLS
         << " wasted)" << endl;

    // print time usages
    cout << "Simulation time:             " << setw(fixWidth) << control.totalSimTime
         << " (Seconds)" << endl;
    cout << " - % Assembling ............." << setw(fixWidth)
         << 100.0 * control.totalAssembleMatTime / control.totalSimTime << " ("
         << control.totalAssembleMatTime << "s)" << endl;
    cout << " - % Linear solver .........." << setw(fixWidth)
         << 100.0 * control.totalLStime / control.totalSimTime << " ("
         << control.totalLStime << "s)" << endl;
    cout << " - % Updating properties ...." << setw(fixWidth)
         << 100.0 * control.totalUpdatePropertyTime / control.totalSimTime << " ("
         << control.totalUpdatePropertyTime << "s)" << endl;
    cout << " - % Scheduled output ......." << setw(fixWidth)
         << 100.0 * output.outputTime / control.totalSimTime << " ("
         << output.outputTime << "s)" << endl;

    cout << "==================================================" << endl;

    output.PrintInfo();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Dec/05/2021      Format file                          */
/*----------------------------------------------------------------------------*/