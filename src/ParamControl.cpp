/*! \file    ParamControl.cpp
 *  \brief   ParamControl class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamControl.hpp"

/// Initialize control parameters with default values.
void ParamControl::Init(string& indir)
{
    dir = indir;
    InitMethod();
    InitTime();
    InitTuning();
}

/// Initialize with default solution method and linear solver.
void ParamControl::InitMethod()
{
    method      = "IMPEC";
    linearSolve = "./csr.fasp";
}

/// Initialize TUNING parameters.
void ParamControl::InitTuning()
{
    tuning.resize(3);

    // Timestepping controls, * means this param is available
    // Limits: timestep and change factor.
    tuning[0].resize(10);
    tuning[0][0] = 1.0;   //* Maximum initial time stepsize of next timestep
    tuning[0][1] = 365.0; //* Maximum length of timesteps after the next
    tuning[0][2] = 0.1;   //* Minimum length of all timesteps
    tuning[0][3] = 3.0;   //* Maximum time stepsize increase factor
    tuning[0][4] = 0.15;  //* Minimum choppable timestep
    tuning[0][5] = 0.3;   //* Factor by which timestep is cut after convergence failure
    tuning[0][6] = 0.1;   // ???
    tuning[0][7] = 1.25;  // Maximum increase factor after a convergence failure
    tuning[0][8] = (method == "IMPEC") ? 0.2 : 1E20; // ???
    tuning[0][9] = -1; // Maximum next time stepsize following a well modification

    // Timestepping controls, * means this param is available
    // Prediction: an ideal maximum change of variables at next time step.
    // So they're used to calculate change factor of time step by predicting linearly.
    tuning[1].resize(13);
    //* dPlim: ideal maximum Pressure change at next time step.
    tuning[1][0] = (method == "IMPEC") ? 200.0 : 300.0;
    //* dSlim: ideal maximum Saturation change at next time step.
    tuning[1][1] = (method == "IMPEC") ? 0.2 : 0.2;
    //* dNlim: ideal maximum relative Ni(moles of components) change at next time step.
    tuning[1][2] = 0.3; // Target material balance error
    //* dVerrlim: ideal maximum relative Verr(error between fluid and pore) change at
    // next time step.
    tuning[1][3] = (method == "IMPEC") ? 1E-3 : 1E-3;

    tuning[1][4] = 10.0; // Maximum time truncation error
    // Maximum non-linear convergence error
    tuning[1][5] = (method == "IMPEC") ? 0.75 : 0.01;
    tuning[1][6] = 1E-6; // Maximum material balance error
    // Maximum linear convergence error
    tuning[1][7]  = (method == "IMPEC") ? 1E-4 : 1E-3;
    tuning[1][8]  = 1E-3;  // Maximum well flow rate convergence error
    tuning[1][9]  = 0.025; // Target Fluid-in-place error for LGR runs
    tuning[1][10] = -1;    // Target surfactant change (Surfactant Model only)
    tuning[1][11] = 0.01;  // Threshold for damping in ion exchange calc. (Multi-Comp.
                           // Brine Model only)
    tuning[1][12] = 1;     // Weighting factor for active tracer updates

    // Nonlinear Solver controls, * means this param is available
    tuning[2].resize(10);
    //* Maximum number of Newton iterations in a timestep
    tuning[2][0] = (method == "IMPEC") ? 1 : 10;
    //* Maximum non-linear convergence error
    tuning[2][1] = 1e-3;
    //* Maximum Pressure change in a Newton iteration
    tuning[2][2] = 200.0;
    //* Maximum Saturation change in a Newton iteration
    tuning[2][3] = 0.2;
    //* Minimum Pressure change in a Newton iteration
    tuning[2][4] = 1;
    //* Minimum Saturation change in a Newton iteration
    tuning[2][5] = 0.01;
    //* Maximum Verr(error between fluid and pore) change in a Newton iteration
    tuning[2][6] = 0.01;

    tuning[2][7] = 1E6; // Maximum saturation change at last Newton iteration
    // Target maximum pressure change in a timestep
    tuning[2][8] = (method == "IMPEC") ? 100 : 1E6;
    tuning[2][9] = -1; // Maximum tolerable pressure change in a timestep
}

/// Initialize solution method.
void ParamControl::InputMETHOD(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    if (vbuf[0] == "FIM") {
        method      = "FIM";
        linearSolve = "./bsr.fasp";
    }

    if (vbuf.size() > 1) linearSolve = vbuf[1];

    cout << "\n---------------------" << endl
         << "METHOD"
         << "\n---------------------" << endl;
    cout << "   " << method << "  " << linearSolve << endl;
}

/// Read TUNING parameters.
void ParamControl::InputTUNING(ifstream& ifs)
{
    assert(criticalTime.size() >= 1);

    TUNING         tmp(tuning);
    USI            d   = criticalTime.size() - 1;
    USI            row = 0;
    vector<string> vbuf;

    while (ReadLine(ifs, vbuf)) {
        DealDefault(vbuf);
        OCP_INT len = vbuf.size();

        for (OCP_INT i = 0; i < len - 1; i++) {
            tmp[row][i] = stod(vbuf[i]);
        }
        if (vbuf[len - 1] != "/") {
            tmp[row][len - 1] = stod(vbuf[len - 1]);
        } else {
            row++;
        }
        if (row == 3) break;
    }

    tuning_T.push_back(TuningPair(d, tmp));
    DisplayTuning();
}

/// Print TUNING parameters.
void ParamControl::DisplayTuning() const
{
    cout << "\n---------------------" << endl
         << "TUNING"
         << "\n---------------------" << endl;
    for (auto v : tuning_T) {
        cout << v.d << endl;
        for (auto v1 : v.Tuning) {
            for (auto v2 : v1) {
                cout << v2 << "   ";
            }
            cout << "/ " << endl;
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/