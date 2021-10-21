/*! \file    OCPMethod.cpp
 *  \brief   OCPMethod class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMethod.hpp"

void OCP_IMPES::SetupParam(const string& dir, const string& file)
{
    solver.SetupParam(dir, file);
}

void OCP_IMPES::AllocateMat(const Reservoir& rs)
{
    solver.AllocateMem(rs.bulk.GetBulkNum() + rs.wellgroup.GetWellNum());
    rs.conn.AllocateMat(solver);
    rs.wellgroup.AllocateMat(solver);
    solver.AllocateColValMem();
}

void OCP_IMPES::Run(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output)
{
    GetWallTime timer;
    timer.Start();

    USI numdates = ctrl.GetNumDates();
    output.PrintInfoSched(rs, ctrl, timer.Stop());
    for (USI d = 0; d < numdates - 1; d++) {
        rs.wellgroup.ApplyControl(d);
        ctrl.ApplyControl(d);
        ctrl.InitTime(d);
        while (ctrl.criticalTime[d + 1] - ctrl.current_time > TINY) {
            GoOneStep(rs, ctrl);
            output.SetVal(rs, ctrl);
        }
        output.PrintInfoSched(rs, ctrl, timer.Stop());
    }
}

void OCP_IMPES::GoOneStep(Reservoir& rs, OCP_Control& ctrl)
{
    OCP_DBL ve        = 0.01;
    double& dt        = ctrl.current_dt;
    int     flagCheck = 0;

    // cout << setprecision(3) << ctrl.current_time << " Days\n";

    // Init wells
    rs.wellgroup.PrepareWell(rs.bulk);

    OCP_DBL cfl = rs.CalCFL01(dt);
    if (cfl > 1) dt /= (cfl + 1);

    while (true) {
        if (dt < MIN_TIME_STEP) OCP_ABORT("Time stepsize is too small!");

        rs.AssembleMat(solver, dt);
        SolveP(rs, ctrl);

        // first check : Pressure check
        flagCheck = rs.CheckP();
        if (flagCheck == 1) {
            dt /= 2;
            continue;
        } else if (flagCheck == 2) {
            dt /= 2;
            continue;
        }

        rs.conn.CalFlux(rs.bulk);
        rs.wellgroup.CalFlux(rs.bulk);

        // second check : CFL check
        cfl = rs.CalCFL01(dt);
        if (cfl > 1) {
            dt /= 2;
            rs.ResetVal();
            cout << "CFL is too big" << endl;
            continue;
        }

        rs.conn.MassConserve(rs.bulk, dt);
        rs.wellgroup.MassConserve(rs.bulk, dt);

        // third check: Ni check
        if (!rs.CheckNi()) {
            dt /= 2;
            rs.ResetVal01();
            cout << "Negative Ni occurs\n";
            continue;
        }

        rs.bulk.FlashNi();
        rs.bulk.CalVporo();

        // fouth check: Volume error check
        if (!rs.CheckVe(ve)) {
            cout << "###WARNING: volume error is too big\n";
            dt /= 2;
            rs.ResetVal02();
            continue;
        }

        rs.bulk.CalKrPc();
        rs.conn.CalFlux(rs.bulk);

        break;
    }

    rs.wellgroup.CalIPRT(rs.bulk, dt);
    ctrl.tstep += 1;
    ctrl.iterNR = 1;
    ctrl.iterNR_total += 1;

    rs.bulk.CalMaxChange();
    ctrl.CalNextTstep(rs);
    rs.bulk.UpdateLastStep();
    rs.conn.UpdateLastStep();
    rs.wellgroup.UpdateLastStep();

    //cout << ctrl.current_time << "\t";
    //for (USI p = 0; p < 3; p++) {
    //    cout << rs.wellgroup.GetWellDg(8, p) << "\t";
    //}
    //cout << "\n";
}

/// First assemble linear, then solve and return solution
void OCP_IMPES::SolveP(Reservoir& rs, OCP_Control& ctrl)
{

#ifdef DEBUG
    solver.CheckVal();
#endif // DEBUG

#ifdef __SOLVER_FASP__

    solver.AssembleMat_Fasp();
    GetWallTime Timer;
    Timer.Start();
    int status = solver.FaspSolve();
    ctrl.timeLS += Timer.Stop() / 1000;

#ifdef DEBUG
    solver.PrintfMatCSR("testA.out", "testb.out");
    solver.PrintfSolution("testx.out");
#endif // DEBUG

    solver.Free_Fasp();

    ctrl.iterLS = status;
    ctrl.iterLS_total += status;

#endif // __SOLVER_FASP__

    rs.GetSolution_IMPES(solver.GetSolution());
    solver.ClearData();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/