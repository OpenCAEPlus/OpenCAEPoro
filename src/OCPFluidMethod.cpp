/*! \file    OCPFluidMethod.cpp
 *  \brief   Definition of solution methods for fluid part in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Nov/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFluidMethod.hpp"

////////////////////////////////////////////
// OCP_IMPEC
////////////////////////////////////////////

void OCP_IMPEC::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for IMPEC
    rs.AllocateIMPEC_IsoT();
    // Allocate Memory of Matrix for IMPEC
    rs.AllocateMatIMPEC(myLS);

    myLS.SetupLinearSolver(SCALARFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

/// Initialize reservoir
void OCP_IMPEC::InitReservoir(Reservoir& rs) const { rs.InitIMPEC(); }

void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    OCP_DBL cfl = rs.CalCFL(dt);
    if (cfl > 1) dt /= (cfl + 1);
}


void OCP_IMPEC::AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const
{
    rs.AssembleMatIMPEC(myLS, dt);
}


void OCP_IMPEC::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    myLS.CheckEquation();
#endif // DEBUG

    myLS.AssembleMatLinearSolver();

#ifdef DEBUG
    // myLS.OutputLinearSystem("testA_IMPEC.out", "testb_IMPEC.out");
#endif // DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = myLS.Solve();
    if (status < 0) {
        status = myLS.GetNumIters();
    }

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

#ifdef DEBUG
    // myLS.OutputSolution("testx.out");
#endif // DEBUG

    rs.GetSolutionIMPEC(myLS.GetSolution());
    myLS.ClearData();
}

OCP_BOOL OCP_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP(OCP_TRUE, OCP_TRUE);
    switch (flagCheck) {
        case 1:
            // Negative Bulk P, well P, or Perforation P
            dt /= 2;
            return OCP_FALSE;
        case 2:
            // Switch Well opt Mode, or close the cross-flow perforation
            dt /= 1;
            return OCP_FALSE;
        default:
            // All right
            break;
    }

    // Calculate Flux between bulks and between bulks and wells
    rs.CalFLuxIMPEC();

    // second check : CFL check
    OCP_DBL cfl = rs.CalCFL(dt);
    if (cfl > 1) {
        dt /= 2;
        rs.ResetVal01IMPEC();
        cout << "CFL is too big" << endl;
        return OCP_FALSE;
    }

    rs.MassConserveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal02IMPEC();
        cout << "Negative Ni occurs\n";
        return OCP_FALSE;
    }

    rs.CalRock();
    rs.CalFlashIMPEC();

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        // cout << ctrl.GetCurTime() << "Days" << "=======" << endl;
        dt /= 2;
        rs.ResetVal03IMPEC();
        return OCP_FALSE;
    }

    rs.CalKrPc();
    rs.CalConnFluxIMPEC();

    return OCP_TRUE;
}

OCP_BOOL OCP_IMPEC::FinishNR(const Reservoir& rs)
{

    return OCP_TRUE;
}


void OCP_IMPEC::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStepIMPEC();
    ctrl.CalNextTstepIMPEC(rs);
    ctrl.UpdateIters();
}

////////////////////////////////////////////
// OCP_FIM
////////////////////////////////////////////

void OCP_FIM::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateFIM_IsoT();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM_IsoT(myLS);
    // Allocate memory for residual of FIM
    resFIM.Setup_IsoT(rs.GetBulkNum(), rs.GetWellNum(), rs.GetComNum());
    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void OCP_FIM::InitReservoir(Reservoir& rs) const { rs.InitFIM(); }

void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResFIM(resFIM, dt, OCP_TRUE);
}

void OCP_FIM::AssembleMat(LinearSystem&    myLS,
                          const Reservoir& rs,
                          const OCP_DBL&   dt) const
{
    rs.AssembleMatFIM(myLS, dt);
    myLS.AssembleRhsCopy(resFIM.res);
}

void OCP_FIM::SolveLinearSystem(LinearSystem& myLS,
                                Reservoir&    rs,
                                OCPControl&   ctrl) const
{
#ifdef DEBUG
    myLS.CheckEquation();
#endif // DEBUG

    myLS.AssembleMatLinearSolver();

    GetWallTime Timer;
    Timer.Start();
    int status = myLS.Solve();
    if (status < 0) {
        status = myLS.GetNumIters();
    }
    // cout << "LS step = " << status << endl;

#ifdef DEBUG
    myLS.OutputLinearSystem("testA_FIM_new.out", "testb_FIM_new.out");
    myLS.OutputSolution("testx_FIM_new.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionFIM(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.GetSolution01FIM(myLS.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    myLS.ClearData();
}

OCP_BOOL OCP_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(OCP_TRUE, OCP_FALSE) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM();
        rs.CalResFIM(resFIM, dt, OCP_TRUE);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    rs.CalFlashFIM();
    rs.CalKrPcFIM();
    rs.CalRock();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt, OCP_FALSE);
    return OCP_TRUE;
}

OCP_BOOL OCP_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();

    if (ctrl.printLevel >= PRINT_MOST) {

        if (OCP_TRUE) {
            vector<OCP_INT> totalPhaseNum(3, 0);
            vector<OCP_INT> ltotalPhaseNum(3, 0);
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (rs.bulk.NRphaseNum[n] == 1) ltotalPhaseNum[0]++;
                if (rs.bulk.NRphaseNum[n] == 2) ltotalPhaseNum[1]++;
                if (rs.bulk.NRphaseNum[n] == 3) ltotalPhaseNum[2]++;
                if (rs.bulk.phaseNum[n] == 1) totalPhaseNum[0]++;
                if (rs.bulk.phaseNum[n] == 2) totalPhaseNum[1]++;
                if (rs.bulk.phaseNum[n] == 3) totalPhaseNum[2]++;
            }
            cout << to_string(totalPhaseNum[0]) + " (" +
                        to_string((totalPhaseNum[0] - ltotalPhaseNum[0]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      "
                 << to_string(totalPhaseNum[1]) + " (" +
                        to_string((totalPhaseNum[1] - ltotalPhaseNum[1]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      "
                 << to_string(totalPhaseNum[2]) + " (" +
                        to_string((totalPhaseNum[2] - ltotalPhaseNum[2]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      ";
        }

        if (OCP_TRUE) {
            OCP_USI eVnum = 0;
            OCP_USI eNnum = 0;
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (resFIM.resRelV[n] < ctrl.ctrlNR.NRtol) eVnum++;
                if (resFIM.resRelN[n] < ctrl.ctrlNR.NRtol) eNnum++;
            }
            cout << "eVnum : " << setprecision(ceil(log(rs.GetBulkNum()) / log(10)))
                 << fixed << setw(10) << (1 - eVnum * 1.0 / rs.GetBulkNum()) * 100.0
                 << "%   eNnum : " << setw(10)
                 << (1 - eNnum * 1.0 / rs.GetBulkNum()) * 100.0 << "%" << endl;
        }

        const USI sp = rs.grid.GetNumDigitIJK();
        USI       tmpI, tmpJ, tmpK;
        string    tmps;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, rs.bulk.index_maxNRdSSP);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "### NR : " + to_string(ctrl.iterNR) + "    Res:    " << setprecision(2)
             << scientific << resFIM.maxRelRes0_v << setw(12) << resFIM.maxRelRes_v
             << setw(12) << resFIM.maxRelRes_mol << setw(11) << resFIM.maxWellRelRes_mol
             << setw(20) << NRdPmax << setw(11) << NRdNmax << setw(11) << NRdSmax
             << setw(40) << sqrt(rs.bulk.NRdSSP) / rs.bulk.numBulk << setw(12)
             << rs.bulk.maxNRdSSP << "  " << rs.bulk.phaseNum[rs.bulk.index_maxNRdSSP]
             << "  " << rs.bulk.NRphaseNum[rs.bulk.index_maxNRdSSP] << "  " << tmps
             << "   " << rs.bulk.ePEC[rs.bulk.index_maxNRdSSP] << endl
             << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, dSn);
        tmps = GetIJKformat(tmpI, tmpJ, tmpK, sp);
        cout << "S: " << tmps << setw(12) << resFIM.resRelV[dSn] << setw(12)
             << resFIM.resRelN[dSn] << setw(12) << rs.bulk.ePEC[dSn] << setw(5)
             << rs.bulk.phaseNum[dSn] << setw(5) << rs.bulk.NRphaseNum[dSn] << setw(12)
             << rs.bulk.dPNR[dSn] << setw(8) << "  |" << setw(12) << rs.bulk.S[dSn * 3]
             << setw(12) << rs.bulk.dSNR[dSn * 3] << setw(12) << rs.bulk.dSNRP[dSn * 3]
             << setw(12) << rs.bulk.S[dSn * 3 + 1] << setw(12)
             << rs.bulk.dSNR[dSn * 3 + 1] << setw(12) << rs.bulk.dSNRP[dSn * 3 + 1]
             << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIM.maxId_v);
        tmps = GetIJKformat(tmpI, tmpJ, tmpK, sp);
        cout << "V: " << tmps << setw(12) << resFIM.resRelV[resFIM.maxId_v] << setw(12)
             << resFIM.resRelN[resFIM.maxId_v] << setw(12) << rs.bulk.ePEC[resFIM.maxId_v]
             << setw(5) << rs.bulk.phaseNum[resFIM.maxId_v] << setw(5)
             << rs.bulk.NRphaseNum[resFIM.maxId_v] << setw(12)
             << rs.bulk.dPNR[resFIM.maxId_v] << setw(8) << "  |" << setw(12)
             << rs.bulk.S[resFIM.maxId_v * 3] << setw(12)
             << rs.bulk.dSNR[resFIM.maxId_v * 3] << setw(12)
             << rs.bulk.dSNRP[resFIM.maxId_v * 3] << setw(12)
             << rs.bulk.S[resFIM.maxId_v * 3 + 1] << setw(12)
             << rs.bulk.dSNR[resFIM.maxId_v * 3 + 1] << setw(12)
             << rs.bulk.dSNRP[resFIM.maxId_v * 3 + 1] << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIM.maxId_mol);
        tmps = GetIJKformat(tmpI, tmpJ, tmpK, sp);
        cout << "N: " << tmps << setw(12) << resFIM.resRelV[resFIM.maxId_mol] << setw(12)
             << resFIM.resRelN[resFIM.maxId_mol] << setw(12)
             << rs.bulk.ePEC[resFIM.maxId_mol] << setw(5)
             << rs.bulk.phaseNum[resFIM.maxId_mol] << setw(5)
             << rs.bulk.NRphaseNum[resFIM.maxId_mol] << setw(12)
             << rs.bulk.dPNR[resFIM.maxId_mol] << setw(8) << "  |" << setw(12)
             << rs.bulk.S[resFIM.maxId_mol * 3] << setw(12)
             << rs.bulk.dSNR[resFIM.maxId_mol * 3] << setw(12)
             << rs.bulk.dSNRP[resFIM.maxId_mol * 3] << setw(12)
             << rs.bulk.S[resFIM.maxId_mol * 3 + 1] << setw(12)
             << rs.bulk.dSNR[resFIM.maxId_mol * 3 + 1] << setw(12)
             << rs.bulk.dSNRP[resFIM.maxId_mol * 3 + 1] << endl;
        // cout << "Res2   " << Dnorm2(resFIM.res.size(), &resFIM.res[0]) /
        // resFIM.res.size() << "   "; cout << "SdP " << Dnorm1(rs.bulk.dPNR.size(),
        // &rs.bulk.dPNR[0]) / rs.bulk.dPNR.size() << "   "; cout << "SdNi " <<
        // Dnorm1(rs.bulk.dNNR.size(), &rs.bulk.dNNR[0]) / rs.bulk.dNNR.size() << "   ";
        cout << endl;
        // rs.ShowRes(resFIM.res);
    }

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM();
        rs.CalResFIM(resFIM, ctrl.current_dt, OCP_TRUE);
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    }

    if (((resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
          resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
          resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol) &&
         resFIM.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
            case 1:
                ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
                rs.ResetFIM();
                rs.CalResFIM(resFIM, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                if (ctrl.printLevel >= PRINT_MORE) {
                    OCP_WARNING("Reset Newton iteration!");
                }
                return OCP_FALSE;
            case 2:
                ctrl.current_dt /= 1;
                rs.ResetFIM();
                rs.CalResFIM(resFIM, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                if (ctrl.printLevel >= PRINT_MORE) {
                    OCP_WARNING("Reset Newton iteration!");
                }
                return OCP_FALSE;
            default:
                return OCP_TRUE;
                break;
        }
    } else {
        return OCP_FALSE;
    }
}

void OCP_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStepFIM();
    ctrl.CalNextTstepFIM(rs);
    ctrl.UpdateIters();
}

////////////////////////////////////////////
// OCP_FIMn
////////////////////////////////////////////


/// Setup FIMn
void OCP_FIMn::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateFIMn_IsoT();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM_IsoT(myLS);
    // Allocate memory for residual of FIM
    resFIMn.Setup_IsoT(rs.GetBulkNum(), rs.GetWellNum(), rs.GetComNum());

    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}


void OCP_FIMn::InitReservoir(Reservoir& rs) const { rs.InitFIMn(); }


/// Prepare for Assembling matrix.
void OCP_FIMn::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResFIM(resFIMn, dt, OCP_TRUE);
}


/// Assemble Matrix
void OCP_FIMn::AssembleMat(LinearSystem&    myLS,
                           const Reservoir& rs,
                           const OCP_DBL&   dt) const
{
    rs.AssembleMatFIMn(myLS, dt);
    myLS.AssembleRhsAccumulate(resFIMn.res);
}

/// Solve the linear system.
void OCP_FIMn::SolveLinearSystem(LinearSystem& myLS,
                                 Reservoir&    rs,
                                 OCPControl&   ctrl) const
{
#ifdef DEBUG
    myLS.CheckEquation();
#endif // DEBUG

    myLS.AssembleMatLinearSolver();

    GetWallTime Timer;
    Timer.Start();
    int status = myLS.Solve();
    if (status < 0) {
        status = myLS.GetNumIters();
    }
    // cout << "LS step = " << status << endl;

#ifdef DEBUG
    myLS.OutputLinearSystem("testA_FIMn_new.out", "testb_FIMn_new.out");
    myLS.OutputSolution("testx_FIMn_new.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionFIMn(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.GetSolution01FIM(myLS.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    myLS.ClearData();
}

/// Update properties of fluids.
OCP_BOOL OCP_FIMn::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(OCP_TRUE, OCP_FALSE) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIMn();
        rs.CalResFIM(resFIMn, dt, OCP_TRUE);
        cout << "Cut time stepsize and repeat!\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    rs.CalFlashFIMn();
    rs.CalKrPcFIM();
    rs.CalRock();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResFIM(resFIMn, dt, OCP_FALSE);

    return OCP_TRUE;
}


/// Finish a Newton-Raphson iteration.
OCP_BOOL OCP_FIMn::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();

    if (ctrl.printLevel >= PRINT_MOST) {

        if (OCP_TRUE) {
            vector<OCP_INT> totalPhaseNum(3, 0);
            vector<OCP_INT> ltotalPhaseNum(3, 0);
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (rs.bulk.NRphaseNum[n] == 1) ltotalPhaseNum[0]++;
                if (rs.bulk.NRphaseNum[n] == 2) ltotalPhaseNum[1]++;
                if (rs.bulk.NRphaseNum[n] == 3) ltotalPhaseNum[2]++;
                if (rs.bulk.phaseNum[n] == 1) totalPhaseNum[0]++;
                if (rs.bulk.phaseNum[n] == 2) totalPhaseNum[1]++;
                if (rs.bulk.phaseNum[n] == 3) totalPhaseNum[2]++;
            }
            cout << to_string(totalPhaseNum[0]) + " (" +
                        to_string((totalPhaseNum[0] - ltotalPhaseNum[0]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      "
                 << to_string(totalPhaseNum[1]) + " (" +
                        to_string((totalPhaseNum[1] - ltotalPhaseNum[1]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      "
                 << to_string(totalPhaseNum[2]) + " (" +
                        to_string((totalPhaseNum[2] - ltotalPhaseNum[2]) * 100.0 /
                                  rs.GetBulkNum()) +
                        "%)      ";
        }

        if (OCP_TRUE) {
            OCP_USI eVnum = 0;
            OCP_USI eNnum = 0;
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (resFIMn.resRelV[n] < ctrl.ctrlNR.NRtol) eVnum++;
                if (resFIMn.resRelN[n] < ctrl.ctrlNR.NRtol) eNnum++;
            }
            cout << "eVnum : " << setprecision(ceil(log(rs.GetBulkNum()) / log(10)))
                 << fixed << setw(10) << (1 - eVnum * 1.0 / rs.GetBulkNum()) * 100.0
                 << "%   eNnum : " << setw(10)
                 << (1 - eNnum * 1.0 / rs.GetBulkNum()) * 100.0 << "%" << endl;
        }

        const USI sp = rs.grid.GetNumDigitIJK();
        USI       tmpI, tmpJ, tmpK;
        string    tmps;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, rs.bulk.index_maxNRdSSP);
        tmps = GetIJKformat(tmpI, tmpJ, tmpK, sp);
        cout << "### NR : " + to_string(ctrl.iterNR) + "    Res:    " << setprecision(2)
             << scientific << resFIMn.maxRelRes0_v << setw(12) << resFIMn.maxRelRes_v
             << setw(12) << resFIMn.maxRelRes_mol << setw(11) << resFIMn.maxWellRelRes_mol
             << setw(20) << NRdPmax << setw(11) << NRdNmax << setw(11) << NRdSmax
             << setw(40) << sqrt(rs.bulk.NRdSSP) / rs.bulk.numBulk << setw(12)
             << rs.bulk.maxNRdSSP << "  " << rs.bulk.phaseNum[rs.bulk.index_maxNRdSSP]
             << "  " << rs.bulk.NRphaseNum[rs.bulk.index_maxNRdSSP] << "  " << tmps
             << "   " << rs.bulk.ePEC[rs.bulk.index_maxNRdSSP] << endl
             << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, dSn);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "S: " << tmps << setw(12) << resFIMn.resRelV[dSn] << setw(12)
             << resFIMn.resRelN[dSn] << setw(12) << rs.bulk.ePEC[dSn] << setw(5)
             << rs.bulk.phaseNum[dSn] << setw(5) << rs.bulk.NRphaseNum[dSn] << setw(12)
             << rs.bulk.dPNR[dSn] << setw(8) << "  |" << setw(12) << rs.bulk.S[dSn * 3]
             << setw(12) << rs.bulk.dSNR[dSn * 3] << setw(12) << rs.bulk.dSNRP[dSn * 3]
             << setw(12) << rs.bulk.S[dSn * 3 + 1] << setw(12)
             << rs.bulk.dSNR[dSn * 3 + 1] << setw(12) << rs.bulk.dSNRP[dSn * 3 + 1]
             << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIMn.maxId_v);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "V: " << tmps << setw(12) << resFIMn.resRelV[resFIMn.maxId_v] << setw(12)
             << resFIMn.resRelN[resFIMn.maxId_v] << setw(12) << rs.bulk.ePEC[resFIMn.maxId_v]
             << setw(5) << rs.bulk.phaseNum[resFIMn.maxId_v] << setw(5)
             << rs.bulk.NRphaseNum[resFIMn.maxId_v] << setw(12)
             << rs.bulk.dPNR[resFIMn.maxId_v] << setw(8) << "  |" << setw(12)
             << rs.bulk.S[resFIMn.maxId_v * 3] << setw(12)
             << rs.bulk.dSNR[resFIMn.maxId_v * 3] << setw(12)
             << rs.bulk.dSNRP[resFIMn.maxId_v * 3] << setw(12)
             << rs.bulk.S[resFIMn.maxId_v * 3 + 1] << setw(12)
             << rs.bulk.dSNR[resFIMn.maxId_v * 3 + 1] << setw(12)
             << rs.bulk.dSNRP[resFIMn.maxId_v * 3 + 1] << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIMn.maxId_mol);
        tmps = GetIJKformat(tmpI, tmpJ, tmpK, sp);
        cout << "N: " << tmps << setw(12) << resFIMn.resRelV[resFIMn.maxId_mol] << setw(12)
             << resFIMn.resRelN[resFIMn.maxId_mol] << setw(12)
             << rs.bulk.ePEC[resFIMn.maxId_mol] << setw(5)
             << rs.bulk.phaseNum[resFIMn.maxId_mol] << setw(5)
             << rs.bulk.NRphaseNum[resFIMn.maxId_mol] << setw(12)
             << rs.bulk.dPNR[resFIMn.maxId_mol] << setw(8) << "  |" << setw(12)
             << rs.bulk.S[resFIMn.maxId_mol * 3] << setw(12)
             << rs.bulk.dSNR[resFIMn.maxId_mol * 3] << setw(12)
             << rs.bulk.dSNRP[resFIMn.maxId_mol * 3] << setw(12)
             << rs.bulk.S[resFIMn.maxId_mol * 3 + 1] << setw(12)
             << rs.bulk.dSNR[resFIMn.maxId_mol * 3 + 1] << setw(12)
             << rs.bulk.dSNRP[resFIMn.maxId_mol * 3 + 1] << endl;
        // cout << "Res2   " << Dnorm2(resFIM.res.size(), &resFIM.res[0]) /
        // resFIM.res.size() << "   "; cout << "SdP " << Dnorm1(rs.bulk.dPNR.size(),
        // &rs.bulk.dPNR[0]) / rs.bulk.dPNR.size() << "   "; cout << "SdNi " <<
        // Dnorm1(rs.bulk.dNNR.size(), &rs.bulk.dNNR[0]) / rs.bulk.dNNR.size() << "   ";
        cout << endl;
        // rs.ShowRes(resFIM.res);
    }

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIMn();
        rs.CalResFIM(resFIMn, ctrl.current_dt, OCP_TRUE);
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    }

    if (((resFIMn.maxRelRes_v <= resFIMn.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
          resFIMn.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
          resFIMn.maxRelRes_mol <= ctrl.ctrlNR.NRtol) &&
         resFIMn.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
            case 1:
                ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
                rs.ResetFIMn();
                rs.CalResFIM(resFIMn, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                if (ctrl.printLevel >= PRINT_MORE) {
                    OCP_WARNING("Reset Newton iteration!");
                }
                return OCP_FALSE;
            case 2:
                ctrl.current_dt /= 1;
                rs.ResetFIMn();
                rs.CalResFIM(resFIMn, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                if (ctrl.printLevel >= PRINT_MORE) {
                    OCP_WARNING("Reset Newton iteration!");
                }
                return OCP_FALSE;
            default:
                return OCP_TRUE;
                break;
        }
    } else {
        return OCP_FALSE;
    }
}

/// Finish a time step.
void OCP_FIMn::FinishStep(Reservoir& rs, OCPControl& ctrl) const
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStepFIMn();
    ctrl.CalNextTstepFIM(rs);
    ctrl.UpdateIters();
}


////////////////////////////////////////////
// OCP_AIMc
////////////////////////////////////////////

void OCP_AIMc::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateAIMc_IsoT();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM_IsoT(myLS);
    // Allocate memory for resiual of FIM
    resAIMc.Setup_IsoT(rs.GetBulkNum(), rs.GetWellNum(), rs.GetComNum());

    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void OCP_AIMc::InitReservoir(Reservoir& rs) const { rs.InitAIMc(); }

void OCP_AIMc::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResAIMc(resAIMc, dt, OCP_TRUE);

    // Set FIM Bulk
    rs.CalCFL(dt);
    rs.SetupWellBulk();
    rs.SetupFIMBulk();
    //  Calculate FIM Bulk properties
    rs.CalFlashAIMcI();
    rs.CalKrPcAIMcI();
    rs.UpdateLastStepFIM();

    // rs.bulk.ShowFIMBulk(OCP_FALSE);
}

void OCP_AIMc::AssembleMat(LinearSystem&    myLS,
                           const Reservoir& rs,
                           const OCP_DBL&   dt) const
{
    rs.AssembleMatAIMc(myLS, dt);
    myLS.AssembleRhsCopy(resAIMc.res);
}

void OCP_AIMc::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    myLS.CheckEquation();
#endif // DEBUG

    myLS.AssembleMatLinearSolver();

    GetWallTime Timer;
    Timer.Start();
    int status = myLS.Solve();
    if (status < 0) {
        status = myLS.GetNumIters();
    }
    // cout << "LS step = " << status << endl;

#ifdef DEBUG
    myLS.OutputLinearSystem("testA_AIMc_new.out", "testb_AIMc_new.out");
    myLS.OutputSolution("testx_AIMc_new.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionAIMc(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    myLS.ClearData();
}

OCP_BOOL OCP_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(OCP_TRUE, OCP_FALSE) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetAIMc();
        rs.CalResAIMc(resAIMc, dt, OCP_TRUE);
        cout << "Cut time stepsize and repeat!\n";
        return OCP_FALSE;
    }

    rs.CalFlashAIMcI();
    rs.CalFlashAIMcEp();   // think more and more
    rs.CalKrPcAIMcI();
    rs.CalRock();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResAIMc(resAIMc, dt, OCP_FALSE);
    return OCP_TRUE;
}

OCP_BOOL OCP_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI       dSn;
    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

#ifdef DEBUG
    cout << "### DEBUG: Residuals = " << setprecision(3) << scientific
         << resAIMc.maxRelRes0_v << "  " << resAIMc.maxRelRes_v << "  "
         << resAIMc.maxRelRes_mol << "  " << NRdPmax << "  " << NRdSmax << endl;
#endif

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetAIMc();
        rs.CalResAIMc(resAIMc, ctrl.current_dt, OCP_TRUE);
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    }

    //cout << scientific << setprecision(16) << resFIM.maxRelRes0_v << "   "
    //    << resFIM.maxRelRes_v << "   "
    //    << resFIM.maxRelRes_mol << "   "
    //    << resFIM.maxWellRelRes_mol << "   "
    //    << NRdPmax << "   "
    //    << NRdSmax << endl;

    if (((resAIMc.maxRelRes_v <= resAIMc.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
          resAIMc.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
          resAIMc.maxRelRes_mol <= ctrl.ctrlNR.NRtol) &&
         resAIMc.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
            case 1:
                ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
                rs.ResetAIMc();
                rs.CalResAIMc(resAIMc, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                return OCP_FALSE;
            case 2:
                // ctrl.current_dt /= 1;
                rs.ResetAIMc();
                rs.CalResAIMc(resAIMc, ctrl.current_dt, OCP_TRUE);
                ctrl.ResetIterNRLS();
                return OCP_FALSE;
            default:
                // Update IMPEC Bulk Properties
                rs.CalFlashAIMcEa();
                rs.CalKrPcAIMcE();
                return OCP_TRUE;
                break;
        }
    } else {
        return OCP_FALSE;
    }
}

/// Finish a time step.
void OCP_AIMc::FinishStep(Reservoir& rs, OCPControl& ctrl) const
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStepAIMc();
    ctrl.CalNextTstepFIM(rs);
    ctrl.UpdateIters();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/