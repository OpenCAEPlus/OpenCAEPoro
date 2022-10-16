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
    rs.AllocateAuxIMPEC();
    // Allocate Memory of Matrix for IMPEC
    rs.AllocateMatIMPEC(myLS);

    myLS.SetupLinearSolver(SCALARFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

/// Init
void OCP_IMPEC::InitReservoir(Reservoir& rs) const
{
    rs.InitIMPEC();
}

void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    OCP_DBL cfl = rs.CalCFL(dt);
    if (cfl > 1) dt /= (cfl + 1);
}

void OCP_IMPEC::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
{
#ifdef _DEBUG
    myLS.CheckEquation();
#endif // DEBUG

    myLS.AssembleMatLinearSolver();

#ifdef _DEBUG
    // myLS.OutputLinearSystem("testA.out", "testb.out");
#endif // _DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = myLS.Solve();
    if (status < 0) {
        status = myLS.GetNumIters();
    }

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
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
            // Switch Well opt Mode, or close the crossflow perforation
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

    rs.MassConseveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal02IMPEC();
        cout << "Negative Ni occurs\n";
        return OCP_FALSE;
    }

    rs.CalVpore();
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
    // rs.allWells.ShowWellStatus(rs.bulk);

    return OCP_TRUE;
}


OCP_BOOL OCP_IMPEC::FinishNR(const Reservoir& rs)
{ 
    //for (USI j = 0; j < rs.bulk.numPhase; j++)
    //    cout << rs.bulk.totalPhaseNum[j] << "   ";
    //cout << endl;
    return OCP_TRUE; 
}


OCP_BOOL OCP_IMPEC::UpdateProperty01(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP();
    switch (flagCheck) {
    case 1:
        // Negative Bulk P, well P, or Perforation P
        dt /= 2;
        return OCP_FALSE;
    case 2:
        // Switch Well opt Mode, or close the crossflow perforation
        dt /= 1;
        return OCP_FALSE;
    default:
        // All right
        break;
    }

    rs.CalFLuxIMPEC();

    // second check : CFL check
    OCP_DBL cfl = rs.CalCFL(dt);
    if (cfl > 1) {
        dt /= 2;
        rs.ResetVal01IMPEC();
        cout << "CFL is too big" << endl;
        return OCP_FALSE;
    }

    rs.MassConseveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal02IMPEC();
        cout << "Negative Ni occurs\n";
        return OCP_FALSE;
    }

    rs.CalVpore();
    rs.CalFlashIMPEC();
    rs.CalKrPc();
    rs.CalConnFluxIMPEC();

    return OCP_TRUE;
}

OCP_BOOL OCP_IMPEC::FinishNR01(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;
    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        // continue NR 
        rs.bulk.ResetNi();
        rs.allWells.CalTrans(rs.bulk);
        rs.allWells.CalFlux(rs.bulk);
        rs.allWells.CalProdWeight(rs.bulk);
        rs.allWells.CaldG(rs.bulk);
        return OCP_FALSE;
    }
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
    rs.AllocateAuxFIM();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM(myLS);
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);

    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void OCP_FIM::InitReservoir(Reservoir& rs) const
{
    rs.InitFIM();
}

void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
}

void OCP_FIM::AssembleMat(LinearSystem& myLS, const Reservoir& rs,
                          const OCP_DBL& dt) const
{
    rs.AssembleMatFIM(myLS, dt);
    myLS.AssembleRhs(resFIM.res);
}

void OCP_FIM::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl) const
{
#ifdef _DEBUG
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
    myLS.OutputLinearSystem("testA_FIM.out", "testb_FIM.out");
    myLS.OutputSolution("testx_FIM.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
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
        rs.ResetFIM(OCP_FALSE);
        rs.CalResFIM(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time step size and repeat! current dt = " << fixed << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    rs.CalFlashDerivFIM();
    rs.CalKrPcDerivFIM();
    rs.CalVpore();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);

    //if (rs.bulk.NRdPmax < 1E-4) {
    //    // correct
    //    cout << "correct" << endl;
    //    rs.bulk.CorrectNi(resFIM.res);
    //    // Update reservoir properties
    //    rs.CalFlashDerivFIM();
    //    rs.CalKrPcDerivFIM();
    //    rs.CalVpore();
    //    rs.CalWellTrans();
    //    rs.CalWellFlux();
    //    rs.CalResFIM(resFIM, dt);
    //}

    return OCP_TRUE;
}

OCP_BOOL OCP_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();
    OCP_DBL NRdSmaxP = rs.GetNRdSmaxP();



    if (ctrl.printLevel > 1) {

        if (OCP_TRUE) {
            vector<OCP_INT> totalPhaseNum(3, 0);
            vector<OCP_INT> ltotalPhaseNum(3, 0);
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (rs.bulk.NRphaseNum[n] == 0)  ltotalPhaseNum[0]++;
                if (rs.bulk.NRphaseNum[n] == 1)  ltotalPhaseNum[1]++;
                if (rs.bulk.NRphaseNum[n] == 2)  ltotalPhaseNum[2]++;
                if (rs.bulk.phaseNum[n] == 0)  totalPhaseNum[0]++;
                if (rs.bulk.phaseNum[n] == 1)  totalPhaseNum[1]++;
                if (rs.bulk.phaseNum[n] == 2)  totalPhaseNum[2]++;
            }
            cout << to_string(totalPhaseNum[0]) + " (" + to_string((totalPhaseNum[0] - ltotalPhaseNum[0]) * 100.0 / rs.GetBulkNum())
                + "%)      " << to_string(totalPhaseNum[1]) + " (" + to_string((totalPhaseNum[1] - ltotalPhaseNum[1]) * 100.0 / rs.GetBulkNum())
                + "%)      " << to_string(totalPhaseNum[2]) + " (" + to_string((totalPhaseNum[2] - ltotalPhaseNum[2]) * 100.0 / rs.GetBulkNum())
                + "%)      " ;
        }

        if (OCP_TRUE) {
            OCP_USI eVnum = 0;
            OCP_USI eNnum = 0;
            for (OCP_USI n = 0; n < rs.GetBulkNum(); n++) {
                if (rs.bulk.eV[n] < ctrl.ctrlNR.NRtol)
                    eVnum++;
                if (rs.bulk.eN[n] < ctrl.ctrlNR.NRtol)
                    eNnum++;
            }
            cout << "eVnum : " << setprecision(ceil(log(rs.GetBulkNum()) / log(10))) << fixed << setw(10) << (1 - eVnum * 1.0 / rs.GetBulkNum()) * 100.0 << "%   eNnum : "
                << setw(10) << (1 - eNnum * 1.0 / rs.GetBulkNum()) * 100.0 << "%" << endl;
        }

        const USI sp = rs.grid.GetNumDigitIJK();
        USI tmpI, tmpJ, tmpK;
        string tmps;
     
        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, rs.bulk.index_maxNRdSSP);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "### NR : " + to_string(ctrl.iterNR) + "    Res:    " << setprecision(2) << scientific << resFIM.maxRelRes0_v << setw(12)
            << resFIM.maxRelRes_v << setw(12) << resFIM.maxRelRes_mol << setw(11) << resFIM.maxWellRelRes_mol
            << setw(20) << NRdPmax << setw(11) << NRdNmax << setw(11) << NRdSmax          
            << setw(40) << sqrt(rs.bulk.NRdSSP) / rs.bulk.numBulk << setw(12) << rs.bulk.maxNRdSSP << "  "
            << rs.bulk.phaseNum[rs.bulk.index_maxNRdSSP] << "  " << rs.bulk.NRphaseNum[rs.bulk.index_maxNRdSSP] << "  "
            << tmps << "   " << rs.bulk.ePEC[rs.bulk.index_maxNRdSSP] << endl << endl;


        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, dSn);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "S: " << tmps << setw(12) << rs.bulk.eV[dSn] << setw(12) << rs.bulk.eN[dSn]
            << setw(12) << rs.bulk.ePEC[dSn] << setw(5) << rs.bulk.phaseNum[dSn] << setw(5) << rs.bulk.NRphaseNum[dSn] << setw(12) << rs.bulk.dPNR[dSn]
            << setw(8) << "  |" << setw(12) << rs.bulk.S[dSn * 3] << setw(12) << rs.bulk.dSNR[dSn * 3] << setw(12) << rs.bulk.dSNRP[dSn * 3]
            << setw(12) << rs.bulk.S[dSn * 3 + 1] << setw(12) << rs.bulk.dSNR[dSn * 3 + 1] << setw(12) << rs.bulk.dSNRP[dSn * 3 + 1] << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIM.maxId_v);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "V: " << tmps << setw(12) << rs.bulk.eV[resFIM.maxId_v] << setw(12) << rs.bulk.eN[resFIM.maxId_v]
            << setw(12) << rs.bulk.ePEC[resFIM.maxId_v] << setw(5) << rs.bulk.phaseNum[resFIM.maxId_v] << setw(5) << rs.bulk.NRphaseNum[resFIM.maxId_v] << setw(12) << rs.bulk.dPNR[resFIM.maxId_v]
            << setw(8) << "  |" << setw(12) << rs.bulk.S[resFIM.maxId_v * 3] << setw(12) << rs.bulk.dSNR[resFIM.maxId_v * 3] << setw(12) << rs.bulk.dSNRP[resFIM.maxId_v * 3]
            << setw(12) << rs.bulk.S[resFIM.maxId_v * 3 + 1] << setw(12) << rs.bulk.dSNR[resFIM.maxId_v * 3 + 1] << setw(12) << rs.bulk.dSNRP[resFIM.maxId_v * 3 + 1] << endl;

        rs.grid.GetIJKBulk(tmpI, tmpJ, tmpK, resFIM.maxId_mol);
        tmps = GetIJKformat(to_string(tmpI), to_string(tmpJ), to_string(tmpK), sp);
        cout << "N: " << tmps << setw(12) << rs.bulk.eV[resFIM.maxId_mol] << setw(12) << rs.bulk.eN[resFIM.maxId_mol]
            << setw(12) << rs.bulk.ePEC[resFIM.maxId_mol] << setw(5) << rs.bulk.phaseNum[resFIM.maxId_mol] << setw(5) << rs.bulk.NRphaseNum[resFIM.maxId_mol] << setw(12) << rs.bulk.dPNR[resFIM.maxId_mol]
            << setw(8) << "  |" << setw(12) << rs.bulk.S[resFIM.maxId_mol * 3] << setw(12) << rs.bulk.dSNR[resFIM.maxId_mol * 3] << setw(12) << rs.bulk.dSNRP[resFIM.maxId_mol * 3]
            << setw(12) << rs.bulk.S[resFIM.maxId_mol * 3 + 1] << setw(12) << rs.bulk.dSNR[resFIM.maxId_mol * 3 + 1] << setw(12) << rs.bulk.dSNRP[resFIM.maxId_mol * 3 + 1] << endl;
        // cout << "Res2   " << Dnorm2(resFIM.res.size(), &resFIM.res[0]) / resFIM.res.size() << "   ";
        // cout << "SdP " << Dnorm1(rs.bulk.dPNR.size(), &rs.bulk.dPNR[0]) / rs.bulk.dPNR.size() << "   ";
        // cout << "SdNi " << Dnorm1(rs.bulk.dNNR.size(), &rs.bulk.dNNR[0]) / rs.bulk.dNNR.size() << "   ";
        cout << endl;
        // rs.ShowRes(resFIM.res);
    }


    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(OCP_FALSE);
        rs.CalResFIM(resFIM, ctrl.current_dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  current dt = " 
            << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    }

    if (((resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol)
        && resFIM.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin && fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
            case 1:
                ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
                rs.ResetFIM(OCP_TRUE);
                rs.CalResFIM(resFIM, ctrl.current_dt);
                resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
                ctrl.ResetIterNRLS();
                cout << "-----" << endl;
                return OCP_FALSE;
            case 2:
                ctrl.current_dt /= 1;
                rs.ResetFIM(OCP_TRUE);
                rs.CalResFIM(resFIM, ctrl.current_dt);
                resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
                ctrl.ResetIterNRLS();
                cout << "-----" << endl;
                return OCP_FALSE;
            default:
                return OCP_TRUE;
                break;
        }       
    } else {
        return OCP_FALSE;
    }
}

void OCP_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl) const
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


void OCP_FIMn::InitReservoir(Reservoir& rs)const 
{ 
    rs.InitFIM_n(); 
}


/// Assemble Matrix
void OCP_FIMn::AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const
{
    rs.AssembleMatFIM_n(myLS, dt);
    myLS.AssembleRhs(resFIM.res);
}

/// Solve the linear system.
void OCP_FIMn::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl) const
{
#ifdef _DEBUG
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
    myLS.OutputLinearSystem("testA.out", "testb.out");
    myLS.OutputSolution("testx.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionFIM_n(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
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
        rs.ResetFIM(OCP_FALSE);
        rs.CalResFIM(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    rs.CalFlashDerivFIM_n();
    rs.CalKrPcDerivFIM();
    rs.CalVpore();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);

    return OCP_TRUE;
}


////////////////////////////////////////////
// OCP_AIMc
////////////////////////////////////////////

void OCP_AIMc::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateAuxAIMc();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM(myLS);
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);

    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());

}

void OCP_AIMc::InitReservoir(Reservoir& rs) const
{
    rs.InitAIMc();
}

void OCP_AIMc::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalWellFlux();
    rs.CalResAIMc(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;

    // Set FIM Bulk
    rs.CalCFL(dt);
    rs.SetupWellBulk();
    rs.SetupFIMBulk();
    //rs.bulk.FIMBulk.resize(rs.bulk.numBulk, 0);
    //for (OCP_USI n = 0; n < rs.bulk.numBulk; n++) {
    //    rs.bulk.FIMBulk[n] = n;
    //    rs.bulk.map_Bulk2FIM[n] = n;
    //}
    //rs.bulk.numFIMBulk = rs.bulk.numBulk;
    // Calculate FIM Bulk properties
	rs.CalFlashDerivAIMc();
	rs.CalKrPcDerivAIMc();
	// rs.bulk.CheckDiff();
	rs.UpdateLastStepFIM();

    rs.bulk.ShowFIMBulk(OCP_FALSE);
}

void OCP_AIMc::AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const
{
    rs.AssembleMatAIMc(myLS, dt);
    myLS.AssembleRhs(resFIM.res);
}

void OCP_AIMc::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
{
#ifdef _DEBUG
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
    myLS.OutputLinearSystem("testA.out", "testb.out");
    myLS.OutputSolution("testx.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionAIMc(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.GetSolution01FIM(myLS.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    myLS.ClearData();
}

OCP_BOOL OCP_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt; 

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(OCP_TRUE, OCP_FALSE) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(OCP_FALSE);
        rs.CalResAIMc(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!\n";
        return OCP_FALSE;
    }

    rs.CalFlashDerivAIMc();
    rs.CalKrPcDerivAIMc();
    rs.CalFlashAIMc();
    // Important, Pj must be updated with current and last Pc for IMPEC Bulk
    rs.UpdatePj();

    rs.CalVpore();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResAIMc(resFIM, dt);
    return OCP_TRUE;
}

OCP_BOOL OCP_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();
    const OCP_DBL NRdSmax = rs.GetNRdSmaxP();

//#ifdef _DEBUG
    cout << "### DEBUG: Residuals = " << setprecision(3) << scientific << resFIM.maxRelRes0_v << "  "
       << resFIM.maxRelRes_v << "  " << resFIM.maxRelRes_mol << "  " << NRdPmax
       << "  " << NRdSmax << endl;
    //for (OCP_USI n = 0; n < resFIM.res.size(); n++) {
    //    cout << resFIM.res[n] << endl;
    //}
//#endif

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(OCP_FALSE);
        rs.CalResAIMc(resFIM, ctrl.current_dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!\n";
        return OCP_FALSE;
    }

    if (((resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol)
        && resFIM.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
        case 1:
            ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
            rs.ResetFIM(OCP_TRUE);
            rs.CalResAIMc(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return OCP_FALSE;
        case 2:
            // ctrl.current_dt /= 1;
            rs.ResetFIM(OCP_TRUE);
            rs.CalResAIMc(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return OCP_FALSE;
        default:
            // Update IMPEC Bulk Properties
            rs.CalFlashAIMc01();
            rs.CalKrPcAIMc();
            return OCP_TRUE;
            break;
        }
    }
    else {
        //// Set FIMBulk
        //rs.CalCFL(ctrl.current_dt);
        //rs.SetupWellBulk();
        //rs.SetupFIMBulk(OCP_TRUE);
        //// Calculate FIM Bulk properties
        //rs.CalFlashDerivAIMc();
        //rs.CalKrPcDerivAIMc();
        //// Show FIM Bulk
        //rs.bulk.ShowFIMBulk(OCP_FALSE);
        return OCP_FALSE;
    }
}


////////////////////////////////////////////
// OCP_AIMs
////////////////////////////////////////////

void OCP_AIMs::Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for AIMt 
    rs.AllocateAuxAIMs();
    // Allocate Memory of Matrix for FIM
    rs.AllocateMatFIM(myLS);
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);

    myLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void OCP_AIMs::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    // for Well
    rs.PrepareWell();
    rs.CalWellFlux();

    // Set FIM Bulk
    rs.SetupWellBulk();
    rs.SetupFIMBulk();
    rs.SetupFIMBulkBoundAIMs();

    rs.bulk.ShowFIMBulk();

    // Calculate Resiual
    rs.CalResAIMs(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;

    // Calculat property of FIM Bulk
    rs.CalFlashDerivAIM(OCP_TRUE);
    rs.CalKrPcDerivAIM(OCP_TRUE);

    // Store particular property of FIM Bulk
    rs.bulk.UpdateLastStepAIM();
}

void OCP_AIMs::AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt)
{
    rs.AssembleMatAIMs(myLS, resFIM.res, dt);
    myLS.AssembleRhs(resFIM.res);
}

void OCP_AIMs::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
{
#ifdef _DEBUG
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
    myLS.OutputLinearSystem("testA.out", "testb.out");
    myLS.OutputSolution("testx.out");
    myLS.CheckSolution();
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionAIMs(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.GetSolution01FIM(myLS.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    myLS.ClearData();
}

OCP_BOOL OCP_AIMs::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(OCP_TRUE, OCP_FALSE) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetValAIM();
        rs.CalResAIMs(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!  --  01\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    rs.CalFlashIMPEC();
    rs.CalKrPc();
    rs.CalFlashDerivAIM(OCP_TRUE);
    rs.CalKrPcDerivAIM(OCP_TRUE);
    rs.CalVpore();
    rs.CalFLuxIMPEC();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalConnFluxIMPEC();
    rs.CalResAIMs(resFIM, dt);
      
    return OCP_TRUE;
}

OCP_BOOL OCP_AIMs::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL NRdPmax = rs.GetNRdPmax();
    OCP_DBL NRdSmax = rs.GetNRdSmaxP();

#ifdef _DEBUG
     cout << "### DEBUG: Residuals = " << setprecision(3) << scientific << resFIM.maxRelRes0_v << "  "
        << resFIM.maxRelRes_v << "  " << resFIM.maxRelRes_mol << "  " << NRdSmax
        << "  " << NRdPmax << endl;
#endif

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetValAIM();
        rs.CalResAIMs(resFIM, ctrl.current_dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time stepsize and repeat!\n";
        return OCP_FALSE;
    }

    if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(OCP_FALSE, OCP_TRUE);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
        case 1:
            ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
            rs.ResetValAIM();
            rs.CalResAIMs(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return OCP_FALSE;
        case 2:
            ctrl.current_dt /= 1;
            rs.ResetValAIM();
            rs.CalResAIMs(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return OCP_FALSE;
        default:
            break;
        }
        // Mass Conserve
        rs.bulk.InFIMNi();
        rs.conn.MassConserveIMPEC(rs.bulk, ctrl.current_dt);
        rs.bulk.OutFIMNi();
        if (!rs.CheckNi()) {
            ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
            rs.ResetValAIM();
            rs.CalResAIMs(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            cout << "Cut time stepsize and repeat!\n";
            return OCP_FALSE;
        }
        if (!rs.CheckVe(0.01)) {
            // cout << ctrl.GetCurTime() << "Days" << "=======" << endl;

            rs.AddFIMBulk();
            rs.SetupFIMBulkBoundAIMs();

            rs.bulk.ShowFIMBulk();

            // Calculate Resiual
            rs.CalResAIMs(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;

            // Calculat property of FIM Bulk
            rs.CalFlashDerivAIM(OCP_TRUE);
            rs.CalKrPcDerivAIM(OCP_TRUE);

            // Store particular property of FIM Bulk
            rs.bulk.UpdateLastStepAIM();



            // ctrl.current_dt /= 2;
            rs.ResetValAIM();
            ctrl.ResetIterNRLS();

            cout << "Cut time stepsize and repeat!  --  02\n";
            return OCP_FALSE;
        }
        OCP_DBL cfl = rs.CalCFLAIM(ctrl.current_dt);
        if (cfl > 1) {
            cout << "CFL is too big" << endl;
            ctrl.current_dt /= 2;
            rs.ResetValAIM();
            ctrl.ResetIterNRLS();           
            return OCP_FALSE;
        }
        return OCP_TRUE;
    }
    else {
        return OCP_FALSE;
    }
}

void OCP_AIMs::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    // rs.GetNTQT(ctrl.current_dt);

    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStepAIM();
    ctrl.CalNextTstepIMPEC(rs);
    ctrl.UpdateIters();

    
}

////////////////////////////////////////////
// OCP_AIMt
////////////////////////////////////////////

void OCP_AIMt::Setup(Reservoir& rs, LinearSystem& myLS, LinearSystem& myAuxLS, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for AIMt 
    rs.AllocateAuxAIMt();
    // Allocate Memory of Matrix for IMPEC
    rs.AllocateMatIMPEC(myLS);
    // Allocate memory for internal matrix structure for local FIM
    rs.AllocateMatAIMt(myAuxLS);
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetMaxFIMBulk() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);

    myLS.SetupLinearSolver(SCALARFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
    myAuxLS.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), "./bsr.fasp");
}

void OCP_AIMt::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    OCP_DBL cfl = rs.CalCFL(dt);
    if (cfl > 1) dt /= (cfl + 1);

    // setup WellbulkId
    rs.SetupWellBulk();
}


OCP_BOOL OCP_AIMt::UpdateProperty(Reservoir& rs, OCPControl& ctrl, LinearSystem& myAuxLS)
{
    OCP_DBL& dt = ctrl.current_dt;

    rs.CalFLuxIMPEC();
    rs.CalCFL(dt);
    rs.MassConseveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal03IMPEC();
        cout << "Negative Ni occurs\n";
        return OCP_FALSE;
    }

    rs.CalVpore();
    rs.CalFlashIMPEC();

    // Perform FIM in local grid
    // Init
    rs.SetupFIMBulk();  

    // cout << "FIM Bulk : " << rs.bulk.numFIMBulk << endl;

    //for (USI i = 0; i < rs.bulk.numFIMBulk; i++) {
    //    cout << rs.bulk.P[rs.bulk.FIMBulk[i]] << "   ";
    //}
    //cout << endl << endl;

    rs.bulk.FlashDerivAIM(OCP_FALSE);
    rs.CalKrPcDerivAIM(OCP_FALSE);
   
    rs.CalResAIMt(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;

    //for (USI i = 0; i < rs.bulk.numFIMBulk; i++) {
    //    cout << rs.bulk.P[rs.bulk.FIMBulk[i]] << "   ";
    //}   
    //cout << endl << endl;
    ctrl.iterNR = 0;
    // cout << ctrl.iterNR << "   " << resFIM.maxRelRes0_v << endl;
    while (OCP_TRUE) {
        rs.AssembleMatAIMt(myAuxLS, dt);
        myAuxLS.AssembleRhs(resFIM.res);
        myAuxLS.AssembleMatLinearSolver();
        int status = myAuxLS.Solve();

        rs.GetSolutionAIMt(myAuxLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
        myAuxLS.ClearData();

        // third check: Ni check
        if (!rs.CheckNi()) {
            dt /= 2;
            rs.ResetVal03IMPEC();
            cout << "Negative Ni occurs\n";
            return OCP_FALSE;
        }
       
        rs.bulk.FlashDerivAIM(OCP_FALSE);
        rs.CalKrPcDerivAIM(OCP_FALSE);
        rs.CalVpore();
        rs.CalWellTrans();
        rs.CalWellFlux();
        rs.CalResAIMt(resFIM, dt);

        //for (USI i = 0; i < rs.bulk.numFIMBulk; i++) {
        //    cout << fixed << rs.bulk.P[rs.bulk.FIMBulk[i]] << "   ";
        //}
        //cout << scientific << resFIM.maxRelRes_v;
        //cout << endl << endl;
        OCP_DBL NRdPmax = rs.GetNRdPmax();
        OCP_DBL NRdSmax = rs.GetNRdSmaxP();
        ctrl.iterNR++;
        //cout << ctrl.iterNR << "   " << resFIM.maxRelRes_v << "   "
        //    << resFIM.maxRelRes_mol << "   " << NRdPmax << "   "
        //    << NRdSmaxP << "   " << endl;
        
        if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
            resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
            resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
            (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {
            break;
        }
        if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
            ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
            rs.ResetVal03IMPEC();
            cout << "Local FIM Failed!" << endl;
            return OCP_FALSE;
        }
    }
    
    // Pressure check
    OCP_INT flagCheck = rs.CheckP(OCP_TRUE, OCP_TRUE);
    switch (flagCheck) {
    case 1:
        // Negative Bulk P, well P, or Perforation P
        dt /= 2;
        return OCP_FALSE;
    case 2:
        // Switch Well opt Mode, or close the crossflow perforation
        dt /= 1;
        return OCP_FALSE;
    default:
        // All right
        break;
    }

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


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/