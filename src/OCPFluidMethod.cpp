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
    myLS.OutputLinearSystem("testA.out", "testb.out");
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
    myLS.OutputSolution("testx.out");
#endif // DEBUG

    rs.GetSolutionIMPEC(myLS.GetSolution());
    myLS.ClearData();
}

bool OCP_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP(true, true);
    switch (flagCheck) {
        case 1:
            // Negative Bulk P, well P, or Perforation P
            dt /= 2;
            return false;
        case 2:
            // Switch Well opt Mode, or close the crossflow perforation
            dt /= 1;
            return false;
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
        return false;
    }

    rs.MassConseveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal02IMPEC();
        cout << "Negative Ni occurs\n";
        return false;
    }

    rs.CalVpore();
    rs.CalFlashIMPEC();

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        // cout << ctrl.GetCurTime() << "Days" << "=======" << endl;
        dt /= 2;
        rs.ResetVal03IMPEC();
        return false;
    }

    rs.CalKrPc();
    rs.CalConnFluxIMPEC();
    // rs.allWells.ShowWellStatus(rs.bulk);

    return true;
}

bool OCP_IMPEC::UpdateProperty01(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP();
    switch (flagCheck) {
    case 1:
        // Negative Bulk P, well P, or Perforation P
        dt /= 2;
        return false;
    case 2:
        // Switch Well opt Mode, or close the crossflow perforation
        dt /= 1;
        return false;
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
        return false;
    }

    rs.MassConseveIMPEC(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal02IMPEC();
        cout << "Negative Ni occurs\n";
        return false;
    }

    rs.CalVpore();
    rs.CalFlashIMPEC();
    rs.CalKrPc();
    rs.CalConnFluxIMPEC();

    return true;
}

bool OCP_IMPEC::FinishNR01(Reservoir& rs, OCPControl& ctrl)
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
        return false;
    }
    return true;
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

void OCP_FIM::SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl)
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

    rs.GetSolutionFIM(myLS.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.GetSolution01FIM(myLS.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    myLS.ClearData();
}

bool OCP_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;


    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(true, false) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(false);
        rs.CalResFIM(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!\n";
        return false;
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

    return true;
}

bool OCP_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();
    const OCP_DBL NRdSmax = rs.GetNRdSmax();

#ifdef DEBUG
    cout << "### DEBUG: Residuals = " << setprecision(3) << scientific << resFIM.maxRelRes0_v << "  "
        << resFIM.maxRelRes_v << "  " << resFIM.maxRelRes_mol << "  " << NRdPmax
        << "  " << NRdSmax << "    ";
    // cout << "Res2   " << Dnorm2(resFIM.res.size(), &resFIM.res[0]) / resFIM.res.size() << "   ";
    // cout << "SdP " << Dnorm1(rs.bulk.dPNR.size(), &rs.bulk.dPNR[0]) / rs.bulk.dPNR.size() << "   ";
    // cout << "SdNi " << Dnorm1(rs.bulk.dNiNR.size(), &rs.bulk.dNiNR[0]) / rs.bulk.dNiNR.size() << "   ";
    cout << endl;
    // rs.ShowRes(resFIM.res);
#endif

    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(false);
        rs.CalResFIM(resFIM, ctrl.current_dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time stepsize and repeat!\n";
        return false;
    }

    if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(false, true);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
            case 1:
                ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
                rs.ResetFIM(true);
                rs.CalResFIM(resFIM, ctrl.current_dt);
                resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
                ctrl.ResetIterNRLS();
                cout << "-----" << endl;
                return false;
            case 2:
                ctrl.current_dt /= 1;
                rs.ResetFIM(true);
                rs.CalResFIM(resFIM, ctrl.current_dt);
                resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
                ctrl.ResetIterNRLS();
                cout << "-----" << endl;
                return false;
            default:
                return true;
                break;
        }       
    } else {
        return false;
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

    rs.bulk.ShowFIMBulk(false);
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

bool OCP_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt; 

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(true, false) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(false);
        rs.CalResAIMc(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!\n";
        return false;
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
    return true;
}

bool OCP_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    const OCP_DBL NRdNmax = rs.GetNRdNmax();
    const OCP_DBL NRdSmax = rs.GetNRdSmax();

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
        rs.ResetFIM(false);
        rs.CalResAIMc(resFIM, ctrl.current_dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        ctrl.ResetIterNRLS();
        cout << "### WARNING: NR not fully converged! Cut time stepsize and repeat!\n";
        return false;
    }

    if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(false, true);
#if DEBUG
        if (flagCheck > 0) {
            cout << ">> Switch well constraint: Case " << flagCheck << endl;
        }
#endif

        switch (flagCheck) {
        case 1:
            ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
            rs.ResetFIM(true);
            rs.CalResAIMc(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return false;
        case 2:
            // ctrl.current_dt /= 1;
            rs.ResetFIM(true);
            rs.CalResAIMc(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return false;
        default:
            // Update IMPEC Bulk Properties
            rs.CalFlashAIMc01();
            rs.CalKrPcAIMc();
            return true;
            break;
        }
    }
    else {
        //// Set FIMBulk
        //rs.CalCFL(ctrl.current_dt);
        //rs.SetupWellBulk();
        //rs.SetupFIMBulk(true);
        //// Calculate FIM Bulk properties
        //rs.CalFlashDerivAIMc();
        //rs.CalKrPcDerivAIMc();
        //// Show FIM Bulk
        //rs.bulk.ShowFIMBulk(false);
        return false;
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
    rs.CalFlashDerivAIM(true);
    rs.CalKrPcDerivAIM(true);

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

bool OCP_AIMs::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check and bulk Pressure check
    if (!rs.CheckNi() || rs.CheckP(true, false) != 0) {
        dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetValAIM();
        rs.CalResAIMs(resFIM, dt);
        resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
        cout << "Cut time stepsize and repeat!  --  01\n";
        return false;
    }

    // Update reservoir properties
    rs.CalFlashIMPEC();
    rs.CalKrPc();
    rs.CalFlashDerivAIM(true);
    rs.CalKrPcDerivAIM(true);
    rs.CalVpore();
    rs.CalFLuxIMPEC();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalConnFluxIMPEC();
    rs.CalResAIMs(resFIM, dt);
      
    return true;
}

bool OCP_AIMs::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL NRdPmax = rs.GetNRdPmax();
    OCP_DBL NRdSmax = rs.GetNRdSmax();

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
        return false;
    }

    if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP(false, true);
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
            return false;
        case 2:
            ctrl.current_dt /= 1;
            rs.ResetValAIM();
            rs.CalResAIMs(resFIM, ctrl.current_dt);
            resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
            ctrl.ResetIterNRLS();
            return false;
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
            return false;
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
            rs.CalFlashDerivAIM(true);
            rs.CalKrPcDerivAIM(true);

            // Store particular property of FIM Bulk
            rs.bulk.UpdateLastStepAIM();



            // ctrl.current_dt /= 2;
            rs.ResetValAIM();
            ctrl.ResetIterNRLS();

            cout << "Cut time stepsize and repeat!  --  02\n";
            return false;
        }
        OCP_DBL cfl = rs.CalCFLAIM(ctrl.current_dt);
        if (cfl > 1) {
            cout << "CFL is too big" << endl;
            ctrl.current_dt /= 2;
            rs.ResetValAIM();
            ctrl.ResetIterNRLS();           
            return false;
        }
        return true;
    }
    else {
        return false;
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


bool OCP_AIMt::UpdateProperty(Reservoir& rs, OCPControl& ctrl, LinearSystem& myAuxLS)
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
        return false;
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

    rs.bulk.FlashDerivAIM(false);
    rs.CalKrPcDerivAIM(false);
   
    rs.CalResAIMt(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;

    //for (USI i = 0; i < rs.bulk.numFIMBulk; i++) {
    //    cout << rs.bulk.P[rs.bulk.FIMBulk[i]] << "   ";
    //}   
    //cout << endl << endl;
    ctrl.iterNR = 0;
    // cout << ctrl.iterNR << "   " << resFIM.maxRelRes0_v << endl;
    while (true) {
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
            return false;
        }
       
        rs.bulk.FlashDerivAIM(false);
        rs.CalKrPcDerivAIM(false);
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
        OCP_DBL NRdSmax = rs.GetNRdSmax();
        ctrl.iterNR++;
        //cout << ctrl.iterNR << "   " << resFIM.maxRelRes_v << "   "
        //    << resFIM.maxRelRes_mol << "   " << NRdPmax << "   "
        //    << NRdSmax << "   " << endl;
        
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
            return false;
        }
    }
    
    // Pressure check
    OCP_INT flagCheck = rs.CheckP(true, true);
    switch (flagCheck) {
    case 1:
        // Negative Bulk P, well P, or Perforation P
        dt /= 2;
        return false;
    case 2:
        // Switch Well opt Mode, or close the crossflow perforation
        dt /= 1;
        return false;
    default:
        // All right
        break;
    }

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        // cout << ctrl.GetCurTime() << "Days" << "=======" << endl;
        dt /= 2;
        rs.ResetVal03IMPEC();
        return false;
    }

    rs.CalKrPc();
    rs.CalConnFluxIMPEC();

    return true;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/