#include "OCPFluidMethod.hpp"

void OCP_IMPEC::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Memory
    rs.AllocateAuxIMPEC();
    rs.AllocateMatIMPEC(ls);
}

void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    OCP_DBL cfl = rs.CalCFL01IMPEC(dt);
    if (cfl > 1) dt /= (cfl + 1);
}

void OCP_IMPEC::SolveLinearSystem(LinearSystem& ls, Reservoir& rs,
                                  OCPControl& ctrl)
{
#ifdef _DEBUG
    ls.CheckVal();
#endif // DEBUG

    ls.AssembleMatLinearSolver();
    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();
    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();
#ifdef DEBUG
    ls.OutputLinearSystem("testA.out", "testb.out");
    ls.OutputSolution("testx.out");
#endif // DEBUG

    rs.GetSolutionIMPEC(ls.GetSolution());
    ls.ClearData();
}


bool OCP_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP();
    switch (flagCheck) {
    case 1:
        cout << "well change" << endl;
        dt /= 2;
        rs.ResetVal00IMPEC();
        return false;
    case 2:
        cout << "well change" << endl;
        dt /= 1;
        rs.ResetVal00IMPEC();
        // rs.ResetWellIMPEC();
        return false;
    default:
        break;
    }

    rs.CalFLuxIMPEC();

    // second check : CFL check
    OCP_DBL cfl = rs.CalCFL01IMPEC(dt);
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

    rs.CalFlashIMPEC();
    rs.CalVpore();

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        dt /= 2;
        rs.ResetVal03IMPEC();
        return false;
    }

    rs.CalKrPc();
    rs.CalConnFluxIMPEC();
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

void OCP_FIM::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateAuxFIM();
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM(ls);
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);
}

void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
}

void OCP_FIM::AssembleMat(LinearSystem& lsolver, const Reservoir& rs,
                          const OCP_DBL& dt) const
{
    rs.AssembleMatFIM(lsolver, dt);
    lsolver.AssembleRhs(resFIM.res);
}

void OCP_FIM::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{

    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");

#ifdef DEBUG
    ls.CheckVal();
#endif // DEBUG

    ls.AssembleMatLinearSolver();

#ifdef DEBUG
    // ls.PrintfMatCSR("testA.out", "testb.out");
#endif // DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();

#ifdef DEBUG
    // ls.PrintfSolution("testx.out");
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    rs.GetSolutionFIM(ls.GetSolution(), ctrl.ctrlNR.NRdPmax, ctrl.ctrlNR.NRdSmax);
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

bool OCP_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // Second check: Ni check.
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetFIM(false);
        rs.CalResFIM(resFIM, dt);
        cout << "Negative Ni occurs\n";
        return false;
    }
    rs.CalFlashDerivFIM();
    rs.CalKrPcDerivFIM();
    rs.CalVpore();
    rs.CalWellTrans();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);
    return true;
}

bool OCP_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
 
    if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {    
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        rs.ResetFIM(false);
        rs.CalResFIM(resFIM, ctrl.current_dt);
        ctrl.ResetIterNR();
        cout << "NR Failed, Cut time Step and reset!\n";
        return false;
    }

    OCP_DBL NRdPmax = rs.GetNRdPmax();
    OCP_DBL NRdSmax = rs.GetNRdSmax();

#ifdef DEBUG
    cout << "### DEBUG: Residuals = " << scientific << resFIM.maxRelRes0_v << "  "
         << resFIM.maxRelRes_v << "  " << resFIM.maxRelRes_mol << "  " << NRdSmax
         << "  " << NRdPmax << endl;
#endif

    if (resFIM.maxRelRes_v <= resFIM.maxRelRes0_v * ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_v <= ctrl.ctrlNR.NRtol ||
        resFIM.maxRelRes_mol <= ctrl.ctrlNR.NRtol ||
        (NRdPmax <= ctrl.ctrlNR.NRdPmin && NRdSmax <= ctrl.ctrlNR.NRdSmin)) {

        OCP_INT flagCheck = rs.CheckP();
        switch (flagCheck) {
        case 1:
            cout << "well change, Repeat --- 1" << endl;
            ctrl.current_dt /= 2;
            rs.ResetFIM(true);
            rs.CalResFIM(resFIM, ctrl.current_dt);
            ctrl.ResetIterNR();
            return false;
        case 2:
            cout << "well change, Repeat --- 2" << endl;
            ctrl.current_dt /= 1;
            rs.ResetFIM(true);
            rs.CalResFIM(resFIM, ctrl.current_dt);
            ctrl.ResetIterNR();
            return false;
        default:
            break;
        }

        return true;
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