#include "OCPFluidMethod.hpp"


void OCP_IMPEC::Setup(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl)
{
    // Allocate Memory
    rs.AllocateRsIMPEC();
    ls.SetupParam(ctrl.GetWorkDir(), ctrl.GetLsFile());
    rs.AllocateMatIMPEC(ls);
    // For Fasp
    ls.AllocateFasp();
}


void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    OCP_DBL cfl = rs.CalCFL01(dt);
    if (cfl > 1) dt /= (cfl + 1);
}


void OCP_IMPEC::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs,
    OCPControl& ctrl)
{
#ifdef DEBUG
    solver.CheckVal();
#endif // DEBUG

#ifdef __SOLVER_FASP__

    lsolver.AssembleMat_Fasp();
    GetWallTime Timer;
    Timer.Start();
    int status = lsolver.FaspSolve();
    ctrl.UpdateTimeLS(Timer.Stop() / 1000);

#ifdef DEBUG
    lsolver.PrintfMatCSR("testA.out", "testb.out");
    lsolver.PrintfSolution("testx.out");
#endif // DEBUG

    ctrl.UpdateIterLS(status);

#endif // __SOLVER_FASP__

    rs.GetSolution_IMPEC(lsolver.GetSolution());
    lsolver.ClearData();
}


bool OCP_IMPEC::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{
    // first check : Pressure check
    OCP_INT flagCheck = rs.CheckP();
    switch (flagCheck) {
    case 1:
        cout << "well change" << endl;
        dt /= 2;
        return false;
    case 2:
        cout << "well change" << endl;
        dt /= 2;
        return false;
    default:
        break;
    }

    rs.CalFLux();

    // second check : CFL check
    OCP_DBL cfl = rs.CalCFL01(dt);
    if (cfl > 1) {
        dt /= 2;
        rs.ResetVal();
        cout << "CFL is too big" << endl;
        return false;
    }

    rs.MassConseve(dt);

    // third check: Ni check
    if (!rs.CheckNi()) {
        dt /= 2;
        rs.ResetVal01();
        cout << "Negative Ni occurs\n";
        return false;
    }

    rs.CalFlash();
    rs.CalVpore();

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        dt /= 2;
        rs.ResetVal02();
        return false;
    }

    rs.CalKrPc();
    rs.CalConnFlux();
    return true;
}


void OCP_IMPEC::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStep();
    ctrl.CalNextTstepIMPEC(rs);
    ctrl.UpdateIters();
}



void OCP_FIM::Setup(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    rs.AllocateRsFIM();
    // Read Ls Params
    ls.SetupParamB(ctrl.GetWorkDir(), ctrl.GetLsFile());
    // Allocate memory for internal matrix structure
    rs.AllocateMatFIM(ls);
    // Allocate memory for BFasp matrix structure
    ls.AllocateBFasp();
    // Allocate memory for resiual of FIM
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    resFIM.res.resize(num);
}


void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResFIM(resFIM, dt);
    resFIM.maxRelRes0_v = resFIM.maxRelRes_v;
}


void OCP_FIM::AssembleMat(LinearSolver& lsolver, const Reservoir& rs,
    const OCP_DBL& dt) const
{
    rs.AssembleMatFIM(lsolver, dt);
    lsolver.AssembleRhs_BFasp(resFIM.res);
}


void OCP_FIM::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    solver.CheckVal();
#endif // DEBUG

#ifdef __SOLVER_FASP__

    lsolver.AssembleMat_BFasp();

#ifdef _DEBUG
    //lsolver.PrintfMatCSR("testA.out", "testb.out");
#endif // DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = lsolver.BFaspSolve();

#ifdef _DEBUG
    //lsolver.PrintfSolution("testx.out");
#endif // DEBUG

    ctrl.UpdateTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);

#endif // __SOLVER_FASP__

    rs.GetSolution_FIM(lsolver.GetSolution());
    lsolver.ClearData();
}


bool OCP_FIM::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{

    // first check : Pressure check.
    OCP_INT flagCheck = rs.CheckP();
    switch (flagCheck) {
    case 1:
        cout << "well change" << endl;
        dt /= 2;
        return false;
    case 2:
        cout << "well change" << endl;
        dt /= 2;
        return false;
    default:
        break;
    }
    // Second check: Ni check.
    if (!rs.CheckNi()) {
        dt /= 2;
        cout << "Negative Ni occurs\n";
        return false;
    }
    rs.CalFlashDeriv();
    rs.CalKrPcDeriv();
    rs.CalVpore();
    rs.CalWellFlux();
    rs.CalResFIM(resFIM, dt);
}


bool OCP_FIM::FinishNR()
{
    if (resFIM.maxRelRes_v < resFIM.maxRelRes0_v * 1E-3 ||
        resFIM.maxRelRes_v < 1E-3 ||
        resFIM.maxRelRes_mol < 1E-3) {
        return true;
    }
    else {
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