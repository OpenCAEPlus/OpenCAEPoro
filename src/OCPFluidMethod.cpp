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


void OCP_FIM::Setup(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl)
{
    // Allocate Memory
    rs.AllocateRsFIM();
    ls.SetupParamB(ctrl.GetWorkDir(), ctrl.GetLsFile());
    rs.AllocateMatFIM(ls);
    // For Block Fasp
    ls.AllocateBFasp();
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    res.resize(num);
}


void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResFIM(res, dt);
}


void OCP_FIM::AssembleMat(LinearSolver& lsolver, const Reservoir& rs,
    const OCP_DBL& dt) const
{
    rs.AssembleMatFIM(lsolver, dt);
    lsolver.AssembleRhs_BFasp(res);
}


void OCP_FIM::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    solver.CheckVal();
#endif // DEBUG

#ifdef __SOLVER_FASP__

    lsolver.AssembleMat_BFasp();

#ifdef DEBUG
    lsolver.PrintfMatCSR("testA.out", "testb.out");
#endif // DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = lsolver.BFaspSolve();

#ifdef DEBUG
    lsolver.PrintfSolution("testx.out");
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
        rs.ResetVal01();
        cout << "Negative Ni occurs\n";
        return false;
    }
    rs.CalFlashDeriv();
    rs.CalKrPcDeriv();
    rs.CalVpore();
    rs.CalResFIM(res, dt);
}


bool OCP_FIM::FinishNR()
{

}


