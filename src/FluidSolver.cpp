#include "FluidSolver.hpp"

// FluidSolver

void FluidSolver::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    switch (method) {
        case IMPEC:
            impec.Prepare(rs, dt);
            break;
        case FIM:
            fim.Prepare(rs, dt);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void FluidSolver::AssembleMat(const Reservoir& rs, const OCP_DBL& dt)
{
    switch (method) {
        case IMPEC:
            rs.AssembleMatIMPEC(FLSolver, dt);
            break;
        case FIM:
            fim.AssembleMat(FLSolver, rs, dt);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void FluidSolver::SolveLinearSystem(Reservoir& rs, OCP_Control& ctrl)
{
    switch (method) {
        case IMPEC:
            impec.SolveLinearSystem(FLSolver, rs, ctrl);
            break;
        case FIM:
            fim.SolveLinearSystem(FLSolver, rs, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


bool FluidSolver::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{
    switch (method) {
        case IMPEC:
            return impec.UpdateProperty(rs, dt);
        case FIM:
            return fim.UpdateProperty(rs, dt);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


bool FluidSolver::FinishNR()
{
    switch (method) {
        case IMPEC:
            return impec.FinishNR();
        case FIM:
            return fim.FinishNR();
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void FluidSolver::FinishStep(Reservoir& rs, OCP_Control& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStep();
    ctrl.CalNextTstep(rs);
    ctrl.UpdateIters();
}


void FluidSolver::SetupMethod(const Reservoir& rs, const OCP_Control& ctrl)
{
    method = ctrl.GetMethod();
    switch (method)
    {
    case IMPEC:
        break;
    case FIM:
        fim.Setup(rs);
        break;
    default:
        OCP_ABORT("Wrong Method!");
        break;
    }
}


void FluidSolver::AllocateMat(const Reservoir& rs)
{
    switch (method) {
        case IMPEC:
            rs.AllocateMatIMPEC(FLSolver);
            break;
        case FIM:
            rs.AllocateMatFIM(FLSolver);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void FluidSolver::SetupParamLS(const string& dir, const string& file)
{
    switch (method) {
        case IMPEC:
            FLSolver.SetupParam(dir, file);
            break;
        case FIM:
            FLSolver.SetupParamB(dir, file);
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void FluidSolver::InitReservoir(Reservoir& rs) const
{
    switch (method) {
        case IMPEC:
            rs.InitIMPEC();
            break;
        case FIM:
            rs.InitFIM();
            break;
        default:
            OCP_ABORT("Wrong method!");
    }
}


void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt) { rs.Prepare(dt); }


void OCP_IMPEC::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs,
                                  OCP_Control& ctrl)
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


void OCP_FIM::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCP_Control& ctrl)
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


void OCP_FIM::Setup(const Reservoir& rs)
{
    OCP_USI num = (rs.GetBulkNum() + rs.GetWellNum()) * (rs.GetComNum() + 1);
    res.resize(num, 0);
    relRes.resize(num, 0);
}


void OCP_FIM::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.PrepareWell();
    rs.CalResFIM(res, dt);
    CalMaxRes(rs);
    maxRes0 = maxRes;
}


void OCP_FIM::AssembleMat(LinearSolver& lsolver, const Reservoir& rs,
                          const OCP_DBL& dt) const
{
    rs.AssembleMatFIM(lsolver, dt);
    lsolver.AssembleRhs_BFasp(res);
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
    CalMaxRes(rs);
}


bool OCP_FIM::FinishNR()
{
    if (maxRes < 1E-3 * maxRes0 || maxRes < 1E-3 || maxRelRes < 1E-3)
        return true;
}


void OCP_FIM::CalMaxRes(const Reservoir& rs)
{
    maxRes = 0;
    maxRelRes = 0;
    OCP_USI num = res.size();

    for (OCP_USI n = 0; n < num; n++) {
        maxRes = max(res[n], maxRes);
        maxRelRes = max(relRes[n], maxRelRes);
    }
}