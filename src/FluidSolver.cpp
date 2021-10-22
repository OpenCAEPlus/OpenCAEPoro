#include "FluidSolver.hpp"

// FluidSolver

void FluidSolver::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    switch (method)
    {
    case IMPEC:
        impes.Prepare(rs, dt);
        break;
    case FIM:
        break;
    default:
        OCP_ABORT("Wrong method!");
        break;
    }
}


void FluidSolver::AssembleMat(const Reservoir& rs, const OCP_DBL& dt)
{
    switch (method)
    {
    case IMPEC:
        impes.AssembleMat(FLSolver, rs, dt);
        break;
    case FIM:
        break;
    default:
        OCP_ABORT("Wrong method!");
        break;
    }
}

void FluidSolver::SolveLinearSystem(Reservoir& rs, OCP_Control& ctrl)
{
    switch (method)
    {
    case IMPEC:
        impes.SolveLinearSystem(FLSolver, rs, ctrl);
        break;
    case FIM:
        break;
    default:
        OCP_ABORT("Wrong method!");
        break;
    }
}

bool FluidSolver::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{
    switch (method)
    {
    case IMPEC:
        return impes.UpdateProperty(rs, dt);
    case FIM:
        break;
    default:
        OCP_ABORT("Wrong method!");
        break;
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

void FluidSolver::AllocateMat(const Reservoir& rs) 
{
    rs.AllocateMat(FLSolver);
}


void FluidSolver::SetupParamLS(const string& dir, const string& file)
{
    FLSolver.SetupParam(dir, file);
}


void OCP_IMPEC::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    rs.Prepare(dt);
}

void OCP_IMPEC::AssembleMat(LinearSolver& lsolver, const Reservoir& rs, const OCP_DBL& dt)
{
    rs.AssembleMatIMPEC(lsolver, dt);
}

void OCP_IMPEC::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCP_Control& ctrl)
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

    lsolver.Free_Fasp();

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

    rs.CalFlashIMPEC();
    rs.CalVpore();

    // fouth check: Volume error check
    if (!rs.CheckVe(0.01)) {
        cout << "###WARNING: volume error is too big\n";
        dt /= 2;
        rs.ResetVal02();
        return false;
    }

    rs.CalKrPc();
    rs.CalConnFlux();
    return true;
}


void OCP_IMPEC::SetupParam(LinearSolver& lsolver, const string& dir, const string& file)
{
    lsolver.SetupParam(dir, file);
}

void OCP_IMPEC::AllocateMat(LinearSolver& lsolver, const Reservoir& rs)
{
    lsolver.AllocateMem(rs.bulk.GetBulkNum() + rs.wellgroup.GetWellNum());
    rs.conn.AllocateMat(lsolver);
    rs.wellgroup.AllocateMat(lsolver);
    lsolver.AllocateColValMem();
}
