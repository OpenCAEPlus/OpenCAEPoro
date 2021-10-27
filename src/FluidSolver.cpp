#include "FluidSolver.hpp"

// FluidSolver

void FluidSolver::Prepare(Reservoir& rs, OCP_DBL& dt)
{
    switch (method)
    {
    case IMPEC:
        rs.Prepare(dt);
        break;
    case FIM:
        rs.Prepare(dt);
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
        rs.AssembleMatIMPEC(FLSolver, dt);
        break;
    case FIM:
        rs.AssembleMatFIM(FLSolver, dt);
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
        fim.SolveLinearSystem(FLSolver, rs, ctrl);
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
        return fim.UpdateProperty(rs, dt);
        break;
    default:
        OCP_ABORT("Wrong method!");
        break;
    }
}

bool FluidSolver::FinishNR()
{
    switch (method)
    {
    case IMPEC:
        return impes.FinishNR();
    case FIM:
        return fim.FinishNR();
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
    switch (method)
    {
    case IMPEC:
        rs.AllocateMatIMPEC(FLSolver);
        break;
    case FIM:
        rs.AllocateMatFIM(FLSolver);
        break;
    default:
        OCP_ABORT("Wrong Method!");
        break;
    }
}


void FluidSolver::SetupParamLS(const string& dir, const string& file)
{
    switch (method)
    {
    case IMPEC:
        FLSolver.SetupParam(dir, file);
        break;
    case FIM:
        FLSolver.SetupParamB(dir, file);
        break;
    default:
        OCP_ABORT("Wrong Mthod!");
        break;
    }
}

void FluidSolver::InitReservoir(Reservoir& rs) const
{
    switch (method)
    {
    case IMPEC:
        rs.InitIMPEC();
        break;
    case FIM:
        rs.InitFIM();
        break;
    default:
        OCP_ABORT("Wrong Mthod!");
        break;
    }
}


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

    // lsolver.Free_Fasp();

    ctrl.UpdateIterLS(status);

#endif // __SOLVER_FASP__

    rs.GetSolution_IMPEC(lsolver.GetSolution());
    lsolver.ClearData();
}

void OCP_FIM::SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCP_Control& ctrl)
{

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




bool OCP_FIM::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{

}

bool OCP_FIM::FinishNR()
{

}