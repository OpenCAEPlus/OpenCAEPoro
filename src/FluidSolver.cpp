#include "FluidSolver.hpp"


void FluidSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    method = ctrl.GetMethod();
    switch (method)
    {
    case IMPEC:
        impec.Setup(rs, FLSolver, ctrl);
        break;
    case FIM:
        fim.Setup(rs, FLSolver, ctrl);
        break;
    default:
        OCP_ABORT("Wrong Method!");
        break;
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


void FluidSolver::SolveLinearSystem(Reservoir& rs, OCPControl& ctrl)
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


void FluidSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    rs.UpdateLastStep();
    ctrl.CalNextTstep(rs);
    ctrl.UpdateIters();
}