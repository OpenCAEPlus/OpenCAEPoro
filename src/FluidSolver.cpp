/*! \file    FluidSolver.cpp
 *  \brief   FluidSolver class definition
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "FluidSolver.hpp"

/// Setup solution methods, including IMPEC and FIM.
void FluidSolver::SetupMethod(Reservoir &rs, const OCPControl &ctrl)
{
    method = ctrl.GetMethod();

    switch (method)
    {
    case AIMt:
        aimt.Setup(rs, FLSolver, auxFLSolver, ctrl);
        FLSolver.SetupLinearSolver(1, ctrl.GetWorkDir(), ctrl.GetLsFile());
        auxFLSolver.SetupLinearSolver(2, ctrl.GetWorkDir(), "./bsr.fasp");
        break;
    case IMPEC:
        impec.Setup(rs, FLSolver, ctrl);
        FLSolver.SetupLinearSolver(1, ctrl.GetWorkDir(), ctrl.GetLsFile());       
        break;
    case FIM_IMPEC:
        fimImpec.Setup(rs, FLSolver, auxFLSolver, ctrl);
        FLSolver.SetupLinearSolver(2, ctrl.GetWorkDir(), ctrl.GetLsFile());
        auxFLSolver.SetupLinearSolver(1, ctrl.GetWorkDir(), "./csr.fasp");
        break;
    case FIM:
    default:
        fim.Setup(rs, FLSolver, ctrl);
        FLSolver.SetupLinearSolver(2, ctrl.GetWorkDir(), ctrl.GetLsFile());        
        break;
    }   
}

/// Setup solution methods, including IMPEC and FIM.
void FluidSolver::InitReservoir(Reservoir &rs) const
{
    switch (method)
    {
    case IMPEC:
    case AIMt:
        rs.InitIMPEC();
        break;
    case FIM_IMPEC:
    case FIM:
        rs.InitFIM();
        break;
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/// Prepare solution methods, including IMPEC and FIM.
void FluidSolver::Prepare(Reservoir &rs, OCP_DBL &dt)
{
    switch (method)
    {
    case IMPEC:
        impec.Prepare(rs, dt);
        break;
    case FIM:
        fim.Prepare(rs, dt);
        break;
    case AIMt:
        aimt.Prepare(rs, dt);
        break;
    case FIM_IMPEC:
        fimImpec.Prepare(rs, dt, auxFLSolver);
        break;
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/// Assemble linear systems for IMPEC and FIM.
void FluidSolver::AssembleMat(const Reservoir &rs, const OCP_DBL &dt)
{
    switch (method)
    {
    case IMPEC:
    case AIMt:
        rs.AssembleMatIMPEC(FLSolver, dt);
        break;
    case FIM:
        fim.AssembleMat(FLSolver, rs, dt);
        break;
    case FIM_IMPEC:
        fimImpec.AssembleMat(FLSolver, rs, dt);
        break;
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/// Solve linear systems for IMPEC and FIM.
void FluidSolver::SolveLinearSystem(Reservoir &rs, OCPControl &ctrl)
{
    switch (method)
    {
    case IMPEC:
    case AIMt:
        impec.SolveLinearSystem(FLSolver, rs, ctrl);
        break;
    case FIM_IMPEC:
    case FIM:
        fim.SolveLinearSystem(FLSolver, rs, ctrl);
        break;
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/// Update physical properties for IMPEC and FIM.
bool FluidSolver::UpdateProperty(Reservoir &rs, OCPControl &ctrl)
{
    switch (method) {
        case IMPEC:
            return impec.UpdateProperty(rs, ctrl);
            // return impec.UpdateProperty01(rs, ctrl);
        case FIM:
            return fim.UpdateProperty(rs, ctrl);
        case AIMt:
            return aimt.UpdateProperty(rs, ctrl, auxFLSolver);
        case FIM_IMPEC:
            return fimImpec.UpdateProperty(rs, ctrl, auxFLSolver);
        default:
            OCP_ABORT("Wrong method type!");
    }
}

/// Finish up Newton-Raphson iteration for IMPEC and FIM.
bool FluidSolver::FinishNR(Reservoir &rs, OCPControl &ctrl)
{
    switch (method) {
        case IMPEC:
        case AIMt:
            return impec.FinishNR();
            // return impec.FinishNR01(rs, ctrl);
        case FIM:
            return fim.FinishNR(rs, ctrl);
        case FIM_IMPEC:
            return fimImpec.FinishNR(rs, ctrl, auxFLSolver);
        default:
            OCP_ABORT("Wrong method type!");
    }
}

/// Finish up time step for IMPEC and FIM.
void FluidSolver::FinishStep(Reservoir &rs, OCPControl &ctrl)
{
    switch (method)
    {
    case IMPEC:
    case AIMt:
        return impec.FinishStep(rs, ctrl);
    case FIM_IMPEC:
    case FIM:
        return fim.FinishStep(rs, ctrl);
    default:
        OCP_ABORT("Wrong method type!");
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/