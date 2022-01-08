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

void FluidSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    method = ctrl.GetMethod();
    
    switch (method) {
        case IMPEC:
            cout << "Calling IMPEC method ..." << endl;
            FLSolver.SetupLinearSolver(1, ctrl.GetWorkDir(), ctrl.GetLsFile());
            impec.Setup(rs, FLSolver, ctrl);
            break;
        case FIM:
            cout << "Calling FIM method ..." << endl;
            FLSolver.SetupLinearSolver(2, ctrl.GetWorkDir(), ctrl.GetLsFile());
            fim.Setup(rs, FLSolver, ctrl);
            break;
        default:
            OCP_ABORT("Wrong method type!");
    }
    FLSolver.AllocateLinearSolver();
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
            OCP_ABORT("Wrong method type!");
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
            OCP_ABORT("Wrong method type!");
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
            OCP_ABORT("Wrong method type!");
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
            OCP_ABORT("Wrong method type!");
    }
}

bool FluidSolver::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case IMPEC:
            return impec.UpdateProperty(rs, ctrl);
        case FIM:
            return fim.UpdateProperty(rs, ctrl);
        default:
            OCP_ABORT("Wrong method type!");
    }
}


bool FluidSolver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case IMPEC:
            return impec.FinishNR();
        case FIM:
            return fim.FinishNR(rs, ctrl);
        default:
            OCP_ABORT("Wrong method type!");
    }
}

void FluidSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    switch (method) {
        case IMPEC:
            return impec.FinishStep(rs, ctrl);
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
/*----------------------------------------------------------------------------*/
