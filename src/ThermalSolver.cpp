/*! \file    ThermalSolver.cpp
 *  \brief   ThermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoro header files
#include "ThermalSolver.hpp"

void ThermalSolver::SetupMethod(Reservoir& rs, const OCPControl& ctrl)
{
    fim.Setup(rs, LSolver, ctrl);
}

void ThermalSolver::InitReservoir(Reservoir& rs) const { fim.InitReservoir(rs); }

void ThermalSolver::Prepare(Reservoir& rs, const OCPControl& ctrl)
{
    fim.Prepare(rs, ctrl);
}

void ThermalSolver::AssembleMat(const Reservoir& rs, OCPControl& ctrl)
{
    const OCP_DBL dt = ctrl.GetCurDt();

    GetWallTime timer;
    timer.Start();

    fim.AssembleMat(LSolver, rs, ctrl.GetCurTime() + dt, dt);

    ctrl.RecordTimeAssembleMat(timer.Stop() / 1000);
}

void ThermalSolver::SolveLinearSystem(Reservoir& rs, OCPControl& ctrl)
{
    fim.SolveLinearSystem(LSolver, rs, ctrl);
}

/// Update properties of fluid.
OCP_BOOL ThermalSolver::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    return fim.UpdateProperty(rs, ctrl);
}

/// Finish the Newton-Raphson iteration.
OCP_BOOL ThermalSolver::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    return fim.FinishNR(rs, ctrl);
}

/// Finish the current time step.
void ThermalSolver::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    fim.FinishStep(rs, ctrl);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/