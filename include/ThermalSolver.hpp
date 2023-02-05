/*! \file    ThermalSolver.hpp
 *  \brief   ThermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __THERMALSOLVER_HEADER__
#define __THERMALSOLVER_HEADER__

// OpenCAEPoro header files
#include "ThermalMethod.hpp"

/// ThermalSolver class for fluid solution method.
class ThermalSolver
{
public:
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);
    void InitReservoir(Reservoir& rs) const;
    void Prepare(Reservoir& rs, const OCPControl& ctrl);
    void AssembleMat(const Reservoir& rs, OCPControl& ctrl);
    /// Solve the linear system in single problem.
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluid.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish the Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish the current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    LinearSystem LSolver;
    LinearSystem auxLSolver;
    T_FIM        fim;
};

#endif /* end if __THERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/