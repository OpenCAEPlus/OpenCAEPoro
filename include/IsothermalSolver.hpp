/*! \file    IsothermalSolver.hpp
 *  \brief   IsothermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ISOTHERMALSOLVER_HEADER__
#define __ISOTHERMALSOLVER_HEADER__

// OpenCAEPoro header files
#include "OCPFluidMethod.hpp"

/// IsothermalSolver class for fluid solution method.
class IsothermalSolver
{
public:
    /// Setup the fluid solver.
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);
    /// Initialize the Reservoir and prepare variables for some method.
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for assembling Mat.
    void Prepare(Reservoir& rs, OCPControl& ctrl);
    /// Assemble Mat.
    void AssembleMat(const Reservoir& rs, OCPControl& ctrl);
    /// Solve the linear system in single problem.
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluid.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish the Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish the current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

private:
    USI          method = FIM;
    LinearSystem LSolver;
    LinearSystem auxLSolver;
    IsoT_IMPEC   impec;
    IsoT_FIM     fim;
    IsoT_FIMn    fim_n;
    IsoT_AIMc    aimc;
};

#endif /* end if __ISOTHERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/