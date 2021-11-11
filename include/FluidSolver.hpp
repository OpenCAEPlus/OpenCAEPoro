/*! \file    Solver.hpp
 *  \brief   Solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLUIDSOLVER_HEADER__
#define __FLUIDSOLVER_HEADER__

// OpenCAEPoro header files
#include "OCPFluidMethod.hpp"

class FluidSolver
{
public:
    /// Setup Method
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);
    /// Initialize the Reservoir and prepare variables for some method.
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for assembling Mat.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble Mat.
    void AssembleMat(const Reservoir& rs, const OCP_DBL& dt);
    /// Solve the linear system in single problem.
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluid.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
    /// Determine if NR iteration finishes.
    bool FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);
    
private:
    USI          method = FIM;
    LinearSolver FLSolver;
    OCP_IMPEC    impec;
    OCP_FIM      fim;
};

#endif /* end if __FLUIDSOLVER_HEADER__ */
