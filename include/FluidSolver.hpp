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
#include "LinearSolver.hpp"
#include "OCPControl.hpp"
#include "OCPOutput.hpp"
#include "Reservoir.hpp"
#include "UtilTiming.hpp"

/// OCP_IMPEC is IMPEC (implict pressure explict saturation) method.
class OCP_IMPEC
{
public:

    void Setup(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl);

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Solve the linear system.
    void SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCPControl& ctrl);

    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);

    /// Determine if NR iteration finishes.
    bool FinishNR() { return true; }
};

/// OCP_FIM is FIM (Fully Implicit Method).
class OCP_FIM
{
public:

    /// Setup FIM
    void Setup(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl);

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Matrix
    void AssembleMat(LinearSolver& lsolver, const Reservoir& rs,
                     const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCPControl& ctrl);

    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);

    /// Determine if NR iteration finishes.
    bool FinishNR();

    /// Calculate maximum Res.
    void CalMaxRes(const Reservoir& rs);

private:
    vector<OCP_DBL> res;
    vector<OCP_DBL> relRes;
    OCP_DBL         maxRes0;
    OCP_DBL         maxRes;
    OCP_DBL         maxRelRes;
};

class FluidSolver
{
public:
    /// Prepare for assembling Mat.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble Mat.
    void AssembleMat(const Reservoir& rs, const OCP_DBL& dt);
    /// Solve the linear system in single problem.
    void SolveLinearSystem(Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluid.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
    /// Determine if NR iteration finishes.
    bool FinishNR();
    /// Finish current time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

    /// Setup Method
    void SetupMethod(Reservoir& rs, const OCPControl& ctrl);

    /// Allocate Mat
    void AllocateMat(const Reservoir& rs);

    /// Setup linear solver params.
    void SetupParamLS(const string& dir, const string& file);

    void InitReservoir(Reservoir& rs) const;

private:
    USI          method = FIM;
    LinearSolver FLSolver;
    OCP_IMPEC    impec;
    OCP_FIM      fim;
};

#endif /* end if __FLUIDSOLVER_HEADER__ */
