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
#include "OCPControl.hpp"
#include "OCPOutput.hpp"
#include "Reservoir.hpp"
#include "LinearSolver.hpp"
#include "UtilTiming.hpp"
#include "LinearSolver.hpp"


 /// OCP_IMPEC is IMPEC (implict pressure explict saturation) method.
class OCP_IMPEC
{
public:

    /// Prepare and assemble the mat.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Mat.
    void AssembleMat(LinearSolver& lsolver, const Reservoir& rs, const OCP_DBL& dt);

    /// Solve the linear system.
    void SolveLinearSystem(LinearSolver& lsolver, Reservoir& rs, OCP_Control& ctrl);

    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);

    /// Setup parameters needed for IMPEC.
    void SetupParam(LinearSolver& lsolver, const string& dir, const string& file);

    /// Allocate maximal possible memory needed by linear solver.
    void AllocateMat(LinearSolver& lsolver, const Reservoir& rs);


};



class FluidSolver
{
public:
    /// Prepare for assembling Mat.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble Mat.
    void AssembleMat(const Reservoir& rs, const OCP_DBL& dt);
    /// Solve the linear system in single problem.
    void SolveLinearSystem(Reservoir& rs, OCP_Control& ctrl);
    /// Update properties of fluid.
    bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
    /// Finish current time step.
    void FinishStep(Reservoir& rs, OCP_Control& ctrl);

    /// Allocate Mat
    void AllocateMat(const Reservoir& rs);

    /// Setup linear solver params.
    void SetupParamLS(const string& dir, const string& file);

private:
    USI             method = IMPEC;
	LinearSolver	FLSolver;
    OCP_IMPEC       impes;

};


#endif /* end if __FLUIDSOLVER_HEADER__ */