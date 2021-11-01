/*! \file    OCPFluidMethod.hpp
 *  \brief   Solution methods in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPFLUIDMETHOD_HEADER__
#define __OCPFLUIDMETHOD_HEADER__

 // OpenCAEPoro header files
#include "LinearSolver.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilTiming.hpp"

 // OpenCAEPoro header files
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

    void FinishStep(Reservoir& rs, OCPControl& ctrl);
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
    bool FinishNR(Reservoir& rs, const OCPControl& ctrl);

    void FinishStep(Reservoir& rs, OCPControl& ctrl);


private:
    /// Resiual for FIM
    ResFIM                   resFIM;
};



#endif /* end if __OCP_METHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/