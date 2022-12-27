/*! \file    OCPFluidMethod.hpp
 *  \brief   Declaration of solution methods for fluid part in OpenCAEPoro
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
#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"

/// OCP_IMPEC is IMPEC (implicit pressure explict saturation) method.
class OCP_IMPEC
{
public:
    /// Setup IMPEC
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);

    /// Init
    void InitReservoir(Reservoir& rs) const;

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl);

    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Determine if NR iteration finishes.
    OCP_BOOL FinishNR(const Reservoir& rs);

    void FinishStep(Reservoir& rs, OCPControl& ctrl);
};

/// OCP_FIM is FIM (Fully Implicit Method).
class OCP_FIM
{
public:
    /// Setup FIM
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);

    /// Init
    void InitReservoir(Reservoir& rs) const;

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl) const;

    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);


protected:
    /// Residual for FIM
    OCPRes resFIM;
};

class OCP_FIMn
{
public:
    /// Setup FIM
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);

    /// Init
    void InitReservoir(Reservoir& rs) const;

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl) const;

    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl) const;

protected:
    /// Residual for FIMn
    OCPRes resFIMn;
};

class OCP_AIMc
{
public:
    /// Setup AIMc
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);

    /// Init
    void InitReservoir(Reservoir& rs) const;

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl);

    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);

    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl) const;

protected:
    /// Residual for AIMc
    OCPRes resAIMc;
};

#endif /* end if __OCPFLUIDMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/