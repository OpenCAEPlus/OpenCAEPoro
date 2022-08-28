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
#include "UtilTiming.hpp"

/// OCP_IMPEC is IMPEC (implict pressure explict saturation) method.
class OCP_IMPEC
{
public:
    /// Setup IMPEC
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);

    /// Init
    void InitReservoir(Reservoir& rs) const;

    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl);

    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    bool UpdateProperty01(Reservoir& rs, OCPControl& ctrl);

    /// Determine if NR iteration finishes.
    bool FinishNR() { return true; }
    bool FinishNR01(Reservoir& rs, OCPControl& ctrl);

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
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Finish a Newton-Raphson iteration.
    bool FinishNR(Reservoir& rs, OCPControl& ctrl);

    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl) const;


protected:
    /// Resiual for FIM
    ResFIM resFIM;
};


class OCP_FIMn : public OCP_FIM
{
public:
    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt) const;

    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl) const;

    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);

};

class OCP_AIMc : public OCP_FIM
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
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);

    /// Finish a Newton-Raphson iteration.
    bool FinishNR(Reservoir& rs, OCPControl& ctrl);
};

/// perform AIM in space, that is, some grids will be implicit, others will be explicit at the same time step
class OCP_AIMs
{
public:
    /// Setup AIMs
    void Setup(Reservoir& rs, LinearSystem& myLS, const OCPControl& ctrl);
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& myLS, const Reservoir& rs, const OCP_DBL& dt);
    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& myLS, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    bool FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

private:
    /// Resiual for AIMs
    ResFIM resFIM;

};

/// perform AIM in time, that is, local FIM will be performed after global IMPEC performs
class OCP_AIMt
{
public:
    /// Setup AIMt
    void Setup(Reservoir& rs, LinearSystem& myLS, LinearSystem& myAuxLS, const OCPControl& ctrl);
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCP_DBL& dt);
    /// Update properties of fluids.
    bool UpdateProperty(Reservoir& rs, OCPControl& ctrl, LinearSystem& myAuxLS);

private:
    /// Resiual for FIM
    ResFIM resFIM;

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