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

class IsothermalMethod
{
public:
    void InitRock(Bulk& bk) const;
    void CalRock(Bulk& bk) const;
};

/// IsoT_IMPEC is IMPEC (implicit pressure explict saturation) method.
class IsoT_IMPEC : virtual public IsothermalMethod
{
public:
    /// Setup IMPEC
    void Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, OCPControl& ctrl);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Determine if NR iteration finishes.
    OCP_BOOL FinishNR(const Reservoir& rs);
    void     FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Perform Flash with Sj and calculate values needed for FIM
    void InitFlash(Bulk& bk) const;
    /// Calculate relative permeability and capillary pressure needed for FIM
    void CalKrPc(Bulk& bk) const;
    /// Pass value needed for FIM from flash to bulk
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;

private:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Allocate memory for linear system
    void
    AllocateLinearSystem(LinearSystem& ls, const Reservoir& rs, const OCPControl& ctrl);
    /// Perform Flash with Ni and calculate values needed for FIM
    void CalFlash(Bulk& bk);
    /// Calculate flux between bulks and wells
    void CalFlux(Reservoir& rs) const;
    /// Calculate flux between bulks
    void CalBulkFlux(Reservoir& rs) const;
    /// Update mole composition of each bulk according to mass conservation for IMPEC
    void MassConserve(Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for bulks
    void
    AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for wells
    void
    AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for injection wells
    void AssembleMatInjWells(LinearSystem&  ls,
                             const Bulk&    bk,
                             const Well&    wl,
                             const OCP_DBL& dt) const;
    /// Assemble linear system for production wells
    void AssembleMatProdWells(LinearSystem&  ls,
                              const Bulk&    bk,
                              const Well&    wl,
                              const OCP_DBL& dt) const;
    /// Update P, Pj, BHP after linear system is solved
    void GetSolution(Reservoir& rs, const vector<OCP_DBL>& u);
    /// Reset variables to last time step
    void ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl);
    void ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl);
    void ResetToLastTimeStep03(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for FIM.
    void UpdateLastTimeStep(Reservoir& rs) const;
};

/// IsoT_FIM is FIM (Fully Implicit Method).
class IsoT_FIM : virtual public IsothermalMethod
{
public:
    /// Setup FIM
    void Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl) const;
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Allocate memory for linear system
    void
    AllocateLinearSystem(LinearSystem& ls, const Reservoir& rs, const OCPControl& ctrl);
    /// Pass value needed for FIM from flash to bulk
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;
    /// Calculate relative permeability and capillary pressure needed for FIM
    void CalKrPc(Bulk& bk) const;
    /// Calculate residual
    void CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& resetRes0) const;
    /// Assemble linear system for wells
    void
    AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void
    AssembleMatWellsNew(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for FIM.
    void UpdateLastTimeStep(Reservoir& rs) const;

private:
    /// Perform Flash with Sj and calculate values needed for FIM
    void InitFlash(Bulk& bk) const;
    /// Perform Flash with Ni and calculate values needed for FIM
    void CalFlash(Bulk& bk);
    /// Assemble linear system for bulks
    void
    AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for injection wells
    void AssembleMatInjWells(LinearSystem&  ls,
                             const Bulk&    bk,
                             const Well&    wl,
                             const OCP_DBL& dt) const;
    /// Assemble linear system for production wells
    void AssembleMatProdWells(LinearSystem&  ls,
                              const Bulk&    bk,
                              const Well&    wl,
                              const OCP_DBL& dt) const;
    /// Assemble linear system for bulks
    void
    AssembleMatBulksNew(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void AssembleMatBulksNewS(LinearSystem&    ls,
                              const Reservoir& rs,
                              const OCP_DBL&   dt) const;
    /// Assemble linear system for injection wells
    void AssembleMatInjWellsNew(LinearSystem&  ls,
                                const Bulk&    bk,
                                const Well&    wl,
                                const OCP_DBL& dt) const;
    /// Assemble linear system for production wells
    void AssembleMatProdWellsNew(LinearSystem&  ls,
                                 const Bulk&    bk,
                                 const Well&    wl,
                                 const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void
    GetSolution(Reservoir& rs, const vector<OCP_DBL>& u, const OCPControl& ctrl) const;
};

class IsoT_FIMn : protected IsoT_FIM
{
public:
    /// Setup FIM
    void Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl) const;
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl) const;

protected:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Perform Flash with Sj and calculate values needed for FIMn
    void InitFlash(Bulk& bk) const;
    /// Perform Flash with Ni and calculate values needed for FIMn
    void CalFlash(Bulk& bk);
    /// Pass value needed for FIMn from flash to bulk
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;
    /// Assemble linear system for bulks
    void
    AssembleMatBulksNew(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for wells
    void
    AssembleMatWellsNew(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Assemble linear system for injection wells
    void AssembleMatInjWellsNew(LinearSystem&  ls,
                                const Bulk&    bk,
                                const Well&    wl,
                                const OCP_DBL& dt) const;
    /// Assemble linear system for production wells
    void AssembleMatProdWellsNew(LinearSystem&  ls,
                                 const Bulk&    bk,
                                 const Well&    wl,
                                 const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void
    GetSolution(Reservoir& rs, const vector<OCP_DBL>& u, const OCPControl& ctrl) const;
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for FIMn
    void UpdateLastTimeStep(Reservoir& rs) const;
};

class IsoT_AIMc : protected IsoT_IMPEC, protected IsoT_FIM
{
public:
    /// Setup AIMc
    void Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    /// Init
    void InitReservoir(Reservoir& rs) const;
    /// Prepare for Assembling matrix.
    void Prepare(Reservoir& rs, const OCP_DBL& dt);
    /// Assemble Matrix
    void AssembleMat(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Solve the linear system.
    void SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl);
    /// Update properties of fluids.
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    /// Finish a Newton-Raphson iteration.
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    /// Finish a time step.
    void FinishStep(Reservoir& rs, OCPControl& ctrl) const;

protected:
    /// Allocate memory for reservoir
    void AllocateReservoir(Reservoir& rs);
    /// Determine which bulk are treated Implicit
    void SetFIMBulk(Reservoir& rs);
    /// Perform flash calculation with Ni for Explicit bulk -- Update partial properties
    void CalFlashEp(Bulk& bk);
    /// Perform flash calculation with Ni for Explicit bulk -- Update all properties
    void CalFlashEa(Bulk& bk);
    /// Perform flash calculation with Ni for Implicit bulk
    void CalFlashI(Bulk& bk);
    /// Pass flash value needed for Explicit bulk -- Update partial properties
    void PassFlashValueEp(Bulk& bk, const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure for Explicit bulk
    void CalKrPcE(Bulk& bk);
    /// Calculate relative permeability and capillary pressure for Implicit bulk
    void CalKrPcI(Bulk& bk);
    /// Assemble linear system for bulks
    void
    AssembleMatBulks(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    /// Update P, Ni, BHP after linear system is solved
    void
    GetSolution(Reservoir& rs, const vector<OCP_DBL>& u, const OCPControl& ctrl) const;
    /// Reset variables to last time step
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
    /// Update values of last step for AIMc.
    void UpdateLastTimeStep(Reservoir& rs) const;
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