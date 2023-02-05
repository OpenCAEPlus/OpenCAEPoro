/*! \file    ThermalMethod.hpp
 *  \brief   ThermalMethod class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __THERMALMETHOD_HEADER__
#define __THERMALMETHOD_HEADER__

#include "LinearSystem.hpp"
#include "OCPControl.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"

class T_FIM
{
public:
    void     Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl);
    void     InitReservoir(Reservoir& rs) const;
    void     Prepare(Reservoir& rs, const OCPControl& ctrl);
    void     AssembleMat(LinearSystem&    ls,
                         const Reservoir& rs,
                         const OCP_DBL&   t,
                         const OCP_DBL&   dt) const;
    void     SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl) const;
    OCP_BOOL UpdateProperty(Reservoir& rs, OCPControl& ctrl);
    OCP_BOOL FinishNR(Reservoir& rs, OCPControl& ctrl);
    void     FinishStep(Reservoir& rs, OCPControl& ctrl);

protected:
    void AllocateReservoir(Reservoir& rs);
    void
    AllocateLinearSystem(LinearSystem& ls, const Reservoir& rs, const OCPControl& ctrl);
    void InitRock(Bulk& bk) const;
    void CalRock(Bulk& bk) const;
    void InitFlash(Bulk& bk) const;
    void CalFlash(Bulk& bk) const;
    void PassFlashValue(Bulk& bk, const OCP_USI& n) const;
    void CalKrPc(Bulk& bk) const;
    void CalThermalConduct(BulkConn& conn, Bulk& bk) const;
    void CalHeatLoss(Bulk& bk, const OCP_DBL& t, const OCP_DBL& dt) const;
    void UpdateLastTimeStep(Reservoir& rs) const;
    void CalRes(Reservoir&      rs,
                const OCP_DBL&  t,
                const OCP_DBL&  dt,
                const OCP_BOOL& resetRes0);
    void AssembleMatBulks(LinearSystem&    ls,
                          const Reservoir& rs,
                          const OCP_DBL&   t,
                          const OCP_DBL&   dt) const;
    void
    AssembleMatWells(LinearSystem& ls, const Reservoir& rs, const OCP_DBL& dt) const;
    void AssembleMatInjWells(LinearSystem&  ls,
                             const Bulk&    bk,
                             const Well&    wl,
                             const OCP_DBL& dt) const;
    void AssembleMatProdWells(LinearSystem&  ls,
                              const Bulk&    bk,
                              const Well&    wl,
                              const OCP_DBL& dt) const;
    void
    GetSolution(Reservoir& rs, const vector<OCP_DBL>& u, const OCPControl& ctrl) const;
    void ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl);
};

#endif /* end if __THERMALMETHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/