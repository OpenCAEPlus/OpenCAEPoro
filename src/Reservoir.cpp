/*! \file    Reservoir.cpp
 *  \brief   Reservoir class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Reservoir.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////

void Reservoir::InputParam(ParamRead& param)
{
    OCP_FUNCNAME;

    grid.InputParam(param.paramRs, param.paramOutput);
    bulk.InputParam(param.paramRs);
    allWells.InputParam(param.paramWell, param.paramOutput);
    optFeatures.InputParam(param.paramRs);
}

void Reservoir::SetupIsoT()
{
    OCP_FUNCNAME;

    grid.SetupIsoT();
    bulk.SetupIsoT(grid);
    conn.SetupIsoT(grid, bulk);
    allWells.Setup(grid, bulk);

    bulk.SetupOptionalFeatures(grid, optFeatures);
}

void Reservoir::SetupT()
{
    grid.SetupT();
    bulk.SetupT(grid);
    conn.SetupIsoT(grid, bulk);
    allWells.Setup(grid, bulk);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    allWells.ApplyControl(i);
    allWells.SetupWellGroup(bulk);
}

void Reservoir::CalMaxChange()
{
    OCP_FUNCNAME;

    bulk.CalMaxChange();
    allWells.CalMaxBHPChange();
}

void Reservoir::CalIPRT(const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Calculate injection / production rate for current step
    allWells.CalIPRT(bulk, dt);
    // Calculate Reinjection fluid for next step
    allWells.CalReInjFluid(bulk);
}

OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt) const
{
    fill(bulk.cfl.begin(), bulk.cfl.end(), 0.0);
    const USI np = bulk.numPhase;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        for (USI j = 0; j < np; j++) {
            const OCP_USI uId = conn.upblock[c * np + j];

            if (bulk.phaseExist[uId * np + j]) {
                bulk.cfl[uId * np + j] += fabs(conn.upblock_Velocity[c * np + j]) * dt;
            }
        }
    }

    for (const auto& wl : allWells.wells) {
        if (wl.IsOpen() && wl.WellType() == PROD) {
            for (USI p = 0; p < wl.PerfNum(); p++) {
                if (wl.PerfState(p) == OPEN) {
                    const OCP_USI k = wl.PerfLocation(p);

                    for (USI j = 0; j < np; j++) {
                        bulk.cfl[k * np + j] += fabs(wl.PerfProdQj_ft3(p, j)) * dt;
                    }
                }
            }
        }
    }

    bulk.maxCFL       = 0;
    const OCP_USI len = bulk.numBulk * np;
    for (OCP_USI n = 0; n < len; n++) {
        if (bulk.phaseExist[n]) {
            bulk.cfl[n] /= bulk.vj[n];
#ifdef DEBUG
            if (!isfinite(bulk.cfl[n])) {
                OCP_ABORT("cfl is nan!");
            }
#endif // DEBUG
            if (bulk.maxCFL < bulk.cfl[n]) bulk.maxCFL = bulk.cfl[n];
        }
    }

    return bulk.maxCFL;
}

void Reservoir::PrintSolFIM(const string& outfile) const
{
    ofstream outu(outfile);
    if (!outu.is_open()) cout << "Can not open " << outfile << endl;
    const OCP_USI nb = bulk.numBulk;
    const OCP_USI nc = bulk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        // Pressure
        outu << bulk.P[n] << "\n";
        // Ni
        for (USI i = 0; i < nc; i++) {
            outu << bulk.Ni[n * nc + i] << "\n";
        }
    }
    // Well Pressure
    for (USI w = 0; w < allWells.numWell; w++) {
        outu << allWells.GetWBHP(w) << "\n";
    }
    outu.close();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/