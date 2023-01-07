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

    grid.InputParam(param.paramRs);
    bulk.InputParam(param.paramRs);
    allWells.InputParam(param.paramWell);
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

OCP_INT Reservoir::CheckP(const OCP_BOOL& bulkCheck, const OCP_BOOL& wellCheck)
{
    OCP_FUNCNAME;

    if (bulkCheck) {
        if (!bulk.CheckP()) {
            // negative Pressure
            return 1;
        }
    }

    if (wellCheck) {
        OCP_INT flag = 0;
        flag         = allWells.CheckP(bulk);
        return flag;
    }

    return 0;
}

OCP_BOOL Reservoir::CheckNi()
{
    OCP_FUNCNAME;

    return bulk.CheckNi();
}

OCP_BOOL Reservoir::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    return bulk.CheckVe(Vlim);
}


OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    conn.CalCFL(bulk, dt);
    allWells.CalCFL(bulk, dt);
    bulk.CalCFL();

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