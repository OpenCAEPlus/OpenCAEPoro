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

void Reservoir::InputParam(ParamRead& param)
{
    grid.InputParam(param.paramRs);
    bulk.InputParam(param.paramRs);
    wellgroup.InputParam(param.paramWell);
}

void Reservoir::Setup()
{
    grid.Setup();
    bulk.Setup(grid);
    conn.Setup(grid, bulk);
    wellgroup.Setup(grid, bulk);
}

void Reservoir::Init()
{
    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBlk(50);
    else if (bulk.GetMixMode() == EoS_PVTW)
        bulk.InitSjPcComp(50);

    bulk.CalVporo();
    bulk.FlashSj();
    bulk.CalKrPc();
    bulk.SetLastStep();
    conn.CalFlux(bulk);
    wellgroup.Init(bulk);
}

OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt)
{
    OCP_DBL cflB = conn.CalCFL(bulk, dt);
    OCP_DBL cflW = wellgroup.CalCFL(bulk, dt);

    cfl = max(cflB, cflW);

    return cfl;
}

// assemble mat
void Reservoir::AssembleMat(Solver<OCP_DBL>& mysolver, const OCP_DBL& dt) const
{
    conn.InitAssembleMat(mysolver);
    conn.AssembleMat_IMPES(mysolver, bulk, dt);
    wellgroup.AssemblaMat_WB_IMPES(mysolver, bulk, dt);
}

void Reservoir::GetSolution_IMPES(const vector<OCP_DBL>& u)
{
    bulk.GetSolIMPES(u);
    wellgroup.GetSol_IMPES(u, bulk.GetBulkNum());
}

OCP_INT Reservoir::CheckP()
{
    if (!bulk.CheckP()) return 1;

    OCP_INT flag = 0;
    flag         = wellgroup.CheckP(bulk);
    return flag;
}

void Reservoir::ResetVal01()
{
    bulk.ResetP();
    bulk.ResetPj();
    conn.CalFlux(bulk);
}

void Reservoir::ResetVal02()
{
    bulk.ResetP();
    bulk.ResetPj();
    conn.CalFlux(bulk);

    bulk.ResetNi();
    bulk.FlashNi();

    bulk.ResetVp();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/