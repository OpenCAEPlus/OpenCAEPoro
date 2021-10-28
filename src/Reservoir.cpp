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

void Reservoir::InitIMPEC()
{
    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBlk(50);
    else if (bulk.GetMixMode() == EoS_PVTW)
        bulk.InitSjPcComp(50);

    bulk.CalVpore();
    bulk.FlashSj();
    bulk.CalKrPc();
    bulk.UpdateLastStep();
    conn.CalFlux(bulk);
    conn.UpdateLastStep();
    wellgroup.Init(bulk);
}

void Reservoir::InitFIM()
{
    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBlk(50);
    else if (bulk.GetMixMode() == EoS_PVTW)
        bulk.InitSjPcComp(50);

    bulk.CalVpore();
    bulk.FlashSj();
    bulk.FlashNiDeriv();
    bulk.CalKrPcDeriv();
    bulk.UpdateLastStep();
    conn.CalFlux(bulk);
    conn.UpdateLastStep();
    wellgroup.Init(bulk);
}


void Reservoir::Prepare(OCP_DBL& dt) {
    wellgroup.PrepareWell(bulk);
    OCP_DBL cfl = CalCFL01(dt);
    if (cfl > 1) dt /= (cfl + 1);
}

OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt)
{
    OCP_DBL cflB = conn.CalCFL(bulk, dt);
    OCP_DBL cflW = wellgroup.CalCFL(bulk, dt);


    // cout << dt << "Days\t" << cflB << "\t" << cflW << endl;

    cfl = max(cflB, cflW);

    return cfl;
}

OCP_DBL Reservoir::CalCFL01(const OCP_DBL& dt) 
{
    bulk.InitCFL();
    
    conn.CalCFL01(bulk, dt);

    //if (cfl > 1) {
    //    cout << "bkCFL:\t" << bulk.cfl[2686] << "\t";
    //}
    wellgroup.CalCFL01(bulk, dt);

    //if (cfl > 1) {
    //    cout << "wellCFL:\t" << bulk.cfl[2686] << endl;
    //}
    cfl = bulk.CalCFL(false);

    return cfl;
}


// allocate memory
void Reservoir::AllocateMatIMPEC(LinearSolver& mySolver) const
{
    mySolver.AllocateRowMem(bulk.GetBulkNum() + wellgroup.GetWellNum(), 1);
    conn.AllocateMat(mySolver);
    wellgroup.AllocateMat(mySolver, bulk.GetBulkNum());
    mySolver.AllocateColMem();
    mySolver.AllocateFasp();
}

void Reservoir::AllocateMatFIM(LinearSolver& mySolver) const
{
    mySolver.AllocateRowMem(bulk.GetBulkNum() + wellgroup.GetWellNum(), bulk.GetComNum() + 1);
    conn.AllocateMat(mySolver);
    wellgroup.AllocateMat(mySolver, bulk.GetBulkNum());
    mySolver.AllocateColMem();
    mySolver.AllocateBFasp();
}


// assemble mat
void Reservoir::AssembleMatIMPEC(LinearSolver& mysolver, const OCP_DBL& dt) const
{
    conn.InitAssembleMat(mysolver);
    conn.AssembleMat_IMPEC(mysolver, bulk, dt);
    wellgroup.AssemblaMat_WB_IMPEC(mysolver, bulk, dt);
}

void Reservoir::AssembleMatFIM(LinearSolver& mysolver, const OCP_DBL& dt) const
{
    conn.InitAssembleMat(mysolver);
    conn.AssembleMat_FIM(mysolver, bulk, dt);
    wellgroup.AssemblaMat_WB_FIM(mysolver, bulk, dt);
}

void Reservoir::GetSolution_IMPEC(const vector<OCP_DBL>& u)
{
    bulk.GetSolIMPEC(u);
    wellgroup.GetSol_IMPEC(u, bulk.GetBulkNum());
}


void Reservoir::GetSolution_FIM(const vector<OCP_DBL>& u)
{
    bulk.GetSolFIM(u);
    wellgroup.GetSol_FIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}

void Reservoir::CalResFIM(vector<OCP_DBL>& res, const OCP_DBL& dt)
{
    res.assign(res.size(), 0);
    conn.CalResFIM(res, bulk, dt);
    wellgroup.CalResFIM(res, bulk, dt);
    Dscalar(res.size(), -1, res.data());
}


OCP_INT Reservoir::CheckP()
{
    if (!bulk.CheckP()) return 1;

    OCP_INT flag = 0;
    flag         = wellgroup.CheckP(bulk);
    return flag;
}

void Reservoir::ResetVal()
{
    bulk.ResetP();
    bulk.ResetPj();
    conn.Reset();
}

void Reservoir::ResetVal01()
{
    bulk.ResetP();
    bulk.ResetPj();
    bulk.ResetNi();
    conn.Reset();
}

void Reservoir::ResetVal02()
{
    bulk.ResetP();
    bulk.ResetPj();
    bulk.ResetNi();
    bulk.ResetFlash();
    bulk.ResetVp();
    conn.Reset();

    // bulk.FlashNi();
    // bulk.CheckDiff();

    // conn.CalFlux(bulk);
    // conn.CheckDiff();
    // Becareful! if recalculate the flash, result may be different because the initial
    // flash was calculated by FlashSj not FlashNi.
     
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/