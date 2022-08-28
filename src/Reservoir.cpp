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

}

void Reservoir::Setup()
{
    OCP_FUNCNAME;

    grid.Setup();
    bulk.Setup(grid);
    conn.Setup(grid, bulk);
    allWells.Setup(grid, bulk);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    allWells.ApplyControl(i);
    allWells.SetupWellGroup(bulk);
}

void Reservoir::PrepareWell()
{
    OCP_FUNCNAME;

    allWells.PrepareWell(bulk);
}

void Reservoir::CalWellFlux()
{
    OCP_FUNCNAME;

    allWells.CalFlux(bulk);
}

void Reservoir::CalWellTrans()
{
    OCP_FUNCNAME;

    allWells.CalTrans(bulk);
}

void Reservoir::CalVpore()
{
    OCP_FUNCNAME;

    bulk.CalVpore();
}

void Reservoir::CalKrPc()
{
    OCP_FUNCNAME;

    bulk.CalKrPc();
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

OCP_INT Reservoir::CheckP(const bool& bulkCheck, const bool& wellCheck)
{
    OCP_FUNCNAME;

    if (bulkCheck) {
        if (!bulk.CheckP()) {
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

bool Reservoir::CheckNi()
{
    OCP_FUNCNAME;

    return bulk.CheckNi();
}

bool Reservoir::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    return bulk.CheckVe(Vlim);
}

void Reservoir::GetNTQT(const OCP_DBL& dt) 
{
    OCP_DBL NT = bulk.CalNT();
    OCP_DBL QT = allWells.CalWellQT() * dt;
    cout << setprecision(8) << NT << "   " << QT << endl;
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateAuxIMPEC()
{
    OCP_FUNCNAME;

    bulk.AllocateAuxIMPEC();
    conn.AllocateAuxIMPEC(bulk.GetPhaseNum());
}

void Reservoir::InitIMPEC()
{
    OCP_FUNCNAME;

    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBo(50);
    else if (bulk.GetMixMode() == EOS_PVTW)
        bulk.InitSjPcComp(50, grid);

    bulk.CalVpore();
    bulk.InitFlash(true);
    bulk.CalKrPc();
    bulk.UpdateLastStepIMPEC();
    conn.CalFluxIMPEC(bulk);
    conn.UpdateLastStep();
    allWells.InitBHP(bulk);
    allWells.UpdateLastBHP();
}


OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    bulk.SetCFL2Zero();
    conn.CalCFL(bulk, dt);
    allWells.CalCFL(bulk, dt);
    cfl = bulk.CalCFL();

    return cfl;
}

void Reservoir::CalFLuxIMPEC()
{
    OCP_FUNCNAME;

    conn.CalFluxIMPEC(bulk);
    allWells.CalFlux(bulk);
}

void Reservoir::CalConnFluxIMPEC()
{
    OCP_FUNCNAME;

    conn.CalFluxIMPEC(bulk);
}

void Reservoir::MassConseveIMPEC(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    conn.MassConserveIMPEC(bulk, dt);
    allWells.MassConserveIMPEC(bulk, dt);
}

void Reservoir::CalFlashIMPEC()
{
    OCP_FUNCNAME;

    bulk.Flash();
}

void Reservoir::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;
    bulk.UpdateLastStepIMPEC();
    conn.UpdateLastStep();
    // useless in IMPEC now
    // allWells.UpdateLastBHP();
    // allWells.UpdateLastDg();
}

void Reservoir::AllocateMatIMPEC(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(bulk.GetBulkNum() + allWells.GetWellNum(), 1);
    conn.AllocateMat(myLS);
    allWells.AllocateMat(myLS, bulk.GetBulkNum());
    myLS.AllocateColMem();
}

void Reservoir::AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMatIMPEC(myLS, bulk, dt);
    allWells.AssemblaMatIMPEC(myLS, bulk, dt);
}

void Reservoir::GetSolutionIMPEC(const vector<OCP_DBL>& u)
{
    OCP_FUNCNAME;

    bulk.GetSolIMPEC(u);
    allWells.GetSolIMPEC(u, bulk.GetBulkNum());
}

void Reservoir::ResetWellIMPEC()
{
    // allWells.ResetDg();
    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CalFlux(bulk);
    allWells.CaldG(bulk);
}


void Reservoir::ResetVal01IMPEC()
{
    OCP_FUNCNAME;
    bulk.ResetPj();
    conn.Reset();
}

void Reservoir::ResetVal02IMPEC()
{
    OCP_FUNCNAME;
    
    bulk.ResetPj();
    bulk.ResetNi();
    conn.Reset();
}

void Reservoir::ResetVal03IMPEC()
{
    OCP_FUNCNAME;
    bulk.ResetphaseNum();
    bulk.ResetminEigenSkip();
    bulk.ResetflagSkip();
    bulk.ResetziSkip();
    bulk.ResetPSkip();
    bulk.ResetKs();

    bulk.ResetPj();
    bulk.ResetNi();
    bulk.ResetNt();
    bulk.ResetFlash();
    bulk.ResetVp();
    conn.Reset();

    // Becareful! if recalculate the flash, result may be different because the initial
    // flash was calculated by InitFlash not Flash.
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateAuxFIM()
{
    OCP_FUNCNAME;

    bulk.AllocateAuxFIM();
    conn.AllocateAuxFIM(bulk.GetPhaseNum());
}

void Reservoir::InitFIM()
{
    OCP_FUNCNAME;

    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBo(50);
    else if (bulk.GetMixMode() == EOS_PVTW)
        bulk.InitSjPcComp(50, grid);

    bulk.CalVpore();
    bulk.InitFlash(false);
    bulk.FlashDeriv();
    bulk.CalKrPcDeriv();
    conn.CalFluxFIM(bulk);
    allWells.InitBHP(bulk);
    UpdateLastStepFIM();
}

void Reservoir::CalFlashDerivFIM()
{
    OCP_FUNCNAME;

    bulk.FlashDeriv();
}

void Reservoir::CalFlashDerivFIM_n()
{
    OCP_FUNCNAME;

    bulk.FlashDeriv_n();
}

void Reservoir::CalKrPcDerivFIM()
{
    OCP_FUNCNAME;

    bulk.CalKrPcDeriv();
}

void Reservoir::UpdateLastStepFIM()
{
    OCP_FUNCNAME;

    bulk.UpdateLastStepFIM();
    allWells.UpdateLastBHP();
}

void Reservoir::AllocateMatFIM(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(bulk.GetBulkNum() + allWells.GetWellNum(),
                        bulk.GetComNum() + 1);
    conn.AllocateMat(myLS);
    allWells.AllocateMat(myLS, bulk.GetBulkNum());
    myLS.AllocateColMem();
}

void Reservoir::AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);

#ifdef OCP_NEW_FIM
    conn.AssembleMat_FIM_new(myLS, bulk, dt);
    allWells.AssemblaMatFIM_new(myLS, bulk, dt);
#else
    conn.AssembleMat_FIM(myLS, bulk, dt);
    allWells.AssemblaMatFIM(myLS, bulk, dt);
#endif // OCP_NEW_FIM
}


void Reservoir::AssembleMatFIM_n(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_FIM_new_n(myLS, bulk, dt);
    allWells.AssemblaMatFIM_new_n(myLS, bulk, dt);
}

void Reservoir::GetSolutionFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
                               const OCP_DBL& dSmax)
{
    bulk.GetSolFIM(u, dPmax, dSmax); 
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}


void Reservoir::GetSolutionFIM_n(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    OCP_FUNCNAME;

    bulk.GetSolFIM_n(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}


// Not useful
void Reservoir::GetSolution01FIM(const vector<OCP_DBL>& u)
{
    bulk.GetSol01FIM(u);
    allWells.GetSol01FIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1, 1);
}

void Reservoir::CalResFIM(ResFIM& resFIM, const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResFIM(resFIM.res, bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resFIM, bulk, dt);
    // Calculate RelRes
    bulk.CalRelResFIM(resFIM);
    Dscalar(resFIM.res.size(), -1, resFIM.res.data());

    // Calculate Res2 and ResMax
    // OCP_DBL resmax = 0;
    // OCP_USI maxId = 0;
    // OCP_DBL res2 = 0;
    // for (OCP_USI i = 0; i < resFIM.res.size(); i++) {
    //    res2 += pow(resFIM.res[i], 2);
    //    if (resmax < fabs(resFIM.res[i])) {
    //        resmax = fabs(resFIM.res[i]);
    //        maxId = i;
    //    }
    //}
    // cout << "Res2  " << pow(res2 / resFIM.res.size(), 0.5) << "     "
    //    << resFIM.res.size() << endl;
    // cout << "ResMax   " << resmax << "   " << maxId << endl;
    // cout << "BHP   \n";
    // for (USI w = 0; w < allWells.numWell; w++) {
    //    cout << allWells.GetWBHP(w) << "   ";
    //    for (USI p = 0; p < allWells.GetWellPerfNum(w); p++) {
    //        cout << allWells.wells[w].GetPerfPre(p) << "   ";
    //    }
    //}
    // cout << endl;
}

void Reservoir::ResetFIM(const bool& flag)
{
    bulk.ResetFIM();
    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CalFlux(bulk);

    if (flag) {
        allWells.CaldG(bulk);
        allWells.CalFlux(bulk);
    }
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

void Reservoir::ShowRes(const vector<OCP_DBL>& res) const
{
    bulk.ShowRes(res);
    allWells.ShowRes(res, bulk);
}


/////////////////////////////////////////////////////////////////////
// AIMt
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateAuxAIMt()
{
    bulk.AllocateAuxIMPEC();
    conn.AllocateAuxIMPEC(bulk.GetPhaseNum());
    conn.AllocateAuxAIMt();
    
    bulk.AllocateWellBulkId(allWells.GetWellPerfNum() * 10);
    bulk.AllocateAuxAIM(0.05);
}

void Reservoir::AllocateMatAIMt(LinearSystem& myLS) const
{
    OCP_FUNCNAME;
    
    myLS.AllocateRowMem(bulk.GetMaxFIMBulk() + allWells.GetWellNum(),
        bulk.GetComNum() + 1);
    myLS.AllocateColMem(10);
}


void Reservoir::CalFlashDerivAIM(const bool& IfAIMs)
{
    bulk.FlashDerivAIM(IfAIMs);
}


void Reservoir::CalKrPcDerivAIM(const bool& IfAIMs)
{
    bulk.CalKrPcDerivAIM(IfAIMs);
}


void Reservoir::CalResAIMt(ResFIM& resFIM, const OCP_DBL& dt)
{
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResAIMt(resFIM.res, bulk, dt);
    // Well to Bulk
    allWells.CalResAIMt(resFIM, bulk, dt);
    // Calculate RelRes
    bulk.CalRelResAIMt(resFIM);
    Dscalar(resFIM.res.size(), -1, resFIM.res.data());
}


void Reservoir::AssembleMatAIMt(LinearSystem& myLS, const OCP_DBL& dt) const
{
    conn.SetupMatSparsityAIMt(myLS, bulk);
    conn.AssembleMat_AIMt(myLS, bulk, dt);
    allWells.AssemblaMatAIMt(myLS, bulk, dt);
}

void Reservoir::GetSolutionAIMt(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    bulk.GetSolAIMt(u, dPmax, dSmax);
    allWells.GetSolAIMt(u, bulk.numFIMBulk, bulk.GetComNum() + 1);
}


void Reservoir::AllocateAuxAIMs()
{
    bulk.AllocateAuxIMPEC();  
    bulk.AllocateWellBulkId(allWells.GetWellPerfNum() * 10);
    bulk.AllocateAuxAIM(1);

    conn.AllocateAuxIMPEC(bulk.GetPhaseNum());
    conn.AllocateAuxAIMt();
}

void Reservoir::CalResAIMs(ResFIM& resFIM, const OCP_DBL& dt)
{
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResAIMs(resFIM.res, bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resFIM, bulk, dt);
    // Calculate RelRes
    bulk.CalRelResAIMs(resFIM);
    Dscalar(resFIM.res.size(), -1, resFIM.res.data());
}

void Reservoir::AssembleMatAIMs(LinearSystem& myLS, vector<OCP_DBL>& res, const OCP_DBL& dt) const
{
    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_AIMs(myLS, res, bulk, dt);
    allWells.AssemblaMatAIMs(myLS, bulk, dt);
}

void Reservoir::GetSolutionAIMs(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    bulk.GetSolAIMs(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}

void Reservoir::ResetValAIM()
{
    ResetVal03IMPEC();    
    bulk.ResetFIMBulk();
    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CalFlux(bulk);

    if (false) {
        allWells.CaldG(bulk);
        allWells.CalFlux(bulk);
    }  
}

OCP_DBL Reservoir::CalCFLAIM(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    bulk.SetCFL2Zero();
    conn.CalCFL(bulk, dt);
    cfl = bulk.CalCFL();

    return cfl;
}

void Reservoir::UpdateLastStepAIM()
{
    OCP_FUNCNAME;
    bulk.UpdateLastStepIMPEC();
    bulk.UpdateLastStepAIM();
    conn.UpdateLastStep();
    allWells.UpdateLastBHP();
    allWells.UpdateLastDg();
}

void Reservoir::AllocateAuxAIMc()
{
    OCP_FUNCNAME;

    bulk.AllocateAuxFIM();
    bulk.AllocateAuxAIMc();
    conn.AllocateAuxAIMc(bulk.GetPhaseNum());
}

void Reservoir::AssembleMatAIMc(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_AIMc01(myLS, bulk, dt);
    allWells.AssemblaMatFIM(myLS, bulk, dt);
}

void Reservoir::CalResAIMc(ResFIM& resFIM, const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResAIMc(resFIM.res, bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resFIM, bulk, dt);
    // Calculate RelRes
    bulk.CalRelResFIM(resFIM);
    Dscalar(resFIM.res.size(), -1, resFIM.res.data());
}

void Reservoir::CalFlashAIMc()
{
    bulk.FlashAIMc();
}

void Reservoir::CalFlashAIMc01()
{
    bulk.FlashAIMc01();
}

void Reservoir::CalKrPcAIMc()
{
    bulk.CalKrPcAIMc();
}


/// Calculate Flash for local FIM, some derivatives are needed
void Reservoir::CalFlashDerivAIMc()
{
    bulk.FlashDerivAIMc();
}


/// Calculate Relative Permeability and Capillary and some derivatives for each Bulk
void Reservoir::CalKrPcDerivAIMc()
{
    bulk.CalKrPcDerivAIMc();
}

void Reservoir::GetSolutionAIMc(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    bulk.GetSolAIMc01(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}

void Reservoir::InitAIMc()
{
    OCP_FUNCNAME;

    if (bulk.GetMixMode() == BLKOIL)
        bulk.InitSjPcBo(50);
    else if (bulk.GetMixMode() == EOS_PVTW)
        bulk.InitSjPcComp(50, grid);

    bulk.CalVpore();
    bulk.InitFlash(true);
    bulk.CalKrPc();
    conn.CalFluxFIM(bulk);
    allWells.InitBHP(bulk);
    UpdateLastStepFIM();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/