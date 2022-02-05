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
    wellgroup.InputParam(param.paramWell);
}

void Reservoir::Setup()
{
    OCP_FUNCNAME;

    grid.Setup();
    bulk.Setup(grid);
    conn.Setup(grid, bulk);
    wellgroup.Setup(grid, bulk);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    wellgroup.ApplyControl(i);
}

void Reservoir::PrepareWell()
{
    OCP_FUNCNAME;

    wellgroup.PrepareWell(bulk);
}

void Reservoir::CalWellFlux()
{
    OCP_FUNCNAME;

    wellgroup.CalFlux(bulk);
}

void Reservoir::CalWellTrans()
{
    OCP_FUNCNAME;

    wellgroup.CalTrans(bulk);
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
    wellgroup.CalMaxBHPChange();
}

void Reservoir::CalIPRT(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    wellgroup.CalIPRT(bulk, dt);
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
        flag         = wellgroup.CheckP(bulk);
        return flag;
    }

    return 0;
}

bool Reservoir::CheckNi() const
{
    OCP_FUNCNAME;

    return bulk.CheckNi();
}

bool Reservoir::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    return bulk.CheckVe(Vlim);
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
        bulk.InitSjPcComp(50);

    bulk.CalVpore();
    bulk.InitFlash(true);
    bulk.CalKrPc();
    bulk.UpdateLastStepIMPEC();
    conn.CalFluxIMPEC(bulk);
    conn.UpdateLastStep();
    wellgroup.InitBHP(bulk);
}

OCP_DBL Reservoir::CalCFLIMPEC(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    OCP_DBL cflB = conn.CalCFLIMPEC(bulk, dt);
    OCP_DBL cflW = wellgroup.CalCFLIMPEC(bulk, dt);
    cfl          = max(cflB, cflW);

    return cfl;
}

OCP_DBL Reservoir::CalCFL01IMPEC(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    bulk.InitCFLIMPEC();
    conn.CalCFL01IMPEC(bulk, dt);
    wellgroup.CalCFL01IMPEC(bulk, dt);
    cfl = bulk.CalCFL01IMPEC();

    return cfl;
}

void Reservoir::CalFLuxIMPEC()
{
    OCP_FUNCNAME;

    conn.CalFluxIMPEC(bulk);
    wellgroup.CalFlux(bulk);
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
    wellgroup.MassConserveIMPEC(bulk, dt);
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
    wellgroup.UpdateLastDg();
}

void Reservoir::AllocateMatIMPEC(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(bulk.GetBulkNum() + wellgroup.GetWellNum(), 1);
    conn.AllocateMat(myLS);
    wellgroup.AllocateMat(myLS, bulk.GetBulkNum());
    myLS.AllocateColMem();
}

void Reservoir::AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMatIMPEC(myLS, bulk, dt);
    wellgroup.AssemblaMatIMPEC(myLS, bulk, dt);
}

void Reservoir::GetSolutionIMPEC(const vector<OCP_DBL>& u)
{
    OCP_FUNCNAME;

    bulk.GetSolIMPEC(u);
    wellgroup.GetSolIMPEC(u, bulk.GetBulkNum());
}

void Reservoir::ResetWellIMPEC()
{
    // wellgroup.ResetDg();
    wellgroup.ResetBHP();
    wellgroup.CalTrans(bulk);
    wellgroup.CalFlux(bulk);
    wellgroup.CaldG(bulk);
}

void Reservoir::ResetVal00IMPEC()
{
    OCP_FUNCNAME;

    bulk.ResetP();
    wellgroup.CalTrans(bulk);
}

void Reservoir::ResetVal01IMPEC()
{
    OCP_FUNCNAME;
    bulk.ResetP();
    bulk.ResetPj();
    conn.Reset();
}

void Reservoir::ResetVal02IMPEC()
{
    OCP_FUNCNAME;

    bulk.ResetP();
    bulk.ResetPj();
    bulk.ResetNi();
    conn.Reset();
}

void Reservoir::ResetVal03IMPEC()
{
    OCP_FUNCNAME;

    bulk.ResetP();
    bulk.ResetPj();
    bulk.ResetNi();
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
        bulk.InitSjPcComp(50);

    bulk.CalVpore();
    bulk.InitFlash();
    bulk.FlashDeriv();
    bulk.CalKrPcDeriv();
    conn.CalFluxFIM(bulk);
    wellgroup.InitBHP(bulk);
    UpdateLastStepFIM();
}

void Reservoir::CalFlashDerivFIM()
{
    OCP_FUNCNAME;

    bulk.FlashDeriv();
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
    conn.UpdateLastUpblockFIM();
    wellgroup.UpdateLastBHP();
    wellgroup.UpdateLastDg();
}

void Reservoir::AllocateMatFIM(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(bulk.GetBulkNum() + wellgroup.GetWellNum(),
                        bulk.GetComNum() + 1);
    conn.AllocateMat(myLS);
    wellgroup.AllocateMat(myLS, bulk.GetBulkNum());
    myLS.AllocateColMem();
}

void Reservoir::AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_FIM(myLS, bulk, dt);
    wellgroup.AssemblaMatFIM(myLS, bulk, dt);
}

void Reservoir::GetSolutionFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
                               const OCP_DBL& dSmax)
{
    OCP_FUNCNAME;

    bulk.GetSolFIM(u, dPmax, dSmax);
    wellgroup.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}

// Not useful
void Reservoir::GetSolution01FIM(const vector<OCP_DBL>& u)
{
    OCP_DBL alpha = bulk.GetSol01FIM(u);
    wellgroup.GetSol01FIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1, alpha);

    cout << alpha << endl;
}

void Reservoir::CalResFIM(ResFIM& resFIM, const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResFIM(resFIM.res, bulk, dt);
    // Well to Bulk
    wellgroup.CalResFIM(resFIM, bulk, dt);
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
    // for (USI w = 0; w < wellgroup.numWell; w++) {
    //    cout << wellgroup.GetWBHP(w) << "   ";
    //    for (USI p = 0; p < wellgroup.GetWellPerfNum(w); p++) {
    //        cout << wellgroup.wellGroup[w].GetPerfPre(p) << "   ";
    //    }
    //}
    // cout << endl;
}

void Reservoir::ResetFIM(const bool& flag)
{
    bulk.ResetFIM();
    conn.ResetUpblockFIM();
    wellgroup.ResetBHP();
    wellgroup.CalTrans(bulk);
    wellgroup.CalFlux(bulk);

    if (flag) {
        wellgroup.CaldG(bulk);
        wellgroup.CalFlux(bulk);
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
    for (USI w = 0; w < wellgroup.numWell; w++) {
        outu << wellgroup.GetWBHP(w) << "\n";
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