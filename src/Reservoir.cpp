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

void Reservoir::SetupIsoT()
{
    OCP_FUNCNAME;

    grid.SetupIsoT();
    bulk.SetupIsoT(grid);
    conn.Setup(grid);
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

void Reservoir::CalRock()
{
    OCP_FUNCNAME;

    bulk.CalRock();
}

void Reservoir::CalKrPc()
{
    OCP_FUNCNAME;

    bulk.CalKrPcIMPEC();
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


/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateIMPEC_IsoT()
{
    OCP_FUNCNAME;

    bulk.AllocateIMPEC_IsoT();
    conn.AllocateIMPEC_IsoT(bulk.GetPhaseNum());
}

void Reservoir::InitIMPEC()
{
    OCP_FUNCNAME;

    bulk.InitRock();
    bulk.InitSjPc(50);
    bulk.CalRock();
    bulk.InitFlashIMPEC();
    bulk.CalKrPcIMPEC();
    bulk.UpdateLastStepIMPEC();

    conn.CalAkd(bulk);
    conn.CalFluxIMPEC(bulk);
    conn.UpdateLastStepIMPEC();

    allWells.InitBHP(bulk);
    allWells.UpdateLastBHP();
}

void Reservoir::AllocateMatIMPEC(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(GetBulkNum() + GetWellNum(), 1);  
    myLS.AllocateColMem(conn.GetNeighborNum(), allWells.GetWell2Bulk());
}

void Reservoir::AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.AssembleMatIMPEC(myLS, bulk, dt);
    allWells.AssemblaMatIMPEC(myLS, bulk, dt);
}

void Reservoir::GetSolutionIMPEC(const vector<OCP_DBL>& u)
{
    OCP_FUNCNAME;

    bulk.GetSolIMPEC(u);
    allWells.GetSolIMPEC(u, bulk.GetBulkNum());
}


void Reservoir::CalFLuxIMPEC()
{
    OCP_FUNCNAME;

    conn.CalFluxIMPEC(bulk);
    allWells.CalFlux(bulk);
}

void Reservoir::MassConserveIMPEC(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    conn.MassConserveIMPEC(bulk, dt);
    allWells.MassConserveIMPEC(bulk, dt);
}


void Reservoir::ResetVal01IMPEC()
{
    OCP_FUNCNAME;
    bulk.ResetVal01IMPEC();
    conn.ResetIMPEC();
}

void Reservoir::ResetVal02IMPEC()
{
    OCP_FUNCNAME;

    bulk.ResetVal02IMPEC();
    conn.ResetIMPEC();
}

void Reservoir::ResetVal03IMPEC()
{
    OCP_FUNCNAME;
    bulk.ResetVal03IMPEC();
    conn.ResetIMPEC();
}


void Reservoir::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;
    bulk.UpdateLastStepIMPEC();
    conn.UpdateLastStepIMPEC();
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateFIM_IsoT()
{
    OCP_FUNCNAME;

    bulk.AllocateFIM_IsoT();
    conn.AllocateFIM_IsoT(bulk.GetPhaseNum());
}

void Reservoir::InitFIM()
{
    OCP_FUNCNAME;

    bulk.InitRock();
    bulk.InitSjPc(50);
    bulk.CalRock();
    bulk.InitFlashFIM();
    bulk.CalKrPcFIM();

    conn.CalAkd(bulk);
 
    allWells.InitBHP(bulk);

    UpdateLastStepFIM();
}

void Reservoir::AllocateMatFIM_IsoT(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(GetBulkNum() + GetWellNum(), GetComNum() + 1);
    myLS.AllocateColMem(conn.GetNeighborNum(), allWells.GetWell2Bulk());
}


void Reservoir::AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

#ifdef OCP_OLD_FIM
    conn.AssembleMat_FIM(myLS, bulk, dt);
    allWells.AssemblaMatFIM(myLS, bulk, dt);
#else
    conn.AssembleMat_FIM_new(myLS, bulk, dt);
    allWells.AssemblaMatFIM_new(myLS, bulk, dt);
#endif // OCP_OLD_FIM
}


void Reservoir::GetSolutionFIM(const vector<OCP_DBL>& u,
    const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    bulk.GetSolFIM(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}


void Reservoir::CalResFIM(OCPRes& resFIM, const OCP_DBL& dt, const OCP_BOOL& resetRes0)
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
    Dscalar(resFIM.res.size(), -1.0, resFIM.res.data());

    if (resetRes0) {
        resFIM.SetInitRes();
    }
}

void Reservoir::ResetFIM()
{
    bulk.ResetFIM();
    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CaldG(bulk);
    allWells.CalFlux(bulk);
}

void Reservoir::UpdateLastStepFIM()
{
    OCP_FUNCNAME;

    bulk.UpdateLastStepFIM();
    allWells.UpdateLastBHP();
}


/////////////////////////////////////////////////////////////////////
// FIMn
/////////////////////////////////////////////////////////////////////


void Reservoir::InitFIMn()
{
    OCP_FUNCNAME;

    bulk.InitRock();
    bulk.InitSjPc(50);
    bulk.CalRock();
    bulk.InitFlashFIMn();
    bulk.CalKrPcFIMn();

    conn.CalAkd(bulk);

    allWells.InitBHP(bulk);

    UpdateLastStepFIM();
}


void Reservoir::UpdateLastStepFIMn()
{
    OCP_FUNCNAME;

    bulk.UpdateLastStepFIMn();
    allWells.UpdateLastBHP();
}

void Reservoir::AssembleMatFIMn(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.AssembleMat_FIM_new_n(myLS, bulk, dt);
    allWells.AssemblaMatFIM_new_n(myLS, bulk, dt);
}


void Reservoir::GetSolutionFIMn(const vector<OCP_DBL>& u,
                                 const OCP_DBL&         dPmax,
                                 const OCP_DBL&         dSmax)
{
    OCP_FUNCNAME;

    bulk.GetSolFIMn(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}


void Reservoir::AllocateFIMn_IsoT()
{
    OCP_FUNCNAME;

    bulk.AllocateFIMn_IsoT();
    conn.AllocateFIM_IsoT(bulk.GetPhaseNum());
}

void Reservoir::ResetFIMn()
{
    bulk.ResetFIMn();
    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CaldG(bulk);
    allWells.CalFlux(bulk);
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


void Reservoir::AllocateAIMc_IsoT()
{
    OCP_FUNCNAME;

    bulk.AllocateFIM_IsoT();
    bulk.AllocateAIMc_IsoT();
    conn.AllocateAIMc_IsoT(bulk.GetPhaseNum());
}

void Reservoir::AssembleMatAIMc(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_AIMc(myLS, bulk, dt);
    allWells.AssemblaMatFIM_new(myLS, bulk, dt);
}

void Reservoir::CalResAIMc(OCPRes& resAIMc, const OCP_DBL& dt, const OCP_BOOL& resetRes0)
{
    OCP_FUNCNAME;
    // Initialize
    resAIMc.SetZero();
    // Bulk to Bulk
    conn.CalResAIMc(resAIMc.res, bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resAIMc, bulk, dt);
    // Calculate RelRes
    bulk.CalRelResFIM(resAIMc);
    Dscalar(resAIMc.res.size(), -1, resAIMc.res.data());

    if (resetRes0) {
        resAIMc.SetInitRes();
    }
}


void Reservoir::GetSolutionAIMc(const vector<OCP_DBL>& u,
                                const OCP_DBL&         dPmax,
                                const OCP_DBL&         dSmax)
{
    bulk.GetSolAIMc(u, dPmax, dSmax);
    allWells.GetSolFIM(u, bulk.GetBulkNum(), bulk.GetComNum() + 1);
}

void Reservoir::InitAIMc()
{
    OCP_FUNCNAME;

    bulk.InitRock();
    bulk.InitSjPc(50);
    bulk.CalRock();
    bulk.InitFlashIMPEC();
    bulk.CalKrPcIMPEC();

    conn.CalAkd(bulk);

    allWells.InitBHP(bulk);

    UpdateLastStepAIMc();
}

void Reservoir::UpdateLastStepAIMc()
{
    OCP_FUNCNAME;

    bulk.UpdateLastStepAIMc();
    allWells.UpdateLastBHP();
}

void Reservoir::ResetAIMc()
{
    OCP_FUNCNAME;

    bulk.ResetAIMc();

    allWells.ResetBHP();
    allWells.CalTrans(bulk);
    allWells.CaldG(bulk);
    allWells.CalFlux(bulk);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/