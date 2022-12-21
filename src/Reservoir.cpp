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
    allWells.InputParam(param.paramWell);
}

void Reservoir::SetupIsoT()
{
    OCP_FUNCNAME;

    grid.SetupIsoT();
    conn.Setup(grid.initInfo);
    allWells.Setup(grid);
}

void Reservoir::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;

    allWells.ApplyControl(i);
    allWells.SetupWellGroup(grid.bulk);
}

void Reservoir::PrepareWell()
{
    OCP_FUNCNAME;

    allWells.PrepareWell(grid.bulk);
}

void Reservoir::CalWellFlux()
{
    OCP_FUNCNAME;

    allWells.CalFlux(grid.bulk);
}

void Reservoir::CalWellTrans()
{
    OCP_FUNCNAME;

    allWells.CalTrans(grid.bulk);
}

void Reservoir::CalRock()
{
    OCP_FUNCNAME;

    grid.bulk.CalRock();
}

void Reservoir::CalKrPc()
{
    OCP_FUNCNAME;

    grid.bulk.CalKrPcIMPEC();
}

void Reservoir::CalMaxChange()
{
    OCP_FUNCNAME;

    grid.bulk.CalMaxChange();
    allWells.CalMaxBHPChange();
}

void Reservoir::CalIPRT(const OCP_DBL& dt)
{
    OCP_FUNCNAME;
    // Calculate injection / production rate for current step
    allWells.CalIPRT(grid.bulk, dt);
    // Calculate Reinjection fluid for next step
    allWells.CalReInjFluid(grid.bulk);
}

OCP_INT Reservoir::CheckP(const OCP_BOOL& bulkCheck, const OCP_BOOL& wellCheck)
{
    OCP_FUNCNAME;

    if (bulkCheck) {
        if (!grid.bulk.CheckP()) {
            // negative Pressure
            return 1;
        }
    }

    if (wellCheck) {
        OCP_INT flag = 0;
        flag         = allWells.CheckP(grid.bulk);
        return flag;
    }

    return 0;
}

OCP_BOOL Reservoir::CheckNi()
{
    OCP_FUNCNAME;

    return grid.bulk.CheckNi();
}

OCP_BOOL Reservoir::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    return grid.bulk.CheckVe(Vlim);
}


OCP_DBL Reservoir::CalCFL(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    conn.CalCFL(grid.bulk, dt);
    allWells.CalCFL(grid.bulk, dt);
    grid.bulk.CalCFL();

    return grid.bulk.maxCFL;
}


/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateIMPEC_IsoT()
{
    OCP_FUNCNAME;

    grid.bulk.AllocateIMPEC_IsoT();
    conn.AllocateIMPEC_IsoT(grid.bulk.GetPhaseNum());
}

void Reservoir::InitIMPEC()
{
    OCP_FUNCNAME;

    grid.bulk.InitRock();
    grid.bulk.InitSjPc(50);
    grid.bulk.CalRock();
    grid.bulk.InitFlashIMPEC();
    grid.bulk.CalKrPcIMPEC();
    grid.bulk.UpdateLastStepIMPEC();

    conn.CalAkd(grid.bulk);
    conn.CalFluxIMPEC(grid.bulk);
    conn.UpdateLastStepIMPEC();

    allWells.InitBHP(grid.bulk);
    allWells.UpdateLastBHP();
}

void Reservoir::AllocateMatIMPEC(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(GetBulkNum() + GetWellNum(), 1);
    conn.AllocateMat(myLS);
    allWells.AllocateMat(myLS, GetBulkNum());
    myLS.AllocateColMem();
}

void Reservoir::AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMatIMPEC(myLS, grid.bulk, dt);
    allWells.AssemblaMatIMPEC(myLS, grid.bulk, dt);
}

void Reservoir::GetSolutionIMPEC(const vector<OCP_DBL>& u)
{
    OCP_FUNCNAME;

    grid.bulk.GetSolIMPEC(u);
    allWells.GetSolIMPEC(u, grid.bulk.GetBulkNum());
}


void Reservoir::CalFLuxIMPEC()
{
    OCP_FUNCNAME;

    conn.CalFluxIMPEC(grid.bulk);
    allWells.CalFlux(grid.bulk);
}

void Reservoir::MassConserveIMPEC(const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    conn.MassConserveIMPEC(grid.bulk, dt);
    allWells.MassConserveIMPEC(grid.bulk, dt);
}


void Reservoir::ResetVal01IMPEC()
{
    OCP_FUNCNAME;
    grid.bulk.ResetVal01IMPEC();
    conn.ResetIMPEC();
}

void Reservoir::ResetVal02IMPEC()
{
    OCP_FUNCNAME;

    grid.bulk.ResetVal02IMPEC();
    conn.ResetIMPEC();
}

void Reservoir::ResetVal03IMPEC()
{
    OCP_FUNCNAME;
    grid.bulk.ResetVal03IMPEC();
    conn.ResetIMPEC();
}


void Reservoir::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;
    grid.bulk.UpdateLastStepIMPEC();
    conn.UpdateLastStepIMPEC();
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Reservoir::AllocateFIM_IsoT()
{
    OCP_FUNCNAME;

    grid.bulk.AllocateFIM_IsoT();
    conn.AllocateFIM_IsoT(grid.bulk.GetPhaseNum());
}

void Reservoir::InitFIM()
{
    OCP_FUNCNAME;

    grid.bulk.InitRock();
    grid.bulk.InitSjPc(50);
    grid.bulk.CalRock();
    grid.bulk.InitFlashFIM();
    grid.bulk.CalKrPcFIM();

    conn.CalAkd(grid.bulk);
 
    allWells.InitBHP(grid.bulk);

    UpdateLastStepFIM();
}

void Reservoir::AllocateMatFIM_IsoT(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.AllocateRowMem(GetBulkNum() + GetWellNum(), GetComNum() + 1);
    conn.AllocateMat(myLS);
    allWells.AllocateMat(myLS, GetBulkNum());
    myLS.AllocateColMem();
}


void Reservoir::AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);

#ifdef OCP_OLD_FIM
    conn.AssembleMat_FIM(myLS, grid.bulk, dt);
    allWells.AssemblaMatFIM(myLS, grid.bulk, dt);
#else
    conn.AssembleMat_FIM_new(myLS, grid.bulk, dt);
    allWells.AssemblaMatFIM_new(myLS, grid.bulk, dt);
#endif // OCP_OLD_FIM
}


void Reservoir::GetSolutionFIM(const vector<OCP_DBL>& u,
    const OCP_DBL& dPmax,
    const OCP_DBL& dSmax)
{
    grid.bulk.GetSolFIM(u, dPmax, dSmax);
    allWells.GetSolFIM(u, grid.bulk.GetBulkNum(), grid.bulk.GetComNum() + 1);
}


void Reservoir::CalResFIM(OCPRes& resFIM, const OCP_DBL& dt, const OCP_BOOL& resetRes0)
{
    OCP_FUNCNAME;
    // Initialize
    resFIM.SetZero();
    // Bulk to Bulk
    conn.CalResFIM(resFIM.res, grid.bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resFIM, grid.bulk, dt);
    // Calculate RelRes
    grid.bulk.CalRelResFIM(resFIM);
    Dscalar(resFIM.res.size(), -1.0, resFIM.res.data());

    if (resetRes0) {
        resFIM.SetInitRes();
    }
}

void Reservoir::ResetFIM()
{
    grid.bulk.ResetFIM();
    allWells.ResetBHP();
    allWells.CalTrans(grid.bulk);
    allWells.CaldG(grid.bulk);
    allWells.CalFlux(grid.bulk);
}

void Reservoir::UpdateLastStepFIM()
{
    OCP_FUNCNAME;

    grid.bulk.UpdateLastStepFIM();
    allWells.UpdateLastBHP();
}


/////////////////////////////////////////////////////////////////////
// FIMn
/////////////////////////////////////////////////////////////////////


void Reservoir::InitFIMn()
{
    OCP_FUNCNAME;

    grid.bulk.InitRock();
    grid.bulk.InitSjPc(50);
    grid.bulk.CalRock();
    grid.bulk.InitFlashFIMn();
    grid.bulk.CalKrPcFIMn();

    conn.CalAkd(grid.bulk);

    allWells.InitBHP(grid.bulk);

    UpdateLastStepFIM();
}


void Reservoir::UpdateLastStepFIMn()
{
    OCP_FUNCNAME;

    grid.bulk.UpdateLastStepFIMn();
    allWells.UpdateLastBHP();
}

void Reservoir::AssembleMatFIMn(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_FIM_new_n(myLS, grid.bulk, dt);
    allWells.AssemblaMatFIM_new_n(myLS, grid.bulk, dt);
}


void Reservoir::GetSolutionFIMn(const vector<OCP_DBL>& u,
                                 const OCP_DBL&         dPmax,
                                 const OCP_DBL&         dSmax)
{
    OCP_FUNCNAME;

    grid.bulk.GetSolFIMn(u, dPmax, dSmax);
    allWells.GetSolFIM(u, grid.bulk.GetBulkNum(), grid.bulk.GetComNum() + 1);
}


void Reservoir::AllocateFIMn_IsoT()
{
    OCP_FUNCNAME;

    grid.bulk.AllocateFIMn_IsoT();
    conn.AllocateFIM_IsoT(grid.bulk.GetPhaseNum());
}

void Reservoir::ResetFIMn()
{
    grid.bulk.ResetFIMn();
    allWells.ResetBHP();
    allWells.CalTrans(grid.bulk);
    allWells.CaldG(grid.bulk);
    allWells.CalFlux(grid.bulk);
}

void Reservoir::PrintSolFIM(const string& outfile) const
{
    ofstream outu(outfile);
    if (!outu.is_open()) cout << "Can not open " << outfile << endl;
    const OCP_USI nb = grid.bulk.numBulk;
    const OCP_USI nc = grid.bulk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        // Pressure
        outu << grid.bulk.P[n] << "\n";
        // Ni
        for (USI i = 0; i < nc; i++) {
            outu << grid.bulk.Ni[n * nc + i] << "\n";
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

    grid.bulk.AllocateFIM_IsoT();
    grid.bulk.AllocateAIMc_IsoT();
    conn.AllocateAIMc_IsoT(grid.bulk.GetPhaseNum());
}

void Reservoir::AssembleMatAIMc(LinearSystem& myLS, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    conn.SetupMatSparsity(myLS);
    conn.AssembleMat_AIMc(myLS, grid.bulk, dt);
    allWells.AssemblaMatFIM_new(myLS, grid.bulk, dt);
}

void Reservoir::CalResAIMc(OCPRes& resAIMc, const OCP_DBL& dt, const OCP_BOOL& resetRes0)
{
    OCP_FUNCNAME;
    // Initialize
    resAIMc.SetZero();
    // Bulk to Bulk
    conn.CalResAIMc(resAIMc.res, grid.bulk, dt);
    // Well to Bulk
    allWells.CalResFIM(resAIMc, grid.bulk, dt);
    // Calculate RelRes
    grid.bulk.CalRelResFIM(resAIMc);
    Dscalar(resAIMc.res.size(), -1, resAIMc.res.data());

    if (resetRes0) {
        resAIMc.SetInitRes();
    }
}


void Reservoir::GetSolutionAIMc(const vector<OCP_DBL>& u,
                                const OCP_DBL&         dPmax,
                                const OCP_DBL&         dSmax)
{
    grid.bulk.GetSolAIMc(u, dPmax, dSmax);
    allWells.GetSolFIM(u, grid.bulk.GetBulkNum(), grid.bulk.GetComNum() + 1);
}

void Reservoir::InitAIMc()
{
    OCP_FUNCNAME;

    grid.bulk.InitRock();
    grid.bulk.InitSjPc(50);
    grid.bulk.CalRock();
    grid.bulk.InitFlashIMPEC();
    grid.bulk.CalKrPcIMPEC();

    conn.CalAkd(grid.bulk);

    allWells.InitBHP(grid.bulk);

    UpdateLastStepAIMc();
}

void Reservoir::UpdateLastStepAIMc()
{
    OCP_FUNCNAME;

    grid.bulk.UpdateLastStepAIMc();
    allWells.UpdateLastBHP();
}

void Reservoir::ResetAIMc()
{
    OCP_FUNCNAME;

    grid.bulk.ResetAIMc();

    allWells.ResetBHP();
    allWells.CalTrans(grid.bulk);
    allWells.CaldG(grid.bulk);
    allWells.CalFlux(grid.bulk);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/