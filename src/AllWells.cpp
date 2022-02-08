/*! \file    AllWells.cpp
 *  \brief   AllWells class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "AllWells.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////

void AllWells::InputParam(const ParamWell& paramWell)
{
    OCP_FUNCNAME;

    USI len = paramWell.solSet.size();
    solvents.resize(len);
    for (USI i = 0; i < len; i++) {
        solvents[i] = paramWell.solSet[i];
    }

    numWell = paramWell.well.size();
    wells.resize(numWell);
    USI         t = paramWell.criticalTime.size();
    vector<USI> wellCriticalTime;
    for (USI w = 0; w < numWell; w++) {
        wells[w].name  = paramWell.well[w].name;
        wells[w].depth = paramWell.well[w].depth;
        wells[w].I     = paramWell.well[w].I - 1;
        wells[w].J     = paramWell.well[w].J - 1;

        wells[w].InputPerfo(paramWell.well[w]);

        // opt
        wells[w].optSet.resize(t);
        USI n = paramWell.well[w].optParam.size();
        wellCriticalTime.clear();
        wellCriticalTime.resize(n + 1);
        for (USI i = 0; i < n; i++) {
            wellCriticalTime[i] = paramWell.well[w].optParam[i].d;
        }
        wellCriticalTime.back() = t;
        for (USI i = 0; i < n; i++) {
            for (USI d = wellCriticalTime[i]; d < wellCriticalTime[i + 1]; d++) {
                wells[w].optSet[d] = WellOpt(paramWell.well[w].optParam[i].opt);
            }
        }
    }
}

void AllWells::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;

    SetupWell(myGrid, myBulk);
    SetupMixture(myBulk);
}

void AllWells::SetupWell(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w].Setup(myGrid, myBulk, solvents);
    }
}

void AllWells::SetupMixture(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    flashCal = myBulk.GetMixture();
}

void AllWells::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;
    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        wells[w].opt = wells[w].optSet[i];
        if (wells[w].WellState()) {
            wells[w].wEId = wId;
            wId++;
        }
    }
}

void AllWells::InitBHP(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w].InitBHP(myBulk);
    }
}

void AllWells::PrepareWell(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            wells[w].CalTrans(myBulk);
            wells[w].CalFlux(myBulk, true);
            wells[w].CalProdWeight(myBulk);
            wells[w].CaldG(myBulk);
            // test
            // wells[w].SmoothdG();
            wells[w].CheckOptMode(myBulk);
        }
    }
}

void AllWells::CalTrans(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalTrans(myBulk);
        }
    }
}

void AllWells::CalFlux(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalFlux(myBulk);
        }
    }
}

void AllWells::CaldG(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CaldG(myBulk);
        }
    }
}

void AllWells::CalIPRT(const Bulk& myBulk, OCP_DBL dt)
{
    OCP_FUNCNAME;

    FGIR = 0;
    FWIR = 0;
    FOPR = 0;
    FGPR = 0;
    FWPR = 0;
    for (USI w = 0; w < numWell; w++) {
        wells[w].WGIR = 0;
        wells[w].WWIR = 0;
        wells[w].WOPR = 0;
        wells[w].WGPR = 0;
        wells[w].WWPR = 0;

        if (wells[w].WellState()) {
            if (wells[w].WellType() == PROD) {
                wells[w].CalProdQj(myBulk, dt);
            } else {
                wells[w].CalInjQi(myBulk, dt);
            }
        }
        FGIR += wells[w].WGIR;
        FWIR += wells[w].WWIR;
        FOPR += wells[w].WOPR;
        FGPR += wells[w].WGPR;
        FWPR += wells[w].WWPR;
    }
    FGIT += FGIR * dt;
    FWIT += FWIR * dt;
    FOPT += FOPR * dt;
    FGPt += FGPR * dt;
    FWPT += FWPR * dt;
}

void AllWells::AllocateMat(LinearSystem& myLS, const USI& bulknum) const
{
    OCP_FUNCNAME;

    USI maxNum = GetMaxWellPerNum() + 1;
    for (USI w = 0; w < numWell; w++) {
        wells[w].AllocateMat(myLS);
        myLS.EnlargeRowCap(bulknum + w, maxNum);
    }
}

void AllWells::ResetBHP()
{
    for (auto& w : wells) {
        if (w.WellState()) {
            w.BHP = w.lBHP;
            w.SetBHP();
        }
    }
}

OCP_INT AllWells::CheckP(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    // 0 : All correct
    // 1   : negative P, cut the timestep and resolve
    // 2.1 : change well mode to BHP, resolve
    // 2.2 : crossflow happens, then close corresponding perf, resolve
    bool flag2 = false;
    bool flag3 = false;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            OCP_INT flag = wells[w].CheckP(myBulk);
#ifdef _DEBUG
            // wells[w].ShowPerfStatus();
#endif // _DEBUG
            switch (flag) {
                case 1:
                    return 1;
                case 2:
                    flag2 = true;
                    break;
                case 3:
                    flag3 = true;
                    break;
                default:
                    break;
            }
        }
    }

    if (flag2 || flag3) return 2;

    return 0;
}

// return the index of Specified well name
USI AllWells::GetIndex(const string& name) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].name == name) {
            return w;
        }
    }
    OCP_ABORT("Well name not found!");
}

USI AllWells::GetMaxWellPerNum() const
{
    OCP_FUNCNAME;

    USI m = 0;
    for (USI w = 0; w < numWell; w++) {
        m = max(m, wells[w].numPerf);
    }
    return m;
}

void AllWells::CalMaxBHPChange()
{
    dPmax = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            dPmax = max(dPmax, fabs(wells[w].BHP - wells[w].lBHP));
        }
    }
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

OCP_DBL AllWells::CalCFLIMPEC(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    OCP_DBL cflw = 0;
    OCP_DBL tmp  = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            tmp = wells[w].CalCFLIMPEC(myBulk, dt);
            if (cflw < tmp) cflw = tmp;
        }
    }
    return cflw;
}

void AllWells::CalCFL01IMPEC(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalCFL01IMPEC(myBulk, dt);
        }
    }
}

void AllWells::MassConserveIMPEC(Bulk& myBulk, OCP_DBL dt)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].MassConserveIMPEC(myBulk, dt);
        }
    }
}

void AllWells::AssemblaMatIMPEC(LinearSystem& myLS, const Bulk& myBulk,
                                 const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
                case INJ:
                    wells[w].AssembleMatINJ_IMPEC(myBulk, myLS, dt);
                    break;
                case PROD:
                    wells[w].AssembleMatPROD_IMPEC(myBulk, myLS, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}

void AllWells::GetSolIMPEC(const vector<OCP_DBL>& u, const OCP_USI& bId)
{
    OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].BHP = u[bId + wId];
            wells[w].UpdatePerfP();
            wId++;
        }
    }
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void AllWells::AssemblaMatFIM(LinearSystem& myLS, const Bulk& myBulk,
                               const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
                case INJ:
                    wells[w].AssembleMatINJ_FIM(myBulk, myLS, dt);
                    break;
                case PROD:
                    wells[w].AssembleMatPROD_FIM(myBulk, myLS, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}

void AllWells::GetSolFIM(const vector<OCP_DBL>& u, const OCP_USI& bId, const USI& len)
{
    OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].BHP += u[(bId + wId) * len];
            wells[w].UpdatePerfP();
            wId++;
        }
    }
}

void AllWells::GetSol01FIM(const vector<OCP_DBL>& u, const OCP_USI& bId,
                            const USI& len, const OCP_DBL& alpha)
{
    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].BHP += u[(bId + wId) * len] * alpha;
            wells[w].UpdatePerfP();
            wId++;
        }
    }
}

void AllWells::CalResFIM(ResFIM& resFIM, const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalResFIM(resFIM, myBulk, dt, wId);
            wId++;
        }
    }

    // cout << "Well  " << resFIM.maxRelRes_v;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/