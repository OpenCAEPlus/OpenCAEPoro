/*! \file    WellGroup.cpp
 *  \brief   WellGroup class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "WellGroup.hpp"


 /////////////////////////////////////////////////////////////////////
 // General
 /////////////////////////////////////////////////////////////////////


void WellGroup::InputParam(const ParamWell& paramWell) { OCP_FUNCNAME;

    numWell = paramWell.well.size();
    wellGroup.resize(numWell);
    USI         t = paramWell.criticalTime.size();
    vector<USI> wellCriticalTime;
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].name  = paramWell.well[w].name;
        wellGroup[w].depth = paramWell.well[w].depth;
        wellGroup[w].I     = paramWell.well[w].I - 1;
        wellGroup[w].J     = paramWell.well[w].J - 1;

        wellGroup[w].InputPerfo(paramWell.well[w]);

        // opt
        wellGroup[w].optSet.resize(t);
        USI n = paramWell.well[w].optParam.size();
        wellCriticalTime.clear();
        wellCriticalTime.resize(n + 1);
        for (USI i = 0; i < n; i++) {
            wellCriticalTime[i] = paramWell.well[w].optParam[i].d;
        }
        wellCriticalTime.back() = t;
        for (USI i = 0; i < n; i++) {
            for (USI d = wellCriticalTime[i]; d < wellCriticalTime[i + 1]; d++) {
                wellGroup[w].optSet[d] = WellOpt(paramWell.well[w].optParam[i].opt);
            }
        }
    }
}


void WellGroup::Setup(const Grid& myGrid, const Bulk& myBulk) { OCP_FUNCNAME;

    SetupWell(myGrid, myBulk);
    SetupMixture(myBulk);
}


void WellGroup::SetupWell(const Grid& myGrid, const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].Setup(myGrid, myBulk);
    }
}


void WellGroup::SetupMixture(const Bulk& myBulk) { OCP_FUNCNAME;

    flashCal = myBulk.GetMixture();
}


void WellGroup::ApplyControl(const USI& i) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].opt = wellGroup[w].optSet[i];
    }
}


void WellGroup::InitBHP(const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].InitBHP(myBulk);
    }
}


void WellGroup::PrepareWell(const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            wellGroup[w].CalTrans(myBulk);
            wellGroup[w].CalFlux(myBulk, true);
            wellGroup[w].CaldG(myBulk);
            // test
            // wellGroup[w].SmoothdG();
            wellGroup[w].CheckOptMode(myBulk);
        }
    }
}


void WellGroup::CalTrans(const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CalTrans(myBulk);
        }
    }
}


void WellGroup::CalFlux(const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CalFlux(myBulk);
        }
    }
}


void WellGroup::CaldG(const Bulk& myBulk) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CaldG(myBulk);
        }
    }
}


void WellGroup::CalIPRT(const Bulk& myBulk, OCP_DBL dt) {OCP_FUNCNAME;

    FGIR = 0;
    FWIR = 0;
    FOPR = 0;
    FGPR = 0;
    FWPR = 0;
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].WGIR = 0;
        wellGroup[w].WWIR = 0;
        wellGroup[w].WOPR = 0;
        wellGroup[w].WGPR = 0;
        wellGroup[w].WWPR = 0;

        if (wellGroup[w].WellState()) {
            if (wellGroup[w].WellType() == PROD) {
                wellGroup[w].CalProdQiBO(myBulk, dt);
            }
            else {
                wellGroup[w].CalInjQiBO(myBulk, dt);
            }
        }
        FGIR += wellGroup[w].WGIR;
        FWIR += wellGroup[w].WWIR;
        FOPR += wellGroup[w].WOPR;
        FGPR += wellGroup[w].WGPR;
        FWPR += wellGroup[w].WWPR;
    }
    FGIT += FGIR * dt;
    FWIT += FWIR * dt;
    FOPT += FOPR * dt;
    FGPt += FGPR * dt;
    FWPT += FWPR * dt;
}


void WellGroup::AllocateMat(LinearSolver& mySolver, const USI& bulknum) const { OCP_FUNCNAME;
    
    USI maxNum = GetMaxWellPerNum() + 1;
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].AllocateMat(mySolver);
        mySolver.RowCapPlus(bulknum + w, maxNum);
    }
}


OCP_INT WellGroup::CheckP(const Bulk& myBulk) { OCP_FUNCNAME;

    // 0 : All correct
    // 1   : negative P, cut the timestep and resolve
    // 2.1 : change well mode to BHP, resolve
    // 2.2 : crossflow happens, then close corresponding perf, resolve
    bool flag2 = false;
    bool flag3 = false;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            OCP_INT flag = wellGroup[w].CheckP(myBulk);
#ifdef _DEBUG
            wellGroup[w].ShowPerfStatus();
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
USI WellGroup::GetIndex(const string& name) const { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].name == name) {
            return w;
        }
    }
    OCP_ABORT("Well name not found!");
}


USI WellGroup::GetMaxWellPerNum() const { OCP_FUNCNAME;

    USI m = 0;
    for (USI w = 0; w < numWell; w++) {
        m = max(m, wellGroup[w].numPerf);
    }
    return m;
}


 /////////////////////////////////////////////////////////////////////
 // IMPEC
 /////////////////////////////////////////////////////////////////////


OCP_DBL WellGroup::CalCFLIMPEC(const Bulk& myBulk, const OCP_DBL& dt) const { OCP_FUNCNAME;

    OCP_DBL cflw = 0;
    OCP_DBL tmp  = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            tmp = wellGroup[w].CalCFLIMPEC(myBulk, dt);
            if (cflw < tmp) cflw = tmp;
        }
    }
    return cflw;
}


void WellGroup::CalCFL01IMPEC(const Bulk& myBulk, const OCP_DBL& dt) const { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CalCFL01IMPEC(myBulk, dt);
        }
    }
}


void WellGroup::MassConserveIMPEC(Bulk& myBulk, OCP_DBL dt) { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].MassConserveIMPEC(myBulk, dt);
        }
    }
}


void WellGroup::AssemblaMatIMPEC(LinearSolver& mySolver, const Bulk& myBulk,
                                     const OCP_DBL& dt) const { OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            switch (wellGroup[w].WellType()) {
                case INJ:
                    wellGroup[w].AssembleMatINJ_IMPEC(myBulk, mySolver, dt);
                    break;
                case PROD:
                    wellGroup[w].AssembleMatPROD_BO_IMPEC(myBulk, mySolver, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}


void WellGroup::GetSolIMPEC(const vector<OCP_DBL>& u, const OCP_USI& bId) { OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].BHP = u[bId + wId];
            wellGroup[w].UpdatePerfP();
            wId++;
        }
    }
}


 /////////////////////////////////////////////////////////////////////
 // FIM
 /////////////////////////////////////////////////////////////////////


void WellGroup::AssemblaMatFIM(LinearSolver& mySolver, const Bulk& myBulk,
    const OCP_DBL& dt) const { OCP_FUNCNAME;
     
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            switch (wellGroup[w].WellType()) {
            case INJ:
                wellGroup[w].AssembleMatINJ_FIM(myBulk, mySolver, dt);
                break;
            case PROD:
                wellGroup[w].AssembleMatPROD_BO_FIM(myBulk, mySolver, dt);
                break;
            default:
                OCP_ABORT("Wrong well type");
            }
        }
    }
}


void WellGroup::GetSolFIM(const vector<OCP_DBL>& u, const OCP_USI& bId, const USI& len) { OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].BHP += u[(bId + wId)*len];
            wellGroup[w].UpdatePerfP();
            wId++;
        }
    }
}


void WellGroup::CalResFIM(ResFIM& resFIM, const Bulk& myBulk, const OCP_DBL& dt) const { OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CalResFIM(resFIM, myBulk, dt, wId);
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