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

void WellGroup::InputParam(const ParamWell& param_Well)
{
    numWell = param_Well.well.size();
    wellGroup.resize(numWell);
    USI         t = param_Well.criticalTime.size();
    vector<USI> wellCriticalTime;
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].name       = param_Well.well[w].name;
        wellGroup[w].depth      = param_Well.well[w].depth;
        wellGroup[w].radius     = param_Well.well[w].diameter / 2;
        wellGroup[w].kh         = param_Well.well[w].kh;
        wellGroup[w].skinFactor = param_Well.well[w].skinFactor;
        wellGroup[w].WI         = param_Well.well[w].WI;
        wellGroup[w].I          = param_Well.well[w].I - 1;
        wellGroup[w].J          = param_Well.well[w].J - 1;
        wellGroup[w].K1         = param_Well.well[w].K1 - 1;
        wellGroup[w].K2         = param_Well.well[w].K2 - 1;

        // opt
        wellGroup[w].optSet.resize(t);
        USI n = param_Well.well[w].optParam.size();
        wellCriticalTime.clear();
        wellCriticalTime.resize(n + 1);
        for (USI i = 0; i < n; i++) {
            wellCriticalTime[i] = param_Well.well[w].optParam[i].d;
        }
        wellCriticalTime.back() = t;
        for (USI i = 0; i < n; i++) {
            for (USI d = wellCriticalTime[i]; d < wellCriticalTime[i + 1]; d++) {
                wellGroup[w].optSet[d] = WellOpt(param_Well.well[w].optParam[i].opt);
            }
        }
    }

    cout << "WellGroup::InputParam" << endl;
}

void WellGroup::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    SetupWell(myGrid, myBulk);
    SetupMixture(myBulk);
}

void WellGroup::SetupWell(const Grid& myGrid, const Bulk& myBulk)
{
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].Setup(myGrid, myBulk);
    }
}

void WellGroup::SetupMixture(const Bulk& myBulk)
{

    flashCal = myBulk.GetMixture();
    cout << "WellGroup::setupMixture" << endl;
}

void WellGroup::Init(const Bulk& myBulk)
{
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].Init(myBulk);
    }
}

void WellGroup::ApplyControl(USI i)
{
    for (USI w = 0; w < numWell; w++) {
        wellGroup[w].opt = wellGroup[w].optSet[i];
    }
}

OCP_DBL WellGroup::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_DBL cflw = 0;
    OCP_DBL tmp  = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            tmp = wellGroup[w].CalCFL(myBulk, dt);
            if (cflw < tmp) cflw = tmp;
        }
    }
    return cflw;
}

void WellGroup::CalFlux(const Bulk& myBulk)
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].CalFlux(myBulk);
        }
    }
}

void WellGroup::MassConserve(Bulk& myBulk, OCP_DBL dt)
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].MassConserve(myBulk, dt);
        }
    }
}

void WellGroup::PrepareWell(const Bulk& myBulk)
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            wellGroup[w].CalTrans(myBulk);
            wellGroup[w].CalFlux(myBulk, true);
            wellGroup[w].CaldG(myBulk);
            // test
            wellGroup[w].SmoothdG();
            wellGroup[w].CheckOptMode(myBulk);
        }
    }
}

void WellGroup::AssemblaMat_WB_IMPES(Solver<OCP_DBL>& mySolver, const Bulk& myBulk,
                                     const OCP_DBL& dt) const
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {

            switch (wellGroup[w].WellType()) {
                case INJ:
                    wellGroup[w].AssembleMat_INJ_IMPES(myBulk, mySolver, dt);
                    break;
                case PROD:
                    wellGroup[w].AssembleMat_PROD_BLK_IMPES(myBulk, mySolver, dt);
                    break;
                default:
                    ERRORcheck("Wrong Well Type in function");
                    exit(0);
            }
        }
    }
}

void WellGroup::GetSol_IMPES(const vector<OCP_DBL>& u, const OCP_USI& bid)
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].WellState()) {
            wellGroup[w].BHP = u[bid + w];
            wellGroup[w].UpdatePerfP();
        }
    }
}

void WellGroup::CalIPRT(const Bulk& myBulk, OCP_DBL dt)
{
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
                wellGroup[w].CalProdQi_Blk(myBulk, dt);
            } else {
                wellGroup[w].CalInjQi_Blk(myBulk, dt);
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

OCP_INT WellGroup::CheckP(const Bulk& myBulk)
{
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
USI WellGroup::GetIndex(const string& name) const
{
    for (USI w = 0; w < numWell; w++) {
        if (wellGroup[w].name == name) {
            return w;
        }
    }
    ERRORcheck("No such well name!");
    exit(0);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/