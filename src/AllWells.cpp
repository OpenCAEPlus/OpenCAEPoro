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
        wells[w].group = paramWell.well[w].group;
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


void AllWells::SetupWellGroup(const Bulk& myBulk)
{
    wellGroup.clear();
    // Field Group, contain all wells
    wellGroup.push_back(WellGroup("Field"));
    for (USI w = 0; w < numWell; w++) {
        wellGroup[0].wId.push_back(w);
        if (wells[w].WellType() == INJ)  
            wellGroup[0].wIdINJ.push_back(w);
        else 
            wellGroup[0].wIdPROD.push_back(w);
    }

    // other subgroups
    USI glen = 1;
    USI g = 1;
    for (USI w = 0; w < numWell; w++) {
        for (g = 1; g < glen; g++) {
            if (wells[w].group == wellGroup[g].name) {
                // existing group
                wellGroup[g].wId.push_back(w);
                if (wells[w].WellType() == INJ)
                    wellGroup[g].wIdINJ.push_back(w);
                else
                    wellGroup[g].wIdPROD.push_back(w);
                break;
            }
        }
        if (g == glen && wells[w].group != "Field") {
            // new group
            wellGroup.push_back(WellGroup(wells[w].group));
            wellGroup[glen].wId.push_back(w);
            if (wells[w].WellType() == INJ)  wellGroup[glen].wIdINJ.push_back(w);
            else wellGroup[glen].wIdPROD.push_back(w);
            glen++;
        }
    } 

    numGroup = wellGroup.size();
    // relation between wellGroup should be completed if "node groups" exist

    // control of group should be update according to input file

    // for test
    if (OCP_FALSE) {
        wellGroup[0].reInj = OCP_TRUE;
        wellGroup[0].saleRate = 1500;
        wellGroup[0].injPhase = GAS;
        wellGroup[0].prodGroup = 0;
    }
    wellGroup[0].zi.resize(myBulk.GetComNum());
    CalReInjFluid(myBulk);

}


void AllWells::SetupMixture(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    flashCal = myBulk.GetMixture();
}

void AllWells::SetupWellBulk(Bulk& myBulk) const 
{
    myBulk.ClearWellBulkId();
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].SetupWellBulk(myBulk);
        }
    }
}

void AllWells::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;
    wellChange = OCP_FALSE;
    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        wells[w].opt = wells[w].optSet[i];
        if (wells[w].WellState()) {
            wells[w].wEId = wId;
            wId++;
        }
        if (i > 0 && wells[w].opt != wells[w].optSet[i - 1])
            wellChange = OCP_TRUE;
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
            wells[w].CalFlux(myBulk, OCP_TRUE);
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
            wells[w].CalFlux(myBulk, OCP_FALSE);
        }
    }
}

void AllWells::CalProdWeight(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalProdWeight(myBulk);
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


void AllWells::CalReInjFluid(const Bulk& myBulk)
{
    for (auto& wG : wellGroup) {
        if (wG.reInj) {
            const USI nc = myBulk.GetComNum();
            // const USI np = myBulk.GetPhaseNum();
            fill(wG.zi.begin(), wG.zi.end(), 0.0);
            for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                Daxpy(nc, 1.0, &wells[prod].qi_lbmol[0], &wG.zi[0]);
            }
            OCP_DBL qt = Dnorm1(nc, &wG.zi[0]);
            if (fabs(qt) < 1E-8) {
                // Recalculate zi
                fill(wG.zi.begin(), wG.zi.end(), 0.0);
                for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                    wells[prod].CalReInjFluid(myBulk, wG.zi);
                }
                qt = Dnorm1(nc, &wG.zi[0]);
            }
            flashCal[0]->Flash(PRESSURE_STD, TEMPERATURE_STD, &wG.zi[0], 0, 0, 0);
            Dcopy(nc, &wG.zi[0], &flashCal[0]->xij[wG.injPhase * nc]);
            wG.xi = flashCal[0]->xi[wG.injPhase];
            wG.factor = wG.xi * flashCal[0]->v[wG.injPhase] / qt;
            // assign to every open injection well in wG
            for (auto& w : wG.wIdINJ) {
                if (wells[w].WellState()) {
                    wells[w].opt.reInj = OCP_TRUE;
                    wells[w].opt.connWell = wG.wIdPROD;
                    wells[w].opt.injPhase = wG.injPhase;
                    wells[w].opt.zi = wG.zi;
                    wells[w].opt.xiINJ = wG.xi;
                    wells[w].opt.factor = wG.factor;
                    wells[w].opt.maxRate = -wG.saleRate * wG.xi * 1000; // Mscf -> ft3 -> lbmol
                }
            }
            //for (USI i = 0; i < nc; i++) {
            //    cout << scientific << wG.zi[i] << "   ";
            //}
            //cout << "xi:  " << wG.xi << "    Factor:   " << wG.factor;
            //cout << endl;
        }
    }
}


void AllWells::AllocateMat(LinearSystem& myLS, const USI& bulknum) const
{
    OCP_FUNCNAME;

    USI maxNum = (GetMaxWellPerNum() + 1) * numWell;
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
            // w.dG = w.ldG;
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
    OCP_BOOL flag2 = OCP_FALSE;
    OCP_BOOL flag3 = OCP_FALSE;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            OCP_INT flag = wells[w].CheckP(myBulk);
//#ifdef _DEBUG
            // wells[w].ShowPerfStatus();
//#endif // _DEBUG
            switch (flag) {
                case 1:
                    return 1;
                case 2:
                    flag2 = OCP_TRUE;
                    break;
                case 3:
                    flag3 = OCP_TRUE;
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

USI AllWells::GetWellPerfNum() const
{
    USI numPerf = 0;
    for (USI w = 0; w < numWell; w++) {
        numPerf += wells[w].numPerf;
    }
    return numPerf;
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

OCP_DBL AllWells::CalWellQT()
{
    QT = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            for (auto& q : wells[w].qi_lbmol) {
                QT += q;
            }
        }
    }
    return QT;
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////


void AllWells::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalCFL(myBulk, dt);
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

    // for Reinjection
    for (auto& wG : wellGroup) {
        if (wG.reInj) {
            for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                if (wells[prod].WellState()) {
                    wells[prod].AssembleMatReinjection_IMPEC(myBulk, myLS, dt, wells, wG.wIdINJ);
                }
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

    // for Reinjection
    for (auto& wG : wellGroup) {
        if (wG.reInj) {
            for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                if (wells[prod].WellState()) {
                    wells[prod].AssembleMatReinjection_FIM(myBulk, myLS, dt, wells, wG.wIdINJ);
                }
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


void AllWells::CalResFIM(ResFIM& resFIM, const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalResFIM(resFIM, myBulk, dt, wId, wells);
            wId++;
        }
    }

    // cout << "Well  " << resFIM.maxRelRes_v;
}

void AllWells::ShowRes(const vector<OCP_DBL>& res, const Bulk& myBulk) const
{
    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            cout << endl;
            wells[w].ShowRes(wId, res, myBulk);
            wId++;
        }
    }
    cout << endl;
}


/////////////////////////////////////////////////////////////////////
// AFIM(new)
/////////////////////////////////////////////////////////////////////


void AllWells::AssemblaMatFIM_new(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
            case INJ:
                wells[w].AssembleMatINJ_FIM_new(myBulk, myLS, dt);
                break;
            case PROD:
                wells[w].AssembleMatPROD_FIM_new(myBulk, myLS, dt);
                break;
            default:
                OCP_ABORT("Wrong well type");
            }
        }
    }
}


void AllWells::AssemblaMatFIM_new_n(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
            case INJ:
                wells[w].AssembleMatINJ_FIM_new_n(myBulk, myLS, dt);
                break;
            case PROD:
                wells[w].AssembleMatPROD_FIM_new_n(myBulk, myLS, dt);
                break;
            default:
                OCP_ABORT("Wrong well type");
            }
        }
    }
}



/////////////////////////////////////////////////////////////////////
// AIMt
/////////////////////////////////////////////////////////////////////

void AllWells::AssemblaMatAIMt(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
            case INJ:
                wells[w].AssembleMatINJ_AIMt(myBulk, myLS, dt);
                break;
            case PROD:
                wells[w].AssembleMatPROD_AIMt(myBulk, myLS, dt);
                break;
            default:
                OCP_ABORT("Wrong well type");
            }
        }
    }
}

void AllWells::CalResAIMt(ResFIM& resFIM, const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI wId = 0;
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {
            wells[w].CalResAIMt(resFIM, myBulk, dt, wId, wells);
            wId++;
        }
    }

    // cout << "Well  " << resFIM.maxRelRes_v;
}

void AllWells::GetSolAIMt(const vector<OCP_DBL>& u, const OCP_USI& bId, const USI& len)
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

void AllWells::AssemblaMatAIMs(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    for (USI w = 0; w < numWell; w++) {
        if (wells[w].WellState()) {

            switch (wells[w].WellType()) {
            case INJ:
                wells[w].AssembleMatINJ_AIMs(myBulk, myLS, dt);
                break;
            case PROD:
                wells[w].AssembleMatPROD_AIMs(myBulk, myLS, dt);
                break;
            default:
                OCP_ABORT("Wrong well type");
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/