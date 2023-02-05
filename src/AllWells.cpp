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

void AllWells::InputParam(const ParamWell& paramWell, const ParamOutput& output_param)
{
    OCP_FUNCNAME;

    for (auto& s : paramWell.solSet) {
        solvents.push_back(s);
    }

    Psurf = paramWell.Psurf;
    Tsurf = paramWell.Tsurf;

    numWell = paramWell.well.size();
    wells.resize(numWell);
    USI         t = paramWell.criticalTime.size();
    vector<USI> wellOptTime;
    for (USI w = 0; w < numWell; w++) {
        wells[w].name  = paramWell.well[w].name;
        wells[w].group = paramWell.well[w].group;
        wells[w].depth = paramWell.well[w].depth;
        wells[w].I     = paramWell.well[w].I - 1;
        wells[w].J     = paramWell.well[w].J - 1;
        wells[w].InputPerfo(paramWell.well[w]);
        wells[w].Psurf         = Psurf;
        wells[w].Tsurf         = Tsurf;
        wells[w].ifUseUnweight = paramWell.well[w].ifUseUnweight;
        // opt
        wells[w].optSet.resize(t);

        // If identical optParam.d occurs, use the last one
        vector<WellOptPair> tmpOptParam;
        const USI           n0 = paramWell.well[w].optParam.size();
        for (USI i = 0; i < n0 - 1; i++) {
            if (paramWell.well[w].optParam[i].d !=
                paramWell.well[w].optParam[i + 1].d) {
                tmpOptParam.push_back(paramWell.well[w].optParam[i]);
            }
        }
        tmpOptParam.push_back(paramWell.well[w].optParam.back());

        const USI n = tmpOptParam.size();
        wellOptTime.clear();
        wellOptTime.resize(n + 1);
        for (USI i = 0; i < n; i++) {
            wellOptTime[i] = tmpOptParam[i].d;
        }
        wellOptTime.back() = t;
        for (USI i = 0; i < n; i++) {
            for (USI d = wellOptTime[i]; d < wellOptTime[i + 1]; d++) {
                wells[w].optSet[d] = WellOpt(tmpOptParam[i].opt);
            }
        }
    }

    // for output
    useVTK = output_param.outVTKParam.useVTK;
}

void AllWells::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;
    SetupWell(myGrid, myBulk);
    SetPolyhedronWell(myGrid);
    SetupMixture(myBulk);
}

void AllWells::SetupWell(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        wells[w].Setup(myGrid, myBulk, solvents);
    }
    SetupConnWell2Bulk(myBulk);
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
    USI g    = 1;
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
            if (wells[w].WellType() == INJ)
                wellGroup[glen].wIdINJ.push_back(w);
            else
                wellGroup[glen].wIdPROD.push_back(w);
            glen++;
        }
    }

    numGroup = wellGroup.size();
    // relation between wellGroup should be completed if "node groups" exist

    // control of group should be update according to input file

    // for test
    if (OCP_FALSE) {
        wellGroup[0].reInj      = OCP_TRUE;
        wellGroup[0].saleRate   = 1500;
        wellGroup[0].reInjPhase = GAS;
        wellGroup[0].prodGroup  = 0;
    }
    wellGroup[0].reInjZi.resize(myBulk.GetComNum());
    CalReInjFluid(myBulk);
}

void AllWells::SetupMixture(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    flashCal = myBulk.GetMixture();
}

void AllWells::SetupWellBulk(Bulk& myBulk) const
{
    for (auto& w : wells) {
        if (w.IsOpen()) {
            for (auto& p : w.perf) {
                myBulk.AddWellBulkId(p.Location());
            }
        }
    }
}

void AllWells::SetupConnWell2Bulk(const Bulk& myBulk)
{
    well2bulk.resize(numWell);
    USI wId = 0;
    for (auto& w : wells) {
        for (auto& p : w.perf) well2bulk[wId].push_back(p.Location());
        wId++;
    }
}

void AllWells::ApplyControl(const USI& i)
{
    OCP_FUNCNAME;
    wellChange = OCP_FALSE;
    USI wId    = 0;
    for (USI w = 0; w < numWell; w++) {
        wells[w].opt = wells[w].optSet[i];
        if (wells[w].IsOpen()) {
            wells[w].wOId = wId;
            wId++;
        }
        if (i > 0 && wells[w].opt != wells[w].optSet[i - 1]) wellChange = OCP_TRUE;
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
        if (wells[w].IsOpen()) {

            wells[w].CalTrans(myBulk);
            wells[w].CaldG(myBulk);
            wells[w].CheckOptMode(myBulk);
            wells[w].CalFlux(myBulk, OCP_TRUE);
        }
    }
}

void AllWells::CalTrans(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].IsOpen()) {
            wells[w].CalTrans(myBulk);
        }
    }
}

void AllWells::CalFlux(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].IsOpen()) {
            wells[w].CalFlux(myBulk, OCP_FALSE);
        }
    }
}

void AllWells::CaldG(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].IsOpen()) {
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

        if (wells[w].IsOpen()) {
            if (wells[w].WellType() == PROD) {
                wells[w].CalProdQj(myBulk, dt);
            } else {
                wells[w].CalInjQj(myBulk, dt);
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
            fill(wG.reInjZi.begin(), wG.reInjZi.end(), 0.0);
            for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                Daxpy(nc, 1.0, &wells[prod].qi_lbmol[0], &wG.reInjZi[0]);
            }
            OCP_DBL qt = Dnorm1(nc, &wG.reInjZi[0]);
            if (fabs(qt) < 1E-8) {
                // Recalculate zi
                fill(wG.reInjZi.begin(), wG.reInjZi.end(), 0.0);
                for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
                    wells[prod].CalReInjFluid(myBulk, wG.reInjZi);
                }
                qt = Dnorm1(nc, &wG.reInjZi[0]);
            }
            flashCal[0]->Flash(Psurf, Tsurf, &wG.reInjZi[0]);
            Dcopy(nc, &wG.reInjZi[0], &flashCal[0]->xij[wG.reInjPhase * nc]);
            wG.reInjXi     = flashCal[0]->xi[wG.reInjPhase];
            wG.reInjFactor = wG.reInjXi * flashCal[0]->vj[wG.reInjPhase] / qt;
            // assign to every open injection well in wG
            for (auto& w : wG.wIdINJ) {
                if (wells[w].IsOpen()) {
                    wells[w].opt.reInj      = OCP_TRUE;
                    wells[w].opt.connWell   = wG.wIdPROD;
                    wells[w].opt.reInjPhase = wG.reInjPhase;
                    wells[w].opt.injZi      = wG.reInjZi;
                    wells[w].opt.factorINJ =
                        wG.reInjXi * 1000; // Blackoil Model prohibited
                    wells[w].opt.reInjFactor = wG.reInjFactor;
                    wells[w].opt.maxRate =
                        -wG.saleRate * wells[w].opt.factorINJ; // Mscf -> ft3 -> lbmol
                }
            }
        }
    }
}

void AllWells::ResetBHP()
{
    for (auto& w : wells) {
        if (w.IsOpen()) {
            w.bhp = w.lbhp;
            w.CorrectBHP();
        }
    }
}

OCP_INT AllWells::CheckP(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    OCP_BOOL flagSwitch = OCP_FALSE;
    OCP_BOOL flagCrossf = OCP_FALSE;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].IsOpen()) {

            OCP_INT flag = wells[w].CheckP(myBulk);

            switch (flag) {
                case WELL_NEGATIVE_PRESSURE:
                    return WELL_NEGATIVE_PRESSURE;

                case WELL_SWITCH_TO_BHPMODE:
                    flagSwitch = OCP_TRUE;
                    break;

                case WELL_CROSSFLOW:
                    flagCrossf = OCP_TRUE;
                    break;

                case WELL_SUCCESS:
                default:
                    break;
            }
        }
    }

    if (flagSwitch || flagCrossf) return WELL_SWITCH_TO_BHPMODE;

    return WELL_SUCCESS;
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
        if (wells[w].IsOpen()) {
            dPmax = max(dPmax, fabs(wells[w].bhp - wells[w].lbhp));
        }
    }
}

void AllWells::SetPolyhedronWell(const Grid& myGrid)
{
    if (!useVTK) return;

    wellVal.resize(numWell);
    polyhedronWell.resize(numWell);
    for (USI w = 0; w < numWell; w++) {
        wells[w].SetPolyhedronWell(myGrid, polyhedronWell[w]);
    }
}

void AllWells::SetWellVal() const
{
    if (!useVTK) return;

    for (USI w = 0; w < numWell; w++) {
        if (wells[w].opt.state) {
            if (wells[w].opt.type == INJ) {
                if (wells[w].opt.optMode == BHP_MODE)
                    wellVal[w] = wells[w].opt.maxBHP;
                else
                    wellVal[w] = wells[w].opt.maxRate;
            } else {
                if (wells[w].opt.optMode == BHP_MODE)
                    wellVal[w] = wells[w].opt.minBHP;
                else
                    wellVal[w] = wells[w].opt.maxRate;
            }
        } else {
            wellVal[w] = 0;
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