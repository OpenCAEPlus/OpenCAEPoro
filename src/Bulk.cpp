/*! \file    Bulk.cpp
 *  \brief   Bulk class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include <algorithm>
#include <cmath>
#include <ctime>

#include "Bulk.hpp"

 /////////////////////////////////////////////////////////////////////
 // General
 /////////////////////////////////////////////////////////////////////


void Bulk::InputParam(ParamReservoir& rs_param) { OCP_FUNCNAME;

    rockPref = rs_param.rock.Pref;
    rockC1   = rs_param.rock.Cr;
    rockC2   = rockC1;

    T        = rs_param.rsTemp;
    blackOil = rs_param.blackOil;
    comps    = rs_param.comps;
    oil      = rs_param.oil;
    gas      = rs_param.gas;
    water    = rs_param.water;
    disGas   = rs_param.disGas;

    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    if (rs_param.PBVD_T.data.size() > 0) {
        EQUIL.PBVD.Setup(rs_param.PBVD_T.data[0]);
    }
    
    if (blackOil) {
        if (water && !oil && !gas) {
            // water
            numPhase = 1;
            numCom   = 1;
            SATmode  = PHASE_W;
            PVTmode  = PHASE_W;
        } else if (water && oil && !gas) {
            // water, dead oil
            numPhase    = 2;
            numCom      = 2;
            EQUIL.DOWC  = rs_param.EQUIL[2];
            EQUIL.PcOWC = rs_param.EQUIL[3];
            SATmode     = PHASE_OW;
            PVTmode     = PHASE_OW;
        } else if (water && oil && gas && !disGas) {
            // water, dead oil, dry gas
            numPhase    = 3;
            numCom      = 3;
            EQUIL.DOWC  = rs_param.EQUIL[2];
            EQUIL.PcOWC = rs_param.EQUIL[3];
            EQUIL.DGOC  = rs_param.EQUIL[4];
            EQUIL.PcGOC = rs_param.EQUIL[5];
            SATmode     = PHASE_DOGW;
            PVTmode     = PHASE_DOGW; // maybe it should be added later
        } else if (water && oil && gas && disGas) {
            // water, live oil, dry gas
            numPhase    = 3;
            numCom      = 3;
            EQUIL.DOWC  = rs_param.EQUIL[2];
            EQUIL.PcOWC = rs_param.EQUIL[3];
            EQUIL.DGOC  = rs_param.EQUIL[4];
            EQUIL.PcGOC = rs_param.EQUIL[5];
            SATmode     = PHASE_ODGW;
            PVTmode     = PHASE_ODGW;
        }
        rs_param.numPhase = numPhase;
        rs_param.numCom   = numCom;
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit(rs_param, SATmode, i));
        if (oil & gas & water) {
            for (USI i = 0; i < rs_param.NTSFUN; i++) {
                flow[i]->Generate_SWPCWG();
            }
        }
        for (USI i = 0; i < rs_param.NTPVT; i++)
            flashCal.push_back(new BOMixture(rs_param, PVTmode, i));
        cout << "Bulk::InputParam" << endl;
    } else if (comps) {
        initZi = rs_param.initZi;
    }
}

void Bulk::Setup(const Grid& myGrid) { OCP_FUNCNAME;

    numBulk = myGrid.activeGridNum;
    dx.resize(numBulk, 0);
    dy.resize(numBulk, 0);
    dz.resize(numBulk, 0);
    depth.resize(numBulk, 0);
    ntg.resize(numBulk, 0);
    rockVpInit.resize(numBulk, 0);
    rockVp.resize(numBulk, 0);
    rockKxInit.resize(numBulk, 0);
    rockKyInit.resize(numBulk, 0);
    rockKzInit.resize(numBulk, 0);
    SATNUM.resize(numBulk, 0);
    PVTNUM.resize(numBulk, 0);

    for (OCP_USI bIdb = 0; bIdb < numBulk; bIdb++) {
        OCP_USI bIdg = myGrid.activeMap_B2G[bIdb];

        dx[bIdb]    = myGrid.dx[bIdg];
        dy[bIdb]    = myGrid.dy[bIdg];
        dz[bIdb]    = myGrid.dz[bIdg];
        depth[bIdb] = myGrid.depth[bIdg];
        ntg[bIdb]   = myGrid.ntg[bIdg];

        rockVpInit[bIdb] = myGrid.v[bIdg] * myGrid.ntg[bIdg] * myGrid.poro[bIdg];
        // Rock_PoroInit[bIdb] = myGrid.poro[bIdg];
        rockKxInit[bIdb] = myGrid.kx[bIdg];
        rockKyInit[bIdb] = myGrid.ky[bIdg];
        rockKzInit[bIdb] = myGrid.kz[bIdg];

        SATNUM[bIdb] = myGrid.SATNUM[bIdg];
        PVTNUM[bIdb] = myGrid.PVTNUM[bIdg];
    }
    rockVp = rockVpInit;
    // Rock_Poro = Rock_PoroInit;
    rockKx = rockKxInit;
    rockKy = rockKyInit;
    rockKz = rockKzInit;

    // physical variable
    P.resize(numBulk);
    Pj.resize(numBulk * numPhase);
    Pc.resize(numBulk * numPhase);
    phaseExist.resize(numBulk * numPhase);
    S.resize(numBulk * numPhase);
    rho.resize(numBulk * numPhase);
    xi.resize(numBulk * numPhase);
    cij.resize(numBulk * numPhase * numCom);
    Ni.resize(numBulk * numCom);
    mu.resize(numBulk * numPhase);
    kr.resize(numBulk * numPhase);
    vj.resize(numBulk * numPhase);
    vf.resize(numBulk);

    phase2Index.resize(3);

    if (blackOil) {
        Pb.resize(numBulk);
        switch (PVTmode) {
            case PHASE_W:
                index2Phase.resize(1);
                index2Phase[0] = WATER;

                phase2Index[WATER] = 0;

                break;
            case PHASE_OW:
                index2Phase.resize(2);
                index2Phase[0] = OIL;
                index2Phase[1] = WATER;

                phase2Index[OIL]   = 0;
                phase2Index[WATER] = 1;
                break;
            case PHASE_OG:
                index2Phase.resize(2);
                index2Phase[0] = OIL;
                index2Phase[1] = GAS;

                phase2Index[OIL] = 0;
                phase2Index[GAS] = 1;
                break;
            case PHASE_GW:
                index2Phase.resize(2);
                index2Phase[0] = GAS;
                index2Phase[1] = WATER;

                phase2Index[GAS]   = 0;
                phase2Index[WATER] = 1;
                break;
            case PHASE_ODGW:
            case PHASE_DOGW:
                index2Phase.resize(3);
                index2Phase[0] = OIL;
                index2Phase[1] = GAS;
                index2Phase[2] = WATER;

                phase2Index[OIL]   = 0;
                phase2Index[GAS]   = 1;
                phase2Index[WATER] = 2;
                break;
            default:
                break;
        }
    }
    else if (comps)
    {
        // initZi.resize(numBulk * numCom);
    }

    // CheckSetup();
}


void Bulk::CheckSetup()  const
{
    CheckInitVpore();
    CheckVpore();
}

void Bulk::CheckInitVpore() const
{
    // Check InitVpore
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (rockVpInit[n] < TINY) {
            OCP_ABORT("bulk volume is too small: bulk(" + std::to_string(n) + ") = " +
                std::to_string(rockVpInit[n]));
        }
    }
}


void Bulk::CheckVpore() const
{
    // Check InitVpore
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (rockVp[n] < TINY) {
            OCP_ABORT("bulk volume is too small: bulk(" + std::to_string(n) + ") = " +
                std::to_string(rockVp[n]));
        }
    }
}



void Bulk::InitSjPcBo(const USI& tabrow) { OCP_FUNCNAME;

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOWC = EQUIL.PcOWC;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGOC = EQUIL.PcGOC;

    OCP_DBL Zmin = 1E8;
    OCP_DBL Zmax = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_DBL temp1 = depth[n] - dz[n] / 2;
        OCP_DBL temp2 = depth[n] + dz[n] / 2;
        Zmin = Zmin < temp1 ? Zmin : temp1;
        Zmax = Zmax > temp2 ? Zmax : temp2;
    }
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    // creater table
    OCPTable         DepthP(tabrow, 4);
    vector<OCP_DBL>& Ztmp = DepthP.GetCol(0);
    vector<OCP_DBL>& Potmp = DepthP.GetCol(1);
    vector<OCP_DBL>& Pgtmp = DepthP.GetCol(2);
    vector<OCP_DBL>& Pwtmp = DepthP.GetCol(3);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++) {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    USI beginId = 0;
    // find the RefId
    if (Dref <= Ztmp[0]) {
        beginId = 0;
    }
    else if (Dref >= Ztmp[tabrow - 1]) {
        beginId = tabrow - 1;
    }
    else {
        beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
            [s = Dref](auto t) { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure:
    OCP_DBL Pbb = Pref;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI     mynum = 10;
    OCP_DBL mydz = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    if (Dref < DOGC) {
        // reference pressure is gas pressure
        if (flow[0]->IsEmpty_SGOF()) {
            OCP_ABORT("SGOF is missing !");
        }

        Pgref = Pref;
        gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--) {
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref = 0;
        Ptmp = Pgref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
            Ptmp += gammaGtmp * mydz;
        }
        Ptmp -= PcGOC;
        for (USI i = 0; i < mynum; i++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, DOGC - i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (!EQUIL.PBVD.IsEmpty()) {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOWC;
        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else if (Dref > DOWC) {
        // reference pressure is water pressure
        Pwref = Pref;
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp = Pwref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp += gammaWtmp * mydz;
        }
        Ptmp += PcOWC;

        for (USI i = 0; i < mynum; i++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, DOWC - i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (!EQUIL.PBVD.IsEmpty()) {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (!flow[0]->IsEmpty_SGOF()) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;

            for (USI i = 0; i < mynum; i++) {
                if (!EQUIL.PBVD.IsEmpty()) {
                    Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
                }
                gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
                Ptmp += gammaOtmp * mydz;
            }
            Ptmp += PcGOC;
            for (USI i = 0; i < mynum; i++) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
                Ptmp -= gammaGtmp * mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }
            for (USI id = beginId; id < tabrow - 1; id++) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }
    }
    else {
        // reference pressure is oil pressure
        Poref = Pref;
        if (!EQUIL.PBVD.IsEmpty()) {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure in tab
        for (USI id = beginId; id > 0; id--) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (!flow[0]->IsEmpty_SGOF()) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;

            for (USI i = 0; i < mynum; i++) {
                if (!EQUIL.PBVD.IsEmpty()) {
                    Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
                }
                gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
                Ptmp += gammaOtmp * mydz;
            }
            Ptmp += PcGOC;
            for (USI i = 0; i < mynum; i++) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
                Ptmp -= gammaGtmp * mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }

            for (USI id = beginId; id < tabrow - 1; id++) {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            if (!EQUIL.PBVD.IsEmpty()) {
                Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOWC;
        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0);
    std::vector<OCP_DBL> cdata(4, 0);
    for (OCP_USI n = 0; n < numBulk; n++) {
        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po = data[1];
        OCP_DBL Pg = data[2];
        OCP_DBL Pw = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw = flow[0]->EvalInv_SWOF(3, Pcow, 0);
        OCP_DBL Sg = 0;
        if (!flow[0]->IsEmpty_SGOF()) {
            Sg = flow[0]->Eval_SGOF(3, Pcgo, 0);
        }
        if (Sw + Sg > 1) {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw = flow[0]->EvalInv_SWPCWG(1, Pcgw, 0);
            Sg = 1 - Sw;
        }

        if (1 - Sw < TINY) {
            // all water
            Po = Pw + flow[0]->Eval_SWOF(0, 1.0, 3);
            // Pg = Po + flow->Eval_SGOF(0, 0.0, 3);
        }
        else if (1 - Sg < TINY) {
            // all gas
            Po = Pg - flow[0]->Eval_SGOF(0, 1.0, 3);
            // Pw = Po - flow->Eval_SWOF(0, 0.0, 3);
        }
        else if (1 - Sw - Sg < TINY) {
            // water and gas
            Po = Pg - flow[0]->Eval_SGOF(0, Sg, 3);
            // Pw = Po - flow->Eval_SWOF(0, Sw, 3);
        }
        P[n] = Po;

        if (depth[n] < DOGC) {
            Pbb = Po;
        }
        else if (!EQUIL.PBVD.IsEmpty()) {
            Pbb = EQUIL.PBVD.Eval(0, depth[n], 1);
        }
        Pb[n] = Pbb;

        // cal Sg and Sw
        Sw = 0;
        Sg = 0;
        USI ncut = 10;

        for (USI k = 0; k < ncut; k++) {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep = depth[n] + dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po = data[1];
            Pg = data[2];
            Pw = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            tmpSw = flow[0]->EvalInv_SWOF(3, Pcow, 0);
            if (!flow[0]->IsEmpty_SGOF()) {
                tmpSg = flow[0]->Eval_SGOF(3, Pcgo, 0);
            }
            if (tmpSw + tmpSg > 1) {
                // should me modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw = flow[0]->EvalInv_SWPCWG(1, Pcgw, 0);
                tmpSg = 1 - tmpSw;
            }
            Sw += tmpSw;
            Sg += tmpSg;
        }
        Sw /= ncut;
        Sg /= ncut;
        S[n * numPhase + numPhase - 1] = Sw;
        if (!flow[0]->IsEmpty_SGOF()) {
            S[n * numPhase + numPhase - 2] = Sg;
        }
    }
}


void Bulk::InitSjPcComp(const USI& tabrow) { OCP_FUNCNAME;

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOWC = EQUIL.PcOWC;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGOC = EQUIL.PcGOC;

    OCP_DBL Zmin = 1E8;
    OCP_DBL Zmax = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_DBL temp1 = depth[n] - dz[n] / 2;
        OCP_DBL temp2 = depth[n] + dz[n] / 2;
        Zmin = Zmin < temp1 ? Zmin : temp1;
        Zmax = Zmax > temp2 ? Zmax : temp2;
    }
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    // creater table
    OCPTable         DepthP(tabrow, 4);
    vector<OCP_DBL>& Ztmp = DepthP.GetCol(0);
    vector<OCP_DBL>& Potmp = DepthP.GetCol(1);
    vector<OCP_DBL>& Pgtmp = DepthP.GetCol(2);
    vector<OCP_DBL>& Pwtmp = DepthP.GetCol(3);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++) {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    USI beginId = 0;
    // find the RefId
    if (Dref <= Ztmp[0]) {
        beginId = 0;
    }
    else if (Dref >= Ztmp[tabrow - 1]) {
        beginId = tabrow - 1;
    }
    else {
        beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
            [s = Dref](auto t) { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure:
    OCP_DBL mytemp = T;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI     mynum = 10;
    OCP_DBL mydz = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    if (Dref < DOGC) {
        // reference pressure is gas pressure
        Pgref = Pref;
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref = 0;
        Ptmp = Pgref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaGtmp * mydz;
        }
        Ptmp -= PcGOC;
        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOWC;
        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else if (Dref > DOWC) {
        // reference pressure is water pressure
        Pwref = Pref;
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp = Pwref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp += gammaWtmp * mydz;
        }
        Ptmp += PcOWC;

        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the gas pressure in Dref by Poref
        Pgref = 0;
        Ptmp = Poref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp += PcGOC;
        for (USI i = 0; i < mynum; i++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaGtmp * mydz;
        }
        Pgref = Ptmp;

        // find the gas pressure in tab
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else {
        // reference pressure is oil pressure
        Poref = Pref;
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure
        for (USI id = beginId; id > 0; id--) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the gas pressure in Dref by Poref
        Pgref = 0;
        Ptmp = Poref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp += PcGOC;
        for (USI i = 0; i < mynum; i++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaGtmp * mydz;
        }
        Pgref = Ptmp;

        // find the gas pressure in tab
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++) {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOWC;
        for (USI i = 0; i < mynum; i++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0);
    std::vector<OCP_DBL> cdata(4, 0);

    for (OCP_USI n = 0; n < numBulk; n++) {
        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po = data[1];
        OCP_DBL Pg = data[2];
        OCP_DBL Pw = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw = flow[0]->EvalInv_SWOF(3, Pcow, 0);
        OCP_DBL Sg = 0;
        if (!flow[0]->IsEmpty_SGOF()) {
            Sg = flow[0]->Eval_SGOF(3, Pcgo, 0);
        }
        if (Sw + Sg > 1) {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw = flow[0]->EvalInv_SWPCWG(1, Pcgw, 0);
            Sg = 1 - Sw;
        }

        if (1 - Sw < TINY) {
            // all water
            Po = Pw + flow[0]->Eval_SWOF(0, 1.0, 3);
            // Pg = Po + flow->Eval_SGOF(0, 0.0, 3);
        }
        else if (1 - Sg < TINY) {
            // all gas
            Po = Pg - flow[0]->Eval_SGOF(0, 1.0, 3);
            // Pw = Po - flow->Eval_SWOF(0, 0.0, 3);
        }
        else if (1 - Sw - Sg < TINY) {
            // water and gas
            Po = Pg - flow[0]->Eval_SGOF(0, Sg, 3);
            // Pw = Po - flow->Eval_SWOF(0, Sw, 3);
        }
        P[n] = Po;

        // cal Sg and Sw
        Sw = 0;
        Sg = 0;
        USI ncut = 10;

        for (USI k = 0; k < ncut; k++) {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep = depth[n] + dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po = data[1];
            Pg = data[2];
            Pw = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            tmpSw = flow[0]->EvalInv_SWOF(3, Pcow, 0);
            if (!flow[0]->IsEmpty_SGOF()) {
                tmpSg = flow[0]->Eval_SGOF(3, Pcgo, 0);
            }
            if (tmpSw + tmpSg > 1) {
                // should me modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw = flow[0]->EvalInv_SWPCWG(1, Pcgw, 0);
                tmpSg = 1 - tmpSw;
            }
            Sw += tmpSw;
            Sg += tmpSg;
        }
        Sw /= ncut;
        Sg /= ncut;

        S[n * numPhase + numPhase - 1] = Sw;
        if (!flow[0]->IsEmpty_SGOF()) {
            S[n * numPhase + numPhase - 2] = Sg;
        }
    }
}

// Flash
void Bulk::FlashSj() { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->Flash_Sj(P[n], Pb[n], T, &S[n * numPhase], rockVp[n],
            initZi.data());
        for (USI i = 0; i < numCom; i++) {
            Ni[n * numCom + i] = flashCal[PVTNUM[n]]->Ni[i];
        }
        PassFlashValue(n);
    }
#ifdef _DEBUG
    CheckSat();
#endif // _DEBUG
}


void Bulk::FlashNi() { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->Flash_Ni(P[n], T, &Ni[n * numCom]);
        PassFlashValue(n);
    }
#ifdef _DEBUG
    CheckSat();
#endif // _DEBUG
}


void Bulk::FlashNiDeriv() { OCP_FUNCNAME;

    dSec_dPri.clear();
    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->Flash_Ni_Deriv(P[n], T, &Ni[n * numCom]);
        PassFlashValueDeriv(n);
    }
#ifdef _DEBUG
    CheckSat();
#endif // _DEBUG
}


void Bulk::PassFlashValue(const OCP_USI& n) { OCP_FUNCNAME;

    OCP_USI bId = n * numPhase;
    USI     pvtnum = PVTNUM[n];
    for (USI j = 0; j < numPhase; j++) {
        phaseExist[bId + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be pass no matter if the phase exists.
        // Because it will be used to calculate relative permeability and capillary
        // pressure at every time step, be sure all saturation are updated at every
        // step.
        S[bId + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bId + j]) { // j -> bId + j   fix bugs.
            rho[bId + j] = flashCal[pvtnum]->rho[j];
            xi[bId + j] = flashCal[pvtnum]->xi[j];
            for (USI i = 0; i < numCom; i++) {
                cij[bId * numCom + j * numCom + i] =
                    flashCal[pvtnum]->cij[j * numCom + i];
            }
            mu[bId + j] = flashCal[pvtnum]->mu[j];
            vj[bId + j] = flashCal[pvtnum]->v[j];
        }
    }
    vf[n] = flashCal[pvtnum]->vf;
    vfp[n] = flashCal[pvtnum]->vfp;
    bId = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bId + i] = flashCal[pvtnum]->vfi[i];
    }
}


void Bulk::PassFlashValueDeriv(const OCP_USI& n) { OCP_FUNCNAME;

    OCP_USI bId = n * numPhase;
    USI     pvtnum = PVTNUM[n];
    for (USI j = 0; j < numPhase; j++) {
        phaseExist[bId + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be pass no matter if the phase exists.
        // Because it will be used to calculate relative permeability and capillary
        // pressure at every time step, be sure all saturation are updated at every
        // step.
        S[bId + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bId + j]) { // j -> bId + j   fix bugs.
            rho[bId + j] = flashCal[pvtnum]->rho[j];
            xi[bId + j] = flashCal[pvtnum]->xi[j];
            mu[bId + j] = flashCal[pvtnum]->mu[j];
            // vj[bId + j] = flashCal[pvtnum]->v[j];

            // Derivatives
            muP[bId + j] = flashCal[pvtnum]->muP[j];
            xiP[bId + j] = flashCal[pvtnum]->xiP[j];
            rhoP[bId + j] = flashCal[pvtnum]->rhoP[j];
            for (USI i = 0; i < numCom; i++) {
                cij[bId * numCom + j * numCom + i] =
                    flashCal[pvtnum]->cij[j * numCom + i];
                muC[bId * numCom + j * numCom + i] =
                    flashCal[pvtnum]->muC[j * numCom + i];
                xiC[bId * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xiC[j * numCom + i];
                rhoC[bId * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhoC[j * numCom + i];
            }
        }
    }
    vf[n] = flashCal[pvtnum]->vf;
    Nt[n] = flashCal[pvtnum]->Nt;
    vfp[n] = flashCal[pvtnum]->vfp;
    bId = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bId + i] = flashCal[pvtnum]->vfi[i];
    }
    dSec_dPri.insert(dSec_dPri.end(), flashCal[pvtnum]->dSec_dPri.begin(),
        flashCal[pvtnum]->dSec_dPri.end());
}

void Bulk::ResetFlash() {  OCP_FUNCNAME;

    phaseExist = lphaseExist;
    S = lS;
    rho = lrho;
    xi = lxi;
    cij = lcij;
    mu = lmu;
    vj = lvj;
    vf = lvf;
    vfp = lvfp;
    vfi = lvfi;
}


// relative permeability and capillary pressure
void Bulk::CalKrPc() { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId]);
        for (USI j = 0; j < numPhase; j++)
            Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
    }
}


void Bulk::CalKrPcDeriv() { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPcDeriv(&S[bId], &kr[bId], &Pc[bId],
            &dKr_dS[bId * numPhase],
            &dPcj_dS[bId * numPhase]);
        for (USI j = 0; j < numPhase; j++)
            Pj[bId + j] = P[n] + Pc[bId + j];
    }
}


void Bulk::CalVpore() { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_DBL dP = rockC1 * (P[n] - rockPref);
        // rockVp[n] = rockVpInit[n] * (1 + dP + dP * dP / 2);
        rockVp[n] = rockVpInit[n] * (1 + dP);
    }
}

OCP_DBL Bulk::CalFPR() const { OCP_FUNCNAME;

    OCP_DBL ptmp = 0;
    OCP_DBL vtmp = 0;
    OCP_DBL tmp = 0;

    if (numPhase == 3) {
        for (OCP_USI n = 0; n < numBulk; n++) {
            tmp = rockVp[n] * (1 - S[n * numPhase + 2]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    }
    else if (numPhase < 3) {
        for (OCP_USI n = 0; n < numBulk; n++) {
            tmp = rockVp[n] * (S[n * numPhase]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    }
    else {
        OCP_ABORT("Number of phases is out of range!");
    }
    return ptmp / vtmp;
}


void Bulk::CalMaxChange() { OCP_FUNCNAME;

    dPmax = 0;
    dNmax = 0;
    dSmax = 0;
    dVmax = 0;
    OCP_DBL tmp = 0;
    OCP_USI id;

    for (OCP_USI n = 0; n < numBulk; n++) {

        // dP
        tmp = fabs(P[n] - lP[n]);
        dPmax = dPmax < tmp ? tmp : dPmax;

        // dS
        for (USI j = 0; j < numPhase; j++) {
            id = n * numPhase + j;
            tmp = fabs(S[id] - lS[id]);
            dSmax = dSmax < tmp ? tmp : dSmax;
        }

        // dN
        for (USI i = 0; i < numCom; i++) {
            id = n * numCom + i;

            tmp = fabs(max(Ni[id], lNi[id]));
            if (tmp > TINY) {
                tmp = fabs(Ni[id] - lNi[id]) / tmp;
                dNmax = dNmax < tmp ? tmp : dNmax;
            }
        }

        tmp = fabs(vf[n] - rockVp[n]) / rockVp[n];
        dVmax = dVmax < tmp ? tmp : dVmax;
    }
}


/// Return true if no negative pressure and false otherwise.
bool Bulk::CheckP() const { OCP_FUNCNAME;

    for (auto p : P) {
        if (p < 0.0) return false;
    }
    return true;
}


/// Return true if no negative Ni and false otherwise.
bool Bulk::CheckNi() const { OCP_FUNCNAME;

    for (auto ni : Ni) {
        if (ni < 0.0) {
            cout << "###WARNING: Negative Ni  " << ni << endl;
            return false;
        }         
    }
    return true;
}


/// Return true if all Ve < Vlim and false otherwise.
bool Bulk::CheckVe(const OCP_DBL& Vlim) const { OCP_FUNCNAME;

    OCP_DBL tmp = 0.0;
    for (OCP_USI n = 0; n < numBulk; n++) {
        tmp = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (tmp > Vlim) {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << tmp
                << " is too big!" << endl;
            return false;
        }
    }
    return true;
}


void Bulk::CheckDiff() { OCP_FUNCNAME;

    OCP_DBL tmp;
    OCP_USI id;
    for (OCP_USI n = 0; n < numBulk; n++) {
        for (USI j = 0; j < numPhase; j++) {
            id = n * numPhase + j;
            tmp = fabs(phaseExist[id] - lphaseExist[id]);
            if (tmp != 0.0) {
                cout << "Difference in phaseExist\t" << tmp << "\n";
            }
            if (lphaseExist[id] || phaseExist[id]) {
                tmp = fabs(S[id] - lS[id]);
                if (tmp != 0.0) {
                    cout << "Difference in S\t" << tmp << "  " << phaseExist[id]
                        << "\n";
                }
                tmp = fabs(xi[id] - lxi[id]);
                if (tmp != 0.0) {
                    cout << "Difference in Xi\t" << tmp << "  " << phaseExist[id]
                        << "\n";
                }
                tmp = fabs(rho[id] - lrho[id]);
                if (tmp != 0.0) {
                    cout << "Difference in rho\t" << tmp << "  " << phaseExist[id]
                        << "\n";
                }
            }
        }
    }
}


void Bulk::CheckSat() const { OCP_FUNCNAME;

    OCP_DBL tmp;
    for (OCP_USI n = 0; n < numBulk; n++) {
        tmp = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (S[n * numPhase + j] < 0) {
                OCP_ABORT("Negative Volume!");
            }
            tmp += S[n * numPhase + j];
        }
        if (fabs(tmp - 1) > TINY) {
            OCP_ABORT("Saturation is greater than 1");
        }
    }
}


USI Bulk::GetMixMode() const { OCP_FUNCNAME;

    if (blackOil)
        return BLKOIL;
    else if (comps)
        return EoS_PVTW;
    else
        return 0; // TODO: Make sure code does not reach here!
}



 /////////////////////////////////////////////////////////////////////
 // IMPEC
 /////////////////////////////////////////////////////////////////////


void Bulk::AllocateAuxIMPEC() { OCP_FUNCNAME;

    vfi.resize(numBulk * numCom);
    vfp.resize(numBulk);
    cfl.resize(numBulk * numPhase);

    // for last step
    lP.resize(numBulk);
    lPj.resize(numBulk * numPhase);
    lPc.resize(numBulk * numPhase);
    lphaseExist.resize(numBulk * numPhase);
    lS.resize(numBulk * numPhase);
    lrho.resize(numBulk * numPhase);
    lxi.resize(numBulk * numPhase);
    lcij.resize(numBulk * numPhase * numCom);
    lNi.resize(numBulk * numCom);
    lmu.resize(numBulk * numPhase);
    lkr.resize(numBulk * numPhase);
    lvj.resize(numBulk * numPhase);
    lvf.resize(numBulk);
    lvfi.resize(numBulk * numCom);
    lvfp.resize(numBulk);
    rockLVp.resize(numBulk);
}


void Bulk::GetSolIMPEC(const vector<OCP_DBL>& u) { OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        P[n] = u[n];
        for (USI j = 0; j < numPhase; j++) {
            OCP_USI id = n * numPhase + j;
            if (phaseExist[id]) {
                Pj[id] = P[n] + Pc[id];
            }
        }
    }
}


OCP_DBL Bulk::CalCFL01IMPEC() const { OCP_FUNCNAME;

    OCP_DBL tmp = 0;
    OCP_USI id;
    for (OCP_USI n = 0; n < numBulk; n++) {
        for (USI j = 0; j < numPhase; j++) {
            id = n * numPhase + j;
            if (phaseExist[id]) {

                if (vj[id] <= 0) continue;  // temp

                cfl[id] /= vj[id];
#ifdef _DEBUG
                if (!isfinite(cfl[id])) {
                    OCP_ABORT("cfl is nan!");
                }
#endif // _DEBUG

                if (tmp < cfl[id]) tmp = cfl[id];
            }
        }
    }
    return tmp;
}


void Bulk::UpdateLastStepIMPEC() { OCP_FUNCNAME;

    lP = P;
    lPj = Pj;
    lPc = Pc;
    lphaseExist = phaseExist;
    lS = S;
    lrho = rho;
    lxi = xi;
    lcij = cij;
    lNi = Ni;
    lmu = mu;
    lkr = kr;
    lvj = vj;
    lvf = vf;
    lvfi = vfi;
    lvfp = vfp;
    rockLVp = rockVp;
}


 /////////////////////////////////////////////////////////////////////
 // FIM
 /////////////////////////////////////////////////////////////////////


void Bulk::AllocateAuxFIM() 
{ 
    OCP_FUNCNAME;

    lP.resize(numBulk);
    lS.resize(numBulk * numPhase);
    lNi.resize(numBulk * numCom);

    Nt.resize(numBulk);
    vfi.resize(numBulk * numCom);
    vfp.resize(numBulk);
    muP.resize(numBulk * numPhase);
    xiP.resize(numBulk * numPhase);
    rhoP.resize(numBulk * numPhase);
    muC.resize(numBulk * numCom * numPhase);
    xiC.resize(numBulk * numCom * numPhase);
    rhoC.resize(numBulk * numCom * numPhase);
    dSec_dPri.resize(numBulk * (numCom + 1) * (numCom + 1) * numPhase);
    dKr_dS.resize(numBulk * numPhase * numPhase);
    dPcj_dS.resize(numBulk * numPhase * numPhase);
}


void Bulk::GetSolFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim, const OCP_DBL& dSmaxlim) {  OCP_FUNCNAME;
    
    NRdSmax = 0;
    NRdPmax = 0;
    OCP_DBL dP;
    USI row = numPhase * (numCom + 1);
    USI col = numCom + 1;
    USI bsize = row * col;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL chopmin = 1;
    OCP_DBL choptmp = 0;
    for (OCP_USI n = 0; n < numBulk; n++) {

        chopmin = 1;

        // compute the chop
        dtmp.assign(row, 0);
        DaAxpby(row, col, 1, dSec_dPri.data() + n * bsize, u.data() + n * col, 1, dtmp.data());

        for (USI j = 0; j < numPhase; j++) {
            if (fabs(dtmp[j]) > dSmaxlim) {
                choptmp = dSmaxlim / fabs(dtmp[j]);
            }
            else if (S[n * numPhase + j] + dtmp[j] < 0) {
                choptmp = 0.9 * S[n * numPhase + j] / fabs(dtmp[j]);
            }
            else {
                choptmp = 1;
            }
            chopmin = min(chopmin, choptmp);
            NRdSmax = max(NRdSmax, choptmp * fabs(dtmp[j]));
        }
        dP = u[n * col];
        choptmp = dPmaxlim / fabs(dP);
        chopmin = min(chopmin, choptmp);
        NRdPmax = max(NRdPmax, fabs(dP));
        P[n] += dP;

        for (USI i = 0; i < numCom; i++) {
            Ni[n * numCom + i] += u[n * col + 1 + i] * chopmin;
        }
    }
}


void Bulk::CalRelResFIM(ResFIM& resFIM) const 
{ 
    OCP_FUNCNAME;

    OCP_USI tmpN;
    OCP_DBL tmp1, tmp2, tmp3, tmp4;
    tmp1 = tmp2 = tmp3 = tmp4 = -3.141592653;
    OCP_DBL tmp0 = 0;
    // CheckVpore();

    OCP_DBL tmp;
    const USI len = numCom + 1;
    for (OCP_USI n = 0; n < numBulk; n++) {
        for (USI i = 0; i < len; i++) {
            tmp = fabs(resFIM.res[n * len + i] / rockVp[n]);
            //if (!isfinite(tmp)) {
            //    OCP_ABORT("tmp is nan!");
            //}
            if (resFIM.maxRelRes_v < tmp) {
                resFIM.maxRelRes_v = tmp;
            }
        }
        for (USI i = 0; i < numCom; i++) {
            tmp = fabs(resFIM.res[n * len + 1 + i] / Nt[n]);
            //if (!isfinite(tmp)) {
            //    OCP_ABORT("tmp is nan!");
            //}
            if (resFIM.maxRelRes_mol < tmp) {
                resFIM.maxRelRes_mol = tmp;
            }
        }
    }
}


void Bulk::ResetFIM()
{
    OCP_FUNCNAME;

    P = lP;
    Ni = lNi;
    FlashNiDeriv();
    CalVpore();
    CalKrPcDeriv();
}



void Bulk::UpdateLastStepFIM() { OCP_FUNCNAME;

    lP = P;
    lS = S;
    lNi = Ni; 
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
