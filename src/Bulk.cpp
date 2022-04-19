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

/// Read parameters from rs_param data strucutre.
void Bulk::InputParam(ParamReservoir &rs_param)
{
    OCP_FUNCNAME;

    rockPref = rs_param.rock.Pref;
    rockC1 = rs_param.rock.Cr;
    rockC2 = rockC1;
    T = rs_param.rsTemp + 460;
    ScalePcow = rs_param.ScalePcow;
    blackOil = rs_param.blackOil;
    comps = rs_param.comps;
    oil = rs_param.oil;
    gas = rs_param.gas;
    water = rs_param.water;
    disGas = rs_param.disGas;
    miscible = rs_param.EoSp.miscible;
    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    

    if (rs_param.PBVD_T.data.size() > 0)
        EQUIL.PBVD.Setup(rs_param.PBVD_T.data[0]);

    if (blackOil)
    {

        if (water && !oil && !gas)
        {
            // water
            numPhase = 1;
            numCom = 1;
            SATmode = PHASE_W;
            PVTmode = PHASE_W;
        }
        else if (water && oil && !gas)
        {
            // water, dead oil
            numPhase = 2;
            numCom = 2;
            EQUIL.DOWC = rs_param.EQUIL[2];
            EQUIL.PcOW = rs_param.EQUIL[3];
            SATmode = PHASE_OW;
            PVTmode = PHASE_OW;
        }
        else if (water && oil && gas && !disGas)
        {
            // water, dead oil, dry gas
            numPhase = 3;
            numCom = 3;
            EQUIL.DOWC = rs_param.EQUIL[2];
            EQUIL.PcOW = rs_param.EQUIL[3];
            EQUIL.DGOC = rs_param.EQUIL[4];
            EQUIL.PcGO = rs_param.EQUIL[5];
            SATmode = PHASE_DOGW;
            PVTmode = PHASE_DOGW; // maybe it should be added later
        }
        else if (water && oil && gas && disGas)
        {
            // water, live oil, dry gas
            numPhase = 3;
            numCom = 3;
            EQUIL.DOWC = rs_param.EQUIL[2];
            EQUIL.PcOW = rs_param.EQUIL[3];
            EQUIL.DGOC = rs_param.EQUIL[4];
            EQUIL.PcGO = rs_param.EQUIL[5];
            PVTmode = PHASE_ODGW;

            if (rs_param.SOF3_T.data.size() > 0)
            {
                SATmode = PHASE_ODGW02;
            }
            else
            {
                SATmode = PHASE_ODGW01;
            }
        }
        numCom_1 = numCom - 1;
        rs_param.numPhase = numPhase;
        rs_param.numCom = numCom;

        switch (PVTmode)
        {
        case PHASE_W:
            OCP_ABORT("Wrong Type!");
            break;
        case PHASE_OW:
            for (USI i = 0; i < rs_param.NTPVT; i++)
                flashCal.push_back(new BOMixture_OW(rs_param, i));
            break;
        case PHASE_DOGW:
            OCP_ABORT("Wrong Type!");
            break;
        case PHASE_ODGW:
            for (USI i = 0; i < rs_param.NTPVT; i++)
                flashCal.push_back(new BOMixture_ODGW(rs_param, i));
            break;
        default:
            OCP_ABORT("Wrong Type!");
            break;
        }

        cout << "Bulk::InputParam --- BLACKOIL" << endl;
    }
    else if (comps)
    {

        // Water exists and is excluded in EoS model NOW!
        oil = true;
        gas = true;
        water = true;

        initZi = rs_param.EoSp.zi; // Easy initialization!
        initZi.push_back(0.0); // include water

        numPhase = rs_param.EoSp.numPhase + 1;
        numCom = rs_param.EoSp.numComp + 1;
        numCom_1 = numCom - 1;
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];

        if (rs_param.SOF3_T.data.size() > 0)
        {
            SATmode = PHASE_ODGW02;
        }
        else
        {
            SATmode = PHASE_ODGW01;
        }

        for (USI i = 0; i < rs_param.NTPVT; i++)
            flashCal.push_back(new MixtureComp(rs_param, i));

        cout << "Bulk::InputParam --- COMPOSITIONAL" << endl;
    }

    if (SATmode == PHASE_ODGW01 && miscible) {
        SATmode = PHASE_ODGW01_MISCIBLE;
    }

    switch (SATmode)
    {
    case PHASE_W:
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit_W(rs_param, i));
        break;
    case PHASE_OW:
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit_OW(rs_param, i));
        break;
    case PHASE_ODGW01:
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit_ODGW01(rs_param, i));
        break;
    case PHASE_ODGW01_MISCIBLE:
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit_ODGW01_Miscible(rs_param, i));
        break;
    case PHASE_ODGW02:
        for (USI i = 0; i < rs_param.NTSFUN; i++)
            flow.push_back(new FlowUnit_ODGW02(rs_param, i));
        break;
    default:
        OCP_ABORT("Wrong Type!");
        break;
    }
}

/// Setup bulk information.
void Bulk::Setup(const Grid &myGrid)
{
    OCP_FUNCNAME;

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

    if (myGrid.SwatInit.size() > 0) {
        SwatInitExist = true;
        SwatInit.resize(numBulk);
    }
    if (ScalePcow) {
        ScaleValuePcow.resize(numBulk, 0);
    }
        
    for (OCP_USI bIdb = 0; bIdb < numBulk; bIdb++)
    {
        OCP_USI bIdg = myGrid.activeMap_B2G[bIdb];

        dx[bIdb] = myGrid.dx[bIdg];
        dy[bIdb] = myGrid.dy[bIdg];
        dz[bIdb] = myGrid.dz[bIdg];
        depth[bIdb] = myGrid.depth[bIdg];
        ntg[bIdb] = myGrid.ntg[bIdg];

        rockVpInit[bIdb] = myGrid.v[bIdg] * myGrid.ntg[bIdg] * myGrid.poro[bIdg];
        rockKxInit[bIdb] = myGrid.kx[bIdg];
        rockKyInit[bIdb] = myGrid.ky[bIdg];
        rockKzInit[bIdb] = myGrid.kz[bIdg];

        SATNUM[bIdb] = myGrid.SATNUM[bIdg];
        PVTNUM[bIdb] = myGrid.PVTNUM[bIdg];

        if (SwatInitExist) {
            SwatInit[bIdb] = myGrid.SwatInit[bIdg];
        }
    }

    rockVp = rockVpInit;
    rockKx = rockKxInit;
    rockKy = rockKyInit;
    rockKz = rockKzInit;

    // physical variables
    P.resize(numBulk);
    Pb.resize(numBulk);
    Pj.resize(numBulk * numPhase);
    Pc.resize(numBulk * numPhase);
    phaseExist.resize(numBulk * numPhase);
    S.resize(numBulk * numPhase);
    rho.resize(numBulk * numPhase);
    xi.resize(numBulk * numPhase);
    xij.resize(numBulk * numPhase * numCom);
    Ni.resize(numBulk * numCom);
    mu.resize(numBulk * numPhase);
    kr.resize(numBulk * numPhase);
    vj.resize(numBulk * numPhase);
    vf.resize(numBulk);
    Nt.resize(numBulk);

    phase2Index.resize(3);

    if (blackOil)
    {
        switch (PVTmode)
        {
        case PHASE_W:
            index2Phase.resize(1);
            index2Phase[0] = WATER;
            phase2Index[WATER] = 0;
            break;
        case PHASE_OW:
            index2Phase.resize(2);
            index2Phase[0] = OIL;
            index2Phase[1] = WATER;
            phase2Index[OIL] = 0;
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
            phase2Index[GAS] = 0;
            phase2Index[WATER] = 1;
            break;
        case PHASE_ODGW:
        case PHASE_DOGW:
            index2Phase.resize(3);
            index2Phase[0] = OIL;
            index2Phase[1] = GAS;
            index2Phase[2] = WATER;
            phase2Index[OIL] = 0;
            phase2Index[GAS] = 1;
            phase2Index[WATER] = 2;
            break;
        default:
            OCP_ABORT("Unknown PVT model!");
        }
    }
    else if (comps)
    {
        phase2Index[OIL] = 0;
        phase2Index[GAS] = 1;
        phase2Index[WATER] = 2;

        phaseNum.resize(numBulk);
        lphaseNum.resize(numBulk);
        minEigenSkip.resize(numBulk);
        flagSkip.resize(numBulk);
        ziSkip.resize(numBulk* numCom);
        PSkip.resize(numBulk);
        lminEigenSkip.resize(numBulk);
        lflagSkip.resize(numBulk);
        lziSkip.resize(numBulk* numCom);
        lPSkip.resize(numBulk);
        Ks.resize(numBulk* (numCom - 1));
        lKs.resize(numBulk* (numCom - 1));

        if (miscible) {
            surTen.resize(numBulk);
            Fk.resize(numBulk);
            Fp.resize(numBulk);
        }
    }

    CalSomeInfo(myGrid);

#if DEBUG
    CheckSetup();
#endif
}

void Bulk::CheckSetup() const
{
    CheckInitVpore();
    CheckVpore();
}

/// Check init pore volume.
void Bulk::CheckInitVpore() const
{
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (rockVpInit[n] < TINY)
        {
            OCP_ABORT("Bulk volume is too small: bulk[" + std::to_string(n) +
                      "] = " + std::to_string(rockVpInit[n]));
        }
    }
}

// Check pore volume.
void Bulk::CheckVpore() const
{
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (rockVp[n] < TINY)
        {
            OCP_ABORT("Bulk volume is too small: bulk[" + std::to_string(n) +
                      "] = " + std::to_string(rockVp[n]));
        }
    }
}

/// Here tabrow is maximum number of depth nodes in table of depth vs pressure.
void Bulk::InitSjPcBo(const USI &tabrow)
{
    OCP_FUNCNAME;

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOW = EQUIL.PcOW;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGO = EQUIL.PcGO;
    OCP_DBL Zmin = 1E8;
    OCP_DBL Zmax = 0;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        OCP_DBL temp1 = depth[n] - dz[n] / 2;
        OCP_DBL temp2 = depth[n] + dz[n] / 2;
        Zmin = Zmin < temp1 ? Zmin : temp1;
        Zmax = Zmax > temp2 ? Zmax : temp2;
    }
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    // creater table
    OCPTable DepthP(tabrow, 4);
    vector<OCP_DBL> &Ztmp = DepthP.GetCol(0);
    vector<OCP_DBL> &Potmp = DepthP.GetCol(1);
    vector<OCP_DBL> &Pgtmp = DepthP.GetCol(2);
    vector<OCP_DBL> &Pwtmp = DepthP.GetCol(3);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++)
    {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    // find the RefId
    USI beginId = 0;
    if (Dref <= Ztmp[0])
    {
        beginId = 0;
    }
    else if (Dref >= Ztmp[tabrow - 1])
    {
        beginId = tabrow - 1;
    }
    else
    {
        beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
                                                 [s = Dref](auto t)
                                                 { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure
    OCP_DBL Pbb = Pref;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI mynum = 10;
    OCP_DBL mydz = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    if (Dref < DOGC)
    {

        // reference pressure is gas pressure
        if (!gas)
            OCP_ABORT("SGOF is missing!");

        Pgref = Pref;
        gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--)
        {
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref = 0;
        Ptmp = Pgref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
            Ptmp += gammaGtmp * mydz;
        }
        Ptmp -= PcGO;
        for (USI i = 0; i < mynum; i++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, DOGC - i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (!EQUIL.PBVD.IsEmpty())
        {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else if (Dref > DOWC)
    {

        // reference pressure is water pressure
        Pwref = Pref;
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp = Pwref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp += gammaWtmp * mydz;
        }
        Ptmp += PcOW;

        for (USI i = 0; i < mynum; i++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, DOWC - i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (!EQUIL.PBVD.IsEmpty())
        {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas)
        {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;

            for (USI i = 0; i < mynum; i++)
            {
                if (!EQUIL.PBVD.IsEmpty())
                {
                    Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
                }
                gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
                Ptmp += gammaOtmp * mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
                Ptmp -= gammaGtmp * mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }
            for (USI id = beginId; id < tabrow - 1; id++)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }
    }
    else
    {

        // reference pressure is oil pressure
        Poref = Pref;
        if (!EQUIL.PBVD.IsEmpty())
        {
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);
        }
        gammaOtmp = flashCal[0]->GammaPhaseO(Poref, Pbb);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure in tab
        for (USI id = beginId; id > 0; id--)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Potmp[id], Pbb);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas)
        {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;

            for (USI i = 0; i < mynum; i++)
            {
                if (!EQUIL.PBVD.IsEmpty())
                {
                    Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
                }
                gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
                Ptmp += gammaOtmp * mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Ptmp);
                Ptmp -= gammaGtmp * mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            gammaGtmp = flashCal[0]->GammaPhaseG(Pgref);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }

            for (USI id = beginId; id < tabrow - 1; id++)
            {
                gammaGtmp = flashCal[0]->GammaPhaseG(Pgtmp[id]);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            if (!EQUIL.PBVD.IsEmpty())
            {
                Pbb = EQUIL.PBVD.Eval(0, Dref + i * mydz, 1);
            }
            gammaOtmp = flashCal[0]->GammaPhaseO(Ptmp, Pbb);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0), cdata(4, 0);
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po = data[1];
        OCP_DBL Pg = data[2];
        OCP_DBL Pw = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw = flow[0]->GetSwByPcow(Pcow);
        OCP_DBL Sg = 0;
        if (gas)
        {
            Sg = flow[0]->GetSgByPcgo(Pcgo);
        }
        if (Sw + Sg > 1)
        {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw = flow[0]->GetSwByPcgw(Pcgw);
            Sg = 1 - Sw;
        }

        if (1 - Sw < TINY)
        {
            // all water
            Po = Pw + flow[0]->GetPcowBySw(1.0);
        }
        else if (1 - Sg < TINY)
        {
            // all gas
            Po = Pg - flow[0]->GetPcgoBySg(1.0);
        }
        else if (1 - Sw - Sg < TINY)
        {
            // water and gas
            Po = Pg - flow[0]->GetPcgoBySg(Sg);
        }
        P[n] = Po;

        if (depth[n] < DOGC)
        {
            Pbb = Po;
        }
        else if (!EQUIL.PBVD.IsEmpty())
        {
            Pbb = EQUIL.PBVD.Eval(0, depth[n], 1);
        }
        Pb[n] = Pbb;

        // cal Sg and Sw
        Sw = 0;
        Sg = 0;
        USI ncut = 10;

        for (USI k = 0; k < ncut; k++)
        {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep = depth[n] + dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po = data[1];
            Pg = data[2];
            Pw = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            tmpSw = flow[0]->GetSwByPcow(Pcow);
            if (gas)
            {
                tmpSg = flow[0]->GetSgByPcgo(Pcgo);
            }
            if (tmpSw + tmpSg > 1)
            {
                // should me modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw = flow[0]->GetSwByPcgw(Pcgw);
                tmpSg = 1 - tmpSw;
            }
            Sw += tmpSw;
            Sg += tmpSg;
        }
        Sw /= ncut;
        Sg /= ncut;
        S[n * numPhase + numPhase - 1] = Sw;
        if (gas) {
            S[n * numPhase + numPhase - 2] = Sg;
        }
    }
}

/// Here tabrow is maximum number of depth nodes in table of depth vs pressure.
void Bulk::InitSjPcComp(const USI &tabrow, const Grid& myGrid)
{
    OCP_FUNCNAME;

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOW = EQUIL.PcOW;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGO = EQUIL.PcGO;
    OCP_DBL Zmin = 1E8;
    OCP_DBL Zmax = 0;

    OCP_DBL swco = flow[0]->GetSwco();
    OCP_DBL tmp;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        OCP_DBL temp1 = depth[n] - dz[n] / 2;
        OCP_DBL temp2 = depth[n] + dz[n] / 2;
        Zmin = Zmin < temp1 ? Zmin : temp1;
        Zmax = Zmax > temp2 ? Zmax : temp2;
    }
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    // creater table
    OCPTable DepthP(tabrow, 4);
    vector<OCP_DBL> &Ztmp = DepthP.GetCol(0);
    vector<OCP_DBL> &Potmp = DepthP.GetCol(1);
    vector<OCP_DBL> &Pwtmp = DepthP.GetCol(2);
    vector<OCP_DBL> &Pgtmp = DepthP.GetCol(3);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++)
    {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    // find the RefId
    USI beginId = 0;
    if (Dref <= Ztmp[0])
    {
        beginId = 0;
    }
    else if (Dref >= Ztmp[tabrow - 1])
    {
        beginId = tabrow - 1;
    }
    else
    {
        beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
                                                 [s = Dref](auto t)
                                                 { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure:
    OCP_DBL mytemp = T;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI mynum = 10;
    OCP_DBL mydz = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    if (Dref < DOGC)
    {

        // reference pressure is gas pressure
        Pgref = Pref;
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref = 0;
        Ptmp = Pgref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaGtmp * mydz;
        }
        Ptmp -= PcGO;
        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else if (Dref > DOWC)
    {

        // reference pressure is water pressure
        Pwref = Pref;
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp = Pwref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp += gammaWtmp * mydz;
        }
        Ptmp += PcOW;

        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaOtmp * mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the gas pressure in Dref by Poref
        Pgref = 0;
        Ptmp = Poref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp += PcGO;
        for (USI i = 0; i < mynum; i++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaGtmp * mydz;
        }
        Pgref = Ptmp;

        // find the gas pressure in tab
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else
    {

        // reference pressure is oil pressure
        Poref = Pref;
        gammaOtmp = flashCal[0]->GammaPhaseOG(Poref, mytemp, &initZi[0]);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure
        for (USI id = beginId; id > 0; id--)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Potmp[id], mytemp, &initZi[0]);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the gas pressure in Dref by Poref
        Pgref = 0;
        Ptmp = Poref;
        mydz = (DOGC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp += PcGO;
        for (USI i = 0; i < mynum; i++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp -= gammaGtmp * mydz;
        }
        Pgref = Ptmp;

        // find the gas pressure in tab
        gammaGtmp = flashCal[0]->GammaPhaseOG(Pgref, mytemp, &initZi[0]);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaGtmp = flashCal[0]->GammaPhaseOG(Pgtmp[id], mytemp, &initZi[0]);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;

        for (USI i = 0; i < mynum; i++)
        {
            gammaOtmp = flashCal[0]->GammaPhaseOG(Ptmp, mytemp, &initZi[0]);
            Ptmp += gammaOtmp * mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Ptmp);
            Ptmp -= gammaWtmp * mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        gammaWtmp = flashCal[0]->GammaPhaseW(Pwref);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++)
        {
            gammaWtmp = flashCal[0]->GammaPhaseW(Pwtmp[id]);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    cout << "Depth    "
         << "Poil    "
         << "Pwat    "
         << "Pgas" << endl;
    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0), cdata(4, 0);

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po = data[1];
        OCP_DBL Pw = data[2];
        OCP_DBL Pg = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw = flow[0]->GetSwByPcow(Pcow);
        OCP_DBL Sg = 0;
        if (gas)
        {
            Sg = flow[0]->GetSgByPcgo(Pcgo);
        }
        if (Sw + Sg > 1)
        {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw = flow[0]->GetSwByPcgw(Pcgw);
            Sg = 1 - Sw;
        }

        if (1 - Sw < TINY)
        {
            // all water
            Po = Pw + flow[0]->GetPcowBySw(1.0);
        }
        else if (1 - Sg < TINY)
        {
            // all gas
            Po = Pg - flow[0]->GetPcgoBySg(1.0);
        }
        else if (1 - Sw - Sg < TINY)
        {
            // water and gas
            Po = Pg - flow[0]->GetPcgoBySg(Sg);
        }
        P[n] = Po;

        // cal Sg and Sw
        Sw = 0;
        Sg = 0;
        USI ncut = 10;
        OCP_DBL avePcow = 0;

        for (USI k = 0; k < ncut; k++)
        {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep = depth[n] + dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po = data[1];
            Pw = data[2];
            Pg = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            avePcow += Pcow;
            tmpSw = flow[0]->GetSwByPcow(Pcow);
            if (gas)
            {
                tmpSg = flow[0]->GetSgByPcgo(Pcgo);
            }
            if (tmpSw + tmpSg > 1)
            {
                // should me modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw = flow[0]->GetSwByPcgw(Pcgw);
                tmpSg = 1 - tmpSw;
            }
            Sw += tmpSw;
            Sg += tmpSg;
        }
        Sw /= ncut;
        Sg /= ncut;
        avePcow /= ncut;

        if (SwatInitExist) {
            if (ScalePcow) {
                ScaleValuePcow[n] = 1.0;
            }
            if (SwatInit[n] <= swco) {
                Sw = swco;               
            }
            else {
                Sw = SwatInit[n];
                if (ScalePcow) {
                    if (avePcow > 0) {
                        tmp = flow[0]->GetPcowBySw(Sw);
                        if (tmp > 0) {
                            ScaleValuePcow[n] = avePcow / tmp;                           
                            /*USI I, J, K;
                            myGrid.GetIJKGrid(I, J, K, myGrid.activeMap_B2G[n]);
                            cout << I << "   " << J << "   " << K << "   " 
                                << fixed << setprecision(3) << "   " << ScaleValuePcow[n] * flow[0]->GetPcowBySw(0) << endl;*/
                        }
                    }
                }
            }           
        }        
        S[n * numPhase + numPhase - 1] = Sw;
        if (gas)
        {
            S[n * numPhase + numPhase - 2] = Sg;
        }
    }
}


/// Use initial saturation in blackoil model and initial Zi in compositional model.
/// It gives initial properties and some derivatives for IMPEC
/// It just gives Ni for FIM
void Bulk::InitFlash(const bool &flag)
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        flashCal[PVTNUM[n]]->InitFlash(P[n], Pb[n], T, &S[n * numPhase], rockVp[n],
                                       initZi.data());
        for (USI i = 0; i < numCom; i++)
        {
            Ni[n * numCom + i] = flashCal[PVTNUM[n]]->Ni[i];
        }
        if (flag)
        {
            PassFlashValue(n);
        }
    }

#ifdef DEBUG
    if (flag)
    {
        CheckSat();
    }
#endif // DEBUG
}

/// Use moles of component and pressure both in blackoil and compositional model.
void Bulk::Flash()
{
    OCP_FUNCNAME;

    if (comps) {
        FlashCOMP();
    }
    else {
        FlashBLKOIL();
    }

#ifdef DEBUG
    CheckSat();
#endif // DEBUG
}


void Bulk::FlashBLKOIL()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], 0, 0, 0);
        PassFlashValue(n);
    }
}


void Bulk::FlashCOMP()
{
    USI ftype;
    OCP_USI bId;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    // cout << endl << "==================================" << endl;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        ftype = 1;
        if (flagSkip[n])
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10)
            {
                ftype = 0;
            }
            // cout << setprecision(2) << scientific << minEig / 10 << "   " << fabs(1 - lP[n] / P[n]) << "   ";
            if (ftype == 1)
            {
                bId = n * numCom;
                Ntw = Nt[n] - Ni[bId + numCom - 1];
                for (USI i = 0; i < numCom - 1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }
        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValue(n);
    }
}


/// Use moles of component and pressure both in blackoil and compositional model.
void Bulk::FlashDeriv()
{
    OCP_FUNCNAME;

    if (comps) {
        FlashDerivCOMP();
    }
    else {
        FlashDerivBLKOIL();
    }

#ifdef DEBUG
    CheckSat();
#endif // DEBUG
}


/// Perform flash calculation with Ni in Black Oil Model
void Bulk::FlashDerivBLKOIL()
{
    // dSec_dPri.clear();
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        flashCal[PVTNUM[n]]->FlashDeriv(P[n], T, &Ni[n * numCom], 0, 0, 0);
        PassFlashValueDeriv(n);
    }
}


/// Perform flash calculation with Ni in Compositional Model
void Bulk::FlashDerivCOMP()
{
    USI ftype;
    OCP_USI bId;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    // dSec_dPri.clear();
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        ftype = 1;
        if (flagSkip[n])
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10)
            {
                ftype = 0;
            }
            // cout << setprecision(2) << scientific << minEig / 10 << "   " << fabs(1 - lP[n] / P[n]) << "   ";
            if (ftype == 1)
            {
                bId = n * numCom;
                Ntw = Nt[n] - Ni[bId + numCom_1];
                for (USI i = 0; i < numCom_1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }

        flashCal[PVTNUM[n]]->FlashDeriv(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValueDeriv(n);
    }
}

void Bulk::PassFlashValue(const OCP_USI &n)
{
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI pvtnum = PVTNUM[n];
    USI nptmp = 0;
    for (USI j = 0; j < numPhase; j++)
    {
        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bIdp + j])
        { // j -> bId + j   fix bugs.
            nptmp++;
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            for (USI i = 0; i < numCom; i++)
            {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
            }
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];
            vj[bIdp + j] = flashCal[pvtnum]->v[j];
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfp[n] = flashCal[pvtnum]->vfp;
    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++)
    {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

    if (comps)
    {
        phaseNum[n] = nptmp - 1; // water is excluded
        if (nptmp == 3)
        {
            // num of hydrocarbon phase equals 2
            // Calculate Ks
            OCP_USI bIdc1 = n * numCom_1;
            for (USI i = 0; i < numCom_1; i++)
            {
                Ks[bIdc1 + i] = flashCal[pvtnum]->xij[i] / flashCal[pvtnum]->xij[numCom + i];
            }
        }

        if (flashCal[pvtnum]->GetFtype() == 0)
        {
            flagSkip[n] = flashCal[pvtnum]->GetFlagSkip();
            if (flagSkip[n])
            {
                minEigenSkip[n] = flashCal[pvtnum]->GetMinEigenSkip();
                for (USI j = 0; j < numPhase - 1; j++)
                {
                    if (phaseExist[bIdp + j])
                    {
                        for (USI i = 0; i < numCom - 1; i++)
                        {
                            ziSkip[bIdc + i] = flashCal[pvtnum]->xij[j * numCom + i];
                        }
                        break;
                    }
                }
                PSkip[n] = P[n];
            }
        }

        if (miscible) {
            surTen[n] = flashCal[pvtnum]->GetSurTen();
        }
    }
}

void Bulk::PassFlashValueAIMc(const OCP_USI& n)
{
    // only var about volume needs, some flash var also
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI pvtnum = PVTNUM[n];


    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfp[n] = flashCal[pvtnum]->vfp;
    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++)
    {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

    if (comps)
    {
        USI nptmp = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (flashCal[pvtnum]->phaseExist[j]) {
                nptmp++;
            }
        }

        phaseNum[n] = nptmp - 1; // water is excluded
        if (nptmp == 3)
        {
            // num of hydrocarbon phase equals 2
            // Calculate Ks
            OCP_USI bIdc1 = n * numCom_1;
            for (USI i = 0; i < numCom_1; i++)
            {
                Ks[bIdc1 + i] = flashCal[pvtnum]->xij[i] / flashCal[pvtnum]->xij[numCom + i];
            }
        }

        if (flashCal[pvtnum]->GetFtype() == 0)
        {
            flagSkip[n] = flashCal[pvtnum]->GetFlagSkip();
            if (flagSkip[n])
            {
                minEigenSkip[n] = flashCal[pvtnum]->GetMinEigenSkip();
                for (USI j = 0; j < numPhase - 1; j++)
                {
                    if (phaseExist[bIdp + j])
                    {
                        for (USI i = 0; i < numCom - 1; i++)
                        {
                            ziSkip[bIdc + i] = flashCal[pvtnum]->xij[j * numCom + i];
                        }
                        break;
                    }
                }
                PSkip[n] = P[n];
            }
        }

        if (miscible) {
            surTen[n] = flashCal[pvtnum]->GetSurTen();
        }
    }
}

void Bulk::PassFlashValueDeriv(const OCP_USI &n)
{
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI pvtnum = PVTNUM[n];
    USI nptmp = 0;
    for (USI j = 0; j < numPhase; j++)
    {
        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bIdp + j])
        { // j -> bId + j fix bugs.
            nptmp++;
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];
            vj[bIdp + j] = flashCal[pvtnum]->v[j];

            // Derivatives
            muP[bIdp + j] = flashCal[pvtnum]->muP[j];
            xiP[bIdp + j] = flashCal[pvtnum]->xiP[j];
            rhoP[bIdp + j] = flashCal[pvtnum]->rhoP[j];
            for (USI i = 0; i < numCom; i++)
            {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
                mux[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->mux[j * numCom + i];
                xix[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xix[j * numCom + i];
                rhox[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhox[j * numCom + i];
            }
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfp[n] = flashCal[pvtnum]->vfp;

    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++)
    {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }
    //dSec_dPri.insert(dSec_dPri.end(), flashCal[pvtnum]->dXsdXp.begin(),
    //                 flashCal[pvtnum]->dXsdXp.end());
    Dcopy(lendSdP, &dSec_dPri[0] + n * lendSdP, &flashCal[pvtnum]->dXsdXp[0]);

    if (comps)
    {
        phaseNum[n] = nptmp - 1; // water is excluded
        if (nptmp == 3)
        {
            // num of hydrocarbon phase equals 2
            // Calculate Ks
            // IMPORTANT!!!
            // Ks will change as long as nums of hydroncarbon phase equals 2, and it will has an effect
            // on phase spliting calculation as a intial value. So you should not expect to obtain
            // the exact same result with identical P, T, Ni if the final mixture contains 2 hydroncarbon phase.
            OCP_USI bIdc1 = n * numCom_1;
            for (USI i = 0; i < numCom_1; i++) {
                Ks[bIdc1 + i] = flashCal[pvtnum]->xij[i] / flashCal[pvtnum]->xij[numCom + i];
            }
        }

        if (flashCal[pvtnum]->GetFtype() == 0)
        {
            flagSkip[n] = flashCal[pvtnum]->GetFlagSkip();
            if (flagSkip[n])
            {
                minEigenSkip[n] = flashCal[pvtnum]->GetMinEigenSkip();
                for (USI j = 0; j < numPhase - 1; j++)
                {
                    if (phaseExist[bIdp + j])
                    {
                        for (USI i = 0; i < numCom - 1; i++)
                        {
                            ziSkip[bIdc + i] = flashCal[pvtnum]->xij[j * numCom + i];
                        }
                        break;
                    }
                }
                PSkip[n] = P[n];
            }
        }
        if (miscible) {
            surTen[n] = flashCal[pvtnum]->GetSurTen();
        }
    }
}

void Bulk::ResetFlash()
{
    OCP_FUNCNAME;

    phaseExist = lphaseExist;
    S = lS;
    rho = lrho;
    xi = lxi;
    xij = lxij;
    mu = lmu;
    vj = lvj;
    vf = lvf;
    vfp = lvfp;
    vfi = lvfi;
}

/// Relative permeability and capillary pressure
void Bulk::CalKrPc()
{
    OCP_FUNCNAME;

    if (!miscible) {
        OCP_DBL tmp = 0;
        for (OCP_USI n = 0; n < numBulk; n++)
        {
            OCP_USI bId = n * numPhase;
            flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId], 0, tmp, tmp);
            for (USI j = 0; j < numPhase; j++)
                Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
        }
    }
    else {
        for (OCP_USI n = 0; n < numBulk; n++)
        {
            OCP_USI bId = n * numPhase;
            flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId], surTen[n], Fk[n], Fp[n]);
            for (USI j = 0; j < numPhase; j++)
                Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
        }
    } 

    if (ScalePcow) {
        // correct
        USI Wid = phase2Index[WATER];
        for (USI n = 0; n < numBulk; n++) {
            Pc[n * numPhase + Wid] *= ScaleValuePcow[n];
            Pj[n * numPhase + Wid] = P[n] + Pc[n * numPhase + Wid];
        }
    }
}

void Bulk::CalKrPcDeriv()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPcDeriv(&S[bId], &kr[bId], &Pc[bId],
                                      &dKr_dS[bId * numPhase],
                                      &dPcj_dS[bId * numPhase]);
        for (USI j = 0; j < numPhase; j++)
            Pj[bId + j] = P[n] + Pc[bId + j];
    }
}

void Bulk::CalVpore()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        OCP_DBL dP = rockC1 * (P[n] - rockPref);
        rockVp[n] = rockVpInit[n] * (1 + dP);
        // rockVp[n] = rockVpInit[n] * (1 + dP + dP * dP / 2);
    }
}

OCP_DBL Bulk::CalFPR() const
{
    OCP_FUNCNAME;

    OCP_DBL ptmp = 0;
    OCP_DBL vtmp = 0;
    OCP_DBL tmp = 0;

    if (numPhase == 3)
    {
        for (OCP_USI n = 0; n < numBulk; n++)
        {
            tmp = rockVp[n] * (1 - S[n * numPhase + 2]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    }
    else if (numPhase < 3)
    {
        for (OCP_USI n = 0; n < numBulk; n++)
        {
            tmp = rockVp[n] * (S[n * numPhase]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    }
    else
    {
        OCP_ABORT("Number of phases is out of range!");
    }
    return ptmp / vtmp;
}

void Bulk::CalMaxChange()
{
    OCP_FUNCNAME;

    dPmax = 0;
    dNmax = 0;
    dSmax = 0;
    dVmax = 0;
    OCP_DBL tmp = 0;
    OCP_USI id;
    OCP_USI ndPmax, ndNmax, ndSmax, ndVmax;

    for (OCP_USI n = 0; n < numBulk; n++)
    {

        // dP
        tmp   = fabs(P[n] - lP[n]);
        if (dPmax < tmp) {
            dPmax = tmp;
            ndPmax = n;
        }

        // dS
        for (USI j = 0; j < numPhase; j++) {
            id    = n * numPhase + j;
            tmp   = fabs(S[id] - lS[id]);
            if (dSmax < tmp) {
                dSmax = tmp;
                ndSmax = n;
            }
        }

        // dN
        for (USI i = 0; i < numCom; i++)
        {
            id = n * numCom + i;

            tmp = fabs(max(Ni[id], lNi[id]));
            if (tmp > TINY) {
                tmp   = fabs(Ni[id] - lNi[id]) / tmp;
                if (dNmax < tmp) {
                    dNmax = tmp;
                    ndNmax = n;
                }
            }
        }

        tmp   = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (dVmax < tmp) {
            dVmax = tmp;
            ndVmax = n;
        }
    }

    //cout << scientific << setprecision(6);
    //cout << dPmax << "   " << dNmax << "   " << dSmax << "   " << setw(20) << dVmax << endl;
    //cout << ndPmax << "   " << ndNmax << "   " << ndSmax << "   " << ndVmax << endl;
}

/// Return true if no negative pressure and false otherwise.
bool Bulk::CheckP() const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (P[n] < 0)
        {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << P[n];
            OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            cout << "P = " << P[n] << endl;
            return false;
        }
    }

    return true;
}

/// Return true if no negative Ni and false otherwise.
bool Bulk::CheckNi() const
{
    OCP_FUNCNAME;

    OCP_USI len = numBulk * numCom;
    for (OCP_USI n = 0; n < len; n++)
    {
        if (Ni[n] < 0.0)
        {
            OCP_USI bId = n / numCom;
            USI cId = n - bId * numCom;
            std::ostringstream NiStringSci;
            NiStringSci << std::scientific << Ni[n];
            OCP_WARNING("Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                        std::to_string(bId) + "] = " + NiStringSci.str());
            return false;
        }
    }
    return true;
}

/// Return true if all Ve < Vlim and false otherwise.
bool Bulk::CheckVe(const OCP_DBL &Vlim) const
{
    OCP_FUNCNAME;

    // bool tmpflag = true;

    OCP_DBL dVe = 0.0;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        dVe = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (dVe > Vlim)
        {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << dVe
                 << " is too big!" << endl;
            // OutputInfo(n);
            return false;
            // tmpflag = false;
        }
    }
    // OutputInfo(39);
    // if (!tmpflag) return false;
    return true;
}

void Bulk::CheckDiff()
{
    OCP_FUNCNAME;

    OCP_DBL tmp;
    OCP_USI id;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        for (USI j = 0; j < numPhase; j++)
        {
            id = n * numPhase + j;
            tmp = fabs(phaseExist[id] - lphaseExist[id]);
            if (tmp >= 1E-10)
            {
                cout << "Difference in phaseExist\t" << tmp << "\n";
            }
            if (lphaseExist[id] || phaseExist[id])
            {
                tmp = fabs(S[id] - lS[id]);
                if (tmp >= 1E-10)
                {
                    cout << scientific << setprecision(16);
                    cout << "Bulk[" << n << "]" << endl;
                    cout << "Saturation" << endl;
                    for (USI j = 0; j < numPhase; j++) {
                        cout << S[n * numPhase + j] << "   " << lS[n * numPhase + j] << endl;
                    }
                    cout << "Pressure" << endl;
                    cout << P[n] << "   " << lP[n] << endl;
                    cout << "Ni" << endl;
                    for (USI i = 0; i < numCom; i++) {
                        cout << Ni[n * numCom + i] << "   " << lNi[n * numCom + i] << endl;
                    }
                    cout << "PhaseNum" << endl;
                    cout << phaseNum[n] << "   " << lphaseNum[n] << endl;
                    cout << "minEigenSkip" << endl;
                    cout << minEigenSkip[n] << "   " << lminEigenSkip[n] << endl;
                    cout << "flagSkip" << endl;
                    cout << flagSkip[n] << "   " << lflagSkip[n] << endl;
                    cout << "PSkip" << endl;
                    cout << PSkip[n] << "   " << lPSkip[n] << endl;
                    cout << "ziSkip" << endl;
                    for (USI i = 0; i < numCom; i++) {
                        cout << ziSkip[n * numCom + i] << "   " << lziSkip[n * numCom + i] << endl;
                    }
                    cout << "Ks" << endl;
                    for (USI i = 0; i < numCom_1; i++) {
                        cout << Ks[n * numCom_1 + i] << "   " << lKs[n * numCom_1 + i] << endl;
                    }

                    cout << "Difference in S\t" << tmp << "  " << phaseExist[id]
                         << "\n";
                }
                tmp = fabs(xi[id] - lxi[id]);
                if (tmp >= 1E-10)
                {
                    cout << "Difference in Xi\t" << tmp << "  " << phaseExist[id]
                         << "\n";
                }
                tmp = fabs(rho[id] - lrho[id]);
                if (tmp >= 1E-10)
                {
                    cout << "Difference in rho\t" << tmp << "  " << phaseExist[id]
                         << "\n";
                }
            }
        }
    }
}

void Bulk::CheckSat() const
{
    OCP_FUNCNAME;

    OCP_DBL tmp;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        tmp = 0.0;
        for (USI j = 0; j < numPhase; j++)
        {
            if (phaseExist[n * numPhase + j])
            {
                if (S[n * numPhase + j] < 0)
                {
                    OCP_ABORT("Negative volume!");
                }
                tmp += S[n * numPhase + j];
            }
        }
        if (fabs(tmp - 1) > TINY)
        {
            OCP_ABORT("Saturation is greater than 1!");
        }
    }
}

USI Bulk::GetMixMode() const
{
    OCP_FUNCNAME;

    if (blackOil)
        return BLKOIL;
    else if (comps)
        return EOS_PVTW;
    else
        OCP_ABORT("Mixture model is not supported!");
}

void Bulk::CalSomeInfo(const Grid &myGrid) const
{
    // test
    OCP_DBL depthMax = 0;
    OCP_USI ndepa = 0;
    OCP_DBL depthMin = 1E8;
    OCP_USI ndepi = 0;
    OCP_DBL dxMax = 0;
    OCP_USI nxa = 0;
    OCP_DBL dxMin = 1E8;
    OCP_USI nxi = 0;
    OCP_DBL dyMax = 0;
    OCP_USI nya = 0;
    OCP_DBL dyMin = 1E8;
    OCP_USI nyi = 0;
    OCP_DBL dzMax = 0;
    OCP_USI nza = 0;
    OCP_DBL dzMin = 1E8;
    OCP_USI nzi = 0;
    OCP_DBL RVMax = 0;
    OCP_USI nRVa = 0;
    OCP_DBL RVMin = 1E8;
    OCP_USI nRVi = 0;
    OCP_DBL RVPMax = 0;
    OCP_USI nRVPa = 0;
    OCP_DBL RVPMin = 1E8;
    OCP_USI nRVPi = 0;
    OCP_DBL PerxMax = 0;
    OCP_USI nPerxa = 0;
    OCP_DBL PerxMin = 1E8;
    OCP_USI nPerxi = 0;
    USI I, J, K;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        // if (!activeMap_G2B[nn].IsAct())
        //     continue;
        // OCP_USI n = activeMap_G2B[nn].GetId();
        if (depthMax < depth[n])
        {
            depthMax = depth[n];
            ndepa = n;
        }
        if (depthMin > depth[n])
        {
            depthMin = depth[n];
            ndepi = n;
        }
        if (dxMax < dx[n])
        {
            dxMax = dx[n];
            nxa = n;
        }
        if (dxMin > dx[n])
        {
            dxMin = dx[n];
            nxi = n;
        }
        if (dyMax < dy[n])
        {
            dyMax = dy[n];
            nya = n;
        }
        if (dyMin > dy[n])
        {
            dyMin = dy[n];
            nyi = n;
        }
        if (dzMax < dz[n])
        {
            dzMax = dz[n];
            nza = n;
        }
        if (dzMin > dz[n])
        {
            dzMin = dz[n];
            nzi = n;
        }
        OCP_DBL tmp = myGrid.v[myGrid.activeMap_B2G[n]];
        if (RVMax < tmp)
        {
            RVMax = tmp;
            nRVa = n;
        }
        if (RVMin > tmp)
        {
            RVMin = tmp;
            nRVi = n;
        }
        tmp = rockVp[n];
        if (RVPMax < tmp)
        {
            RVPMax = tmp;
            nRVPa = n;
        }
        if (RVPMin > tmp)
        {
            RVPMin = tmp;
            nRVPi = n;
        }
        tmp = rockKx[n];
        if (PerxMax < tmp)
        {
            PerxMax = tmp;
            nPerxa = n;
        }
        if (PerxMin > tmp && fabs(tmp) > 1E-8)
        {
            PerxMin = tmp;
            nPerxi = n;
        }
    }

    cout << "---------------------" << endl
         << "BULK" << endl
         << "---------------------" << endl;
    myGrid.GetIJKGrid(I, J, K, ndepa);
    cout << "  Depthmax = " << depthMax << " feet  (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, ndepi);
    cout << "  Depthmin = " << depthMin << " feet  (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, nxa);
    cout << "  DXmax    = " << dxMax << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nxi);
    cout << "  DXmin    = " << dxMin << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nya);
    cout << "  DYmax    = " << dyMax << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nyi);
    cout << "  DYmin    = " << dyMin << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nza);
    cout << "  DZmax    = " << dzMax << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nzi);
    cout << "  DZmin    = " << dzMin << " feet  (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nRVa);
    cout << "  RVmax    = " << RVMax / CONV1 << " rb   (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, nRVi);
    cout << "  RVmin    = " << RVMin / CONV1 << " rb   (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, nRVPa);
    cout << "  RVmax    = " << RVPMax / CONV1 << " rb   (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, nRVPi);
    cout << "  RVmin    = " << RVPMin / CONV1 << " rb   (" << I << ", " << J << ", " << K
         << ")" << endl;
    myGrid.GetIJKGrid(I, J, K, nPerxa);
    cout << "  Perxmax  = " << PerxMax << "   (" << I << ", " << J << ", " << K << ")"
         << endl;
    myGrid.GetIJKGrid(I, J, K, nPerxi);
    cout << "  Perxmin  = " << scientific << PerxMin << "   (" << I << ", " << J << ", "
         << K << ")" << endl;
}


/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateAuxIMPEC()
{
    OCP_FUNCNAME;

    // IMPEC var
    vfi.resize(numBulk * numCom);
    vfp.resize(numBulk);
    cfl.resize(numBulk * numPhase);

    // Last Step
    lP.resize(numBulk);
    lPj.resize(numBulk * numPhase);
    lPc.resize(numBulk * numPhase);
    lphaseExist.resize(numBulk * numPhase);
    lS.resize(numBulk * numPhase);
    lrho.resize(numBulk * numPhase);
    lxi.resize(numBulk * numPhase);
    lxij.resize(numBulk * numPhase * numCom);
    lNi.resize(numBulk * numCom);
    lmu.resize(numBulk * numPhase);
    lkr.resize(numBulk * numPhase);
    lvj.resize(numBulk * numPhase);
    lvf.resize(numBulk);
    lNt.resize(numBulk);
    lvfi.resize(numBulk * numCom);
    lvfp.resize(numBulk);
    rockLVp.resize(numBulk);
}

void Bulk::GetSolIMPEC(const vector<OCP_DBL> &u)
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        P[n] = u[n];
        for (USI j = 0; j < numPhase; j++)
        {
            OCP_USI id = n * numPhase + j;
            if (phaseExist[id])
                Pj[id] = P[n] + Pc[id];
        }
    }
}

OCP_DBL Bulk::CalCFL() const
{
    OCP_FUNCNAME;

    OCP_DBL tmp = 0;
    OCP_USI id;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        for (USI j = 0; j < numPhase; j++)
        {
            id = n * numPhase + j;
            if (phaseExist[id])
            {
                if (vj[id] <= 0)
                    continue; // temp
                cfl[id] /= vj[id];

#ifdef DEBUG
                if (!isfinite(cfl[id]))
                {
                    OCP_ABORT("cfl is nan!");
                }
#endif // DEBUG

                if (tmp < cfl[id])
                    tmp = cfl[id];
            }
        }
    }
    return tmp;
}

void Bulk::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;
    lphaseNum = phaseNum;
    lminEigenSkip = minEigenSkip;
    lflagSkip = flagSkip;
    lziSkip = ziSkip;
    lPSkip = PSkip;
    lKs = Ks;

    lP = P;
    lPj = Pj;
    lPc = Pc;
    lphaseExist = phaseExist;
    lS = S;
    lrho = rho;
    lxi = xi;
    lxij = xij;
    lNi = Ni;
    lmu = mu;
    lkr = kr;
    lvj = vj;
    lvf = vf;
    lvfi = vfi;
    lvfp = vfp;
    rockLVp = rockVp;
    lNt = Nt;
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateAuxFIM()
{
    OCP_FUNCNAME;

    // FIM var
    vfi.resize(numBulk * numCom);
    vfp.resize(numBulk);
    muP.resize(numBulk * numPhase);
    xiP.resize(numBulk * numPhase);
    rhoP.resize(numBulk * numPhase);
    mux.resize(numBulk * numCom * numPhase);
    xix.resize(numBulk * numCom * numPhase);
    rhox.resize(numBulk * numCom * numPhase);
    lendSdP = (numCom + 1) * (numCom + 1) * numPhase;
    dSec_dPri.resize(numBulk * lendSdP);
    dKr_dS.resize(numBulk * numPhase * numPhase);
    dPcj_dS.resize(numBulk * numPhase * numPhase);

    // Last Step
    lP.resize(numBulk);
    lPj.resize(numBulk * numPhase);
    lPc.resize(numBulk * numPhase);
    lphaseExist.resize(numBulk * numPhase);
    lS.resize(numBulk * numPhase);
    lrho.resize(numBulk * numPhase);
    lxi.resize(numBulk * numPhase);
    lxij.resize(numBulk * numPhase * numCom);
    lNi.resize(numBulk * numCom);
    lmu.resize(numBulk * numPhase);
    lkr.resize(numBulk * numPhase);
    lvj.resize(numBulk * numPhase);
    lvf.resize(numBulk);
    lNt.resize(numBulk);
    lvfi.resize(numBulk * numCom);
    lvfp.resize(numBulk);
    rockLVp.resize(numBulk);
    lmuP.resize(numBulk * numPhase);
    lxiP.resize(numBulk * numPhase);
    lrhoP.resize(numBulk * numPhase);
    lmux.resize(numBulk * numCom * numPhase);
    lxix.resize(numBulk * numCom * numPhase);
    lrhox.resize(numBulk * numCom * numPhase);
    ldSec_dPri.resize(numBulk * lendSdP);
    ldKr_dS.resize(numBulk * numPhase * numPhase);
    ldPcj_dS.resize(numBulk * numPhase * numPhase);
}

void Bulk::GetSolFIM(const vector<OCP_DBL> &u, const OCP_DBL &dPmaxlim,
                     const OCP_DBL &dSmaxlim)
{
    OCP_FUNCNAME;

    NRdSmax = 0;
    NRdPmax = 0;
    OCP_DBL dP;
    USI row = numPhase * (numCom + 1);
    USI col = numCom + 1;
    USI bsize = row * col;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL chopmin = 1;
    OCP_DBL choptmp = 0;

    for (OCP_USI n = 0; n < numBulk; n++)
    {

        chopmin = 1;

        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);
        DaAxpby(row, col, 1, dSec_dPri.data() + n * bsize, u.data() + n * col, 1,
                dtmp.data());

        for (USI j = 0; j < numPhase; j++)
        {

            choptmp = 1;
            if (fabs(dtmp[j]) > dSmaxlim)
            {
                choptmp = dSmaxlim / fabs(dtmp[j]);
            }
            else if (S[n * numPhase + j] + dtmp[j] < 0.0)
            {
                choptmp = 0.9 * S[n * numPhase + j] / fabs(dtmp[j]);
            }

            chopmin = min(chopmin, choptmp);
            NRdSmax = max(NRdSmax, choptmp * fabs(dtmp[j]));
        }
        dP = u[n * col];
        choptmp = dPmaxlim / fabs(dP);
        chopmin = min(chopmin, choptmp);
        NRdPmax = max(NRdPmax, fabs(dP));
        P[n] += dP; // seems better

        //// Correct chopmin
        // for (USI i = 0; i < numCom; i++) {
        //    if (Ni[n * numCom + i] + u[n * col + 1 + i] < 0) {
        //        chopmin = 0.9 * min(chopmin, fabs(Ni[n * numCom + i] / u[n * col + 1 +
        //        i]));

        //        //if (chopmin < 0 || !isfinite(chopmin)) {
        //        //    OCP_ABORT("Wrong Chop!");
        //        //}
        //    }
        //}

        for (USI i = 0; i < numCom; i++)
        {
            Ni[n * numCom + i] += u[n * col + 1 + i] * chopmin;

            // if (Ni[n * numCom + i] < 0) {
            //    cout << Ni[n * numCom + i] << "  " << u[n * col + 1 + i] * chopmin <<
            //    "   " << chopmin << endl;
            //}
        }
    }
}

OCP_DBL Bulk::GetSol01FIM(const vector<OCP_DBL> &u)
{
    OCP_DBL tmp;
    OCP_DBL alpha = 1;
    USI len = numCom + 1;
    OCP_DBL Ni0, dNi;

    for (OCP_USI n = 0; n < numBulk; n++)
    {

        for (USI i = 0; i < numCom; i++)
        {
            Ni0 = Ni[n * numCom + i];
            dNi = u[n * len + 1 + i];
            if (Ni0 <= 0 && dNi <= 0)
            {
                continue;
            }
            if (Ni0 + dNi <= 0 && dNi < 0)
            {
                tmp = 0.9 * fabs(Ni0 / dNi);
                alpha = min(tmp, alpha);
            }
        }
    }

    // Newton step
    for (OCP_USI n = 0; n < numBulk; n++)
    {

        P[n] += alpha * u[n * len];

        for (USI i = 0; i < numCom; i++)
        {
            if (Ni[n * numCom + i] <= 0 && u[n * len + 1 + i] <= 0)
            {
                continue;
            }

            tmp = Ni[n * numCom + i];
            Ni[n * numCom + i] += alpha * u[n * len + 1 + i];
        }
    }

    return alpha;
}

void Bulk::CalRelResFIM(ResFIM &resFIM) const
{
    OCP_FUNCNAME;

    //OCP_USI tmpid01 = -1;
    //OCP_USI tmpid02 = -1;
    OCP_DBL tmp;

    const USI len = numCom + 1;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        for (USI i = 0; i < len; i++)
        {
            tmp = fabs(resFIM.res[n * len + i] / rockVp[n]);
            if (resFIM.maxRelRes_v < tmp)
            {
                resFIM.maxRelRes_v = tmp;
                // tmpid01            = n;
            }
        }
        for (USI i = 1; i < len; i++)
        {
            tmp = fabs(resFIM.res[n * len + i] / Nt[n]);
            if (resFIM.maxRelRes_mol < tmp)
            {
                resFIM.maxRelRes_mol = tmp;
                // tmpid02              = n;
            }
        }
    }
     //cout << scientific;
     //if (tmpid01 < numBulk) {
     //    cout << "maxRelRes_v: " << tmpid01 << "   " << S[tmpid01 * numPhase] << "   "
     //        << S[tmpid01 * numPhase + 1] << "   " << S[tmpid01 * numPhase + 2] <<
     //        "   " << resFIM.maxRelRes_v << endl;
     //}
     //if (tmpid02 < numBulk) {
     //    cout << "maxRelRes_mol: " << tmpid02 << "   " << S[tmpid02 * numPhase] << "   "
     //        << S[tmpid02 * numPhase + 1] << "   " << S[tmpid02 * numPhase + 2] <<
     //        "   " << resFIM.maxRelRes_mol << endl;
     //}
}

void Bulk::ResetFIM()
{
    OCP_FUNCNAME;

    phaseNum = lphaseNum;
    minEigenSkip = lminEigenSkip;
    flagSkip = lflagSkip;
    ziSkip = lziSkip;
    PSkip = lPSkip;
    Ks = lKs;

    P = lP;
    Pj = lPj;
    Pc = lPc;
    phaseExist = lphaseExist;
    S = lS;
    rho = lrho;
    xi = lxi;
    xij = lxij;
    Ni = lNi;
    mu = lmu;
    kr = lkr;
    vj = lvj;
    vf = lvf;
    Nt = lNt;
    vfi = lvfi;
    vfp = lvfp;
    rockVp = rockLVp;
    muP = lmuP;
    xiP = lxiP;
    rhoP = lrhoP;
    mux = lmux;
    xix = lxix;
    rhox = lrhox;
    dSec_dPri = ldSec_dPri;
    dKr_dS = ldKr_dS;
    dPcj_dS = ldPcj_dS;
}

void Bulk::UpdateLastStepFIM()
{
    OCP_FUNCNAME;
    lphaseNum = phaseNum;
    lminEigenSkip = minEigenSkip;
    lflagSkip = flagSkip;
    lziSkip = ziSkip;
    lPSkip = PSkip;
    lKs = Ks;

    lP = P;
    lPj = Pj;
    lPc = Pc;
    lphaseExist = phaseExist;
    lS = S;
    lrho = rho;
    lxi = xi;
    lxij = xij;
    lNi = Ni;
    lmu = mu;
    lkr = kr;
    lvj = vj;
    lvf = vf;
    lNt = Nt;
    lvfi = vfi;
    lvfp = vfp;
    rockLVp = rockVp;
    lmuP = muP;
    lxiP = xiP;
    lrhoP = rhoP;
    lmux = mux;
    lxix = xix;
    lrhox = rhox;
    ldSec_dPri = dSec_dPri;
    ldKr_dS = dKr_dS;
    ldPcj_dS = dPcj_dS;
}

void Bulk::OutputInfo(const OCP_USI &n) const
{
    OCP_USI bIdC = n * numCom;
    OCP_USI bIdP = n * numPhase;
    OCP_USI bIdPC = bIdP * numCom;

    cout << "------------------------------" << endl;
    cout << "Bulk[" << n << "]" << endl;
    cout << fixed << setprecision(3);
    for (USI i = 0; i < numCom; i++)
    {
        cout << Ni[bIdC + i] << "   ";
    }
    cout << endl
         << P[n] << "   " << T;
    cout << endl;

    cout << fixed << setprecision(8);
    if (phaseExist[bIdP + 0])
    {
        for (USI i = 0; i < numCom; i++)
        {
            cout << xij[bIdPC + i] << "   ";
        }
    }
    else
    {
        for (USI i = 0; i < numCom; i++)
        {
            cout << 0.000000 << "   ";
        }
    }
    cout << phaseExist[bIdP + 0] << "   ";
    cout << xi[bIdP + 0] << "   ";
    cout << S[bIdP + 0] << "   ";
    cout << endl;

    if (phaseExist[bIdP + 1])
    {
        for (USI i = 0; i < numCom; i++)
        {
            cout << xij[bIdPC + numCom + i] << "   ";
        }
    }
    else
    {
        for (USI i = 0; i < numCom; i++)
        {
            cout << 0.000000 << "   ";
        }
    }

    cout << phaseExist[bIdP + 1] << "   ";
    cout << xi[bIdP + 1] << "   ";
    cout << S[bIdP + 1] << "   ";
    cout << endl;

    cout << fixed << setprecision(3);

    cout << vf[n] << "   " << rockVp[n] << "   ";
    cout << fixed << setprecision(12);
    cout << fabs(vf[n] - rockVp[n]) / rockVp[n] << endl;
    cout << fixed << setprecision(3);
    cout << vj[bIdP] << "   " << vj[bIdP + 1] << "   " << vj[bIdP + 2] << endl;
    cout << "------------------------------" << endl;
}

/////////////////////////////////////////////////////////////////////
// For AIMt
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateAuxAIM(const OCP_DBL& ratio)
{
    OCP_FUNCNAME;

    maxNumFIMBulk = numBulk * ratio;
    if (maxNumFIMBulk < wellBulkId.capacity()) {
        maxNumFIMBulk = wellBulkId.capacity();
    }
    FIMBulk.reserve(maxNumFIMBulk);
    map_Bulk2FIM.resize(numBulk, -1);

    muP.resize(maxNumFIMBulk * numPhase);
    xiP.resize(maxNumFIMBulk * numPhase);
    rhoP.resize(maxNumFIMBulk * numPhase);
    mux.resize(maxNumFIMBulk * numCom * numPhase);
    xix.resize(maxNumFIMBulk * numCom * numPhase);
    rhox.resize(maxNumFIMBulk * numCom * numPhase);
    lendSdP = (numCom + 1) * (numCom + 1) * numPhase;
    dSec_dPri.resize(maxNumFIMBulk * lendSdP);
    dKr_dS.resize(maxNumFIMBulk * numPhase * numPhase);
    dPcj_dS.resize(maxNumFIMBulk * numPhase * numPhase);

    lmuP.resize(maxNumFIMBulk * numPhase);
    lxiP.resize(maxNumFIMBulk * numPhase);
    lrhoP.resize(maxNumFIMBulk * numPhase);
    lmux.resize(maxNumFIMBulk * numCom * numPhase);
    lxix.resize(maxNumFIMBulk * numCom * numPhase);
    lrhox.resize(maxNumFIMBulk * numCom * numPhase);
    ldSec_dPri.resize(maxNumFIMBulk * lendSdP);
    ldKr_dS.resize(maxNumFIMBulk * numPhase * numPhase);
    ldPcj_dS.resize(maxNumFIMBulk * numPhase * numPhase);

    FIMNi.resize(maxNumFIMBulk * numCom);
}

void Bulk::FlashDerivAIM(const bool& IfAIMs)
{
    OCP_FUNCNAME;

    USI ftype;
    OCP_USI n, bIdc;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    dSec_dPri.clear();

    USI len;
    if (IfAIMs) {
        len = FIMBulk.size();
    }
    else {
        len = numFIMBulk;
    }

    for (OCP_USI k = 0; k < len; k++)
    {
        n = FIMBulk[k];
        bIdc = n * numCom;     
        ftype = 1;
        if (flagSkip[n]) 
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10) {
                ftype = 0;
            }           
            if (ftype == 1) {
                
                Ntw = Nt[n] - Ni[bIdc + numCom_1];
                for (USI i = 0; i < numCom_1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bIdc + i] / Ntw - ziSkip[bIdc + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }

        flashCal[PVTNUM[n]]->FlashDeriv(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValueDerivAIM(k);
    }

#ifdef DEBUG
    CheckSat();
#endif // DEBUG
}

void Bulk::PassFlashValueDerivAIM(const OCP_USI& fn)
{
    OCP_FUNCNAME;

    OCP_USI n = FIMBulk[fn];
    OCP_USI bIdp = n * numPhase;
    OCP_USI fnp = fn * numPhase;
    USI pvtnum = PVTNUM[n];
    USI nptmp = 0;
    for (USI j = 0; j < numPhase; j++)
    {
        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bIdp + j])
        { // j -> bId + j fix bugs.
            nptmp++;
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];
            vj[bIdp + j] = flashCal[pvtnum]->v[j];

            // Derivatives
            muP[fnp + j] = flashCal[pvtnum]->muP[j];
            xiP[fnp + j] = flashCal[pvtnum]->xiP[j];
            rhoP[fnp + j] = flashCal[pvtnum]->rhoP[j];
            for (USI i = 0; i < numCom; i++)
            {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
                mux[fnp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->mux[j * numCom + i];
                xix[fnp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xix[j * numCom + i];
                rhox[fnp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhox[j * numCom + i];
            }
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfp[n] = flashCal[pvtnum]->vfp;

    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++)
    {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }
    dSec_dPri.insert(dSec_dPri.end(), flashCal[pvtnum]->dXsdXp.begin(),
        flashCal[pvtnum]->dXsdXp.end());

    if (comps)
    {
        phaseNum[n] = nptmp - 1; // water is excluded
        if (nptmp == 3)
        {
            // num of hydrocarbon phase equals 2
            // Calculate Ks
            OCP_USI bIdc1 = n * numCom_1;
            for (USI i = 0; i < numCom_1; i++)
            {
                Ks[bIdc1 + i] = flashCal[pvtnum]->xij[i] / flashCal[pvtnum]->xij[numCom + i];
            }
        }

        if (flashCal[pvtnum]->GetFtype() == 0)
        {
            flagSkip[n] = flashCal[pvtnum]->GetFlagSkip();
            if (flagSkip[n])
            {
                minEigenSkip[n] = flashCal[pvtnum]->GetMinEigenSkip();
                for (USI j = 0; j < numPhase - 1; j++)
                {
                    if (phaseExist[bIdp + j])
                    {
                        for (USI i = 0; i < numCom - 1; i++)
                        {
                            ziSkip[bIdc + i] = flashCal[pvtnum]->xij[j * numCom + i];
                        }
                        break;
                    }
                }
                PSkip[n] = P[n];
            }
        }
        if (miscible) {
            surTen[n] = flashCal[pvtnum]->GetSurTen();
        }
    }
}

void Bulk::CalKrPcDerivAIM(const bool& IfAIMs)
{
    OCP_FUNCNAME;

    USI len;
    if (IfAIMs) {
        len = FIMBulk.size();
    }
    else {
        len = numFIMBulk;
    }

    for (OCP_USI fn = 0; fn < len; fn++)
    {
        OCP_USI n = FIMBulk[fn];
        OCP_USI fnp = fn * numPhase;
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPcDeriv(&S[bId], &kr[bId], &Pc[bId],
            &dKr_dS[fnp * numPhase],
            &dPcj_dS[fnp * numPhase]);
        for (USI j = 0; j < numPhase; j++)
            Pj[bId + j] = P[n] + Pc[bId + j];
    }
}

void Bulk::CalRelResAIMt(ResFIM& resFIM) const
{
    OCP_FUNCNAME;

    // OCP_USI tmpid01 = -1;
    // OCP_USI tmpid02 = -1;
    OCP_DBL tmp;

    const USI len = numCom + 1;
    for (OCP_USI fn = 0; fn < numFIMBulk; fn++) 
    {
        OCP_USI n = FIMBulk[fn];

        for (USI i = 0; i < len; i++)
        {
            tmp = fabs(resFIM.res[fn * len + i] / rockVp[n]);
            if (resFIM.maxRelRes_v < tmp)
            {
                resFIM.maxRelRes_v = tmp;
                // tmpid01            = n;
            }
        }
        for (USI i = 1; i < len; i++)
        {
            tmp = fabs(resFIM.res[fn * len + i] / Nt[n]);
            if (resFIM.maxRelRes_mol < tmp)
            {
                resFIM.maxRelRes_mol = tmp;
                // tmpid02              = n;
            }
        }
    }
}

void Bulk::GetSolAIMt(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
    const OCP_DBL& dSmaxlim)
{
    OCP_FUNCNAME;

    NRdSmax = 0;
    NRdPmax = 0;
    OCP_DBL dP;
    USI row = numPhase * (numCom + 1);
    USI col = numCom + 1;
    USI bsize = row * col;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL chopmin = 1;
    OCP_DBL choptmp = 0;

    for (OCP_USI fn = 0; fn < numFIMBulk; fn++)
    {
        OCP_USI n = FIMBulk[fn];

        chopmin = 1;
        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);
        DaAxpby(row, col, 1, dSec_dPri.data() + fn * bsize, u.data() + fn * col, 1,
            dtmp.data());

        for (USI j = 0; j < numPhase; j++)
        {

            choptmp = 1;
            if (fabs(dtmp[j]) > dSmaxlim)
            {
                choptmp = dSmaxlim / fabs(dtmp[j]);
            }
            else if (S[n * numPhase + j] + dtmp[j] < 0.0)
            {
                choptmp = 0.9 * S[n * numPhase + j] / fabs(dtmp[j]);
            }

            chopmin = min(chopmin, choptmp);
            NRdSmax = max(NRdSmax, choptmp * fabs(dtmp[j]));
        }
        dP = u[fn * col];
        choptmp = dPmaxlim / fabs(dP);
        chopmin = min(chopmin, choptmp);
        dP *= chopmin;
        NRdPmax = max(NRdPmax, fabs(dP));
        P[n] += dP; // seems better

        //// Correct chopmin
        // for (USI i = 0; i < numCom; i++) {
        //    if (Ni[n * numCom + i] + u[n * col + 1 + i] < 0) {
        //        chopmin = 0.9 * min(chopmin, fabs(Ni[n * numCom + i] / u[n * col + 1 +
        //        i]));

        //        //if (chopmin < 0 || !isfinite(chopmin)) {
        //        //    OCP_ABORT("Wrong Chop!");
        //        //}
        //    }
        //}

        for (USI i = 0; i < numCom; i++)
        {
            Ni[n * numCom + i] += u[fn * col + 1 + i] * chopmin;

            // if (Ni[n * numCom + i] < 0) {
            //    cout << Ni[n * numCom + i] << "  " << u[n * col + 1 + i] * chopmin <<
            //    "   " << chopmin << endl;
            //}
        }
    }
}

void Bulk::CalRelResAIMs(ResFIM& resFIM) const
{
    OCP_FUNCNAME;

    // OCP_USI tmpid01 = -1;
    // OCP_USI tmpid02 = -1;
    OCP_DBL tmp;

    const USI len = numCom + 1;
    for (OCP_USI fn = 0; fn < numFIMBulk; fn++)
    {
        OCP_USI n = FIMBulk[fn];

        for (USI i = 0; i < len; i++)
        {
            tmp = fabs(resFIM.res[n * len + i] / rockVp[n]);
            if (resFIM.maxRelRes_v < tmp)
            {
                resFIM.maxRelRes_v = tmp;
                // tmpid01            = n;
            }
        }
        for (USI i = 1; i < len; i++)
        {
            tmp = fabs(resFIM.res[n * len + i] / Nt[n]);
            if (resFIM.maxRelRes_mol < tmp)
            {
                resFIM.maxRelRes_mol = tmp;
                // tmpid02              = n;
            }
        }
    }
}


void Bulk::GetSolAIMs(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
    const OCP_DBL& dSmaxlim)
{
    OCP_FUNCNAME;

    NRdSmax = 0;
    NRdPmax = 0;
    OCP_DBL dP;
    USI row = numPhase * (numCom + 1);
    USI col = numCom + 1;
    USI bsize = row * col;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL chopmin = 1;
    OCP_DBL choptmp = 0;
    OCP_INT bIde;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        bIde = map_Bulk2FIM[n];
        
        if (bIde < 0 || bIde >= numFIMBulk) {
            // IMPEC Bulk
            // Method 1
            P[n] += u[n * col];
            // Method 2
            // P[n] = u[n * col];
            if (P[n] < 0) {
                cout << "";
            }
        }
        else {
            // FIM Bulk
            chopmin = 1;
            // compute the chop
            fill(dtmp.begin(), dtmp.end(), 0.0);
            DaAxpby(row, col, 1, dSec_dPri.data() + bIde * bsize, u.data() + bIde * col, 1,
                dtmp.data());

            for (USI j = 0; j < numPhase; j++)
            {
                choptmp = 1;
                if (fabs(dtmp[j]) > dSmaxlim)
                {
                    choptmp = dSmaxlim / fabs(dtmp[j]);
                }
                else if (S[n * numPhase + j] + dtmp[j] < 0.0)
                {
                    choptmp = 0.9 * S[n * numPhase + j] / fabs(dtmp[j]);
                }

                chopmin = min(chopmin, choptmp);
                NRdSmax = max(NRdSmax, choptmp * fabs(dtmp[j]));
            }
            dP = u[n * col];
            choptmp = dPmaxlim / fabs(dP);
            chopmin = min(chopmin, choptmp);
            dP *= chopmin;
            NRdPmax = max(NRdPmax, fabs(dP));
            P[n] += dP; // seems better

            //// Correct chopmin
            // for (USI i = 0; i < numCom; i++) {
            //    if (Ni[n * numCom + i] + u[n * col + 1 + i] < 0) {
            //        chopmin = 0.9 * min(chopmin, fabs(Ni[n * numCom + i] / u[n * col + 1 +
            //        i]));

            //        //if (chopmin < 0 || !isfinite(chopmin)) {
            //        //    OCP_ABORT("Wrong Chop!");
            //        //}
            //    }
            //}

            for (USI i = 0; i < numCom; i++)
            {
                Ni[n * numCom + i] += u[n * col + 1 + i] * chopmin;

                // if (Ni[n * numCom + i] < 0) {
                //    cout << Ni[n * numCom + i] << "  " << u[n * col + 1 + i] * chopmin <<
                //    "   " << chopmin << endl;
                //}
            }
        }
    }
}


void Bulk::UpdateLastStepAIM()
{
    lmuP = muP;
    lxiP = xiP;
    lrhoP = rhoP;
    lmux = mux;
    lxix = xix;
    lrhox = rhox;
    ldPcj_dS = dPcj_dS;
    ldKr_dS = dKr_dS;
    ldSec_dPri = dSec_dPri;
}


void Bulk::ResetFIMBulk()
{
    muP = lmuP;
    xiP = lxiP;
    rhoP = lrhoP;
    mux = lmux;
    xix = lxix;
    rhox = lrhox;
    dPcj_dS = ldPcj_dS;
    dKr_dS = ldKr_dS;
    dSec_dPri = ldSec_dPri;
}


void Bulk::ShowFIMBulk(const bool& flag) const
{
    cout << numFIMBulk << "   " << fixed << setprecision(3)
        << numFIMBulk * 100.0 / numBulk << "%" << endl;
    if (flag) {
        for (USI n = 0; n < numFIMBulk; n++) {
            cout << setw(6) << FIMBulk[n] << "   ";
            if ((n + 1) % 10 == 0) {
                cout << endl;
            }
        }
        cout << endl;
    } 
    if (false) {
        for (USI n = 0; n < numBulk; n++) {
            cout << setw(6) << map_Bulk2FIM[n] << "   ";
            if ((n + 1) % 10 == 0) {
                cout << endl;
            }
        }
        cout << endl;
    }
}


bool Bulk::CheckNiFIMBulk() const
{
    OCP_FUNCNAME;

    OCP_USI len = numBulk * numCom;
    for (OCP_USI n = 0; n < len; n++)
    {
        if (Ni[n] < 0.0)
        {
            OCP_USI bId = n / numCom;
            USI cId = n - bId * numCom;
            std::ostringstream NiStringSci;
            NiStringSci << std::scientific << Ni[n];
            OCP_WARNING("Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                std::to_string(bId) + "] = " + NiStringSci.str());
            return false;
        }
    }
    return true;
}


void Bulk::InFIMNi()
{
    OCP_USI bIdb, bIdf;
    for (OCP_USI n = 0; n < numFIMBulk; n++) {
        bIdf = n * numCom;
        bIdb = FIMBulk[n] * numCom;
        for (USI i = 0; i < numCom; i++) {
            FIMNi[bIdf + i] = Ni[bIdb + i];
        }
    }
}


void Bulk::OutFIMNi()
{
    OCP_USI bIdb, bIdf;
    for (OCP_USI n = 0; n < numFIMBulk; n++) {
        bIdf = n * numCom;
        bIdb = FIMBulk[n] * numCom;
        for (USI i = 0; i < numCom; i++) {
            Ni[bIdb + i] = FIMNi[bIdf + i];           
        }
    }
}

void Bulk::AllocateAuxAIMc()
{
    cfl.resize(numBulk * numPhase);
    map_Bulk2FIM.resize(numBulk, -1);
}

void Bulk::FlashAIMc()
{
    if (comps) {
        FlashCOMPAIMc();
    }
    else {
        FlashBLKOILAIMc();
    }
}

void Bulk::FlashBLKOILAIMc()
{
    for (OCP_USI n = 0; n < numBulk; n++) 
    {
        if (map_Bulk2FIM[n] > -1) {
            // FIM bulk
            continue;
        }

        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], 0, 0, 0);
        PassFlashValueAIMc(n);
    }
}

void Bulk::FlashCOMPAIMc()
{
    USI ftype;
    OCP_USI bId;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    // cout << endl << "==================================" << endl;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (map_Bulk2FIM[n] > -1) {
            // FIM bulk
            continue;
        }

        ftype = 1;
        if (flagSkip[n])
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10)
            {
                ftype = 0;
            }
            // cout << setprecision(2) << scientific << minEig / 10 << "   " << fabs(1 - lP[n] / P[n]) << "   ";
            if (ftype == 1)
            {
                bId = n * numCom;
                Ntw = Nt[n] - Ni[bId + numCom - 1];
                for (USI i = 0; i < numCom - 1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }
        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValueAIMc(n);
    }
}

void Bulk::FlashAIMc01()
{
    if (comps) {
        FlashCOMPAIMc01();
    }
    else {
        FlashBLKOILAIMc01();
    }
}

void Bulk::FlashBLKOILAIMc01()
{
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (map_Bulk2FIM[n] > -1) {
            // FIM bulk
            continue;
        }

        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], 0, 0, 0);
        PassFlashValue(n);
    }
}

void Bulk::FlashCOMPAIMc01()
{
    USI ftype;
    OCP_USI bId;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    // cout << endl << "==================================" << endl;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (map_Bulk2FIM[n] > -1) {
            // FIM bulk
            continue;
        }

        ftype = 1;
        if (flagSkip[n])
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10)
            {
                ftype = 0;
            }
            // cout << setprecision(2) << scientific << minEig / 10 << "   " << fabs(1 - lP[n] / P[n]) << "   ";
            if (ftype == 1)
            {
                bId = n * numCom;
                Ntw = Nt[n] - Ni[bId + numCom - 1];
                for (USI i = 0; i < numCom - 1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }
        flashCal[PVTNUM[n]]->Flash(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValue(n);
    }
}

void Bulk::FlashDerivAIMc()
{
    if (comps) {
        FlashDerivCOMPAIMc();
    }
    else {
        FlashDerivBLKOILAIMc();
    }
}

void Bulk::FlashDerivBLKOILAIMc()
{
    // dSec_dPri.clear();
    for (auto& n : FIMBulk) {
        flashCal[PVTNUM[n]]->FlashDeriv(P[n], T, &Ni[n * numCom], 0, 0, 0);
        PassFlashValueDeriv(n);
    }
}

void Bulk::FlashDerivCOMPAIMc()
{
    USI ftype;
    OCP_USI bId;
    OCP_DBL Ntw;
    OCP_DBL minEig;
    // dSec_dPri.clear();
    for (auto& n : FIMBulk)
    {
        ftype = 1;
        if (flagSkip[n])
        {
            minEig = minEigenSkip[n];
            if (fabs(1 - PSkip[n] / P[n]) >= minEig / 10)
            {
                ftype = 0;
            }
            // cout << setprecision(2) << scientific << minEig / 10 << "   " << fabs(1 - lP[n] / P[n]) << "   ";
            if (ftype == 1)
            {
                bId = n * numCom;
                Ntw = Nt[n] - Ni[bId + numCom - 1];
                for (USI i = 0; i < numCom - 1; i++)
                {
                    // cout << fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) << "   ";
                    if (fabs(Ni[bId + i] / Ntw - ziSkip[bId + i]) >= minEig / 10)
                    {
                        ftype = 0;
                        break;
                    }
                }
            }
            // cout << n << endl;
        }
        else {
            ftype = 0;
        }

        flashCal[PVTNUM[n]]->FlashDeriv(P[n], T, &Ni[n * numCom], ftype, phaseNum[n], &Ks[n * numCom_1]);
        PassFlashValueDeriv(n);
        // lphaseNum[n];
    }
}

void Bulk::CalKrPcAIMc()
{
    OCP_DBL tmp = 0;
    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (map_Bulk2FIM[n] > -1) {
            // FIM bulk
            continue;
        }
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId], 0, tmp, tmp);
        for (USI j = 0; j < numPhase; j++)
            Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
    }
}

void Bulk::CalKrPcDerivAIMc()
{
    for (auto& n : FIMBulk)
    {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPcDeriv(&S[bId], &kr[bId], &Pc[bId],
            &dKr_dS[bId * numPhase],
            &dPcj_dS[bId * numPhase]);
        for (USI j = 0; j < numPhase; j++)
            Pj[bId + j] = P[n] + Pc[bId + j];
    }
}

void Bulk::GetSolAIMc(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
    const OCP_DBL& dSmaxlim)
{
    NRdSmax = 0;
    NRdPmax = 0;
    OCP_DBL dP;
    USI row = numPhase * (numCom + 1);
    USI col = numCom + 1;
    USI bsize = row * col;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL chopmin = 1;
    OCP_DBL choptmp = 0;

    for (OCP_USI n = 0; n < numBulk; n++)
    {
        if (map_Bulk2FIM[n] < 0) {
            // IMPEC Bulk
            P[n] += u[n * col];
            for (USI i = 0; i < numCom; i++)
            {
                Ni[n * numCom + i] += u[n * col + 1 + i];

                // if (Ni[n * numCom + i] < 0) {
                //    cout << Ni[n * numCom + i] << "  " << u[n * col + 1 + i] * chopmin <<
                //    "   " << chopmin << endl;
                //}
            }
            continue;
        }
        
        // FIM Bulk

        chopmin = 1;

        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);
        DaAxpby(row, col, 1, dSec_dPri.data() + n * bsize, u.data() + n * col, 1,
            dtmp.data());

        for (USI j = 0; j < numPhase; j++)
        {

            choptmp = 1;
            if (fabs(dtmp[j]) > dSmaxlim)
            {
                choptmp = dSmaxlim / fabs(dtmp[j]);
            }
            else if (S[n * numPhase + j] + dtmp[j] < 0.0)
            {
                choptmp = 0.9 * S[n * numPhase + j] / fabs(dtmp[j]);
            }

            chopmin = min(chopmin, choptmp);
            NRdSmax = max(NRdSmax, choptmp * fabs(dtmp[j]));
        }
        dP = u[n * col];
        choptmp = dPmaxlim / fabs(dP);
        chopmin = min(chopmin, choptmp);
        NRdPmax = max(NRdPmax, fabs(dP));
        P[n] += dP; // seems better

        //// Correct chopmin
        // for (USI i = 0; i < numCom; i++) {
        //    if (Ni[n * numCom + i] + u[n * col + 1 + i] < 0) {
        //        chopmin = 0.9 * min(chopmin, fabs(Ni[n * numCom + i] / u[n * col + 1 +
        //        i]));

        //        //if (chopmin < 0 || !isfinite(chopmin)) {
        //        //    OCP_ABORT("Wrong Chop!");
        //        //}
        //    }
        //}

        for (USI i = 0; i < numCom; i++)
        {
            Ni[n * numCom + i] += u[n * col + 1 + i] * chopmin;

            // if (Ni[n * numCom + i] < 0) {
            //    cout << Ni[n * numCom + i] << "  " << u[n * col + 1 + i] * chopmin <<
            //    "   " << chopmin << endl;
            //}
        }
    }
}

void Bulk::UpdatePj()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        for (USI j = 0; j < numPhase; j++) {
            Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
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
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/