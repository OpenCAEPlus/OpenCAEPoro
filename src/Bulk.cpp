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

// OpenCAEPoro header files
#include "Bulk.hpp"

/////////////////////////////////////////////////////////////////////
// HeatLoss
/////////////////////////////////////////////////////////////////////

void HeatLoss::InputParam(const HLoss& loss)
{
    ifHLoss = loss.ifHLoss;
    if (ifHLoss) {
        obC = loss.obC;
        obK = loss.obK;
        ubC = loss.ubC;
        ubK = loss.ubK;
    }
}

void HeatLoss::Setup(const OCP_USI& nb)
{
    if (ifHLoss) {
        obD     = obK / obC;
        ubD     = ubK / ubC;
        numBulk = nb;
        I.resize(numBulk);
        p.resize(numBulk);
        pT.resize(numBulk);
        lI.resize(numBulk);
        lp.resize(numBulk);
        lpT.resize(numBulk);
    }
}

void HeatLoss::CalHeatLoss(const vector<USI>&     location,
                           const vector<OCP_DBL>& T,
                           const vector<OCP_DBL>& lT,
                           const vector<OCP_DBL>& initT,
                           const OCP_DBL&         t,
                           const OCP_DBL&         dt)
{
    if (ifHLoss) {
        OCP_DBL lambda, d, dT, theta;
        OCP_DBL tmp, q;
        for (OCP_USI n = 0; n < numBulk; n++) {
            if (location[n] > 0) {
                // overburden or underburden
                lambda = location[n] == 1 ? obD : ubD;

                dT    = T[n] - lT[n];
                theta = T[n] - initT[n];
                d     = sqrt(lambda * t) / 2;
                tmp   = 3 * pow(d, 2) + lambda * dt;
                p[n]  = (theta * (lambda * dt / d) + lI[n] -
                        dT * (pow(d, 3) / (lambda * dt))) /
                       tmp;
                pT[n] = (lambda * dt / d - pow(d, 3) / (lambda * dt)) / tmp;
                q     = (2 * p[n] * d - theta + pow(d, 2) * dT / (lambda * dt)) /
                    (2 * pow(d, 2));
                I[n] = theta * d + p[n] * pow(d, 2) + 2 * q * pow(d, 3);
            }
        }
    }
}

void HeatLoss::ResetToLastTimeStep()
{
    I  = lI;
    p  = lp;
    pT = lpT;
}

void HeatLoss::UpdateLastTimeStep()
{
    lI  = I;
    lp  = p;
    lpT = pT;
}

/////////////////////////////////////////////////////////////////////
// Input Param and Setup
/////////////////////////////////////////////////////////////////////

/// Read parameters from rs_param data structure.
void Bulk::InputParam(const ParamReservoir& rs_param)
{
    OCP_FUNCNAME;

    // Common input
    ifBlackOil = rs_param.blackOil;
    ifComps    = rs_param.comps;
    ifThermal  = rs_param.thermal;

    if (ifThermal) {
        ifBlackOil = OCP_FALSE;
        ifComps    = OCP_FALSE;
    }

    NTPVT  = rs_param.NTPVT;
    NTSFUN = rs_param.NTSFUN;
    NTROCC = rs_param.NTROOC;

    if (rs_param.PBVD_T.data.size() > 0) EQUIL.PBVD.Setup(rs_param.PBVD_T.data[0]);

    if (ifBlackOil) {
        // Isothermal blackoil model
        InputParamBLKOIL(rs_param);
    } else if (ifComps) {
        // Isothermal compositional model
        InputParamCOMPS(rs_param);
    } else if (ifThermal) {
        // ifThermal model
        InputParamTHERMAL(rs_param);
    }

    InputSatFunc(rs_param);
}

void Bulk::InputParamBLKOIL(const ParamReservoir& rs_param)
{
    oil    = rs_param.oil;
    gas    = rs_param.gas;
    water  = rs_param.water;
    disGas = rs_param.disGas;

    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];

    if (water && !oil && !gas) {
        // water
        numPhase = 1;
        numCom   = 1;
        SATmode  = PHASE_W;
        PVTmodeB = PHASE_W;
    } else if (water && oil && !gas) {
        // water, dead oil
        numPhase   = 2;
        numCom     = 2;
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        SATmode    = PHASE_OW;
        PVTmodeB   = PHASE_OW;
    } else if (water && oil && gas && !disGas) {
        // water, dead oil, dry gas
        numPhase   = 3;
        numCom     = 3;
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];
        SATmode    = PHASE_DOGW;
        PVTmodeB   = PHASE_DOGW; // maybe it should be added later
    } else if (water && oil && gas && disGas) {
        // water, live oil, dry gas
        numPhase = 3;
        numCom   = 3;

        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];
        PVTmodeB   = PHASE_ODGW;

        if (rs_param.SOF3_T.data.size() > 0) {
            SATmode = PHASE_ODGW02;
        } else {
            SATmode = PHASE_ODGW01;
        }
    }
    numComH = numCom - 1;

    // PVT mode
    switch (PVTmodeB) {
        case PHASE_W:
            OCP_ABORT("Wrong Type!");
            break;
        case PHASE_OW:
            for (USI i = 0; i < NTPVT; i++)
                flashCal.push_back(new BOMixture_OW(rs_param, i));
            break;
        case PHASE_DOGW:
            OCP_ABORT("Wrong Type!");
            break;
        case PHASE_ODGW:
            for (USI i = 0; i < NTPVT; i++)
                flashCal.push_back(new BOMixture_ODGW(rs_param, i));
            break;
        default:
            OCP_ABORT("Wrong Type!");
            break;
    }

    phase2Index.resize(3);
    switch (flashCal[0]->GetMixtureType()) {
        case BLKOIL_W:
            phase2Index[WATER] = 0;
            break;
        case BLKOIL_OW:
            phase2Index[OIL]   = 0;
            phase2Index[WATER] = 1;
            break;
        case BLKOIL_OG:
            phase2Index[OIL] = 0;
            phase2Index[GAS] = 1;
            break;
        case BLKOIL_DOGW:
        case BLKOIL_ODGW:
            phase2Index[OIL]   = 0;
            phase2Index[GAS]   = 1;
            phase2Index[WATER] = 2;
            break;
        default:
            OCP_ABORT("WRONG Mixture Type!");
    }

    InputRockFunc(rs_param);
    cout << endl << "BLACKOIL model is selected" << endl;
}

void Bulk::InputParamCOMPS(const ParamReservoir& rs_param)
{

    // Water exists and is excluded in EoS model NOW!
    oil      = OCP_TRUE;
    gas      = OCP_TRUE;
    water    = OCP_TRUE;
    ifUseEoS = OCP_TRUE;

    numPhase   = rs_param.comsParam.numPhase + 1;
    numCom     = rs_param.comsParam.numCom + 1;
    numComH    = numCom - 1;
    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    EQUIL.DOWC = rs_param.EQUIL[2];
    EQUIL.PcOW = rs_param.EQUIL[3];
    EQUIL.DGOC = rs_param.EQUIL[4];
    EQUIL.PcGO = rs_param.EQUIL[5];

    // Init Zi
    for (auto& v : rs_param.ZMFVD_T.data) {
        initZi_Tab.push_back(OCPTable(v));
    }

    // Init T
    // Use RTEMP
    rsTemp = rs_param.rsTemp;
    vector<vector<OCP_DBL>> temp;
    temp.resize(2);
    // add depth
    temp[0].push_back(0);
    temp[0].push_back(1E8);
    // add temperature
    temp[1].push_back(rsTemp);
    temp[1].push_back(rsTemp);
    initT_Tab.push_back(OCPTable(temp));

    // Saturation mode
    if (rs_param.SOF3_T.data.size() > 0) {
        SATmode = PHASE_ODGW02;
    } else {
        SATmode = PHASE_ODGW01;
        if (rs_param.comsParam.miscible) {
            SATmode = PHASE_ODGW01_MISCIBLE;
        }
    }

    // PVT mode
    for (USI i = 0; i < NTPVT; i++) flashCal.push_back(new MixtureComp(rs_param, i));

    // phase index
    phase2Index.resize(3);
    phase2Index[OIL]   = 0;
    phase2Index[GAS]   = 1;
    phase2Index[WATER] = 2;

    InputRockFunc(rs_param);
    cout << endl << "COMPOSITIONAL model is selected" << endl;
}

void Bulk::InputParamTHERMAL(const ParamReservoir& rs_param)
{
    oil   = rs_param.oil;
    gas   = rs_param.gas;
    water = rs_param.water;

    // Init T
    rsTemp = rs_param.rsTemp;
    for (auto& v : rs_param.TEMPVD_T.data) {
        initT_Tab.push_back(OCPTable(v));
    }
    if (initT_Tab.size() == 0) {
        // Use RTEMP
        vector<vector<OCP_DBL>> temp;
        temp.resize(2);
        // add depth
        temp[0].push_back(0);
        temp[0].push_back(1E8);
        // add temperature
        temp[1].push_back(rsTemp);
        temp[1].push_back(rsTemp);
        initT_Tab.push_back(OCPTable(temp));
    }
    // ifThermal conductivity
    if (oil) {
        thconp.push_back(rs_param.thcono);
    }
    if (gas) {
        thconp.push_back(rs_param.thcong);
    }
    if (water) {
        thconp.push_back(rs_param.thconw);
    }

    // Only Now
    SATmode    = PHASE_OW;
    numPhase   = 2;
    numCom     = 2;
    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];
    EQUIL.DOWC = rs_param.EQUIL[2];
    EQUIL.PcOW = rs_param.EQUIL[3];
    EQUIL.DGOC = rs_param.EQUIL[4];
    EQUIL.PcGO = rs_param.EQUIL[5];

    // PVT mode
    for (USI i = 0; i < NTPVT; i++)
        flashCal.push_back(new MixtureThermal_K01(rs_param, i));

    // phase index
    phase2Index.resize(3);
    phase2Index[OIL]   = 0;
    phase2Index[WATER] = 1;

    InputRockFuncT(rs_param);
}

void Bulk::InputSatFunc(const ParamReservoir& rs_param)
{
    // Setup Saturation function
    satcm.resize(NTSFUN);
    switch (SATmode) {
        case PHASE_W:
            for (USI i = 0; i < NTSFUN; i++)
                flow.push_back(new FlowUnit_W(rs_param, i));
            break;
        case PHASE_OW:
            for (USI i = 0; i < NTSFUN; i++)
                flow.push_back(new FlowUnit_OW(rs_param, i));
            break;
        case PHASE_ODGW01:
            for (USI i = 0; i < NTSFUN; i++) {
                flow.push_back(new FlowUnit_ODGW01(rs_param, i));
                satcm[i] = flow[i]->GetScm();
            }
            break;
        case PHASE_ODGW01_MISCIBLE:
            for (USI i = 0; i < NTSFUN; i++) {
                flow.push_back(new FlowUnit_ODGW01_Miscible(rs_param, i));
                satcm[i] = flow[i]->GetScm();
            }
            break;
        case PHASE_ODGW02:
            for (USI i = 0; i < NTSFUN; i++)
                flow.push_back(new FlowUnit_ODGW02(rs_param, i));
            break;
        default:
            OCP_ABORT("Wrong Type!");
            break;
    }
}

void Bulk::InputRockFunc(const ParamReservoir& rs_param)
{
    for (USI i = 0; i < NTROCC; i++) {
        rock.push_back(new Rock_Linear(rs_param.rockSet[i]));
    }
}

void Bulk::InputRockFuncT(const ParamReservoir& rs_param)
{
    for (USI i = 0; i < NTROCC; i++) {
        if (rs_param.rockSet[i].type == "LINEAR") {
            rock.push_back(new RockT_Linear(rs_param.rockSet[i]));
        } else {
            rock.push_back(new RockT_Exp(rs_param.rockSet[i]));
        }
    }

    // Heat Loss
    hLoss.InputParam(rs_param.hLoss);
}

/// Setup bulk information.
void Bulk::SetupIsoT(const Grid& myGrid)
{
    OCP_FUNCNAME;

    numBulk = myGrid.activeGridNum;
    AllocateGridRockIsoT(myGrid);
    AllocateRegion(myGrid);
    AllocateError();
}

/// Allocate memory for fluid grid for ifThermal model
void Bulk::SetupT(const Grid& myGrid)
{
    numBulk = myGrid.activeGridNum;
    AllocateGridRockT(myGrid);
    AllocateRegion(myGrid);
    SetupBulkType(myGrid);
    // Setup Heat Loss
    hLoss.Setup(numBulk);
}

void Bulk::SetupOptionalFeatures(const Grid& myGrid, OptionalFeatures& optFeatures)
{
    for (USI i = 0; i < NTPVT; i++) {
        flashCal[i]->SetupOptionalFeatures(optFeatures, numBulk);
    }
    for (USI i = 0; i < NTSFUN; i++) {
        flow[i]->SetupOptionalFeatures(myGrid, optFeatures);
    }
}

/////////////////////////////////////////////////////////////////////
// Initial Properties
/////////////////////////////////////////////////////////////////////

void Bulk::InitPTSw(const USI& tabrow)
{
    OCP_FUNCNAME;

    initT.resize(numBulk);

    OCP_DBL Dref = EQUIL.Dref;
    OCP_DBL Pref = EQUIL.Pref;
    OCP_DBL DOWC = EQUIL.DOWC;
    OCP_DBL PcOW = EQUIL.PcOW;
    OCP_DBL DOGC = EQUIL.DGOC;
    OCP_DBL PcGO = EQUIL.PcGO;
    OCP_DBL Zmin = 1E8;
    OCP_DBL Zmax = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_DBL temp1 = depth[n] - dz[n] / 2;
        OCP_DBL temp2 = depth[n] + dz[n] / 2;
        Zmin          = Zmin < temp1 ? Zmin : temp1;
        Zmax          = Zmax > temp2 ? Zmax : temp2;
    }
    OCP_DBL tabdz = (Zmax - Zmin) / (tabrow - 1);

    // create table
    OCPTable         DepthP(tabrow, 4);
    vector<OCP_DBL>& Ztmp  = DepthP.GetCol(0);
    vector<OCP_DBL>& Potmp = DepthP.GetCol(1);
    vector<OCP_DBL>& Pgtmp = DepthP.GetCol(2);
    vector<OCP_DBL>& Pwtmp = DepthP.GetCol(3);

    vector<OCP_DBL> tmpInitZi(numCom, 0);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++) {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    OCP_DBL myTemp = rsTemp;

    // find the RefId
    USI beginId = 0;
    if (Dref <= Ztmp[0]) {
        beginId = 0;
    } else if (Dref >= Ztmp[tabrow - 1]) {
        beginId = tabrow - 1;
    } else {
        beginId =
            distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(),
                                           [s = Dref](auto& t) { return t > s; }));
        beginId--;
    }

    // begin calculating oil pressure:
    OCP_DBL Pbb = Pref;
    OCP_DBL gammaOtmp, gammaWtmp, gammaGtmp;
    OCP_DBL Ptmp;
    USI     mynum = 10;
    OCP_DBL mydz  = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    const OCP_BOOL initZi_flag = initZi_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const OCP_BOOL initT_flag  = initT_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const OCP_BOOL PBVD_flag   = EQUIL.PBVD.IsEmpty() ? OCP_FALSE : OCP_TRUE;

    if (Dref < DOGC) {
        // reference pressure is gas pressure
        Pgref = Pref;
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaGtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Pgref, Pbb, myTemp, tmpInitZi.data(), GAS);
        Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), GAS);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), GAS);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref       = 0;
        Ptmp        = Pgref;
        mydz        = (DOGC - Dref) / mynum;
        OCP_DBL myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaGtmp = GRAVITY_FACTOR *
                        flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), GAS);
            Ptmp += gammaGtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcGO;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), OIL);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Poref, Pbb, myTemp, tmpInitZi.data(), OIL);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp  = Poref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Poref, Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi.data(), WATER);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    } else if (Dref > DOWC) {
        OCP_DBL myz;
        // reference pressure is water pressure
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        Pwref     = Pref;
        gammaWtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi.data(), WATER);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp  = Pwref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Ptmp += gammaWtmp * mydz;
            myz += mydz;
        }
        Ptmp += PcOW;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), OIL);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Poref, Pbb, myTemp, tmpInitZi.data(), OIL);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp  = Poref;
            mydz  = (DOGC - Dref) / mynum;
            myz   = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp =
                    GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), OIL);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp =
                    GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), GAS);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp      = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgref, Pbb, myTemp,
                                                                    tmpInitZi.data(), GAS);
            Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi.data(), GAS);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }
            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi.data(), GAS);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

    } else {
        OCP_DBL myz;
        // reference pressure is oil pressure
        Poref = Pref;
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Poref, Pbb, myTemp, tmpInitZi.data(), OIL);
        Pbegin         = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp  = Poref;
            mydz  = (DOGC - Dref) / mynum;
            myz   = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp =
                    GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), OIL);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp =
                    GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), GAS);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp      = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgref, Pbb, myTemp,
                                                                    tmpInitZi.data(), GAS);
            Pbegin         = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi.data(), GAS);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }

            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp =
                    GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp,
                                                           tmpInitZi.data(), GAS);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp  = Poref;
        mydz  = (DOWC - Dref) / mynum;
        myz   = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag) Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR *
                        flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, tmpInitZi.data(), OIL);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag) initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag) myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR *
                    flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, tmpInitZi.data(), WATER);
        Pbegin         = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag) initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag) myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp,
                                                               tmpInitZi.data(), WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0), cdata(4, 0);
    // whether capillary between water and oil is considered
    vector<OCP_BOOL> FlagPcow(NTSFUN, OCP_TRUE);
    for (USI i = 0; i < NTSFUN; i++) {
        if (fabs(flow[i]->GetPcowBySw(0.0 - TINY)) < TINY &&
            fabs(flow[i]->GetPcowBySw(1.0 + TINY) < TINY)) {
            FlagPcow[i] = OCP_FALSE;
        }
    }

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (initZi_flag) {
            initZi_Tab[0].Eval_All0(depth[n], tmpInitZi);
            for (USI i = 0; i < numComH; i++) {
                Ni[n * numCom + i] = tmpInitZi[i];
            }
        }
        if (initT_flag) {
            myTemp   = initT_Tab[0].Eval(0, depth[n], 1);
            initT[n] = myTemp;
            T[n]     = myTemp;
        }

        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po   = data[1];
        OCP_DBL Pg   = data[2];
        OCP_DBL Pw   = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw   = flow[SATNUM[n]]->GetSwByPcow(Pcow);
        OCP_DBL Sg   = 0;
        if (gas) {
            Sg = flow[SATNUM[n]]->GetSgByPcgo(Pcgo);
        }
        if (Sw + Sg > 1) {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw           = flow[SATNUM[n]]->GetSwByPcgw(Pcgw);
            Sg           = 1 - Sw;
        }

        if (1 - Sw < TINY) {
            // all water
            Po = Pw + flow[SATNUM[n]]->GetPcowBySw(1.0);
        } else if (1 - Sg < TINY) {
            // all gas
            Po = Pg - flow[SATNUM[n]]->GetPcgoBySg(1.0);
        } else if (1 - Sw - Sg < TINY) {
            // water and gas
            Po = Pg - flow[SATNUM[n]]->GetPcgoBySg(Sg);
        }
        P[n] = Po;

        if (depth[n] < DOGC) {
            Pbb = Po;
        } else if (PBVD_flag) {
            Pbb = EQUIL.PBVD.Eval(0, depth[n], 1);
        }
        Pb[n] = Pbb;

        // cal Sw
        OCP_DBL swco = flow[SATNUM[n]]->GetSwco();
        if (!FlagPcow[SATNUM[n]]) {
            S[n * numPhase + numPhase - 1] = swco;
            continue;
        }

        Sw              = 0;
        Sg              = 0;
        USI     ncut    = 10;
        OCP_DBL avePcow = 0;

        for (USI k = 0; k < ncut; k++) {
            OCP_DBL tmpSw = 0;
            OCP_DBL tmpSg = 0;
            OCP_DBL dep   = depth[n] + dz[n] / ncut * (k - (ncut - 1) / 2.0);
            DepthP.Eval_All(0, dep, data, cdata);
            Po   = data[1];
            Pg   = data[2];
            Pw   = data[3];
            Pcow = Po - Pw;
            Pcgo = Pg - Po;
            avePcow += Pcow;
            tmpSw = flow[SATNUM[n]]->GetSwByPcow(Pcow);
            if (gas) {
                tmpSg = flow[SATNUM[n]]->GetSgByPcgo(Pcgo);
            }
            if (tmpSw + tmpSg > 1) {
                // should be modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw        = flow[SATNUM[n]]->GetSwByPcgw(Pcgw);
                tmpSg        = 1 - tmpSw;
            }
            Sw += tmpSw;
            // Sg += tmpSg;
        }
        Sw /= ncut;
        // Sg /= ncut;
        avePcow /= ncut;

        flow[SATNUM[n]]->SetupScale(n, Sw, avePcow);
        S[n * numPhase + numPhase - 1] = Sw;
    }
}

/////////////////////////////////////////////////////////////////////
// Optional Features
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Region
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateRegion(const Grid& myGrid)
{
    SATNUM.resize(numBulk, 0);
    PVTNUM.resize(numBulk, 0);
    ROCKNUM.resize(numBulk, 0);

    // Pass initial grid value
    for (OCP_USI bIda = 0; bIda < numBulk; bIda++) {
        OCP_USI bId = myGrid.map_Act2All[bIda];

        SATNUM[bIda]  = myGrid.SATNUM[bId];
        PVTNUM[bIda]  = myGrid.PVTNUM[bId];
        ROCKNUM[bIda] = myGrid.ROCKNUM[bId];
    }
}

void Bulk::SetupBulkType(const Grid& myGrid)
{
    bType.resize(numBulk, 0);

    OCP_USI count = 0;
    for (OCP_USI n = 0; n < myGrid.numGrid; n++) {
        if (myGrid.map_All2Act[n].IsAct()) {
            if (myGrid.map_All2Flu[n].IsAct()) {
                bType[count]++;
            }
            count++;
        }
    }
}

/////////////////////////////////////////////////////////////////////
// Basic PVT Model Information
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Basic Grid and Basic Rock Information
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateGridRockIsoT(const Grid& myGrid)
{
    dx.resize(numBulk, 0);
    dy.resize(numBulk, 0);
    dz.resize(numBulk, 0);
    v.resize(numBulk, 0);
    depth.resize(numBulk, 0);

    ntg.resize(numBulk, 0);
    poroInit.resize(numBulk, 0);
    rockKx.resize(numBulk, 0);
    rockKy.resize(numBulk, 0);
    rockKz.resize(numBulk, 0);

    for (OCP_USI bIda = 0; bIda < numBulk; bIda++) {
        OCP_USI bId = myGrid.map_Act2All[bIda];

        dx[bIda]    = myGrid.dx[bId];
        dy[bIda]    = myGrid.dy[bId];
        dz[bIda]    = myGrid.dz[bId];
        v[bIda]     = myGrid.v[bId];
        depth[bIda] = myGrid.depth[bId];

        ntg[bIda]      = myGrid.ntg[bId];
        poroInit[bIda] = myGrid.poro[bId];
        rockKx[bIda]   = myGrid.kx[bId];
        rockKy[bIda]   = myGrid.ky[bId];
        rockKz[bIda]   = myGrid.kz[bId];
    }
}

void Bulk::AllocateGridRockT(const Grid& myGrid)
{

    dx.resize(numBulk, 0);
    dy.resize(numBulk, 0);
    dz.resize(numBulk, 0);
    v.resize(numBulk, 0);
    depth.resize(numBulk, 0);

    ntg.resize(numBulk, 0);
    poroInit.resize(numBulk, 0);
    rockKx.resize(numBulk, 0);
    rockKy.resize(numBulk, 0);
    rockKz.resize(numBulk, 0);
    thconr.resize(numBulk, 0);

    bLocation.resize(numBulk, 0);

    for (OCP_USI bIda = 0; bIda < numBulk; bIda++) {
        OCP_USI bId = myGrid.map_Act2All[bIda];

        dx[bIda]    = myGrid.dx[bId];
        dy[bIda]    = myGrid.dy[bId];
        dz[bIda]    = myGrid.dz[bId];
        v[bIda]     = myGrid.v[bId];
        depth[bIda] = myGrid.depth[bId];

        ntg[bIda]      = myGrid.ntg[bId];
        poroInit[bIda] = myGrid.poro[bId];
        rockKx[bIda]   = myGrid.kx[bId];
        rockKy[bIda]   = myGrid.ky[bId];
        rockKz[bIda]   = myGrid.kz[bId];
        thconr[bIda]   = myGrid.thconr[bId];

        bLocation[bIda] = myGrid.gLocation[bId];
    }
}

/////////////////////////////////////////////////////////////////////
// Basic Fluid Information
/////////////////////////////////////////////////////////////////////

OCP_DBL Bulk::CalFPR() const
{
    OCP_FUNCNAME;

    OCP_DBL ptmp = 0;
    OCP_DBL vtmp = 0;
    OCP_DBL tmp  = 0;

    if (numPhase == 3) {
        for (OCP_USI n = 0; n < numBulk; n++) {
            tmp = rockVp[n] * (1 - S[n * numPhase + 2]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    } else if (numPhase < 3) {
        for (OCP_USI n = 0; n < numBulk; n++) {
            tmp = rockVp[n] * (S[n * numPhase]);
            ptmp += P[n] * tmp;
            vtmp += tmp;
        }
    } else {
        OCP_ABORT("Number of phases is out of range!");
    }
    return ptmp / vtmp;
}

OCP_DBL Bulk::CalFTR() const
{
    OCP_FUNCNAME;

    OCP_DBL Ttmp = 0;
    OCP_DBL vtmp = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        Ttmp += T[n] * v[n];
        vtmp += v[n];
    }

    return Ttmp / vtmp;
}

/////////////////////////////////////////////////////////////////////
// Important Indicator Variable and Check
/////////////////////////////////////////////////////////////////////

OCP_DBL Bulk::CalNRdSmax(OCP_USI& index)
{
    NRdSmax     = 0;
    OCP_USI len = numBulk * numPhase;
    for (USI n = 0; n < len; n++) {
        if (fabs(NRdSmax) < fabs(dSNR[n])) {
            NRdSmax = dSNR[n];
            index   = n;
        }
    }
    index /= numPhase;
    return NRdSmax;
}

/// Return OCP_TRUE if no negative pressure and OCP_FALSE otherwise.
OCP_INT Bulk::CheckP() const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (P[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << P[n];
            OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            cout << "P = " << P[n] << endl;
            return BULK_NEGATIVE_PRESSURE;
        }
    }

    return BULK_SUCCESS;
}

OCP_INT Bulk::CheckT() const
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (T[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << T[n];
            OCP_WARNING("Negative pressure: T[" + std::to_string(n) +
                        "] = " + PStringSci.str());
            cout << "T = " << T[n] << endl;
            return BULK_NEGATIVE_TEMPERATURE;
        }
    }

    return BULK_SUCCESS;
}

/// Return OCP_TRUE if no negative Ni and OCP_FALSE otherwise.
OCP_INT Bulk::CheckNi()
{
    OCP_FUNCNAME;

    OCP_USI len = numBulk * numCom;
    for (OCP_USI n = 0; n < len; n++) {
        if (Ni[n] < 0) {
            OCP_USI bId = n / numCom;
            if (Ni[n] > -1E-3 * Nt[bId] && OCP_FALSE) {
                Ni[n] = 1E-8 * Nt[bId];
            } else {
                USI                cId = n - bId * numCom;
                std::ostringstream NiStringSci;
                NiStringSci << std::scientific << Ni[n];
                OCP_WARNING("Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                            std::to_string(bId) + "] = " + NiStringSci.str() + ",  " +
                            "dNi = " + std::to_string(dNNR[n]));

                return BULK_NEGATIVE_COMPONENTS_MOLES;
            }
        }
    }
    return BULK_SUCCESS;
}

/// Return OCP_TRUE if all Ve < Vlim and OCP_FALSE otherwise.
OCP_INT Bulk::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    OCP_DBL dVe = 0.0;
    for (OCP_USI n = 0; n < numBulk; n++) {
        dVe = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (dVe > Vlim) {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << dVe
                 << " is too big!" << endl;
            return BULK_OUTRANGED_VOLUME_ERROR;
        }
    }
    return BULK_SUCCESS;
}

OCP_INT Bulk::CheckCFL(const OCP_DBL& cflLim) const
{
    if (maxCFL > cflLim)
        return BULK_OUTRANGED_CFL;
    else
        return BULK_SUCCESS;
}

void Bulk::CalMaxChange()
{
    OCP_FUNCNAME;

    dPmax       = 0;
    dTmax       = 0;
    dNmax       = 0;
    dSmax       = 0;
    eVmax       = 0;
    OCP_DBL tmp = 0;
    OCP_USI id;

    for (OCP_USI n = 0; n < numBulk; n++) {

        // dP
        tmp = fabs(P[n] - lP[n]);
        if (dPmax < tmp) {
            dPmax = tmp;
        }

        // dT
        tmp = fabs(T[n] - lT[n]);
        if (dTmax < tmp) {
            dTmax = tmp;
        }

        // dS
        for (USI j = 0; j < numPhase; j++) {
            id  = n * numPhase + j;
            tmp = fabs(S[id] - lS[id]);
            if (dSmax < tmp) {
                dSmax = tmp;
            }
        }

        // dN
        for (USI i = 0; i < numCom; i++) {
            id  = n * numCom + i;
            tmp = fabs(max(Ni[id], lNi[id]));
            if (tmp > TINY) {
                tmp = fabs(Ni[id] - lNi[id]) / tmp;
                if (dNmax < tmp) {
                    dNmax = tmp;
                }
            }
        }

        // Ve
        tmp = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (eVmax < tmp) {
            eVmax = tmp;
        }
    }
}

/////////////////////////////////////////////////////////////////////
// Error
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateError()
{
    if (ifUseEoS) {
        ePEC.resize(numBulk);
    }
}

/////////////////////////////////////////////////////////////////////
// For AIMc
/////////////////////////////////////////////////////////////////////

void Bulk::ShowFIMBulk(const OCP_BOOL& flag) const
{

    if (flag) {
        USI iter = 0;
        for (USI n = 0; n < numBulk; n++) {
            if (bulkTypeAIM.IfFIMbulk(n)) {
                // FIM bulk
                cout << setw(6) << n << "   ";
                iter++;
            }

            if ((iter + 1) % 10 == 0) {
                cout << endl;
                iter = 0;
            }
        }
        cout << endl;
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