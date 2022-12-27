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
// Input Param and Setup
/////////////////////////////////////////////////////////////////////

/// Read parameters from rs_param data structure.
void Bulk::InputParam(const ParamReservoir& rs_param)
{
    OCP_FUNCNAME;

    // Common input  
    ifBlackOil   = rs_param.blackOil;
    ifComps      = rs_param.comps;
    ifThermal    = rs_param.thermal;
    
    NTPVT      = rs_param.NTPVT;
    NTSFUN     = rs_param.NTSFUN;
    NTROCC     = rs_param.NTROOC;

    ifScalePcow = rs_param.ScalePcow;

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
    oil = rs_param.oil;
    gas = rs_param.gas;
    water = rs_param.water;
    disGas = rs_param.disGas;

    EQUIL.Dref = rs_param.EQUIL[0];
    EQUIL.Pref = rs_param.EQUIL[1];

    if (water && !oil && !gas) {
        // water
        numPhase = 1;
        numCom = 1;
        SATmode = PHASE_W;
        PVTmodeB = PHASE_W;
    }
    else if (water && oil && !gas) {
        // water, dead oil
        numPhase = 2;
        numCom = 2;
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        SATmode = PHASE_OW;
        PVTmodeB = PHASE_OW;
    }
    else if (water && oil && gas && !disGas) {
        // water, dead oil, dry gas
        numPhase = 3;
        numCom = 3;
        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];
        SATmode = PHASE_DOGW;
        PVTmodeB = PHASE_DOGW; // maybe it should be added later
    }
    else if (water && oil && gas && disGas) {
        // water, live oil, dry gas
        numPhase = 3;
        numCom = 3;

        EQUIL.DOWC = rs_param.EQUIL[2];
        EQUIL.PcOW = rs_param.EQUIL[3];
        EQUIL.DGOC = rs_param.EQUIL[4];
        EQUIL.PcGO = rs_param.EQUIL[5];
        PVTmodeB = PHASE_ODGW;

        if (rs_param.SOF3_T.data.size() > 0) {
            SATmode = PHASE_ODGW02;
        }
        else {
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
    switch (flashCal[0]->GetMixtureType())
    {
    case BLKOIL_W:
        phase2Index[WATER] = 0;
        break;
    case BLKOIL_OW:
        phase2Index[OIL] = 0;
        phase2Index[WATER] = 1;
        break;
    case BLKOIL_OG:
        phase2Index[OIL] = 0;
        phase2Index[GAS] = 1;
        break;
    case BLKOIL_DOGW:
    case BLKOIL_ODGW:
        phase2Index[OIL] = 0;
        phase2Index[GAS] = 1;
        phase2Index[WATER] = 2;
        break;
    default:
        OCP_ABORT("WRONG Mixture Type!");
    }

    InputRockFunc(rs_param);
    cout << "BLACKOIL model" << endl;
}


void Bulk::InputParamCOMPS(const ParamReservoir& rs_param)
{

    // Water exists and is excluded in EoS model NOW!
    oil = OCP_TRUE;
    gas = OCP_TRUE;
    water = OCP_TRUE;
    ifUseEoS = OCP_TRUE;

    numPhase = rs_param.comsParam.numPhase + 1;
    numCom = rs_param.comsParam.numCom + 1;
    numComH = numCom - 1;
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
    RTemp = rs_param.rsTemp;
    vector<vector<OCP_DBL>> temp;
    temp.resize(2);
    // add depth
    temp[0].push_back(0);
    temp[0].push_back(1E8);
    // add temperature
    temp[1].push_back(RTemp);
    temp[1].push_back(RTemp);
    initT_Tab.push_back(OCPTable(temp));

    // Saturation mode
    if (rs_param.SOF3_T.data.size() > 0) {
        SATmode = PHASE_ODGW02;
    }
    else {
        SATmode = PHASE_ODGW01;
        if (rs_param.comsParam.miscible){ SATmode = PHASE_ODGW01_MISCIBLE; }
    }

    // PVT mode
    for (USI i = 0; i < NTPVT; i++)
        flashCal.push_back(new MixtureComp(rs_param, i));

    // phase index
    phase2Index.resize(3);
    phase2Index[OIL]   = 0;
    phase2Index[GAS]   = 1;
    phase2Index[WATER] = 2;

    InputRockFunc(rs_param);
    cout << "COMPOSITIONAL model" << endl;
}


void Bulk::InputParamTHERMAL(const ParamReservoir& rs_param)
{
    oil = rs_param.oil;
    gas = rs_param.gas;
    water = rs_param.water;

    // Init T
    RTemp = rs_param.rsTemp;
    for (auto& v : rs_param.TEMPVD_T.data) {
        initT_Tab.push_back(OCPTable(v));
    }
    if (initT_Tab.size() == 0) {
        // Use RTEMP
        vector<vector<OCP_DBL>> temp;
        temp.resize(2);
        // add depth
        temp[0].push_back(0);       temp[0].push_back(1E8);
        // add temperature
        temp[1].push_back(RTemp);   temp[1].push_back(RTemp);
        initT_Tab.push_back(OCPTable(temp));
    }
    // ifThermal conductivity
    if (oil) { thconp.push_back(rs_param.thcono); }
    if (gas) { thconp.push_back(rs_param.thcong); }
    if (water) { thconp.push_back(rs_param.thconw); }


    // PVT mode
    for (USI i = 0; i < NTPVT; i++)
        flashCal.push_back(new MixtureThermal_K01(rs_param, i));
        

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
        }
        else {
            rock.push_back(new RockT_Exp(rs_param.rockSet[i]));
        }
    }
}


/// Setup bulk information.
void Bulk::SetupIsoT(const Grid& myGrid)
{
    OCP_FUNCNAME;

    numBulk = myGrid.activeGridNum;
    AllocateGridRockIsoT(myGrid);
    AllocateRegion(myGrid);
    AllocateSwatInit(myGrid);
    AllocateScalePcow();
    AllocateError();
    // CalSomeInfo(myGrid);

}

/// Allocate memory for fluid grid for ifThermal model
void Bulk::SetupT(const Grid& myGrid)
{
    numBulk = myGrid.fluidGridNum;
    AllocateGridRockIsoT(myGrid);
    AllocateRegion(myGrid);
    AllocateSwatInit(myGrid);
    AllocateScalePcow();
}

void Bulk::SetupOptionalFeatures(OptionalFeatures& optFeatures)
{
    for (USI i = 0; i < NTPVT; i++) {
        flashCal[i]->SetupOptionalFeatures(optFeatures, numBulk);
    }
    for (USI i = 0; i < NTSFUN; i++) {
        flow[i]->SetupOptionalFeatures(optFeatures, numBulk);
    }
}


/////////////////////////////////////////////////////////////////////
// Initial Properties
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateSwatInit(const Grid& myGrid)
{
    if (myGrid.SwatInit.size() > 0) {
        SwatInitExist = OCP_TRUE;
        SwatInit.resize(numBulk);

        for (OCP_USI bIda = 0; bIda < numBulk; bIda++) {
            OCP_USI bId = myGrid.map_Act2All[bIda];
            SwatInit[bIda] = myGrid.SwatInit[bId];
        }
    }
}


void Bulk::InitSjPc(const USI& tabrow)
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

    vector<OCP_DBL> tmpInitZi(numCom, 0);

    // cal Tab_Ztmp
    Ztmp[0] = Zmin;
    for (USI i = 1; i < tabrow; i++) {
        Ztmp[i] = Ztmp[i - 1] + tabdz;
    }

    OCP_DBL myTemp = RTemp;

    // find the RefId
    USI beginId = 0;
    if (Dref <= Ztmp[0]) {
        beginId = 0;
    }
    else if (Dref >= Ztmp[tabrow - 1]) {
        beginId = tabrow - 1;
    }
    else {
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
    OCP_DBL mydz = 0;
    OCP_DBL Poref, Pgref, Pwref;
    OCP_DBL Pbegin = 0;

    const OCP_BOOL initZi_flag = initZi_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const OCP_BOOL initT_flag = initT_Tab.size() > 0 ? OCP_TRUE : OCP_FALSE;
    const OCP_BOOL PBVD_flag = EQUIL.PBVD.IsEmpty() ? OCP_FALSE : OCP_TRUE;

    if (Dref < DOGC) {
        // reference pressure is gas pressure
        Pgref = Pref;
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag)
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgref, Pbb, myTemp, &tmpInitZi[0], GAS);
        Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
        Pgtmp[beginId] = Pbegin;

        // find the gas pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
            Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
            Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pgref
        Poref = 0;
        Ptmp = Pgref;
        mydz = (DOGC - Dref) / mynum;
        OCP_DBL myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], GAS);
            Ptmp += gammaGtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcGO;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], OIL);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag)
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Poref, Pbb, myTemp, &tmpInitZi[0], OIL);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;
        myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Poref, Pbb, myTemp, &tmpInitZi[0], OIL);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], WATER);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, &tmpInitZi[0], WATER);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }
    else if (Dref > DOWC) {
        OCP_DBL myz;
        // reference pressure is water pressure
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);

        Pwref = Pref;
        gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, &tmpInitZi[0], WATER);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        // find the water pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        // find the oil pressure in Dref by Pwref
        Poref = 0;
        Ptmp = Pwref;
        mydz = (DOWC - Dref) / mynum;
        myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], WATER);
            Ptmp += gammaWtmp * mydz;
            myz += mydz;
        }
        Ptmp += PcOW;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], OIL);
            Ptmp -= gammaOtmp * mydz;
            myz -= mydz;
        }
        Poref = Ptmp;

        // find the oil pressure in tab
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag)
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Poref, Pbb, myTemp, &tmpInitZi[0], OIL);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;
            myz = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], OIL);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], GAS);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgref, Pbb, myTemp, &tmpInitZi[0], GAS);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }
            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

    }
    else {
        OCP_DBL myz;
        // reference pressure is oil pressure
        Poref = Pref;
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);
        if (PBVD_flag)
            Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

        gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Poref, Pbb, myTemp, &tmpInitZi[0], OIL);
        Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
        Potmp[beginId] = Pbegin;

        // find the oil pressure
        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
        }
        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Potmp[id], Pbb, myTemp, &tmpInitZi[0], OIL);
            Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
        }

        if (gas) {
            // find the gas pressure in Dref by Poref
            Pgref = 0;
            Ptmp = Poref;
            mydz = (DOGC - Dref) / mynum;
            myz = Dref;

            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], OIL);
                Ptmp += gammaOtmp * mydz;
                myz += mydz;
            }
            Ptmp += PcGO;
            for (USI i = 0; i < mynum; i++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(myz, tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, myz, 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, myz, 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], GAS);
                Ptmp -= gammaGtmp * mydz;
                myz -= mydz;
            }
            Pgref = Ptmp;

            // find the gas pressure in tab
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Dref, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, Dref, 1);

            gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgref, Pbb, myTemp, &tmpInitZi[0], GAS);
            Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
            Pgtmp[beginId] = Pbegin;

            for (USI id = beginId; id > 0; id--) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
                Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
            }

            for (USI id = beginId; id < tabrow - 1; id++) {
                if (initZi_flag)
                    initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
                if (initT_flag)
                    myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);
                if (PBVD_flag)
                    Pbb = EQUIL.PBVD.Eval(0, Ztmp[id], 1);

                gammaGtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pgtmp[id], Pbb, myTemp, &tmpInitZi[0], GAS);
                Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
            }
        }

        // find the water pressure in Dref by Poref
        Pwref = 0;
        Ptmp = Poref;
        mydz = (DOWC - Dref) / mynum;
        myz = Dref;

        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);
            if (PBVD_flag)
                Pbb = EQUIL.PBVD.Eval(0, myz, 1);

            gammaOtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], OIL);
            Ptmp += gammaOtmp * mydz;
            myz += mydz;
        }
        Ptmp -= PcOW;
        for (USI i = 0; i < mynum; i++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(myz, tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, myz, 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Ptmp, Pbb, myTemp, &tmpInitZi[0], WATER);
            Ptmp -= gammaWtmp * mydz;
            myz -= mydz;
        }
        Pwref = Ptmp;

        // find the water pressure in tab
        if (initZi_flag)
            initZi_Tab[0].Eval_All0(Dref, tmpInitZi);
        if (initT_flag)
            myTemp = initT_Tab[0].Eval(0, Dref, 1);

        gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwref, Pbb, myTemp, &tmpInitZi[0], WATER);
        Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
        Pwtmp[beginId] = Pbegin;

        for (USI id = beginId; id > 0; id--) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
        }

        for (USI id = beginId; id < tabrow - 1; id++) {
            if (initZi_flag)
                initZi_Tab[0].Eval_All0(Ztmp[id], tmpInitZi);
            if (initT_flag)
                myTemp = initT_Tab[0].Eval(0, Ztmp[id], 1);

            gammaWtmp = GRAVITY_FACTOR * flashCal[0]->RhoPhase(Pwtmp[id], Pbb, myTemp, &tmpInitZi[0], WATER);
            Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
        }
    }

    DepthP.Display();

    // calculate Pc from DepthP to calculate Sj
    std::vector<OCP_DBL> data(4, 0), cdata(4, 0);
    // if capillary between water and oil is considered
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
            myTemp = initT_Tab[0].Eval(0, depth[n], 1);
            T[n] = myTemp;
        }

        DepthP.Eval_All(0, depth[n], data, cdata);
        OCP_DBL Po = data[1];
        OCP_DBL Pg = data[2];
        OCP_DBL Pw = data[3];
        OCP_DBL Pcgo = Pg - Po;
        OCP_DBL Pcow = Po - Pw;
        OCP_DBL Sw = flow[SATNUM[n]]->GetSwByPcow(Pcow);
        OCP_DBL Sg = 0;
        if (gas) {
            Sg = flow[SATNUM[n]]->GetSgByPcgo(Pcgo);
        }
        if (Sw + Sg > 1) {
            // should me modified
            OCP_DBL Pcgw = Pcow + Pcgo;
            Sw = flow[SATNUM[n]]->GetSwByPcgw(Pcgw);
            Sg = 1 - Sw;
        }

        if (1 - Sw < TINY) {
            // all water
            Po = Pw + flow[SATNUM[n]]->GetPcowBySw(1.0);
        }
        else if (1 - Sg < TINY) {
            // all gas
            Po = Pg - flow[SATNUM[n]]->GetPcgoBySg(1.0);
        }
        else if (1 - Sw - Sg < TINY) {
            // water and gas
            Po = Pg - flow[SATNUM[n]]->GetPcgoBySg(Sg);
        }
        P[n] = Po;

        if (depth[n] < DOGC) {
            Pbb = Po;
        }
        else if (PBVD_flag) {
            Pbb = EQUIL.PBVD.Eval(0, depth[n], 1);
        }
        Pb[n] = Pbb;

        // cal Sw
        OCP_DBL swco = flow[SATNUM[n]]->GetSwco();
        if (!FlagPcow[SATNUM[n]]) {
            S[n * numPhase + numPhase - 1] = swco;
            continue;
        }

        Sw = 0;
        Sg = 0;
        USI     ncut = 10;
        OCP_DBL avePcow = 0;

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
            avePcow += Pcow;
            tmpSw = flow[SATNUM[n]]->GetSwByPcow(Pcow);
            if (gas) {
                tmpSg = flow[SATNUM[n]]->GetSgByPcgo(Pcgo);
            }
            if (tmpSw + tmpSg > 1) {
                // should be modified
                OCP_DBL Pcgw = Pcow + Pcgo;
                tmpSw = flow[SATNUM[n]]->GetSwByPcgw(Pcgw);
                tmpSg = 1 - tmpSw;
            }
            Sw += tmpSw;
            // Sg += tmpSg;
        }
        Sw /= ncut;
        // Sg /= ncut;
        avePcow /= ncut;

        if (SwatInitExist) {
            if (SwatInit[n] <= swco) {
                Sw = swco;
            }
            else {
                Sw = SwatInit[n];
                if (ifScalePcow) {
                    if (avePcow > 0) {
                        OCP_DBL tmp = flow[SATNUM[n]]->GetPcowBySw(Sw);
                        if (tmp > 0) {
                            ScaleValuePcow[n] = avePcow / tmp;
                        }
                    }
                }
            }
        }
        S[n * numPhase + numPhase - 1] = Sw;
    }
}


/////////////////////////////////////////////////////////////////////
// Optional Features
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateScalePcow()
{
    if (ifScalePcow) {
        ScaleValuePcow.resize(numBulk, 1.0);
    }
}

void Bulk::ScalePcow()
{
    if (ifScalePcow) {
        // correct
        const USI Wid = phase2Index[WATER];
        for (USI n = 0; n < numBulk; n++) {
            Pc[n * numPhase + Wid] *= ScaleValuePcow[n];
            Pj[n * numPhase + Wid] = P[n] + Pc[n * numPhase + Wid];
        }
    }
}


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

        SATNUM[bIda] = myGrid.SATNUM[bId];
        PVTNUM[bIda] = myGrid.PVTNUM[bId];
        ROCKNUM[bIda] = myGrid.ROCKNUM[bId];
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

        dx[bIda] = myGrid.dx[bId];
        dy[bIda] = myGrid.dy[bId];
        dz[bIda] = myGrid.dz[bId];
        v[bIda] = myGrid.v[bId];
        depth[bIda] = myGrid.depth[bId];

        ntg[bIda] = myGrid.ntg[bId];
        poroInit[bIda] = myGrid.poro[bId];
        rockKx[bIda] = myGrid.kx[bId];
        rockKy[bIda] = myGrid.ky[bId];
        rockKz[bIda] = myGrid.kz[bId];
    }
}


void Bulk::InitRock()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        rockVntg[n] = v[n] * ntg[n];
        poro[n] = poroInit[n];
        rockVp[n] = rockVntg[n] * poro[n];
    }
}


void Bulk::CalRock()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        rock[ROCKNUM[n]]->CalPoro(P[n], poroInit[n], poro[n], poroP[n]);
        rockVp[n] = rockVntg[n] * poro[n];
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


/////////////////////////////////////////////////////////////////////
// Important Indicator Variable and Check
/////////////////////////////////////////////////////////////////////

OCP_DBL Bulk::CalNRdSmax(OCP_USI& index)
{
    NRdSmax = 0;
    OCP_USI len = numBulk * numPhase;
    for (USI n = 0; n < len; n++) {
        if (fabs(NRdSmax) < fabs(dSNR[n])) {
            NRdSmax = dSNR[n];
            index = n;
        }
    }
    index /= numPhase;
    return NRdSmax;
}

/// Return OCP_TRUE if no negative pressure and OCP_FALSE otherwise.
OCP_BOOL Bulk::CheckP() const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (P[n] < 0) {
            std::ostringstream PStringSci;
            PStringSci << std::scientific << P[n];
            OCP_WARNING("Negative pressure: P[" + std::to_string(n) +
                "] = " + PStringSci.str());
            cout << "P = " << P[n] << endl;
            return OCP_FALSE;
        }
    }

    return OCP_TRUE;
}


/// Return OCP_TRUE if no negative Ni and OCP_FALSE otherwise.
OCP_BOOL Bulk::CheckNi()
{
    OCP_FUNCNAME;

    OCP_USI len = numBulk * numCom;
    for (OCP_USI n = 0; n < len; n++) {
        if (Ni[n] < 0) {
            OCP_USI bId = n / numCom;
            if (Ni[n] > -1E-3 * Nt[bId] && OCP_FALSE) {
                Ni[n] = 1E-8 * Nt[bId];
                // Ni[n] = 0;
            }
            else {
                USI                cId = n - bId * numCom;
                std::ostringstream NiStringSci;
                NiStringSci << std::scientific << Ni[n];
                OCP_WARNING("Negative Ni: Ni[" + std::to_string(cId) + "] in Bulk[" +
                    std::to_string(bId) + "] = " + NiStringSci.str() + ",  " +
                    "dNi = " + std::to_string(dNNR[n]));
#ifdef DEBUG
                for (USI i = 0; i < numCom; i++) {
                    cout << Ni[bId * numCom + i] << "   ";
                }
                cout << endl;
                for (USI j = 0; j < numPhase - 1; j++) {
                    cout << setprecision(9);
                    if (phaseExist[bId * numPhase + j]) {
                        cout << j << "   ";
                        for (USI i = 0; i < numComH; i++) {
                            cout << xij[bId * numPhase * numCom + j * numCom + i]
                                << "   ";
                        }
                        cout << endl;
                    }
                }
#endif // DEBUG
                return OCP_FALSE;
            }
        }
    }
    return OCP_TRUE;
}


/// Return OCP_TRUE if all Ve < Vlim and OCP_FALSE otherwise.
OCP_BOOL Bulk::CheckVe(const OCP_DBL& Vlim) const
{
    OCP_FUNCNAME;

    OCP_DBL dVe = 0.0;
    for (OCP_USI n = 0; n < numBulk; n++) {
        dVe = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (dVe > Vlim) {
            cout << "Volume error at Bulk[" << n << "] = " << setprecision(6) << dVe
                << " is too big!" << endl;
            return OCP_FALSE;

        }
    }
    return OCP_TRUE;
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

    for (OCP_USI n = 0; n < numBulk; n++) {

        // dP
        tmp = fabs(P[n] - lP[n]);
        if (dPmax < tmp) {
            dPmax = tmp;
        }

        // dS
        for (USI j = 0; j < numPhase; j++) {
            id = n * numPhase + j;
            tmp = fabs(S[id] - lS[id]);
            if (dSmax < tmp) {
                dSmax = tmp;
            }
        }

        // dN
        for (USI i = 0; i < numCom; i++) {
            id = n * numCom + i;
            tmp = fabs(max(Ni[id], lNi[id]));
            if (tmp > TINY) {
                tmp = fabs(Ni[id] - lNi[id]) / tmp;
                if (dNmax < tmp) {
                    dNmax = tmp;
                }
            }
        }

        tmp = fabs(vf[n] - rockVp[n]) / rockVp[n];
        if (dVmax < tmp) {
            dVmax = tmp;
        }
    }
}


OCP_DBL Bulk::CalCFL() const
{
    OCP_FUNCNAME;

    maxCFL = 0;
    const OCP_USI len = numBulk * numPhase;
    for (OCP_USI n = 0; n < len; n++) {
        if (phaseExist[n]) {
            cfl[n] /= vj[n];
#ifdef DEBUG
            if (!isfinite(cfl[n])) {
                OCP_ABORT("cfl is nan!");
            }
#endif // DEBUG
            if (maxCFL < cfl[n]) maxCFL = cfl[n];
        }
    }

    return maxCFL;
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
// IMPEC
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateIMPEC_IsoT()
{
    OCP_FUNCNAME;

    // Rock
    rockVntg.resize(numBulk);
    poro.resize(numBulk);
    rockVp.resize(numBulk);

    lporo.resize(numBulk);
    lrockVp.resize(numBulk);

    // derivatives
    poroP.resize(numBulk);
    lporoP.resize(numBulk);


    // Fluid
    phaseNum.resize(numBulk);
    Nt.resize(numBulk);
    Ni.resize(numBulk * numCom);
    vf.resize(numBulk);
    T.resize(numBulk);
    P.resize(numBulk);
    Pb.resize(numBulk);
    Pj.resize(numBulk * numPhase);
    Pc.resize(numBulk * numPhase);
    phaseExist.resize(numBulk * numPhase);
    S.resize(numBulk * numPhase);
    vj.resize(numBulk * numPhase);
    xij.resize(numBulk * numPhase * numCom);
    rho.resize(numBulk * numPhase);
    xi.resize(numBulk * numPhase);
    mu.resize(numBulk * numPhase);
    kr.resize(numBulk * numPhase);

    lphaseNum.resize(numBulk);
    lNt.resize(numBulk);
    lNi.resize(numBulk * numCom);
    lvf.resize(numBulk);
    lP.resize(numBulk);
    lPj.resize(numBulk * numPhase);
    lPc.resize(numBulk * numPhase);
    lphaseExist.resize(numBulk * numPhase);
    lS.resize(numBulk * numPhase);
    vj.resize(numBulk * numPhase);
    lxij.resize(numBulk * numPhase * numCom);
    lrho.resize(numBulk * numPhase);
    lxi.resize(numBulk * numPhase);
    lmu.resize(numBulk * numPhase);
    lkr.resize(numBulk * numPhase);

    // derivatives  
    vfP.resize(numBulk);
    vfi.resize(numBulk * numCom);

    lvfP.resize(numBulk);
    lvfi.resize(numBulk * numCom);

    // others
    cfl.resize(numBulk * numPhase);
}


void Bulk::InitFlashIMPEC()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->InitFlashIMPEC(P[n], Pb[n], T[n], &S[n * numPhase], rockVp[n],
            Ni.data() + n * numCom, n);
        for (USI i = 0; i < numCom; i++) {
            Ni[n * numCom + i] = flashCal[PVTNUM[n]]->Ni[i];
        }
        PassFlashValueIMPEC(n);
    }
}


void Bulk::CalFlashIMPEC()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {

        flashCal[PVTNUM[n]]->FlashIMPEC(P[n], T[n], &Ni[n * numCom], phaseNum[n],
            &xij[n * numPhase * numCom], n);
        PassFlashValueIMPEC(n);
    }
}


void Bulk::PassFlashValueIMPEC(const OCP_USI& n)
{
    OCP_FUNCNAME;

    const OCP_USI bIdp = n * numPhase;
    const USI     pvtnum = PVTNUM[n];
    phaseNum[n] = 0;
    for (USI j = 0; j < numPhase; j++) {
        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        if (phaseExist[bIdp + j]) {
            phaseNum[n]++;
            for (USI i = 0; i < numCom; i++) {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
            }
            vj[bIdp + j] = flashCal[pvtnum]->vj[j];
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfP[n] = flashCal[pvtnum]->vfP;
    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }
}


void Bulk::CalKrPcIMPEC()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId], n);
        for (USI j = 0; j < numPhase; j++)
            Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
    }

    ScalePcow();
}


void Bulk::GetSolIMPEC(const vector<OCP_DBL>& u)
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        P[n] = u[n];
        for (USI j = 0; j < numPhase; j++) {
            OCP_USI id = n * numPhase + j;
            Pj[id] = P[n] + Pc[id];
        }
    }
}


void Bulk::ResetVal01IMPEC()
{
    Pj = lPj;
}


void Bulk::ResetVal02IMPEC()
{
    Pj = lPj;
    Ni = lNi;
}


void Bulk::ResetVal03IMPEC()
{
    // Rock
    rockVp      = lrockVp;
    poro        = lporo;
    poroP       = lporoP;

    // Fluid
    phaseNum    = lphaseNum;
    Nt          = lNt;
    Ni          = lNi;
    vf          = lvf;
    Pj          = lPj;
    phaseExist  = lphaseExist;
    S           = lS;
    vj          = lvj;
    xij         = lxij;
    rho         = lrho;
    xi          = lxi;   
    mu          = lmu; 

    // derivatives
    vfP         = lvfP;
    vfi         = lvfi;

}


void Bulk::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;

    // Rock
    lporo       = poro;
    lporoP      = poroP;
    lrockVp     = rockVp;

    // Fluid
    lphaseNum   = phaseNum;
    lNt         = Nt;
    lNi         = Ni;
    lvf         = vf;
    lP          = P;
    lPj         = Pj;
    lPc         = Pc;
    lphaseExist = phaseExist;
    lS          = S;
    lvj         = vj;
    lxij        = xij;
    lrho        = rho;
    lxi         = xi;   
    lmu         = mu;
    lkr         = kr;
    
    // derivatives
    lvfP        = vfP;
    lvfi        = vfi;

}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateFIM_IsoT()
{
    OCP_FUNCNAME;

    // Rock
    rockVntg.resize(numBulk);
    poro.resize(numBulk);
    rockVp.resize(numBulk);

    lporo.resize(numBulk);
    lrockVp.resize(numBulk);
    
    // derivatives
    poroP.resize(numBulk);
    lporoP.resize(numBulk);
    
    // Fluid
    phaseNum.resize(numBulk);
    Nt.resize(numBulk);
    Ni.resize(numBulk * numCom);
    vf.resize(numBulk);
    T.resize(numBulk);
    P.resize(numBulk);
    Pb.resize(numBulk);
    Pj.resize(numBulk * numPhase);
    Pc.resize(numBulk * numPhase);
    phaseExist.resize(numBulk * numPhase);
    S.resize(numBulk * numPhase);
    xij.resize(numBulk * numPhase * numCom);
    rho.resize(numBulk * numPhase);
    xi.resize(numBulk * numPhase);
    mu.resize(numBulk * numPhase);
    kr.resize(numBulk * numPhase);

    lphaseNum.resize(numBulk);
    lNt.resize(numBulk);
    lNi.resize(numBulk * numCom);
    lvf.resize(numBulk);
    lP.resize(numBulk);
    lPj.resize(numBulk * numPhase);
    lPc.resize(numBulk * numPhase);
    lphaseExist.resize(numBulk * numPhase);
    lS.resize(numBulk * numPhase);
    lxij.resize(numBulk * numPhase * numCom);
    lrho.resize(numBulk * numPhase);
    lxi.resize(numBulk * numPhase);
    lmu.resize(numBulk * numPhase);
    lkr.resize(numBulk * numPhase);

    // derivatives  
    vfP.resize(numBulk);
    vfi.resize(numBulk * numCom);
    rhoP.resize(numBulk * numPhase);
    rhox.resize(numBulk * numCom * numPhase);
    xiP.resize(numBulk * numPhase);
    xix.resize(numBulk * numCom * numPhase);   
    muP.resize(numBulk * numPhase);
    mux.resize(numBulk * numCom * numPhase);
    dPcj_dS.resize(numBulk * numPhase * numPhase);
    dKr_dS.resize(numBulk * numPhase * numPhase);
  
    lvfP.resize(numBulk);
    lvfi.resize(numBulk * numCom);
    lrhoP.resize(numBulk * numPhase);
    lrhox.resize(numBulk * numCom * numPhase);
    lxiP.resize(numBulk * numPhase);
    lxix.resize(numBulk * numCom * numPhase);
    lmuP.resize(numBulk * numPhase);
    lmux.resize(numBulk * numCom * numPhase);
    ldPcj_dS.resize(numBulk * numPhase * numPhase);
    ldKr_dS.resize(numBulk * numPhase * numPhase);  

    // FIM-Specified
    maxLendSdP = (numCom + 1) * (numCom + 1) * numPhase;
    dSec_dPri.resize(numBulk * maxLendSdP);
    bRowSizedSdP.resize(numBulk);
    pSderExist.resize(numBulk * numPhase);
    pVnumCom.resize(numBulk * numPhase);

    ldSec_dPri.resize(numBulk * maxLendSdP);
    lbRowSizedSdP.resize(numBulk);
    lpSderExist.resize(numBulk * numPhase);
    lpVnumCom.resize(numBulk * numPhase);
 
    // NR
    NRstep.resize(numBulk);
    NRphaseNum.resize(numBulk);
    dSNR.resize(numBulk * numPhase);
    dSNRP.resize(numBulk * numPhase);
    dNNR.resize(numBulk * numCom);
    dPNR.resize(numBulk);
}


void Bulk::InitFlashFIM()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        flashCal[PVTNUM[n]]->InitFlashFIM(P[n], Pb[n], T[n], &S[n * numPhase],
            rockVp[n], Ni.data() + n * numCom, n);
        for (USI i = 0; i < numCom; i++) {
            Ni[n * numCom + i] = flashCal[PVTNUM[n]]->Ni[i];
        }
        PassFlashValueFIM(n);
    }
}


/// Use moles of component and pressure both in blackoil and compositional model.
void Bulk::CalFlashFIM()
{
    OCP_FUNCNAME;

    NRdSSP = 0;
    maxNRdSSP = 0;
    index_maxNRdSSP = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {

        flashCal[PVTNUM[n]]->FlashFIM(P[n], T[n], &Ni[n * numCom], &S[n * numPhase], phaseNum[n],
            &xij[n * numPhase * numCom], n);
        PassFlashValueFIM(n);
    }
}


void Bulk::PassFlashValueFIM(const OCP_USI& n)
{
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI     pvtnum = PVTNUM[n];
    USI     len = 0;
    phaseNum[n] = 0;

    for (USI j = 0; j < numPhase; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        dSNR[bIdp + j] = S[bIdp + j] - dSNR[bIdp + j];
        if (phaseExist[bIdp + j]) {
            NRdSSP +=
                (dSNR[bIdp + j] - dSNRP[bIdp + j]) * (dSNR[bIdp + j] - dSNRP[bIdp + j]);
            if (fabs(maxNRdSSP) < fabs(dSNR[bIdp + j] - dSNRP[bIdp + j])) {
                maxNRdSSP = dSNR[bIdp + j] - dSNRP[bIdp + j];
                index_maxNRdSSP = n;
            }
            // cout << n << "   " << scientific << setprecision(6) << dSNR[bIdp + j] <<
            // "   " << dSNRP[bIdp + j] << endl;
        }

        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        pSderExist[bIdp + j] = flashCal[pvtnum]->pSderExist[j];
        pVnumCom[bIdp + j] = flashCal[pvtnum]->pVnumCom[j];
        if (pSderExist[bIdp + j]) len++;
        len += pVnumCom[bIdp + j];
        if (phaseExist[bIdp + j]) { // j -> bId + j fix bugs.
            phaseNum[n]++;
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];

            // Derivatives
            rhoP[bIdp + j] = flashCal[pvtnum]->rhoP[j];
            xiP[bIdp + j] = flashCal[pvtnum]->xiP[j];
            muP[bIdp + j] = flashCal[pvtnum]->muP[j];
                       
            for (USI i = 0; i < numCom; i++) {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
                rhox[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhox[j * numCom + i];
                xix[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xix[j * numCom + i];
                mux[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->mux[j * numCom + i];               
            }
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfP[n] = flashCal[pvtnum]->vfP;

    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

#ifdef OCP_OLD_FIM
    Dcopy(maxLendSdP, &dSec_dPri[n * maxLendSdP], &flashCal[pvtnum]->dXsdXp[0]);
#else
    bRowSizedSdP[n] = len;
    len *= (numCom + 1);
    Dcopy(len, &dSec_dPri[n * maxLendSdP], &flashCal[pvtnum]->dXsdXp[0]);
#endif // OCP_OLD_FIM

}


void Bulk::CalKrPcFIM()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        OCP_USI bId = n * numPhase;
        flow[SATNUM[n]]->CalKrPcDeriv(&S[bId], &kr[bId], &Pc[bId],
            &dKr_dS[bId * numPhase],
            &dPcj_dS[bId * numPhase], n);
        for (USI j = 0; j < numPhase; j++) Pj[bId + j] = P[n] + Pc[bId + j];
    }

    ScalePcow();
}


void Bulk::GetSolFIM(const vector<OCP_DBL>& u,
                     const OCP_DBL&         dPmaxlim,
                     const OCP_DBL&         dSmaxlim)
{
    OCP_FUNCNAME;

    OCP_DBL         dP;
    USI             row = numPhase * (numCom + 1);
    USI             col = numCom + 1;
    OCP_USI         n_np_j;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    dSNR       = S;
    NRphaseNum = phaseNum;
    NRdSmaxP   = 0;
    NRdPmax    = 0;
    NRdNmax    = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

        chopmin = 1;
        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);

#ifdef OCP_OLD_FIM
        DaAxpby(numPhase, col, 1, &dSec_dPri[n * maxLendSdP], u.data() + n * col, 1,
                dtmp.data());
        const OCP_BOOL newFIM = OCP_FALSE;
#else
        DaAxpby(bRowSizedSdP[n], col, 1, &dSec_dPri[n * maxLendSdP], u.data() + n * col,
                1, dtmp.data());
        const OCP_BOOL newFIM = OCP_TRUE;
#endif // OCP_OLD_FIM

        USI js = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (!pSderExist[n * numPhase + j] && newFIM) {
                continue;
            }
            n_np_j = n * numPhase + j;

            choptmp = 1;
            if (fabs(dtmp[js]) > dSmaxlim) {
                choptmp = dSmaxlim / fabs(dtmp[js]);
            } else if (S[n_np_j] + dtmp[js] < 0.0) {
                choptmp = 0.9 * S[n_np_j] / fabs(dtmp[js]);
            }

            // if (fabs(S[n_np_j] - scm[j]) > TINY &&
            //     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
            //     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));

            chopmin = min(chopmin, choptmp);
            js++;
        }

        // dS
        js = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (!pSderExist[n * numPhase + j] && newFIM) {
                dSNRP[n * numPhase + j] = 0;
                continue;
            }
            dSNRP[n * numPhase + j] = chopmin * dtmp[js];
            S[n * numPhase + j] += dSNRP[n * numPhase + j];
            if (fabs(NRdSmaxP) < fabs(dSNRP[n * numPhase + j]))
                NRdSmaxP = dSNRP[n * numPhase + j];
            js++;
        }

        // dxij   ---- Compositional model only
        if (IfUseEoS()) {
            if (phaseNum[n] >= 3) {
                // num of Hydroncarbon phase >= 2
                OCP_USI  bId     = 0;
                for (USI j = 0; j < 2; j++) {
                    bId = n * numPhase * numCom + j * numCom;
                    for (USI i = 0; i < numComH; i++) {
                        xij[bId + i] += chopmin * dtmp[js];
                        js++;
                    }
                }
            }
        }

        // dP
        dP = u[n * col];
        // choptmp = dPmaxlim / fabs(dP);
        // chopmin = min(chopmin, choptmp);
        if (fabs(NRdPmax) < fabs(dP)) NRdPmax = dP;
        P[n] += dP; // seems better
        dPNR[n] = dP;

        // dNi
        NRstep[n] = chopmin;
        for (USI i = 0; i < numCom; i++) {
            dNNR[n * numCom + i] = u[n * col + 1 + i] * chopmin;
            if (fabs(NRdNmax) < fabs(dNNR[n * numCom + i]) / Nt[n])
                NRdNmax = dNNR[n * numCom + i] / Nt[n];

            Ni[n * numCom + i] += dNNR[n * numCom + i];
        }
        // cout << scientific << setprecision(6) << dP << "   " << n << endl;
    }
}


void Bulk::CalRelResFIM(OCPRes& resFIM) const
{
    OCP_FUNCNAME;

    OCP_DBL tmp;

    const USI len = numCom + 1;
    for (OCP_USI n = 0; n < numBulk; n++) {

        for (USI i = 0; i < len; i++) {
            tmp = fabs(resFIM.res[n * len + i] / rockVp[n]);
            if (resFIM.maxRelRes_v < tmp) {
                resFIM.maxRelRes_v = tmp;
                resFIM.maxId_v     = n;
            }
            resFIM.resRelV[n] += tmp * tmp;
        }
        resFIM.resRelV[n] = sqrt(resFIM.resRelV[n]);

        for (USI i = 1; i < len; i++) {
            tmp = fabs(resFIM.res[n * len + i] / Nt[n]);
            if (resFIM.maxRelRes_mol < tmp) {
                resFIM.maxRelRes_mol = tmp;
                resFIM.maxId_mol     = n;
            }
            resFIM.resRelN[n] += tmp * tmp;
        }
        resFIM.resRelN[n] = sqrt(resFIM.resRelN[n]);
    }
}


void Bulk::ResetFIM()
{
    OCP_FUNCNAME;

    // Rock
    poro         = lporo;
    poroP        = lporoP;
    rockVp       = lrockVp;

    // Fluid
    phaseNum     = lphaseNum;
    Nt           = lNt;
    Ni           = lNi;
    vf           = lvf;
    P            = lP;
    Pj           = lPj;
    Pc           = lPc;
    phaseExist   = lphaseExist;
    S            = lS;
    xij          = lxij;
    rho          = lrho;
    xi           = lxi;    
    mu           = lmu;
    kr           = lkr;
    
    // derivatives
    vfP          = lvfP;
    vfi          = lvfi;
    rhoP         = lrhoP;
    rhox         = lrhox;
    xiP          = lxiP;
    xix          = lxix;
    muP          = lmuP;  
    mux          = lmux;
    dPcj_dS      = ldPcj_dS;
    dKr_dS       = ldKr_dS;

    // FIM-Specified
    bRowSizedSdP = lbRowSizedSdP;
    dSec_dPri    = ldSec_dPri;   
    pSderExist   = lpSderExist;
    pVnumCom     = lpVnumCom;

}


void Bulk::UpdateLastStepFIM()
{
    OCP_FUNCNAME;

    // Rock
    lporo         = poro;
    lporoP        = poroP;
    lrockVp       = rockVp;

    // Fluid
    lphaseNum     = phaseNum;
    lNt           = Nt;
    lNi           = Ni;
    lvf           = vf;
    lP            = P;
    lPj           = Pj;
    lPc           = Pc;
    lphaseExist   = phaseExist;
    lS            = S;
    lxij          = xij;
    lrho          = rho;
    lxi           = xi;    
    lmu           = mu;
    lkr           = kr;
    
    // derivatives
    lvfP          = vfP;
    lvfi          = vfi;
    lrhoP         = rhoP;
    lrhox         = rhox;
    lxiP          = xiP;
    lxix          = xix;
    lmuP          = muP;   
    lmux          = mux;
    ldPcj_dS      = dPcj_dS;
    ldKr_dS       = dKr_dS;

    // FIM-Specified
    lbRowSizedSdP = bRowSizedSdP;
    ldSec_dPri    = dSec_dPri;   
    lpSderExist   = pSderExist;
    lpVnumCom     = pVnumCom;
    
}


/////////////////////////////////////////////////////////////////////
// FIMn
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateFIMn_IsoT()
{
    AllocateFIM_IsoT();

    nj.resize(numBulk * numPhase);
    res_n.resize(numBulk * (numPhase + numPhase * numCom));
    resPc.resize(numBulk);
    resIndex.resize(numBulk + 1, 0);

    lnj.resize(numBulk * numPhase);
    lres_n.resize(numBulk * (numPhase + numPhase * numCom));
    lresPc.resize(numBulk);
    lresIndex.resize(numBulk + 1, 0);
}


void Bulk::InitFlashFIMn()
{
    OCP_FUNCNAME;

    if (ifComps) {
        for (OCP_USI n = 0; n < numBulk; n++) {
            flashCal[PVTNUM[n]]->InitFlashFIMn(P[n], Pb[n], T[n], &S[n * numPhase],
                rockVp[n], Ni.data() + n * numCom, n);
            for (USI i = 0; i < numCom; i++) {
                Ni[n * numCom + i] = flashCal[PVTNUM[n]]->Ni[i];
            }
            PassFlashValueFIMn(n);
        }
    }
    else {
        OCP_ABORT("Not Completed in BLKOIL MODEL!");
    }
}


void Bulk::CalFlashFIMn()
{
    OCP_FUNCNAME;

    if (ifComps) {
        CalFlashFIMn_COMP();
    }
    else {
        CalFlashFIMn_BLKOIL();
    }
}



void Bulk::CalFlashFIMn_BLKOIL() {}



void Bulk::CalFlashFIMn_COMP()
{
    vector<USI> flagB(numPhase, 0);
    NRdSSP = 0;
    maxNRdSSP = 0;
    index_maxNRdSSP = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {

        for (USI j = 0; j < numPhase; j++) flagB[j] = phaseExist[n * numPhase + j];

        flashCal[PVTNUM[n]]->FlashFIMn(
            P[n], T[n], &Ni[n * numCom], &S[n * numPhase], &xij[n * numPhase * numCom],
            &nj[n * numPhase], &flagB[0], phaseNum[n], n);

        PassFlashValueFIMn(n);
    }
}



void Bulk::PassFlashValueFIMn(const OCP_USI& n)
{
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI     pvtnum = PVTNUM[n];
    USI     len = 0;
    phaseNum[n] = 0;

    for (USI j = 0; j < numPhase; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        dSNR[bIdp + j] = S[bIdp + j] - dSNR[bIdp + j];
        if (phaseExist[bIdp + j]) {
            NRdSSP +=
                (dSNR[bIdp + j] - dSNRP[bIdp + j]) * (dSNR[bIdp + j] - dSNRP[bIdp + j]);
            if (fabs(maxNRdSSP) < fabs(dSNR[bIdp + j] - dSNRP[bIdp + j])) {
                maxNRdSSP = dSNR[bIdp + j] - dSNRP[bIdp + j];
                index_maxNRdSSP = n;
            }
        }

        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        pSderExist[bIdp + j] = flashCal[pvtnum]->pSderExist[j];
        pVnumCom[bIdp + j] = flashCal[pvtnum]->pVnumCom[j];
        if (pSderExist[bIdp + j]) len++;
        len += pVnumCom[bIdp + j];
        if (phaseExist[bIdp + j]) { // j -> bId + j fix bugs.
            phaseNum[n]++;
            nj[bIdp + j] = flashCal[pvtnum]->nj[j];
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];

            // Derivatives
            rhoP[bIdp + j] = flashCal[pvtnum]->rhoP[j];
            xiP[bIdp + j] = flashCal[pvtnum]->xiP[j];
            muP[bIdp + j] = flashCal[pvtnum]->muP[j];           
            for (USI i = 0; i < numCom; i++) {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
                rhox[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhox[j * numCom + i];
                xix[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xix[j * numCom + i];
                mux[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->mux[j * numCom + i];                           
            }
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfP[n] = flashCal[pvtnum]->vfP;

    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

    resIndex[n + 1] = resIndex[n] + len;
    Dcopy(len, &res_n[0] + resIndex[n], &flashCal[pvtnum]->res[0]);
    len *= (numCom + 1);
    Dcopy(len, &dSec_dPri[n * maxLendSdP], &flashCal[pvtnum]->dXsdXp[0]);

    resPc[n] = flashCal[pvtnum]->resPc;

}


void Bulk::GetSolFIMn(const vector<OCP_DBL>& u,
    const OCP_DBL& dPmaxlim,
    const OCP_DBL& dSmaxlim)
{
    // For saturations changes:
    // 1. maximum changes must be less than dSmaxlim,
    // 2. if phases become mobile/immobile, then set it to crtical point,
    // 3. After saturations are determined, then scale the nij to conserve Volume
    // equations.

    OCP_FUNCNAME;

    dSNR = S;
    NRphaseNum = phaseNum;
    fill(dSNRP.begin(), dSNRP.end(), 0.0);

    NRdSmaxP = 0;
    NRdPmax = 0;
    NRdNmax = 0;

    USI row = numPhase * (numCom + 1);
    USI ncol = numCom + 1;
    USI n_np_j;

    vector<OCP_DBL> dtmp(row, 0);
    vector<OCP_DBL> tmpNij(numPhase * numCom, 0);

    OCP_DBL dSmax;
    OCP_DBL dP;
    OCP_DBL chop = 1;

    // vector<OCP_DBL> tmpxi(numPhase, 0);
    // vector<OCP_DBL> tmpVf(numPhase, 0);

    for (OCP_USI n = 0; n < numBulk; n++) {
        const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

        dP = u[n * ncol];
        dPNR[n] = dP;
        P[n] += dP; // seems better
        if (fabs(NRdPmax) < fabs(dP)) NRdPmax = dP;

        // rockVp[n] = rockVntg[n] * (1 + rockC1 * (P[n] - rockPref));
        // cout << scientific << setprecision(6) << dP << "   " << n << endl;

        dSmax = 0;
        chop = 1;

        const USI cNp = phaseNum[n];
        const USI len = resIndex[n + 1] - resIndex[n];
        Dcopy(len, &dtmp[0], &res_n[resIndex[n]]);
        DaAxpby(len, ncol, 1.0, &dSec_dPri[n * maxLendSdP], u.data() + n * ncol, 1.0,
            dtmp.data());

        // Calculate dSmax
        USI js = 0;
        for (USI j = 0; j < numPhase; j++) {
            n_np_j = n * numPhase + j;
            if (pSderExist[n_np_j]) {
                dSmax = max(fabs(dtmp[js]), dSmax);
                js++;
            }
        }

        // Calculate chop
        if (dSmax > dSmaxlim) {
            chop = min(chop, dSmaxlim / dSmax);
            dSmax = dSmaxlim;
        }
        js = 0;
        for (USI j = 0; j < numPhase; j++) {
            n_np_j = n * numPhase + j;
            if (pSderExist[n_np_j]) {
                if (fabs(S[n_np_j] - scm[j]) > TINY &&
                    (S[n_np_j] - scm[j]) / (chop * dtmp[js]) < 0)
                    chop *= min(1.0, -((S[n_np_j] - scm[j]) / (chop * dtmp[js])));
                if (S[n_np_j] + chop * dtmp[js] < 0)
                    chop *= min(1.0, -S[n_np_j] / (chop * dtmp[js]));
                js++;
            }
        }

        // Calculate S, phaseExist, xij, nj
        fill(tmpNij.begin(), tmpNij.end(), 0.0);
        // fill(Ni.begin() + n * numCom, Ni.begin() + n * numCom + numCom, 0.0);

        js = 0;
        USI jx = cNp;
        for (USI j = 0; j < numPhase; j++) {
            n_np_j = n * numPhase + j;
            if (pSderExist[n_np_j]) {
                dSNRP[n_np_j] = chop * dtmp[js];
                if (fabs(NRdSmaxP) < fabs(dSNRP[n_np_j])) NRdSmaxP = dSNRP[n_np_j];
                js++;

                // S[n_np_j] += chop * dtmp[js];
                // if (S[n_np_j] < TINY) {
                //     pSderExist[n_np_j] = OCP_FALSE;
                // }
                // js++;
                Daxpy(numCom, nj[n_np_j], &xij[n_np_j * numCom], &tmpNij[j * numCom]);
                for (USI i = 0; i < pVnumCom[j]; i++) {
                    tmpNij[j * numCom + i] += chop * dtmp[jx + i];
                }
                jx += pVnumCom[j];
                nj[n_np_j] = Dnorm1(numCom, &tmpNij[j * numCom]);
                for (USI i = 0; i < numCom; i++) {
                    xij[n_np_j * numCom + i] = tmpNij[j * numCom + i] / nj[n_np_j];
                    // Ni[n * numCom + i] += tmpNij[j * numCom + i];
                }
            }
        }

        // Vf correction
        /*        OCP_DBL dVmin = 100 * rockVp[n];
                USI index = 0;
                for (USI j = 0; j < numPhase; j++) {
                    if (phaseExist[n * numPhase + j]) {
                        tmpxi[j] = flashCal[PVTNUM[n]]->XiPhase(P[n], T, &xij[(n *
           numPhase + j) * numCom]); tmpVf[j] = nj[n * numPhase + j] / (S[n * numPhase +
           j] * tmpxi[j]); if (fabs(tmpVf[j] - rockVp[n]) < dVmin) { dVmin =
           fabs(tmpVf[j] - rockVp[n]); index = j;
                        }
                    }
                } */
                // for (USI j = 0; j < numPhase; j++) {
                //     if (phaseExist[n * numPhase + j]) {
                //         nj[n * numPhase + j] *= (tmpVf[index] / tmpVf[j]);
                //         for (USI i = 0; i < numCom; i++) {
                //             Ni[n * numCom + i] += nj[n * numPhase + j] * xij[(n * numPhase +
                //             j) * numCom + i];
                //         }
                //     }
                // }

        NRstep[n] = chop;
        for (USI i = 0; i < numCom; i++) {
            dNNR[n * numCom + i] = u[n * ncol + 1 + i] * chop;
            if (fabs(NRdNmax) < fabs(dNNR[n * numCom + i]) / Nt[n])
                NRdNmax = dNNR[n * numCom + i] / Nt[n];
            Ni[n * numCom + i] += dNNR[n * numCom + i];
        }
    }
}


void Bulk::ResetFIMn()
{
    ResetFIM();
    nj          = lnj;
    res_n       = lres_n;
    resPc       = lresPc;
    resIndex    = lresIndex;
}


void Bulk::UpdateLastStepFIMn()
{
    UpdateLastStepFIM();
    lnj         = nj;
    lres_n      = res_n;
    lresPc      = resPc;
    lresIndex   = resIndex;
}


/////////////////////////////////////////////////////////////////////
// For AIMc
/////////////////////////////////////////////////////////////////////

void Bulk::AllocateAIMc_IsoT()
{
    AllocateFIM_IsoT();

    vj.resize(numBulk * numPhase);
    xijNR.resize(numBulk * numPhase * numCom);
    cfl.resize(numBulk * numPhase);
    bulkTypeAIM.Setup(numBulk);

    lvj.resize(numBulk * numPhase);
}

void Bulk::CalFlashAIMcEp()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            flashCal[PVTNUM[n]]->FlashIMPEC(P[n], T[n], &Ni[n * numCom], phaseNum[n],
                &xijNR[n * numPhase * numCom], n);
            PassFlashValueAIMcEp(n);
        }
    }
}


void Bulk::CalFlashAIMcEa()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            flashCal[PVTNUM[n]]->FlashIMPEC(P[n], T[n], &Ni[n * numCom], phaseNum[n],
                &xij[n * numPhase * numCom], n);
            PassFlashValueAIMcEa(n);
        }       
    }
}


void Bulk::CalFlashAIMcI()
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk

            flashCal[PVTNUM[n]]->FlashFIM(P[n], T[n], &Ni[n * numCom], &S[n * numPhase], phaseNum[n],
                &xij[n * numPhase * numCom], n);
            PassFlashValueAIMcI(n);
        }
    }
}


void Bulk::PassFlashValueAIMcEp(const OCP_USI& n)
{
    // only var about volume needs, some flash var also
    OCP_FUNCNAME;

    USI     pvtnum = PVTNUM[n];

    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfP[n] = flashCal[pvtnum]->vfP;
    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

    OCP_USI bIdp = n * numPhase;
    phaseNum[n] = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (flashCal[pvtnum]->phaseExist[j]) {
            phaseNum[n]++;

            // IMPORTANT -- need for next Flash
            // But xij in nonlinear equations has been modified
            for (USI i = 0; i < numCom; i++) {
                xijNR[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
            }
        }
    }
}


void Bulk::PassFlashValueAIMcI(const OCP_USI& n)
{
    OCP_FUNCNAME;

    OCP_USI bIdp = n * numPhase;
    USI     pvtnum = PVTNUM[n];
    USI     len = 0;
    phaseNum[n] = 0;

    for (USI j = 0; j < numPhase; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        S[bIdp + j] = flashCal[pvtnum]->S[j];
        dSNR[bIdp + j] = S[bIdp + j] - dSNR[bIdp + j];
        if (phaseExist[bIdp + j]) {
            NRdSSP +=
                (dSNR[bIdp + j] - dSNRP[bIdp + j]) * (dSNR[bIdp + j] - dSNRP[bIdp + j]);
            if (fabs(maxNRdSSP) < fabs(dSNR[bIdp + j] - dSNRP[bIdp + j])) {
                maxNRdSSP = dSNR[bIdp + j] - dSNRP[bIdp + j];
                index_maxNRdSSP = n;
            }
            // cout << n << "   " << scientific << setprecision(6) << dSNR[bIdp + j] <<
            // "   " << dSNRP[bIdp + j] << endl;
        }

        phaseExist[bIdp + j] = flashCal[pvtnum]->phaseExist[j];
        pSderExist[bIdp + j] = flashCal[pvtnum]->pSderExist[j];
        pVnumCom[bIdp + j] = flashCal[pvtnum]->pVnumCom[j];
        if (pSderExist[bIdp + j]) len++;
        len += pVnumCom[bIdp + j];
        if (phaseExist[bIdp + j]) { // j -> bId + j fix bugs.
            phaseNum[n]++;
            vj[bIdp + j] = flashCal[pvtnum]->vj[j];
            rho[bIdp + j] = flashCal[pvtnum]->rho[j];
            xi[bIdp + j] = flashCal[pvtnum]->xi[j];
            mu[bIdp + j] = flashCal[pvtnum]->mu[j];

            // Derivatives
            rhoP[bIdp + j] = flashCal[pvtnum]->rhoP[j];
            xiP[bIdp + j] = flashCal[pvtnum]->xiP[j];
            muP[bIdp + j] = flashCal[pvtnum]->muP[j];

            for (USI i = 0; i < numCom; i++) {
                xij[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xij[j * numCom + i];
                rhox[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->rhox[j * numCom + i];
                xix[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->xix[j * numCom + i];
                mux[bIdp * numCom + j * numCom + i] =
                    flashCal[pvtnum]->mux[j * numCom + i];
            }
        }
    }
    Nt[n] = flashCal[pvtnum]->Nt;
    vf[n] = flashCal[pvtnum]->vf;
    vfP[n] = flashCal[pvtnum]->vfP;

    OCP_USI bIdc = n * numCom;
    for (USI i = 0; i < numCom; i++) {
        vfi[bIdc + i] = flashCal[pvtnum]->vfi[i];
    }

#ifdef OCP_OLD_FIM
    Dcopy(maxLendSdP, &dSec_dPri[n * maxLendSdP], &flashCal[pvtnum]->dXsdXp[0]);
#else
    bRowSizedSdP[n] = len;
    len *= (numCom + 1);
    Dcopy(len, &dSec_dPri[n * maxLendSdP], &flashCal[pvtnum]->dXsdXp[0]);
#endif // OCP_OLD_FIM

}


void Bulk::CalKrPcAIMcE()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk
            OCP_USI bId = n * numPhase;
            flow[SATNUM[n]]->CalKrPc(&S[bId], &kr[bId], &Pc[bId], n);
            for (USI j = 0; j < numPhase; j++)
                Pj[n * numPhase + j] = P[n] + Pc[n * numPhase + j];
        }           
    }

    if (ifScalePcow) {
        // correct
        const USI Wid = phase2Index[WATER];
        for (USI n = 0; n < numBulk; n++) {
            if (bulkTypeAIM.IfIMPECbulk(n)) {
                // Explicit bulk
                Pc[n * numPhase + Wid] *= ScaleValuePcow[n];
                Pj[n * numPhase + Wid] = P[n] + Pc[n * numPhase + Wid];
            }           
        }
    }
}


void Bulk::CalKrPcAIMcI()
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk
            OCP_USI bId = n * numPhase;
            flow[SATNUM[n]]->CalKrPcDeriv(
                &S[bId], &kr[bId], &Pc[bId], &dKr_dS[bId * numPhase],
                &dPcj_dS[bId * numPhase], n);
            for (USI j = 0; j < numPhase; j++) Pj[bId + j] = P[n] + Pc[bId + j];
        }
    }

    if (ifScalePcow) {
        // correct
        const USI Wid = phase2Index[WATER];
        for (OCP_USI n = 0; n < numBulk; n++) {
            if (bulkTypeAIM.IfFIMbulk(n)) {
                // Implicit bulk
                Pc[n * numPhase + Wid] *= ScaleValuePcow[n];
                Pj[n * numPhase + Wid] = P[n] + Pc[n * numPhase + Wid];
            }
        }
    }
}

void Bulk::GetSolAIMc(const vector<OCP_DBL>& u,
                      const OCP_DBL&         dPmaxlim,
                      const OCP_DBL&         dSmaxlim)
{
    dSNR     = S;
    NRdSmaxP = 0;
    NRdPmax  = 0;
    OCP_DBL         dP;
    USI             row = numPhase * (numCom + 1);
    USI             col = numCom + 1;
    OCP_USI         n_np_j;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    for (OCP_USI n = 0; n < numBulk; n++) {
        if (bulkTypeAIM.IfIMPECbulk(n)) {
            // IMPEC Bulk
            // Pressure
            dP      = u[n * col];
            NRdPmax = max(NRdPmax, fabs(dP));
            P[n] += dP; // seems better
            dPNR[n]   = dP;
            NRstep[n] = 1;
            // Ni
            for (USI i = 0; i < numCom; i++) {
                dNNR[n * numCom + i] = u[n * col + 1 + i];
                Ni[n * numCom + i] += dNNR[n * numCom + i];
            }
            // Pj
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                Pj[id] = P[n] + Pc[id];
            }
            continue;
        }

        chopmin = 1;
        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);
        DaAxpby(bRowSizedSdP[n], col, 1, &dSec_dPri[n * maxLendSdP], u.data() + n * col,
                1, dtmp.data());

        USI js = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (!pSderExist[n * numPhase + j]) {
                continue;
            }
            n_np_j = n * numPhase + j;

            choptmp = 1;
            if (fabs(dtmp[js]) > dSmaxlim) {
                choptmp = dSmaxlim / fabs(dtmp[js]);
            } else if (S[n_np_j] + dtmp[js] < 0.0) {
                choptmp = 0.9 * S[n_np_j] / fabs(dtmp[js]);
            }

            // if (fabs(S[n_np_j] - scm[j]) > TINY &&
            //     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
            //     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));

            chopmin = min(chopmin, choptmp);
            js++;
        }

        // dS
        js = 0;
        for (USI j = 0; j < numPhase; j++) {
            if (!pSderExist[n * numPhase + j]) {
                dSNRP[n * numPhase + j] = 0;
                continue;
            }
            dSNRP[n * numPhase + j] = chopmin * dtmp[js];
            if (fabs(NRdSmaxP) < fabs(dSNRP[n * numPhase + j]))
                NRdSmaxP = dSNRP[n * numPhase + j];
            js++;
        }

        // dxij   ---- Compositional model only
        if (IfUseEoS()) {
            if (phaseNum[n] >= 3) {
                OCP_USI  bId     = 0;
                for (USI j = 0; j < 2; j++) {
                    bId = n * numPhase * numCom + j * numCom;
                    for (USI i = 0; i < numComH; i++) {
                        xij[bId + i] += chopmin * dtmp[js];
                        js++;
                    }
                }
            }
        }

        // dP
        dP = u[n * col];
        if (fabs(NRdPmax) < fabs(dP)) NRdPmax = dP;
        P[n] += dP; // seems better
        dPNR[n] = dP;

        // dNi
        NRstep[n] = chopmin;
        for (USI i = 0; i < numCom; i++) {
            dNNR[n * numCom + i] = u[n * col + 1 + i] * chopmin;
            if (fabs(NRdNmax) < fabs(dNNR[n * numCom + i]) / Nt[n])
                NRdNmax = dNNR[n * numCom + i] / Nt[n];

            Ni[n * numCom + i] += dNNR[n * numCom + i];
        }
    }
}


void Bulk::ResetAIMc()
{
    ResetFIM();
    vj = lvj;
    xijNR = lxij;
}

void Bulk::UpdateLastStepAIMc()
{
    UpdateLastStepFIM();
    lvj = vj;
    xijNR = xij;
}

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