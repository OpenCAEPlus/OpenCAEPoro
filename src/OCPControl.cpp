/*! \file    OCPControl.cpp
 *  \brief   OCPControl class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPControl.hpp"

ControlTime::ControlTime(const vector<OCP_DBL>& src)
{
    timeInit    = src[0];
    timeMax     = src[1];
    timeMin     = src[2];
    maxIncreFac = src[3];
    minChopFac  = src[4];
    cutFacNR    = src[5];
}

ControlPreTime::ControlPreTime(const vector<OCP_DBL>& src)
{
    dPlim = src[0];
    dSlim = src[1];
    dNlim = src[2];
    dVlim = src[3];
}

ControlNR::ControlNR(const vector<OCP_DBL>& src)
{
    maxNRiter = src[0];
    NRtol     = src[1];
    NRdPmax   = src[2];
    NRdSmax   = src[3];
    NRdPmin   = src[4];
    NRdSmin   = src[5];
    Verrmax   = src[6];
}

void FastControl::ReadParam(const USI& argc, const char* optset[])
{
    activity = OCP_FALSE;
    timeInit = timeMax = timeMin = -1.0;

    std::stringstream buffer;
    string             tmp;
    string            key;
    string             value;
    for (USI n = 2; n < argc; n++) {
        buffer << optset[n];
        buffer >> tmp;

        OCP_INT pos = tmp.find_last_of('=');

        if (pos == string::npos) OCP_ABORT("Unknown Usage! See -h");

        key         = tmp.substr(0, pos);
        value       = tmp.substr(pos + 1, tmp.size() - pos);

        switch (Map_Str2Int(&key[0],key.size())) 
        {

            case Map_Str2Int("method",6):
                if (value == "FIM") {
                    method = FIM;
                } else if (value == "FIMn") {
                    method = FIMn;
                } else if (value == "IMPEC") {
                    method = IMPEC;
                } else if (value == "AIMc") {
                    method = AIMc;
                }else {
                    OCP_ABORT("Wrong method param in command line!");
                }
                activity = OCP_TRUE;
                if (method == FIM || method == FIMn || method == AIMc) {
                    if (timeInit <= 0) timeInit = 0.1;
                    if (timeMax <= 0) timeMax = 10.0;
                    if (timeMin <= 0) timeMin = 0.1;
                } else {
                    if (timeInit <= 0) timeInit = 0.1;
                    if (timeMax <= 0) timeMax = 1.0;
                    if (timeMin <= 0) timeMin = 0.1;
                }               
                break;

            case Map_Str2Int("dtInit", 6):
                timeInit = stod(value);
                break;

            case Map_Str2Int("dtMin", 5):
                timeMin = stod(value);
                break;

            case Map_Str2Int("dtMax", 5):
                timeMax = stod(value);
                break;

            case Map_Str2Int("pl", 2):
                printLevel = stoi(value);
                break;

            default:
                OCP_ABORT("Unknown Options: " + key + "   See -h");
                break;
        }

        buffer.clear();
    }






    if (argc >= 6 && OCP_FALSE) {
        activity = OCP_TRUE;
        if (string(optset[2]) == "FIM") {
            method = FIM;
        } else if (string(optset[2]) == "FIMn") {
            method = FIMn;
        } else if (string(optset[2]) == "IMPEC") {
            method = IMPEC;
        } else if (string(optset[2]) == "AIMc") {
            method = AIMc;
        }else {
            OCP_ABORT("Wrong method param in command line!");
        }
        timeInit = stod(optset[3]);
        timeMax  = stod(optset[4]);
        timeMin  = stod(optset[5]);
        if (argc >= 7) {
            printLevel = stoi(optset[6]);
        }
    }
}

void OCPControl::InputParam(const ParamControl& CtrlParam)
{
    workDir = CtrlParam.dir;
    if (CtrlParam.method == "IMPEC") {
        method = IMPEC;
    } else if (CtrlParam.method == "FIM") {
        method = FIM;
    } else if (CtrlParam.method == "FIMn") {
        method = FIMn;
    } else if (CtrlParam.method == "AIMc") {
        method = AIMc;
    }else {
        OCP_ABORT("Wrong method specified!");
    }

    linearsolveFile = CtrlParam.linearSolve;
    criticalTime = CtrlParam.criticalTime;
    
    USI t = CtrlParam.criticalTime.size();
    ctrlTimeSet.resize(t);
    ctrlPreTimeSet.resize(t);
    ctrlNRSet.resize(t);

    USI         n = CtrlParam.tuning_T.size();
    vector<USI> ctrlCriticalTime(n + 1);
    for (USI i = 0; i < n; i++) {
        ctrlCriticalTime[i] = CtrlParam.tuning_T[i].d;
    }
    ctrlCriticalTime.back() = t;
    for (USI i = 0; i < n; i++) {
        for (USI d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
            ctrlTimeSet[d]    = ControlTime(CtrlParam.tuning_T[i].Tuning[0]);
            ctrlPreTimeSet[d] = ControlPreTime(CtrlParam.tuning_T[i].Tuning[1]);
            ctrlNRSet[d]      = ControlNR(CtrlParam.tuning_T[i].Tuning[2]);
        }
    }
}

void OCPControl::ApplyControl(const USI& i, const Reservoir& rs)
{
    ctrlTime    = ctrlTimeSet[i];
    ctrlPreTime = ctrlPreTimeSet[i];
    ctrlNR      = ctrlNRSet[i];
    end_time    = criticalTime[i + 1];
    wellChange = rs.allWells.GetWellChange();
    InitTime(i);
}

void OCPControl::InitTime(const USI& i)
{
    OCP_DBL dt = criticalTime[i + 1] - current_time;
    if (dt <= 0) OCP_ABORT("Non-positive time stepsize!");

    static OCP_BOOL firstflag = OCP_TRUE;
    if (wellChange || firstflag) {
        current_dt = min(dt, ctrlTime.timeInit);
        firstflag = OCP_FALSE;
    }
    else {
        current_dt = min(dt, init_dt);
    }
    
}

void OCPControl::SetupFastControl(const USI& argc, const char* optset[])
{
    ctrlFast.ReadParam(argc, optset);
    if (ctrlFast.activity) {
        method = ctrlFast.method;
        switch (method) {
            case IMPEC:
                linearsolveFile = "./csr.fasp";
                break;
            case AIMc:
            case FIM:
            case FIMn:
                linearsolveFile = "./bsr.fasp";
                break;
            default:
                OCP_ABORT("Wrong method in command line!");
                break;
        }
        USI n = ctrlTimeSet.size();
        for (USI i = 0; i < n; i++) {
            ctrlTimeSet[i].timeInit = ctrlFast.timeInit;
            ctrlTimeSet[i].timeMax  = ctrlFast.timeMax;
            ctrlTimeSet[i].timeMin  = ctrlFast.timeMin;
        }
        printLevel = ctrlFast.printLevel;
    }
}

void OCPControl::CalNextTstepIMPEC(const Reservoir& reservoir)
{
    last_dt = current_dt;
    current_time += current_dt;

    OCP_DBL c1, c2, c3, c4, c;
    c1 = c2 = c3 = c4 = 10;

    const OCP_DBL dPmax = reservoir.bulk.GetdPmax();
    const OCP_DBL dNmax = reservoir.bulk.GetdNmax();
    const OCP_DBL dSmax = reservoir.bulk.GetdSmax();
    const OCP_DBL dVmax = reservoir.bulk.GetdVmax();

    //cout << dPmax << "   " << dSmax << endl;

    if (dPmax > TINY) c1 = ctrlPreTime.dPlim / dPmax;
    if (dSmax > TINY) c2 = ctrlPreTime.dSlim / dSmax;
    if (dNmax > TINY) c3 = ctrlPreTime.dNlim / dNmax;
    if (dVmax > TINY) c4 = ctrlPreTime.dVlim / dVmax;

    c = min(min(c1, c2), min(c3, c4));
    c = max(ctrlTime.minChopFac, c);
    c = min(ctrlTime.maxIncreFac, c);

    current_dt *= c;

    if (current_dt > ctrlTime.timeMax) current_dt = ctrlTime.timeMax;
    if (current_dt < ctrlTime.timeMin) current_dt = ctrlTime.timeMin;

    init_dt = current_dt;

    OCP_DBL dt = end_time - current_time;
    if (current_dt > dt) current_dt = dt;


    // cout << "FactorT: " << c << "   Next dt: " << dt << endl;
}

void OCPControl::CalNextTstepFIM(const Reservoir& reservoir)
{
    last_dt = current_dt;
    current_time += current_dt;

    OCP_DBL c1, c2, c;
    c1 = c2 = 10;

    const OCP_DBL dPmaxB = reservoir.bulk.GetdPmax();
    const OCP_DBL dPmaxW = reservoir.allWells.GetdBHPmax();
    const OCP_DBL dPmax  = max(dPmaxB, dPmaxW);
    const OCP_DBL dSmax = reservoir.bulk.GetdSmax();

    if (dPmax > TINY) c1 = ctrlPreTime.dPlim / dPmax;
    if (dSmax > TINY) c2 = ctrlPreTime.dSlim / dSmax;

    OCP_DBL c3 = 1.5;

    if (iterNR < 30) { // temp
        c3 = 2;
    } else if (iterNR > 8) {
        c3 = 0.5;
    }

    c = min(min(c1, c2), c3);
    c = max(ctrlTime.minChopFac, c);
    c = min(ctrlTime.maxIncreFac, c);

    current_dt *= c;
    if (current_dt > ctrlTime.timeMax) current_dt = ctrlTime.timeMax;
    if (current_dt < ctrlTime.timeMin) current_dt = ctrlTime.timeMin;

    init_dt = current_dt;

    const OCP_DBL dt = end_time - current_time;
    if (current_dt > dt) current_dt = dt;

}

void OCPControl::UpdateIters()
{
    numTstep += 1;
    iterNR_total += iterNR;
    iterLS_total += iterLS;
    iterNR = 0;
    iterLS = 0;
}

void OCPControl::ResetIterNRLS()
{
    wastedIterNR += iterNR;
    iterNR = 0;
    wastedIterLS += iterLS;
    iterLS = 0;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/