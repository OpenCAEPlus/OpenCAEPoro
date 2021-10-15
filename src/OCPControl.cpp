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
    timeInit       = src[0];
    timeMax        = src[1];
    timeMin        = src[2];
    timeMinChop    = src[3];
    timeMaxIncr    = src[4];
    timeMinCut     = src[5];
    timeCut_F      = src[6];
    timeMaxIncre_F = src[7];
}

ControlError::ControlError(const vector<OCP_DBL>& src)
{
    errorNL_T = src[1];
    errorMB_T = src[2];
    errorLS_T = src[3];
    errorNL_M = src[5];
    errorMB_M = src[6];
    errorLS_M = src[7];
}

ControlIter::ControlIter(const vector<OCP_DBL>& src)
{
    iterMax_NT  = src[0];
    iterMin_NT  = src[1];
    iterMax_NTL = src[2];
    iterMin_NTL = src[3];
    dPreNT_M    = src[6];
    dSatNT_M    = src[7];
    dPreNT_T    = src[8];
    dpre_M      = src[9];
    if (dpre_M < 0) dpre_M = dPreNT_M;
}

void OCP_Control::InputParam(const ParamControl& CtrlParam)
{
    workDir = CtrlParam.dir;
    if (CtrlParam.method == "IMPES") {
        method = IMPES;
    } else if (CtrlParam.method == "FIM") {
        method = FIM;
    } else {
        ERRORcheck("Wrong Method !");
        exit(0);
    }
    solveFile    = CtrlParam.linearSolve;
    criticalTime = CtrlParam.criticalTime;

    USI t = CtrlParam.criticalTime.size();
    ctrlTimeSet.resize(t);
    ctrlErrorSet.resize(t);
    ctrlIterSet.resize(t);

    USI         n = CtrlParam.tuning_T.size();
    vector<USI> ctrlCriticalTime(n + 1);
    for (USI i = 0; i < n; i++) {
        ctrlCriticalTime[i] = CtrlParam.tuning_T[i].d;
    }
    ctrlCriticalTime.back() = t;
    for (USI i = 0; i < n; i++) {
        for (USI d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
            ctrlTimeSet[d]  = ControlTime(CtrlParam.tuning_T[i].Tuning[0]);
            ctrlErrorSet[d] = ControlError(CtrlParam.tuning_T[i].Tuning[1]);
            ctrlIterSet[d]  = ControlIter(CtrlParam.tuning_T[i].Tuning[2]);
        }
    }

    cout << "OCP_Control::input" << endl;
}

void OCP_Control::ApplyControl(const USI& i)
{
    ctrlTime  = ctrlTimeSet[i];
    ctrlError = ctrlErrorSet[i];
    ctrlIter  = ctrlIterSet[i];

    end_time = criticalTime[i + 1];
}

void OCP_Control::InitTime(const USI& i)
{
    OCP_DBL dt = criticalTime[i + 1] - current_time;
    if (dt < 0) {
        ERRORcheck("Wrong Time Step");
        exit(0);
    }
    current_dt = min(dt, ctrlTime.timeInit);
}

void OCP_Control::SetNextTstep(const Reservoir& reservoir)
{
    current_time += current_dt;

    OCP_DBL dPmax = reservoir.bulk.GetdPmax();
    OCP_DBL dNmax = reservoir.bulk.GetdNmax();
    OCP_DBL dSmax = reservoir.bulk.GetdSmax();
    OCP_DBL dVmax = reservoir.bulk.GetdVmax();

    OCP_DBL dPlim = 300;
    OCP_DBL dNlim = 0.3;
    OCP_DBL dSlim = 0.3;
    OCP_DBL dVlim = 0.001;

    OCP_DBL c1 = dPmax / dPlim;
    OCP_DBL c2 = dNmax / dNlim;
    OCP_DBL c3 = dSmax / dSlim;
    OCP_DBL c4 = dVmax / dVlim;

    OCP_DBL c = max(max(c1, c2), max(c3, c4));
    c         = max(1.0 / 3, c);
    c         = min(10.0 / 3, c);

    current_dt /= c;

    if (current_dt > ctrlTime.timeMax) current_dt = ctrlTime.timeMax;
    if (current_dt < ctrlTime.timeMin) current_dt = ctrlTime.timeMin;

    OCP_DBL dt = end_time - current_time;
    if (current_dt > dt) current_dt = dt;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/