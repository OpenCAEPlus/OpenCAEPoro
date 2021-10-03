#pragma once
#include "OpenCAEPoro_consts.hpp"
#include "ParamControl.hpp"
#include "Reservoir.hpp"
#include <vector>

using namespace std;

class ControlTime
{
public:
    ControlTime() = default;
    ControlTime(vector<OCP_DBL>& src);
    OCP_DBL TimeInit;       // Maximum init step length of next timestep
    OCP_DBL TimeMax;        // Maximum time step during running
    OCP_DBL TimeMin;        // Minmum time step during running
    OCP_DBL TimeMinChop;    // Minimum choppable timestep
    OCP_DBL TimeMaxIncr;    // Maximum timestep increase factor
    OCP_DBL TimeMinCut;     // Minimum timestep cutback factor
    OCP_DBL TimeCut_F;      // Minimum timestep cutback factor after convergence failure
    OCP_DBL TimeMaxIncre_F; // Maximum increase factor after a convergence failure
};

class ControlError
{
public:
    ControlError() = default;
    ControlError(vector<OCP_DBL>& src);
    OCP_DBL ErrorNL_T; // Target non-linear convergence error
    OCP_DBL ErrorMB_T; // Target material balance error
    OCP_DBL ErrorLS_T; // Target linear convergence error
    OCP_DBL ErrorNL_M; // Maximum non-linear convergence error
    OCP_DBL ErrorMB_M; // Maximum material balance error
    OCP_DBL ErrorLS_M; // Maximum linear convergence error
};

class ControlIter
{
public:
    ControlIter() = default;
    ControlIter(vector<OCP_DBL>& src);
    int    ItMax_NT;  // Maximum number of Newton iterations in a timestep
    int    ItMin_NT;  // Minimum number of Newton iterations in a timestep
    int    ItMax_NTL; // Maximum number of linear iterations in a Newton iteration
    int    ItMin_NTL; // Minimum number of linear iterations in a Newton iteration
    OCP_DBL DPreNT_M;  // Maximum pressure change at last Newton iteration
    OCP_DBL DSatNT_M;  // Maximum saturation change at last Newton iteration
    OCP_DBL DPreNT_T;  // Target pressure change at last Newton iteration
    OCP_DBL Dpre_M;    // Target maximum pressure change in a timestep
};

class OCP_Control
{
    friend class OpenCAEPoro;
    friend class OCP_IMPES;

public:
    void inputParam(ParamControl& CtrlParam);
    void ApplyControl(int i);

    void initTime(int i);

    unsigned int getNumDates() { return CriticalTime.size(); }
    OCP_DBL getCurTime() const { return Current_time; }
    int    getLSiter() const { return LS_iter; }
    int    getNRiter() const { return NR_iter; }

    void setNextTstep(Reservoir& reservoir);

private:
    string         Dir;
    int            Method;
    string         SolveFile;
    vector<OCP_DBL> CriticalTime;
    OCP_DBL         Current_dt;
    OCP_DBL         Current_time{0};
    OCP_DBL         End_time;
    OCP_DBL         TotalTime{ 0 };
    int            Tstep{0};
    int            LS_iter{0};
    int            LS_iter_total{0};
    int            NR_iter{0};
    int            NR_iter_total{0};

    OCP_DBL LS_time{0};

    ControlTime         CtrlTime;
    vector<ControlTime> CtrlTimeSet;

    ControlError         CtrlError;
    vector<ControlError> CtrlErrorSet;

    ControlIter         CtrlIter;
    vector<ControlIter> CtrlIterSet;
};
