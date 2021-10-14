/*! \file    OCP_Control.hpp
 *  \brief   OCP_Control class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_CONTROL_HEADER__
#define __OCP_CONTROL_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "OpenCAEPoro_consts.hpp"
#include "ParamControl.hpp"
#include "Reservoir.hpp"


using namespace std;

/// ControlTime contains params used to control the size of time step.
/// the most commonly used params are the first three. the others are under development.
class ControlTime
{
public:
    ControlTime() = default;
    ControlTime(vector<OCP_DBL>& src);
    OCP_DBL TimeInit;       ///< Maximum Init step length of next timestep
    OCP_DBL TimeMax;        ///< Maximum time step during running
    OCP_DBL TimeMin;        ///< Minmum time step during running
    OCP_DBL TimeMinChop;    ///< Minimum choppable timestep
    OCP_DBL TimeMaxIncr;    ///< Maximum timestep increase factor
    OCP_DBL TimeMinCut;     ///< Minimum timestep cutback factor
    OCP_DBL TimeCut_F;      ///< Minimum timestep cutback factor after convergence failure
    OCP_DBL TimeMaxIncre_F; ///< Maximum increase factor after a convergence failure
};

/// ControlError contains params used to control the convergence error or material balance error during the simulation.
class ControlError
{
public:
    ControlError() = default;
    ControlError(vector<OCP_DBL>& src);
    OCP_DBL ErrorNL_T; ///< Target non-linear convergence error
    OCP_DBL ErrorMB_T; ///< Target material balance error
    OCP_DBL ErrorLS_T; ///< Target linear convergence error
    OCP_DBL ErrorNL_M; ///< Maximum non-linear convergence error
    OCP_DBL ErrorMB_M; ///< Maximum material balance error
    OCP_DBL ErrorLS_M; ///< Maximum linear convergence error
};

/// ControlIter contains params used to control the number of Newton iterations or linear iterations.
/// whether a iteration is satisfied will also determined here.
class ControlIter
{
public:
    ControlIter() = default;
    ControlIter(const vector<OCP_DBL>& src);
    USI    ItMax_NT;  // Maximum number of Newton iterations in a timestep
    USI    ItMin_NT;  // Minimum number of Newton iterations in a timestep
    USI    ItMax_NTL; // Maximum number of linear iterations in a Newton iteration
    USI    ItMin_NTL; // Minimum number of linear iterations in a Newton iteration
    OCP_DBL DPreNT_M;  // Maximum pressure change at last Newton iteration
    OCP_DBL DSatNT_M;  // Maximum saturation change at last Newton iteration
    OCP_DBL DPreNT_T;  // Target pressure change at last Newton iteration
    OCP_DBL Dpre_M;    // Target maximum pressure change in a timestep
};

/// OCP_Control is responsible for all of the controller except well controler.
/// these controler includes time controler, error controler and iteration controler,
/// all of which could change at different critical time point, it's up to users.
/// which discrete method will be used is determined here.
class OCP_Control
{
    friend class OpenCAEPoro;
    friend class OCP_IMPES;

public:
    void InputParam(ParamControl& CtrlParam);
    void ApplyControl(const USI& i);

    void InitTime(const USI& i);

    USI getNumDates() const { return criticalTime.size(); }
    OCP_DBL getCurTime() const { return Current_time; }
    USI    getLSiter() const { return LS_iter; }
    USI    getNRiter() const { return NR_iter; }

    void setNextTstep(Reservoir& reservoir);

private:
    USI            Method;
    string         Dir;
    string         solveFile;
    vector<OCP_DBL> criticalTime;
    OCP_DBL         Current_dt;
    OCP_DBL         Current_time{0};
    OCP_DBL         End_time;
    OCP_DBL         TotalTime{ 0 };
    USI            Tstep{0};
    USI            LS_iter{0};
    USI            LS_iter_total{0};
    USI            NR_iter{0};
    USI            NR_iter_total{0};

    OCP_DBL LS_time{0};

    ControlTime         ctrlTime;
    vector<ControlTime> ctrlTimeSet;

    ControlError         ctrlError;
    vector<ControlError> ctrlErrorSet;

    ControlIter         ctrlIter;
    vector<ControlIter> ctrlIterSet;
};

#endif  /* end if __OCP_Control_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/
