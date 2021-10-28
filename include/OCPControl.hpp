/*! \file    OCPControl.hpp
 *  \brief   OCPControl class declaration
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
#include "OCPConst.hpp"
#include "ParamControl.hpp"
#include "Reservoir.hpp"

using namespace std;

/// Params for controlling time stepsized in time marching.
class ControlTime
{
public:
    ControlTime() = default;
    ControlTime(const vector<OCP_DBL>& src);

public:
    // Note: Most commonly used params are the first three.
    OCP_DBL timeInit;    ///< Maximum Init step length of next timestep
    OCP_DBL timeMax;     ///< Maximum time step during running
    OCP_DBL timeMin;     ///< Minmum time step during running
    OCP_DBL timeMinChop; ///< Minimum choppable timestep
    OCP_DBL timeMaxIncr; ///< Maximum timestep increase factor
    OCP_DBL timeMinCut;  ///< Minimum timestep cutback factor
    OCP_DBL timeCut_F;   ///< Minimum timestep cutback factor after convergence failure
    OCP_DBL timeMaxIncre_F; ///< Maximum increase factor after a convergence failure
};

/// Params for controlling convergence and material balance error checks.
class ControlError
{
public:
    ControlError() = default;
    ControlError(const vector<OCP_DBL>& src);

public:
    OCP_DBL errorNL_T; ///< Target non-linear convergence error
    OCP_DBL errorMB_T; ///< Target material balance error
    OCP_DBL errorLS_T; ///< Target linear convergence error
    OCP_DBL errorNL_M; ///< Maximum non-linear convergence error
    OCP_DBL errorMB_M; ///< Maximum material balance error
    OCP_DBL errorLS_M; ///< Maximum linear convergence error
};

/// Params for controlling Newton iterations and linear iterations.
class ControlIter
{
public:
    ControlIter() = default;
    ControlIter(const vector<OCP_DBL>& src);

public:
    // Note: Important for convergence of solution methods
    USI     iterMax_NT;  // Maximum number of Newton iterations in a timestep
    USI     iterMin_NT;  // Minimum number of Newton iterations in a timestep
    USI     iterMax_NTL; // Maximum number of linear iterations in a Newton iteration
    USI     iterMin_NTL; // Minimum number of linear iterations in a Newton iteration
    OCP_DBL dPreNT_M;    // Maximum pressure change at last Newton iteration
    OCP_DBL dSatNT_M;    // Maximum saturation change at last Newton iteration
    OCP_DBL dPreNT_T;    // Target pressure change at last Newton iteration
    OCP_DBL dpre_M;      // Target maximum pressure change in a timestep
};

/// All control parameters except for well controlers.
//  Note: Which discrete method will be used is determined here!
class OCP_Control
{
    friend class OpenCAEPoro;
    friend class OCP_IMPEC;
    friend class OCP_Output;
    friend class DetailInfo;

public:
    /// Input parameters for control.
    void InputParam(const ParamControl& CtrlParam);
    /// Apply control for time step i.
    void ApplyControl(const USI& i);
    /// Initialize time step i.
    void InitTime(const USI& i);

    /// Return the method
    USI GetMethod()const { return method; }

    /// Return number of TSTEPs.
    USI GetNumTSteps() const { return criticalTime.size(); }
    /// Return the current time.
    OCP_DBL GetCurTime() const { return current_time; }
    /// Return current dt.
    OCP_DBL& GetCurDt() { return current_dt; }
    /// Return last dt
    OCP_DBL GetLastCurDt() const { return lcurrent_dt; }
    /// Return the number of linear solver iterations in one time step.
    USI GetLSiter() const { return iterLS; }
    /// Return the number of Newton iterations in one time step.
    USI GetNRiter() const { return iterNR; }

    /// Update num of iterations.
    void UpdateIters();
    /// Update num of linear solver steps.
    void UpdateIterLS(const USI& num)
    {
        iterLS = num;
        iterLS_total += num;
    }
    /// Update time used for linear solver.
    void UpdateTimeLS(const OCP_DBL& t) { timeLS += t; }
    /// Record the total time of simulation.
    void RecordTotalTime(const OCP_DBL& t) { totalTime = t; }
    /// Calculate the next time step according to max change of some variables.
    void CalNextTstep(const Reservoir& reservoir);
    /// Determine whether the critical time point has been reached.
    bool IsCriticalTime(const USI& d)
    {
        return ((criticalTime[d] - current_time) < TINY);
    }

private: // TODO: Add doxygen!
    USI             method;
    string          workDir;
    string          solveFile;
    vector<OCP_DBL> criticalTime;
    OCP_DBL         current_dt;
    OCP_DBL         lcurrent_dt;
    OCP_DBL         current_time{0};
    OCP_DBL         end_time;
    OCP_DBL         totalTime{0};
    USI             tstep{0};
    USI             iterLS{0};
    USI             iterLS_total{0};
    USI             iterNR{0};
    USI             iterNR_total{0};

    OCP_DBL timeLS{0};

    // Includes time controler, error controler, and iteration controler, all of which
    // could change at different critical time step.
    ControlTime         ctrlTime;
    vector<ControlTime> ctrlTimeSet;

    ControlError         ctrlError;
    vector<ControlError> ctrlErrorSet;

    ControlIter         ctrlIter;
    vector<ControlIter> ctrlIterSet;
};

#endif /* end if __OCP_Control_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/
