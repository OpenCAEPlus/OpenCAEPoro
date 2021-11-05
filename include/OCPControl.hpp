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
    OCP_DBL maxIncreFac;  ///< Maximum timestep increase factor
    OCP_DBL minChopFac;   ///< Minimum choppable timestep
    OCP_DBL cutFacNR;     ///< Factor by which timestep is cut after convergence failure
};

/// Params for controlling convergence and material balance error checks.
class ControlPreTime
{
public:
    ControlPreTime() = default;
    ControlPreTime(const vector<OCP_DBL>& src);

public:
    // Used to calculate timestep factor
    OCP_DBL         dPlim; ///< Ideal maximum Pressure change at next time step.
    OCP_DBL         dSlim; ///< Ideal maximum Saturation change at next time step.
    OCP_DBL         dNlim; ///< Ideal maximum relative Ni(moles of components) change at next time step.
    OCP_DBL         dVlim; ///< Ideal maximum relative Verr(error between fluid and pore) change at next time step.
};

/// Params for controlling Newton iterations and linear iterations.
class ControlNR
{
public:
    ControlNR() = default;
    ControlNR(const vector<OCP_DBL>& src);

public:
    // Note: Important for convergence of solution methods
    USI             maxNRiter; ///< Maximum number of Newton iterations in a timestep
    OCP_DBL         NRtol;     ///< Maximum non-linear convergence error
    OCP_DBL         NRdPmax;   ///< Maximum Pressure change in a Newton iteration
    OCP_DBL         NRdSmax;   ///< Maximum Saturation change in a Newton iteration
    OCP_DBL         NRdPmin;   ///< Minimum Pressure change in a Newton iteration
    OCP_DBL         NRdSmin;   ///< Minimum Saturation change in a Newton iteration
    OCP_DBL         Verrmax;   ///< Maximum Verr(error between fluidand pore) change in a Newton iteration
};

/// All control parameters except for well controlers.
//  Note: Which discrete method will be used is determined here!
class OCPControl
{
    friend class OpenCAEPoro;
    friend class OCPOutput;
    friend class DetailInfo;
    friend class Reservoir;

    friend class OCP_FIM;
    friend class OCP_IMPEC;
    friend class Solver;  // temp

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
    OCP_DBL GetLastCurDt() const { return last_dt; }
    /// Return the number of linear solver iterations in one time step.
    USI GetLSiter() const { return iterLS; }
    USI GetLSiterT() const { return iterLS_total; }
    /// Return the number of Newton iterations in one time step.
    USI GetNRiter() const { return iterNR; }
    USI GetNRiterT() const { return iterNR_total; }

    /// Update num of iterations.
    void UpdateIters();
    /// Update num of linear solver steps.
    void UpdateIterLS(const USI& num)
    {
        iterLS = num;
        iterLS_total += num;
    }
    void UpdateIterNR(){ iterNR++; }
    /// Update time used for linear solver.
    void UpdateTimeLS(const OCP_DBL& t) { totalLStime += t; }
    /// Record the total time of simulation.
    void RecordTotalTime(const OCP_DBL& t) { totalSimTime = t; }
    /// Calculate the next time step according to max change of some variables.
    void CalNextTstepIMPEC(const Reservoir& reservoir);
    void CalNextTstepFIM(const Reservoir& reservoir);
    /// Determine whether the critical time point has been reached.
    bool IsCriticalTime(const USI& d) { return ((criticalTime[d] - current_time) < TINY); }

    string GetWorkDir() const { return workDir; }
    string GetLsFile() const { return lsFile; }

private:
    USI             method;   ///< Discrete method
    string          workDir;  ///< Current work directory
    string          lsFile;   ///< File name of linear Solver
    vector<OCP_DBL> criticalTime; ///< Set of Critical time by user

    OCP_DBL         current_dt; ///< Current time step
    OCP_DBL         last_dt; ///< last time step.
    OCP_DBL         current_time{0}; ///< Current time.
    OCP_DBL         end_time; ///< Next Critical time

    // Record
    OCP_DBL         totalSimTime{0}; ///< Total simulation time
    OCP_DBL         totalLStime{ 0 }; ///< Total linear solver time
    USI             numTstep{0};     ///< Num of time step
    USI             iterLS{0};       ///< Current iterations of Linear Solve
    USI             iterLS_total{0}; ///< Total iterations of Linear Solve
    USI             iterNR{0};       ///< Current iterations of NR
    USI             iterNR_total{0}; ///< Total iterations of NR

    // Includes time controler, error controler, and iteration controler, all of which
    // could change at different critical time step.
    ControlTime         ctrlTime;
    vector<ControlTime> ctrlTimeSet;

    ControlPreTime         ctrlPreTime;
    vector<ControlPreTime> ctrlPreTimeSet;

    ControlNR         ctrlNR;
    vector<ControlNR> ctrlNRSet;
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
