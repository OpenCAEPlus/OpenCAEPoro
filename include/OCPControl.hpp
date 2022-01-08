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

/// Params for choosing time stepsize in time marching.
class ControlTime
{
public:
    ControlTime() = default;
    ControlTime(const vector<OCP_DBL>& src);

public:
    // Note: Most commonly used params are the first three
    OCP_DBL timeInit;    ///< Max init step length of next timestep
    OCP_DBL timeMax;     ///< Max time step during running
    OCP_DBL timeMin;     ///< Min time step during running
    OCP_DBL maxIncreFac; ///< Max timestep increase factor
    OCP_DBL minChopFac;  ///< Min choppable timestep
    OCP_DBL cutFacNR;    ///< Factor by which timestep is cut after convergence failure
};

/// Params for convergence and material balance error checks.
class ControlPreTime
{
public:
    ControlPreTime() = default;
    ControlPreTime(const vector<OCP_DBL>& src);

public:
    // Limits for changes at next time step
    OCP_DBL dPlim; ///< Ideal max Pressure change
    OCP_DBL dSlim; ///< Ideal max Saturation change
    OCP_DBL dNlim; ///< Ideal max relative Ni (moles of components) change
    OCP_DBL dVlim; ///< Ideal max relative Verr (pore - fluid) change
};

/// Params for Newton iterations and linear iterations.
class ControlNR
{
public:
    ControlNR() = default;
    ControlNR(const vector<OCP_DBL>& src);

public:
    // Note: Important for convergence of solution methods
    USI     maxNRiter; ///< Maximum number of Newton iterations in a timestep
    OCP_DBL NRtol;     ///< Maximum non-linear convergence error
    OCP_DBL NRdPmax;   ///< Maximum Pressure change in a Newton iteration
    OCP_DBL NRdSmax;   ///< Maximum Saturation change in a Newton iteration
    OCP_DBL NRdPmin;   ///< Minimum Pressure change in a Newton iteration
    OCP_DBL NRdSmin;   ///< Minimum Saturation change in a Newton iteration
    OCP_DBL Verrmax;   ///< Maximum Verr(error between fluidand pore) change in a Newton
                       ///< iteration
};

/// Store shortcut instructions from the command line
class FastControl
{
public:
    void ReadParam(const USI& argc, const char* optset[]);

public:
    bool    activity{false};
    USI     method;        ///< IMPEC or FIM
    OCP_DBL timeInit;      ///< Maximum Init step length of next timestep
    OCP_DBL timeMax;       ///< Maximum time step during running
    OCP_DBL timeMin;       ///< Minmum time step during running
    USI     printLevel{0}; ///< Decide the depth for printfing
};

/// All control parameters except for well controlers.
//  Note: Which solution method will be used is determined here!
class OCPControl
{
    friend class OpenCAEPoro;
    friend class OCPOutput;
    friend class DetailInfo;
    friend class Reservoir;

    friend class OCP_FIM;
    friend class OCP_IMPEC;
    friend class Solver; // temp

public:
    /// Input parameters for control.
    void InputParam(const ParamControl& CtrlParam);

    /// Apply control for time step i.
    void ApplyControl(const USI& i);

    /// Initialize time step i.
    void InitTime(const USI& i);

    /// Setup fast Control.
    void SetupFastControl(const USI& argc, const char* optset[]);

    /// Return type of the solution method.
    USI GetMethod() const { return method; }

    /// Return number of TSTEPs.
    USI GetNumTSteps() const { return criticalTime.size(); }

    /// Return the current time.
    OCP_DBL GetCurTime() const { return current_time; }

    /// Return current dt.
    OCP_DBL& GetCurDt() { return current_dt; }

    /// Return last dt.
    OCP_DBL GetLastDt() const { return last_dt; }

    /// Return the number of linear iterations in one time step.
    USI GetLSiter() const { return iterLS; }

    /// Return the total number of linear iterations.
    USI GetLSiterT() const { return iterLS_total; }

    /// Return the number of Newton iterations in one time step.
    USI GetNRiter() const { return iterNR; }

    /// Return the total nubmer of Newton iterations.
    USI GetNRiterT() const { return iterNR_total; }

    /// Update the number of iterations.
    void UpdateIters();

    /// Update the number of linear iterations.
    void UpdateIterLS(const USI& num) { iterLS += num; }

    /// Update the number of Newton iterations.
    void UpdateIterNR() { iterNR++; }

    /// Reset the number of iterations.
    void ResetIterNRLS();

    /// Update time used for linear solver.
    void UpdateTimeLS(const OCP_DBL& t) { totalLStime += t; }

    /// Record the total time of simulation.
    void RecordTotalTime(const OCP_DBL& t) { totalSimTime = t; }

    /// Calculate the next time step according to max change of some variables.
    void CalNextTstepIMPEC(const Reservoir& reservoir);
    void CalNextTstepFIM(const Reservoir& reservoir);

    /// Determine whether the critical time point has been reached.
    bool IsCriticalTime(const USI& d)
    {
        return ((criticalTime[d] - current_time) < TINY);
    }

    /// Return work dir name.
    string GetWorkDir() const { return workDir; }

    /// Return linear solver file name.
    string GetLsFile() const { return lsFile; }

private:
    USI    method;  ///< Discrete method
    string workDir; ///< Current work directory
    string lsFile;  ///< File name of linear Solver

    vector<OCP_DBL> criticalTime; ///< Set of Critical time by user

    // Record time information
    OCP_DBL current_dt;      ///< Current time step
    OCP_DBL last_dt;         ///< last time step
    OCP_DBL current_time{0}; ///< Current time
    OCP_DBL end_time;        ///< Next Critical time
    OCP_DBL totalSimTime{0}; ///< Total simulation time
    OCP_DBL totalLStime{0};  ///< Total linear solver time

    // Record iteration information
    USI numTstep{0};     ///< Number of time step
    USI iterLS{0};       ///< Current iterations of linear solver
    USI iterLS_total{0}; ///< Total iterations of linear solver
    USI iterNR{0};       ///< Current number of Newton iterations
    USI iterNR_total{0}; ///< Total number of Newton iterations
    USI wastedIterNR{0}; ///< Number of wasted Newton iterations
    USI wastedIterLS{0}; ///< Number of wasted linear iterations

    // Includes time controler, error controler, and iteration controler, all of which
    // could change at different critical time step.
    ControlTime            ctrlTime;
    vector<ControlTime>    ctrlTimeSet;
    ControlPreTime         ctrlPreTime;
    vector<ControlPreTime> ctrlPreTimeSet;
    ControlNR              ctrlNR;
    vector<ControlNR>      ctrlNRSet;
    FastControl            ctrlFast; ///< TODO: What is it for?
};

#endif /* end if __OCP_Control_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/