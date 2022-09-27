/*! \file    OCP.hpp
 *  \brief   Main header file for OpenCAEPoro simulator
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_HEADER__
#define __OCP_HEADER__

// OpenCAEPoro header files
#include "OCPControl.hpp"
#include "OCPOutput.hpp"
#include "ParamRead.hpp"
#include "Reservoir.hpp"
#include "Solver.hpp"
#include "UtilTiming.hpp"

#define OCPVersion "0.1.8" ///< Software version tag used for git

/// Top-level data structure for the OpenCAEPoro simulator.
class OpenCAEPoro
{
public:
    /// Output OpenCAEPoro version information.
    void PrintVersion() const
    {
        std::cout << "=========================================" << std::endl
                  << "OpenCAEPoro Version-" << OCPVersion << std::endl
                  << "=========================================" << std::endl
                  << std::endl;
    };

    /// Provide at least InputFileName for the input data
    void PrintUsage(string cmdname) const
    {
        std::cout << "Usage: " << std::endl
                  << "  " << cmdname << " <InputFileName> [Optional Method Parameters]"
                  << std::endl
                  << std::endl;

        std::cout << "The simplest usage is as follow, and everything will be honor of the input file:" << std::endl
                  << "  " << cmdname
                  << " examples/spe1a/spe1a.data  %% Solve SPE1a in default setting"
                  << std::endl
                  << std::endl;

        std::cout << "Also, you could add some options followed by the input file, but attention, these "
                  << std::endl 
                  << "options will override the corresponding contents in the input file and take effect "
                  << std::endl
                  << "all the time, available options are as follows: " 
                  << std::endl
                  << "method:   determine which method will be used."
                  << std::endl
                  << "dtInit:   determine the initial time step."
                  << std::endl
                  << "dtMax:    determine the maximun time step."
                  << std::endl
                  << "dtMin:    determine the minimum time step."
                  << std::endl
                  << "pl:       determine the print level."
                  << std::endl << std::endl;
        std::cout << "for examples" << std::endl;
        std::cout << "  " << cmdname
                  << " examples/spe1a/spe1a.data method=FIM dtInit=1 dtMax=10 dtMin=0.1 pl=1  %% Solve SPE1a using FIM"
                  << std::endl;
    }

    /// Read input parameters to an internal structure.
    void InputParam(ParamRead &param);

    /// Setup reservoir based on an internal structure.
    void SetupSimulator(ParamRead &param, const USI &argc, const char *optset[]);

    /// Initialize or get intitial status of reserovir.
    void InitReservoir();

    /// Run dynamic simulation.
    void RunSimulation();

    /// Output necessary information for post-processing.
    void OutputResults() const;

private:
    /// The core properties of a reservoir.
    Reservoir reservoir;

    /// Contains discrete methods and linear system solver.
    Solver solver;

    /// Control class handles algorithm params and time steping.
    OCPControl control;

    /// Output class handles output level of the program.
    OCPOutput output;
};

#endif /* end if __OCP_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      New tag info                         */
/*  Chensong Zhang      Sep/21/2022      Add PrintUsage                       */
/*----------------------------------------------------------------------------*/