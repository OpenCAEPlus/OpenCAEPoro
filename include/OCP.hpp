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

#define OCPVersion "0.2.0" ///< Software version tag used for git

/// Top-level data structure for the OpenCAEPoro simulator.
class OpenCAEPoro
{
public:
    /// Output OpenCAEPoro version information.
    void PrintVersion() const
    {
        cout << "=========================================" << endl
             << "OpenCAEPoro Version-" << OCPVersion << endl
             << "=========================================" << endl
             << endl;
    };

    /// Provide at least InputFileName for the input data
    void PrintUsage(string cmdname) const
    {
        cout << "Usage: " << endl
             << "  " << cmdname << " <InputFileName> [Optional Method Parameters]"
             << endl
             << endl;

        cout << "The simplest usage is as follows, where parameters are read from file:"
             << endl
             << "  " << cmdname
             << " examples/spe1a/spe1a.data  %% Solve SPE1a in default setting" << endl
             << endl;

        cout << "You can also pass some command-line options followed by the "
                "input file:"
             << endl
             << "  method:   determine which method will be used " << endl
             << "  dtInit:   determine the initial time step     " << endl
             << "  dtMax:    determine the maximun time step     " << endl
             << "  dtMin:    determine the minimum time step     " << endl
             << "  pl:       determine the print level           " << endl
             << endl;

        cout << "For example: solve SPE1a using FIM" << endl
             << "  " << cmdname
             << " examples/spe1a/spe1a.data method=FIM dtInit=1 dtMax=10 dtMin=0.1 pl=1"
             << endl;

        cout << endl
             << "Attention: " << endl
             << "  - These options override those from the input file;" << endl
             << "  - Only if the method option is specified, other options will work;"
             << endl
             << "  - If (dtInit,dtMax,dtMin) are not set, default values will be used."
             << endl;
    }

    /// Read input parameters to an internal structure.
    void InputParam(ParamRead& param);

    /// Setup reservoir based on an internal structure.
    void SetupSimulator(ParamRead& param, const USI& argc, const char* optset[]);

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