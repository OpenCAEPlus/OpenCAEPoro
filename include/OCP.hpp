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

#define OCPVersion "0.5.0" ///< Software version tag used for git

/// Top-level data structure for the OpenCAEPoro simulator.
class OpenCAEPoro
{
public:
    /// Output OpenCAEPoro version information.
    void PrintVersion() const
    {
        cout << "=========================================" << endl
             << "OpenCAEPoro Version-" << OCPVersion << endl
             << "=========================================" << endl;
    };

    /// Provide at least InputFileName for the input data
    void PrintUsage(string cmdname) const
    {
        cout << "Usage: " << endl
             << "  " << cmdname << " <InputFileName> [<options>]" << endl
             << endl;

        cout << "A simple example is to solve SPE1 Case A in default setting" << endl
             << "  " << cmdname << " examples/spe1a/spe1a.data" << endl;

        cout << endl
             << "Another example is to solve the same problem using FIM" << endl
             << "  " << cmdname
             << " examples/spe1a/spe1a.data method=FIM dtInit=1 dtMax=10 dtMin=0.1"
             << endl
             << endl;

        cout << "You can pass optional cmd arguments after the input file:" << endl
             << "     method = solution method to use " << endl
             << "     dtInit = initial time stepsize  " << endl
             << "      dtMax = maximum time stepsize  " << endl
             << "      dtMin = minimum time stepsize  " << endl
             << "    verbose = print level on screen  " << endl
             << endl;

        cout << "Attention: " << endl
             << "  - Only if `method' is set, other options will take effect;" << endl
             << "  - These cmd options will override those in the input file;" << endl
             << "  - If (dtInit,dtMax,dtMin) are not set, default values will be used."
             << endl
             << endl;
    }

    /// Read Param from an input file.
    void ReadInputFile(const string& filename);

    /// Read input parameters to an internal structure.
    void InputParam(ParamRead& param);

    /// Setup reservoir based on an internal structure.
    void SetupSimulator(const USI& argc, const char* options[]);

    /// Initialize or get initial status of reservoir.
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

    /// Control class handles algorithm params and time stepping.
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