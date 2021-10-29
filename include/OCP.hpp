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
#include "OCPMethod.hpp"
#include "OCPOutput.hpp"
#include "ParamRead.hpp"
#include "Reservoir.hpp"
#include "Solver.hpp"
#include "UtilTiming.hpp"

#define OCPVersion "0.1.1"  ///< Software version number

/// Top level data structure in the OpenCAEPoro simulator
class OpenCAEPoro
{
public:
    /// Output OpenCAEPoro version information.
    void PrintVersion() const
    {
        std::cout << "OpenCAEPoro Version-" << OCPVersion << std::endl
                  << "=========================" << std::endl;
    };

    /// Read input parameters to an internal structure.
    void InputParam(ParamRead& param);

    /// Setup reservoir based on an internal structure.
    void SetupSimulator(ParamRead& param);

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
/*----------------------------------------------------------------------------*/