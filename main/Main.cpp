/*! \file    Main.cpp
 *  \brief   An example to demonstrate main steps of the OCP simulator
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// Standard header files
#include <cstdio>
#include <iostream>
#include <string>

// OpenCAEPoro header files
#include "OCP.hpp"
#include "ParamRead.hpp"

using namespace std;

/// The main() function performs simulation in five steps.
///   -- (1) Read input data;
///   -- (2) Setup static information;
///   -- (3) Initialize reservoir status;
///   -- (4) Run dynamic simulation;
///   -- (5) Output summary and other results.
int main(int argc, const char* argv[])
{
    // Step 0. Print simulator version information.
    OpenCAEPoro simulator;
    if (argc < 2) {
        simulator.PrintUsage(argv[0]);
        return OCP_ERROR_NUM_INPUT; // Need at least one parameter: input data
    } else {
        simulator.PrintVersion();
        if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
            simulator.PrintUsage(argv[0]);
            return OCP_SUCCESS;
        }
    }

    { // Read and process input parameters
        ParamRead rp;

        // Step 1. Read params from an input file to internal params data structure.
        // Remark: The keywords are almost compatible with Ecl100/300; see Keywords.md.
        rp.ReadInputFile(argv[1]);

        // Step 2. Set params using command-line and params read from input file.
        // Remark: It sets up static info, such as active grids and their connections.
        // Remark: Memory allocation for linear systems will also be done at this time.
        simulator.SetupSimulator(rp, argc, argv);
    }

    // Step 3. Initialize the reservoir, which finishes the first step in iterations.
    // Examples: Initial pressure, saturations, moles of components, initial guess of
    // well pressure.
    simulator.InitReservoir();

    // Step 4. Run dynamic simulation using methods like IMPEC and FIM.
    simulator.RunSimulation();

    // Step 5. Output the results according to control params.
    // Remark: It will generate a summary file in your input data directory.
    simulator.OutputResults();

    return OCP_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Jan/07/2022      Update documentation                 */
/*----------------------------------------------------------------------------*/
