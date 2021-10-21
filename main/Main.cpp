/*! \file    Main.cpp
 *  \brief   Main for an example to Display main steps in our simulator
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

// The main function divided reservoir simulation into five steps:
//     (1) Read param from an input file
// --> (2) Setup static information with input parameters
// --> (3) Initialize the reservoir
// --> (4) Run dynamic simulation
// --> (5) Output the results
int main(int argc, const char* argv[])
{
    if (argc != 2) {
        cout << "Wrong number of arguments. Usage: " << argv[0] << " <InputFileName>"
             << endl;
        return(-1);
    }
    OpenCAEPoro simulator;
    simulator.PrintVersion();

    // Step 1. Read params from an input file to internal params data structure.
    // Note: The keywords are almost compatible with Ecl simulator.
    ParamRead rp;
    rp.ReadInputFile(argv[1]);

    // Step 2. Read param from internal params data structure to each modules, and
    // Setup static information, such as active grids, and connections between them.
    // Note: Memory allocation for linear systems will also be done at this time.
    simulator.SetupReservoir(rp);

    // Step 3. Initialize the reservoir, which finishs the first step in iterations.
    // For example: initial pressure, saturation, and moles of components will be
    // calculated. Initial guess of well pressure will also be given here.
    simulator.InitReservoir();

    // Step 4. Run dynamic simulation using methods like IMPEC and FIM. It's a
    // combination of functions of various modules which you could make changes
    // whenever necessary.
    simulator.RunSimulation();

    // Step 5. Output the results according to params users give. It will generate
    // a summary file in your input data directory.
    simulator.OutputResults();

    return 0;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/