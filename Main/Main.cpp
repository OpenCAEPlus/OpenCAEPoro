/*! \file    Main.cpp
 *  \brief   Main for an example to display main steps in our simulator
 *  \author  Shizhe Li
 *  \date    Oct/08/2021
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
#include "OpenCAEPoro.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "Timing.hxx"

using namespace std;

// The main function divided reservoir simulation into five steps:
//     (1) Read param from an input file 
// --> (2) Setup static information with input parameters
// --> (3) Initialize the reservoir
// --> (4) Run dynamic simulation
// --> (5) Output the results
int main(int argc, const char* argv[])
{
    if (argc == 1) {
        cout << "Input file is missing. Usage: ./OpenCAEPoro <filename>" << endl;
        exit(0);
    }

    // Step 1. Read params from an input file to internal params data structure.
    // Note: The keywords are almost compatible with Ecl simulator.
    string    myfile = argv[1];
    ParamRead rp;
    rp.readInputFile(myfile);

    OpenCAEPoro simulator;

    // Step 2. Read param from internal params data structure to each modules, and
    // setup static information, such as active grids, and connections between them.
    // Note: Memory allocation for linear systems will also be done at this time.
    simulator.setup(rp);

    // Step 3. Initialize the reservoir, which finishs the first step in iterations.
    // For example: initial pressure, saturation, and moles of components will be
    // calculated. Initial well pressure will also be given here ---- it's a simple
    // guess now.
    simulator.init();

    // Step 4. Run dynamic simulation using methods like IMPES and FIM. It's a
    // combination of functions of various modules which you could make changes
    // whenever necessary.
    simulator.run();

    // Step 5. Output the results according to params users give. It will generate
    // a summary file in your input data directory.
    simulator.out();

    return 0;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/