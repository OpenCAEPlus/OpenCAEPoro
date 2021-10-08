#include "OpenCAEPoro.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "Timing.hxx"
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;

// the main function, also an example to display main steps in our simulator.
// it basicly divided into five parts.
// read param from input file --> setup static information with input param. 
// --> initialize --> core simulations --> output the results.
int main(int argc, const char* argv[])
{


    if (argc == 1) {
        cout << "Input file is missing. Usage: ./OpenCAEPoro <filename>" << endl;
        exit(0);
    }

    // first step.
    // read params from input file to internal params data structure. 
    // the format of keyword is almost compatible with Eclipse.
    // the data structure is independent of main program.
    string    myfile = argv[1];
    ParamRead rp;
    rp.readInputFile(myfile);

    OpenCAEPoro simulator;

    // second step.
    // read param from internal params data structure to each modules. and then 
    // setup static information, such as active grids, and connections between them,
    // memory allocating for linear system will also be done at this time.
    simulator.setup(rp);

    // third step.
    // initialize the reservoir, which finishs the first step in iterations.
    // for example, initial pressure, saturation and moles of components will be calculated.
    // initial well pressure will also be given here ---- it's a simple guess now.
    simulator.init();

    // fouth step
    // core simulations, you can choose method such as IMPES. actually, it's a combination
    // of functions of each main modules, which means you could do some change if need, even 
    // create new method with these functions.
    simulator.run();

    /// fifth step
    /// output the results, which ups to input param users give. it will generate a summary file
    /// in your data directory.
    simulator.out();

    return 0;
}
