#include "OpenCAEPoro.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "Timing.hxx"
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;


int main(int argc, const char* argv[])
{


    if (argc == 1) {
        cout << "Input file is missing. Usage: ./OpenCAEPoro <filename>" << endl;
        exit(0);
    }

    string    myfile = argv[1];
    ParamRead rp;
    rp.readInputFile(myfile);

    OpenCAEPoro simulator;
    simulator.setup(rp);

    simulator.init();

    simulator.run();

    simulator.out();

    return 0;
}
