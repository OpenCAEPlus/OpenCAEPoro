#include <iostream>
#include <cstdio>
#include <string>
#include "Timing.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "OpenCAEPoro.hpp"

using namespace std;



int main(int argc, const char* argv[])
{

	if (argc == 1) {
		cout << "Input file is missing. Usage: ./OpenCAEPoro <filename>" << endl;
		exit(0);
	}


	string myfile = argv[1];
	ParamRead rp;
	rp.readInputFile(myfile);

	OpenCAEPoro simulator;
	simulator.inputParam(rp);
	simulator.setup();
	simulator.allocateMat();

	simulator.init();
	simulator.run();


	return 0;
}
