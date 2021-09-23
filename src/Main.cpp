#include <iostream>
#include <cstdio>
#include <string>
#include "Timing.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "OpenCAEPoro.hpp"

using namespace std;

int& test(int& a) {
	return a;
}

int main()
{
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64) 
	string myfile{"D:\\Lsz\\PennSim\\input.txt"};
#else
	string myfile{"/mnt/d/Lsz/PennSim/input.txt"};
#endif
	
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
