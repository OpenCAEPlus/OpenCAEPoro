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

	string myfile{"D:\\Lsz\\PennSim\\input.txt"};
	ParamRead rp;
	rp.getDirAndName(myfile);
	rp.init();
	rp.readFile(myfile);
	rp.checkParam();
	cout << "Done !" << endl;

	OpenCAEPoro simulator;
	simulator.inputParam(rp);
	simulator.setup();


	return 0;
}
