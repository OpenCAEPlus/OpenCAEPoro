#include <iostream>
#include <cstdio>
#include <string>
#include "Connection_BB.hpp"
#include "WellGroup.hpp"
#include "Timing.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamRead.hpp"
#include "ReadTool.hpp"

using namespace std;


int main()
{

	vector<string>  vbuf;
	vbuf.push_back("INJE*");
	vbuf.push_back("3*");
	vbuf.push_back("PROD*");
	vbuf.push_back("1*");
	vbuf.push_back("*");

	DealDefault(vbuf);

	string myfile{"D:\\Lsz\\PennSim\\input.txt"};
	ParamRead rp;
	rp.getDirAndName(myfile);
	rp.init();
	rp.readFile(myfile);
	rp.checkParam();
	cout << "Done !" << endl;

	return 0;
}
