#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include "ParamReservoir.hpp"
#include "ParamWell.hpp"
#include "ParamControl.hpp"
#include "ParamOutput.hpp"
#include "ReadTool.hpp"

using namespace std;

class ParamRead
{
public:
	string					FileDir;
	string					FileName;

	ParamReservoir			Rs_param;
	ParamWell				Well_param;
	ParamControl			Control_param;
	ParamOutput				Output_param;


	void getDirAndName(string file);
	void readFile(string file);

	// init
	void init();

	// INCLUDE
	void inputINCLUDE(ifstream& ifs);

	// check
	void checkParam();
};
