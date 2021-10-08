#pragma once

// Standard header files
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoro header files
#include "ParamReservoir.hpp"
#include "ParamWell.hpp"
#include "ParamControl.hpp"
#include "ParamOutput.hpp"
#include "ReadTool.hpp"

using namespace std;

/// ParamRead is pre-processing part in our simulator, which is responsible for inputting params
/// from files supplied by users, it's almost compatible with Eclipse but has own rules for easy to use. 
/// it is extensible friendly.
class ParamRead
{
public:
	string					File;
	string					FileDir;
	string					FileName;

	ParamReservoir			Rs_param;
	ParamWell				Well_param;
	ParamControl			Control_param;
	ParamOutput				Output_param;

	void readInputFile(string& file);
	void getDirAndName();
	void readFile(string file);

	// init
	void init();
	// check
	void checkParam();

	// INCLUDE
	void inputINCLUDE(ifstream& ifs);
};
