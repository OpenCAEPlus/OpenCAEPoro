#pragma once
#include <fstream>
#include <vector>
#include "ReadTool.hpp"

const int IMPES					= 0;
const int FIM					= 1;
const int param_BLKOIL			= 0;
const int param_EoS				= 1;

class ParamControl
{
public:

	int										Model;   // must be given
	int										Method;
	string									LinearSolve;
	std::vector<std::vector<double>>		Tuning;
	
	void init();
	void initMethod();
	void initTuning();
	void inputMETHOD(ifstream& ifs);
	void inputTUNING(ifstream& ifs);
	void showTuning();
};

