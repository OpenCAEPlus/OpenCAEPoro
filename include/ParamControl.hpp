#pragma once

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "ReadTool.hpp"
#include "OpenCAEPoro_consts.hpp"

typedef vector<vector<OCP_DBL>>				TUNING;

class TuningPair
{
public:
	TuningPair(int t, TUNING& tmp) :d(t), Tuning(tmp) {};
	int				d;
	TUNING			Tuning;
};

class ParamControl
{
public:

	string									Dir;
	string									Method;
	string									LinearSolve;
	vector<TuningPair>						Tuning_T;
	TUNING									Tuning;
	vector<OCP_DBL>							CriticalTime;
	
	void init(string& dir);
	void initTime() { CriticalTime.push_back(0); };
	void initMethod();
	void initTuning();
	void inputMETHOD(ifstream& ifs);
	void inputTUNING(ifstream& ifs);
	void showTuning();
};



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/