#pragma once

// Standard header files
#include <iostream>
#include <iomanip>

// OpenCAEPoro header files
#include "Reservoir.hpp"
#include "ParamOutput.hpp"
#include "OpenCAEControl.hpp"

using namespace std;

class CAEIJK
{
public:
	CAEIJK() = default;
	CAEIJK(int i, int j, int k) :I(i), J(j), K(k) {};
	CAEIJK(COOIJK& src) { I = src.I; J = src.J; K = src.K; };
	CAEIJK& operator= (COOIJK& src) { I = src.I; J = src.J; K = src.K; return*this; }
	int			I, J, K;
};

class OCPTypeA
{
public:
	OCPTypeA& operator= (Type_A_o& src) { activity = src.activity; obj = src.obj; return *this; }
	bool				activity{ false };
	vector<string>		obj;
	vector<int>			index;
};

/// 
class OCPTypeB
{
public:
	OCPTypeB& operator= (Type_B_o& src) { activity = src.activity; obj.assign(src.obj.begin(), src.obj.end()); return *this; }
	bool				activity{ false };
	vector<CAEIJK>		obj;
	vector<int>			index;
};

/// SumPair is an auxiliary structure storing summary data to output.
class SumPair
{
public:
	SumPair(string item, string obj, string unit) :Item(item), Obj(obj), Unit(unit) {};
	string				Item;
	string				Obj;
	string				Unit;
	vector<OCP_DBL>		val;
};


/// Summary manages the output for summary file, it contains the most interested information
/// in each time step. usually these data will be convert to chart for analysing following.
class Summary
{
public:
	void inputParam(OutputSummary& summary_param);
	void setup(Reservoir& reservoir);
	void setVal(const Reservoir& reservoir, const OCP_Control& ctrl);
	void printInfo(string& dir);


private:

	vector<SumPair>			Sumdata;

	bool		FPR{ false };
	bool		FOPR{ false };
	bool		FOPT{ false };
	bool		FGPR{ false };
	bool		FGPt{ false };
	bool		FWPR{ false };
	bool		FWPT{ false };
	bool		FGIR{ false };
	bool		FGIT{ false };
	bool		FWIR{ false };
	bool		FWIT{ false };

	OCPTypeA		WOPR;
	OCPTypeA		WOPT;
	OCPTypeA		WGPR;
	OCPTypeA		WGPT;
	OCPTypeA		WWPR;
	OCPTypeA		WWPT;
	OCPTypeA		WGIR;
	OCPTypeA		WGIT;
	OCPTypeA		WWIR;
	OCPTypeA		WWIT;
	OCPTypeA		WBHP;

	OCPTypeB		BPR;
};


/// OCP_Output manages different kinds of ways to output. the most commonly used is summary file.
/// which usually give the information of bulks and wells in each timestep, such as average bulks pressure,
/// oil production rate of wells. if other information at critical time is interested in, you can chose
/// the PRT file(to do). also, some infomation will be printed on the screen at the critical time to make sure
/// the program is at the right way. 
class OCP_Output
{
	friend class OpenCAEPoro;
public:

	void inputParam(ParamOutput& Output_param);
	void setup(Reservoir& reservoir, string& dir);
	void setVal(const Reservoir& reservoir, const OCP_Control& ctrl);
	void printInfo();

private:
	string		Dir;
	Summary		summary;

};
