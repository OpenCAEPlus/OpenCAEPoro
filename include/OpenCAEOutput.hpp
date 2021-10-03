#pragma once
#include <iostream>
#include <iomanip>
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

class CAETypeA
{
public:
	CAETypeA& operator= (Type_A_o& src) { activity = src.activity; obj = src.obj; return *this; }
	bool				activity{ false };
	vector<string>		obj;
	vector<int>			index;
};

class CAETypeB
{
public:
	CAETypeB& operator= (Type_B_o& src) { activity = src.activity; obj.assign(src.obj.begin(), src.obj.end()); return *this; }
	bool				activity{ false };
	vector<CAEIJK>		obj;
	vector<int>			index;
};

class SumPair
{
public:
	SumPair(string item, string obj, string unit) :Item(item), Obj(obj), Unit(unit) {};
	string				Item;
	string				Obj;
	string				Unit;
	vector<OCP_DBL>		val;
};


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

	CAETypeA		WOPR;
	CAETypeA		WOPT;
	CAETypeA		WGPR;
	CAETypeA		WGPT;
	CAETypeA		WWPR;
	CAETypeA		WWPT;
	CAETypeA		WGIR;
	CAETypeA		WGIT;
	CAETypeA		WWIR;
	CAETypeA		WWIT;
	CAETypeA		WBHP;

	CAETypeB		BPR;
};



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
