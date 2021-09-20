#pragma once
#include "Reservoir.hpp"
#include "ParamOutput.hpp"

class CAEIJK
{
public:
	CAEIJK() = default;
	CAEIJK(int i, int j, int k) :I(i), J(j), K(k) {};
	CAEIJK(COOIJK& src) { I = src.I; J = src.J; K = src.K; };
	CAEIJK& operator= (COOIJK& src) { I = src.I; J = src.J; K = src.K; }
	int			I, J, K;
};

class CAETypeA
{
public:
	CAETypeA& operator= (Type_A_o& src) { activity = src.activity; obj = src.obj; }
	bool				activity{ false };
	vector<string>		obj;
};

class CAETypeB
{
public:
	CAETypeB& operator= (Type_B_o& src) { activity = src.activity; obj.assign(src.obj.begin(), src.obj.end()); }
	bool				activity{ false };
	vector<CAEIJK>		obj;
};

class SumPair
{
public:
	SumPair(string item, string obj, string unit) :Item(item), Obj(obj), Unit(unit) {};
	string				Item;
	string				Obj;
	string				Unit;
	vector<double>		val;
};


class Summary
{
public:
	void inputParam(OutputSummary& summary_param);
	void setup(Reservoir& reservoir);

private:

	vector<SumPair>			Sumdata;

	bool		FPR{ false };
	bool		FOPR{ false };
	bool		FOPT{ false };
	bool		FGPR{ false };
	bool		FGPT{ false };
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



class CAEOutput
{
public:

	void inputParam(ParamOutput& Output_param);
	void setup(Reservoir& reservoir);

private:
	Summary		summary;

};
