#pragma once
#include <fstream>
#include <vector>
#include "ReadTool.hpp"

class COOIJK
{
public:
	COOIJK(int i, int j, int k) :I(i), J(j), K(k) {};
	int		I;
	int		J;
	int		K;
};

class Type_B
{
public:
	bool				status{ false };
	vector<COOIJK>		obj;
};

class Type_A
{
public:
	bool				status{ false };
	vector<string>		obj;
};


class ParamOutput
{
public:
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

	Type_A		WOPR;
	Type_A		WOPT;
	Type_A		WGPR;
	Type_A		WGPT;
	Type_A		WWPR;
	Type_A		WWPT;
	Type_A		WGIR;
	Type_A		WGIT;
	Type_A		WWIR;
	Type_A		WWIT;
	Type_A		WBHP;

	Type_B		BPR;

	// Method

	void inputSUMMARY(ifstream& ifs);
	void inputType_A(ifstream& ifs, Type_A& obj);
	void inputType_B(ifstream& ifs, Type_B& obj);
};
