#pragma once

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "ReadTool.hpp"

class COOIJK
{
public:
	COOIJK() = default;
	COOIJK(int i, int j, int k) :I(i), J(j), K(k) {};
	int		I;
	int		J;
	int		K;
};

class Type_B_o
{
public:
	bool				activity{ false };
	vector<COOIJK>		obj;
};

class Type_A_o
{
public:
	bool				activity{ false };
	vector<string>		obj;
};

class OutputSummary 
{

public:
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

	Type_A_o		WOPR;
	Type_A_o		WOPT;
	Type_A_o		WGPR;
	Type_A_o		WGPT;
	Type_A_o		WWPR;
	Type_A_o		WWPT;
	Type_A_o		WGIR;
	Type_A_o		WGIT;
	Type_A_o		WWIR;
	Type_A_o		WWIT;
	Type_A_o		WBHP;

	Type_B_o		BPR;
};

class ParamOutput
{	
public:

	OutputSummary	Summary;
	// Method

	void inputSUMMARY(ifstream& ifs);
	void inputType_A(ifstream& ifs, Type_A_o& obj);
	void inputType_B(ifstream& ifs, Type_B_o& obj);

	
};
