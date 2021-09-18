#pragma once
#include <fstream>
#include <vector>
#include "ReadTool.hpp"



class WellOptParam
{
public:
	WellOptParam(string type, vector<string>& vbuf);
	// WCONINJE & WCONPROD
	string				Type;
	string				FluidType;
	string				State;
	string				OptMode;

	double				MaxRate;		
	double				MaxBHP;	
	double				MinBHP;
};

class WellOptPair
{
public:
	WellOptPair(int i, string type, vector<string>& vbuf) :d(i), Opt(type, vbuf) {};
	int					d;
	WellOptParam		Opt;
};

class WellParam
{
public:
	WellParam(vector<string>& info);
	// static infomation
	// WELSPECS
	string		Name;
	string		Group{ "FEILD" };
	string		Direction{ "z" };
	int			I;
	int			J;
	double		Dref{ -1.0 };
	// COMPDAT
	int			I_perf;
	int			J_perf;
	int			K1;
	int			K2;
	double		Trans{ -1.0 };
	double		Diameter{ 1.0 };
	double		Kh{ -1.0 };
	double		SkinFactor{ 0.0 };

	// dynamic infomation
	vector<WellOptPair>		OptParam;

};

class ParamWell
{
public:

	int								WellNum;
	std::vector<WellParam>			well;
	std::vector<double>				CriticalTime;

	void init() { initTime(); };
	void initTime() { CriticalTime.push_back(0); };
	void inputWELSPECS(ifstream& ifs);
	void inputCOMPDAT(ifstream& ifs);
	void inputWCONINJE(ifstream& ifs);
	void inputWCONPROD(ifstream& ifs);
	void inputTSTEP(ifstream& ifs);
	void inputWELTARG(ifstream& ifs);

	// check
	void checkParam();
	void checkPerf();
};
