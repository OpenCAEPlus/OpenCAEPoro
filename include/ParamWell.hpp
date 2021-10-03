#pragma once
#include <fstream>
#include <vector>
#include "ReadTool.hpp"
#include "OpenCAEPoro_consts.hpp"



class WellOptParam
{
public:
	WellOptParam(string type, vector<string>& vbuf);
	// WCONINJE & WCONPROD
	string				Type;
	string				FluidType;   // inj
	string				State;
	string				OptMode;

	OCP_DBL				MaxRate;		
	OCP_DBL				MaxBHP;	
	OCP_DBL				MinBHP;

	vector<OCP_DBL>		Zi;
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
	OCP_DBL		Dref{ -1.0 };
	// COMPDAT
	int			I_perf;
	int			J_perf;
	int			K1;
	int			K2;
	OCP_DBL		WI{ -1.0 };		// connection factor
	OCP_DBL		Diameter{ 1.0 };
	OCP_DBL		Kh{ -1.0 };
	OCP_DBL		SkinFactor{ 0.0 };

	// dynamic infomation
	vector<WellOptPair>		OptParam;

};

class ParamWell
{
public:

	std::vector<WellParam>			well;
	std::vector<OCP_DBL>				CriticalTime;

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
