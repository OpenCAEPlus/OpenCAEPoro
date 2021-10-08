#pragma once

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "ReadTool.hpp"
#include "OpenCAEPoro_consts.hpp"

using namespace std;

class TableSet
{
public:
	void showTab() 
	{
		cout << Name << "\n";
		for (auto v : data) {
			int len = v[0].size();
			for (int i = 0; i < len; i++) {
				for (int j = 0; j < colNum; j++) {
					cout << v[j][i] << "\t";
				}
				cout << "\n";
			}
			cout << "----------------\n" ;
		}
	}
	string								Name;
	int									colNum;
	vector<vector<vector<OCP_DBL>>>		data;
};

class DIMENS
{
public:
	int Nx;
	int Ny;
	int Nz;
};

class ROCK
{
public:
	OCP_DBL	Pref;
	OCP_DBL	Cr;
};

template<typename T>
class Type_A_r
{
public:
	bool				activity{ false };
	vector<T>			data;
};

class ParamReservoir
{

public:
	// Grid   
	// Cartesian
	DIMENS								Dimens;
	int									Num;
	std::vector<OCP_DBL>					Tops;
	std::vector<OCP_DBL>					Dx;
	std::vector<OCP_DBL>					Dy;
	std::vector<OCP_DBL>					Dz;

	OCP_DBL								RTEMP;

	// Rock
	std::vector<OCP_DBL>					Ntg;
	std::vector<OCP_DBL>					Poro;
	std::vector<OCP_DBL>					PermX;
	std::vector<OCP_DBL>					PermY;
	std::vector<OCP_DBL>					PermZ;
	ROCK								Rock;

	// Restart
	std::vector<OCP_DBL>					Pressure;
	std::vector<OCP_DBL>					Ni;

	// phase property
	Type_A_r<OCP_DBL>					Density;
	Type_A_r<OCP_DBL>					Gravity;

	// Model and Phase
	bool								BLACKOIL{ false };
	bool								COMPS{ false };
	bool								OIL{ false };
	bool								GAS{ false };
	bool								WATER{ false };
	bool								DISGAS{ false };

	// Eos
	std::vector<OCP_DBL>					InitZi;

	// SAT Region & PVT Region
	int									NTSFUN{ 1 };	// SAT num
	int									NTPVT{ 1 };		// PVT num
	Type_A_r<OCP_DBL>					SATNUM;
	Type_A_r<OCP_DBL>					PVTNUM;


	// Saturation table & buble point pressure 
	TableSet							SWOF_T;
	TableSet							SGOF_T;
	TableSet							PBVD_T;
	std::vector<OCP_DBL>					EQUIL;
	
	// PVT property
	int									Np;		// num of phase
	int									Nc;		// num of comp, used for Eos or Restart
	TableSet							PVCO_T;
	TableSet							PVDO_T;
	TableSet							PVDG_T;
	TableSet							PVTW_T;

	// internal method
	vector<OCP_DBL>* FindPtr(string& varName);

	TableSet* FindPtr_T(string& varName);

	// init size of var
	void init();
	void initTab();

	template<typename T>
	void setVal(vector<T>& obj, T val, vector<int>& index);

	template<typename T>
	void copyVal(vector<T>& obj, vector<T>& src, vector<int>& index);

	void multiplyVal(vector<OCP_DBL>& obj, OCP_DBL val, vector<int>& index);

	// COMPS
	void inputCOMPS(ifstream& ifs);

	// DIMENS
	void inputDIMENS(ifstream& ifs);
	void outputDIMENS();

	// RTEMP
	void inputRTEMP(ifstream& ifs);

	// EQUALS
	void inputEQUALS(ifstream& ifs);

	// EQUALS   ----   supplement
	void inputGRID(ifstream& ifs, string& keyword);

	// COPY
	void inputCOPY(ifstream& ifs);

	// MULTIPLY
	void inputMULTIPLY(ifstream& ifs);

	// Table
	void inputTABLE(ifstream& ifs, string& tabName);

	// Rock
	void inputROCK(ifstream& ifs);

	// Gravity
	void inputGRAVITY(ifstream& ifs);

	// Density
	void inputDENSITY(ifstream& ifs);

	// EQUAL
	void inputEQUIL(ifstream& ifs);

	// SATNUM & PVTNUM  -- Region
	void inputTABDIMS(ifstream& ifs);
	void inputRegion(ifstream& ifs, string keyword);


	// check
	void checkParam();
	void checkGrid();
	void checkEQUIL();
	void checkDenGra();
	void checkPhase();
	void checkPhaseTab();
	void checkRegion();
	void checkEqlRegion();

};


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/