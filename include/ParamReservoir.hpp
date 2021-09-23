#pragma once
#include <fstream>
#include <vector>
#include "ReadTool.hpp"

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
	vector<vector<vector<double>>>		data;
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
	double	Pref;
	double	Cr;
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
	std::vector<double>					Tops;
	std::vector<double>					Dx;
	std::vector<double>					Dy;
	std::vector<double>					Dz;

	double								RTEMP;

	// Rock
	std::vector<double>					Ntg;
	std::vector<double>					Poro;
	std::vector<double>					PermX;
	std::vector<double>					PermY;
	std::vector<double>					PermZ;
	ROCK								Rock;

	// Restart
	std::vector<double>					Pressure;
	std::vector<double>					Ni;

	// phase property
	Type_A_r<double>					Density;
	Type_A_r<double>					Gravity;

	// Model and Phase
	bool								BLACKOIL{ false };
	bool								COMPS{ false };
	bool								OIL{ false };
	bool								GAS{ false };
	bool								WATER{ false };
	bool								DISGAS{ false };

	// Eos
	std::vector<double>					InitZi;

	// SAT Region & PVT Region
	int									NTSFUN{ 1 };	// SAT num
	int									NTPVT{ 1 };		// PVT num
	Type_A_r<double>					SATNUM;
	Type_A_r<double>					PVTNUM;


	// Saturation table & buble point pressure 
	TableSet							SWOF_T;
	TableSet							SGOF_T;
	TableSet							PBVD_T;
	std::vector<double>					EQUIL;
	
	// PVT property
	int									Np;		// num of phase
	int									Nc;		// num of comp, used for Eos or Restart
	TableSet							PVCO_T;
	TableSet							PVDO_T;
	TableSet							PVDG_T;
	TableSet							PVTW_T;

	// internal method
	vector<double>* FindPtr(string& varName);

	TableSet* FindPtr_T(string& varName);

	// init size of var
	void init();
	void initTab();

	template<typename T>
	void setVal(vector<T>& obj, T val, vector<int>& index);

	template<typename T>
	void copyVal(vector<T>& obj, vector<T>& src, vector<int>& index);

	void multiplyVal(vector<double>& obj, double val, vector<int>& index);

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
