#pragma once
#include <fstream>
#include<vector>
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

	// Rock
	std::vector<double>					Ntg;
	std::vector<double>					Poro;
	std::vector<double>					PermX;
	std::vector<double>					PermY;
	std::vector<double>					PermZ;
	double								Pref;
	double								C1;
	double								C2;

	// Saturation table & buble point pressure 
	TableSet							SWOF_T;
	TableSet							SGOF_T;
	TableSet							PBVD_T;
	std::vector<double>					EQUIL;
	
	// phase property
	std::vector<double>					Density;
	std::vector<double>					Gravity;

	// PVT property
	TableSet							PVCO_T;
	TableSet							PVDG_T;
	TableSet							PVTW_T;

	// method
	vector<double>* FindPtr(string& varName);
	TableSet* FindPtr_T(string& varName);

	// init size of var
	void initVar();
	void initTab();
	void setVal(vector<double>& obj, double val, vector<int>& index);
	void copyVal(vector<double>& obj, vector<double>& src, vector<int>& index);
	void multiplyVal(vector<double>& obj, double val, vector<int>& index);

	// DIMENS
	void inputDIMENS(ifstream& ifs);
	void outputDIMENS();

	// EUQALS
	void inputEQUALS(ifstream& ifs);

	// COPY
	void inputCOPY(ifstream& ifs);

	// MULTIPLY
	void inputMULTIPLY(ifstream& ifs);

	// Table
	void inputTABLE(ifstream& ifs, string tabName);
};
