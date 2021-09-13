#pragma once
#include <fstream>
#include<vector>
#include "ReadTool.hpp"

using namespace std;

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
	std::vector<std::vector<double>>	SWOF;
	std::vector<std::vector<double>>	SGOF;
	std::vector<double>					EQUIL;
	std::vector<std::vector<double>>	PBVD;
	
	// phase property
	std::vector<double>					Density;
	std::vector<double>					Gravity;

	// PVT property
	std::vector<std::vector<double>>	PVCO;
	std::vector<std::vector<double>>	PVDO;
	std::vector<std::vector<double>>	PVTW;

	// method
	double* FindPtr(string& varName);

	// init size of var
	void initVar();
	void setVal(double* obj, double val, vector<int>& index);
	void copyVal(double* obj, double* src, vector<int>& index);

	// DIMENS
	void inputDIMENS(ifstream& ifs);
	void outputDIMENS();

	// EUQALS
	void inputEQUALS(ifstream& ifs);

	// COPY
	void inputCOPY(ifstream& ifs);
};
