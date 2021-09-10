#pragma once
#include<vector>
#include <iostream>
#include "Grid.h"
#include "Mixture.h"
#include "FlowUnit.h"

using namespace std;
// Bulk contains the infomation of each bulk
// variables are ordered according to the time consuming of program


class Bulk
{
	friend class	Connection_BB;
	friend class	Well;

public:
	Bulk() = default;

	void init(const Grid& myGrid);

	void initSjPc();

	// Flash
	void flash_Sj();
	void flash_Ni();
	void passValue(int n);

	// relative permeability and capillary pressure
	void calKrPc();
	// Rock
	void calporo();

private:
	int						Num;			// num of active bulk
	// Physical infomation
	int						Np;				// num of phase
	int						Nc;				// num of component

	std::vector<double>		T;				// temperature : Num
	std::vector<double>		Pbb;			// buble point pressere: Num
	std::vector<double>		P;				// pressure: Num
	std::vector<double>		Pc;				// capillary pressure of phase: Np*Num
	std::vector<bool>		PhaseExist;		// existence of phase 
	std::vector<double>		Ni;				// molar of ith component in bulk: NC*Num
	std::vector<double>		S;				// saturation of phase j
	std::vector<double>		Xi;				// molar density of phase: Np*Num
	std::vector<double>		Cij;			// Nij / Nj : Np*Nc*Num 
	std::vector<double>		Rho;			// mass density of phase: Np*Num
	std::vector<double>		Mu;				// viscosity of phase: Np*Num
	std::vector<double>		Kr;				// relative permeability of phase: Np*Num


	std::vector<double>		Vf;				// total fluid volume
	std::vector<double>		Vfi;			// dVt / dNi
	std::vector<double>		Vfp;			// dVt / dP


	Mixture*				Flashcal;
	FlowUnit*				Flow;

	
	
	// Bulk rock infomation
	std::vector<double>		Dx;					// dx
	std::vector<double>		Dy;					// dy
	std::vector<double>		Dz;					// dz
	std::vector<double>		Depth;				// depth: Num
	std::vector<double>		Ntg;				// Ntg: Num
	std::vector<double>		Rock_V;				// Vgrid * ntg
	std::vector<double>		Rock_Poro;			// current porosity
	std::vector<double>		Rock_PoroInit;		// initial porosity
	double					Rock_Pref;
	double					Rock_C1;
	double					Rock_C2;
	std::vector<double>		Rock_KxInit;
	std::vector<double>		Rock_Kx;
	std::vector<double>		Rock_KyInit;
	std::vector<double>		Rock_Ky;
	std::vector<double>		Rock_KzInit;
	std::vector<double>		Rock_Kz;

	// Reservoir information
	std::vector<double>		EQUIL;
};
