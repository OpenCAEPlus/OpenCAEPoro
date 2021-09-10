#pragma once
#include <vector>
#include <iostream>

class Mixture
{
	friend class Bulk;
public:

	Mixture() = default;
	virtual ~Mixture() = 0;

	virtual void Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin) = 0;
	virtual void Flash_Ni(const double Pin, const double Tin, const double* Niin) = 0;

	virtual void getProp() {};

	// return rho
	virtual double calRhoO(double Pin, double Pbbin) = 0;
	virtual double calRhoG(double Pin) = 0;
	virtual double calRhoW(double Pin) = 0;

protected:

	int								Np;			// num of phase
	int								Nc;			// num of component
	double							P;			// Pressure
	double							T;			// Temperature

	
	std::vector<double>				Ni;			// molar of component : Nc
	std::vector<bool>				PhaseExist; // existence of phase : Np
	std::vector<double>				Cij;		// Nij / Nj : Np*Nc
	std::vector<double>				S;			// saturation of phase : Np
	std::vector<double>				Xi;			// molar density of phase: Np
	std::vector<double>				Rho;		// mass density of phase : Np
	std::vector<double>				Mu;			// viscosity of phase: Np
	
	double							Vf;			// volume of fluids
	double							Vfp;		// 
	std::vector<double>				Vfi;		// dVf / dNi   : Nc

};
