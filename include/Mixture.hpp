#pragma once
#include <vector>
#include <iostream>

// total fluid, contains all phase
class Mixture
{
	friend class Bulk;
public:

	Mixture() = default;
	virtual ~Mixture() = 0;

	// return type
	int getType() { return MixtureType; }
	// black oil
	virtual bool empty_PVDG() = 0;

	virtual void Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin) = 0;
	virtual void Flash_Ni(const double Pin, const double Tin, const double* Niin) = 0;

	virtual void getProp() {};

	// return rho
	virtual double gammaPhaseO(double Pin, double Pbbin) = 0;
	virtual double gammaPhaseW(double Pin) = 0;
	virtual double gammaPhaseG(double Pin) = 0;
	virtual double gammaPhaseOG(double Pin, double Tin, double* Ziin) = 0;

protected:

	int								MixtureType;

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
