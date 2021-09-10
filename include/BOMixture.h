#pragma once
#include "Mixture.h"
#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.h"

class BOMixture : public Mixture
{
public:
	BOMixture() = default;

	void Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin) override;
	void BOFlash_Sj_W(const double Pin, const double* Sjin, double Vpore, const double* Ziin);
	void BOFlash_Sj_OW(const double Pin, const double* Sjin, double Vpore, const double* Ziin);
	void BOFlash_Sj_OGW(const double Pin, const double Pbbin, const double* Sjin, double Vpore, const double* Ziin);

	void Flash_Ni(const double Pin, const double Tin, const double* Niin) override;
	void BOFlash_Ni_W(const double Pin, const double* Niin);
	void BOFlash_Ni_OW(const double Pin, const double* Niin);
	void BOFlash_Ni_OGW(const double Pin, const double* Niin);

	// return rho
	double calRhoO(double Pin, double Pbbin) override;
	double calRhoG(double Pin) override;
	double calRhoW(double Pin) override;
	double calRhoO_OW(double Pin);
	double calRhoO_OGW(double Pin, double Pbbin);

private:

	int									Mode;
	ReservoirTable<double>				PVCO;
	ReservoirTable<double>				PVDG;
	ReservoirTable<double>				PVTW;
	ReservoirTable<double>				PVDO;


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR
	double								Std_RhoO;
	double								Std_RhoG;
	double								Std_RhoW;

};
