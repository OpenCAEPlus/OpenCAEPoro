#pragma once
#include "Mixture.hpp"
#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.hpp"

class BOMixture : public Mixture
{
public:
	BOMixture() = default;

	bool empty_PVDG() override{ return PVDG.isempty(); }

	void Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin) override;
	void BOFlash_Sj_W(const double Pin, const double* Sjin, double Vpore);
	void BOFlash_Sj_OW(const double Pin, const double* Sjin, double Vpore);
	void BOFlash_Sj_OGW(const double Pin, const double Pbbin, const double* Sjin, double Vpore);

	void Flash_Ni(const double Pin, const double Tin, const double* Niin) override;
	void BOFlash_Ni_W(const double Pin, const double* Niin);
	void BOFlash_Ni_OW(const double Pin, const double* Niin);
	void BOFlash_Ni_OGW(const double Pin, const double* Niin);

	// return rho
	double gammaPhaseO(double Pin, double Pbbin) override;
	double gammaPhaseG(double Pin) override;
	double gammaPhaseW(double Pin) override;
	double gammaPhaseO_OW(double Pin);
	double gammaPhaseO_OGW(double Pin, double Pbbin);

private:

	int									Mode;
	ReservoirTable<double>				PVCO;
	ReservoirTable<double>				PVDG;
	ReservoirTable<double>				PVTW;
	ReservoirTable<double>				PVDO;


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR
	double								Std_RhoO, Std_GammaO;
	double								Std_RhoG, Std_GammaG;
	double								Std_RhoW, Std_GammaW;

};
