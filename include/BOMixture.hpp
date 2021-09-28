#pragma once
#include "Mixture.hpp"
#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamReservoir.hpp"

class BOMixture : public Mixture
{
public:
	BOMixture() = default;
	BOMixture(ParamReservoir& rs_param, int PVTmode, int i);

	bool empty_PVDG() override{ return PVDG.isempty(); }

	// Flash
	void Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin) override;
	void BOFlash_Sj_W(const double Pin, const double* Sjin, double Vpore);
	void BOFlash_Sj_OW(const double Pin, const double* Sjin, double Vpore);
	void BOFlash_Sj_OGW(const double Pin, const double Pbbin, const double* Sjin, double Vpore);

	void Flash_Ni(const double Pin, const double Tin, const double* Niin) override;
	void BOFlash_Ni_W(const double Pin, const double* Niin);
	void BOFlash_Ni_OW(const double Pin, const double* Niin);
	void BOFlash_Ni_OGW(const double Pin, const double* Niin);

	// return Xi  molar density
	double xiPhase(double Pin, double T, double* Ziin) override;
	double xiPhase_OGW(double Pin, double* Ziin);

	// return rho
	double rhoPhase(double Pin, double T, double* Ziin) override;
	double rhoPhase_OGW(double Pin, double* Ziin);


	// return gamma
	double gammaPhaseO(double Pin, double Pbbin) override;
	double gammaPhaseG(double Pin) override;
	double gammaPhaseW(double Pin) override;
	double gammaPhaseO_OW(double Pin);
	double gammaPhaseO_OGW(double Pin, double Pbbin);
	double gammaPhaseOG(double Pin, double Tin, double* Ziin) override { ERRORcheck("should not be used in BLKOIL"); exit(0); };

private:

	int									Mode;
	ReservoirTable<double>				PVCO;
	ReservoirTable<double>				PVDG;
	ReservoirTable<double>				PVTW;
    ReservoirTable<double>              PVDO;

	// Auxiliary parameters for Table interpolation
	int									len{ 0 };
	vector<double>						data;
	vector<double>						cdata;


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR
	double								Std_RhoO, Std_GammaO;
	double								Std_RhoG, Std_GammaG;
	double								Std_RhoW, Std_GammaW;

};
