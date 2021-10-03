#pragma once
#include "Mixture.hpp"
#include "ReservoirTable.hxx"
#include "ParamReservoir.hpp"

class BOMixture : public Mixture
{
public:
	BOMixture() = default;
	BOMixture(ParamReservoir& rs_param, int PVTmode, int i);

	bool empty_PVDG() override{ return PVDG.isempty(); }

	// Flash
	void Flash_Sj(const OCP_DBL Pin, const OCP_DBL Pbbin, const OCP_DBL Tin, const OCP_DBL* Sjin, OCP_DBL Vpore, const OCP_DBL* Ziin) override;
	void BOFlash_Sj_W(const OCP_DBL Pin, const OCP_DBL* Sjin, OCP_DBL Vpore);
	void BOFlash_Sj_OW(const OCP_DBL Pin, const OCP_DBL* Sjin, OCP_DBL Vpore);
	void BOFlash_Sj_OGW(const OCP_DBL Pin, const OCP_DBL Pbbin, const OCP_DBL* Sjin, OCP_DBL Vpore);

	void Flash_Ni(const OCP_DBL Pin, const OCP_DBL Tin, const OCP_DBL* Niin) override;
	void BOFlash_Ni_W(const OCP_DBL Pin, const OCP_DBL* Niin);
	void BOFlash_Ni_OW(const OCP_DBL Pin, const OCP_DBL* Niin);
	void BOFlash_Ni_OGW(const OCP_DBL Pin, const OCP_DBL* Niin);

	// return Xi  molar density
	OCP_DBL xiPhase(const OCP_DBL Pin, const OCP_DBL T, const OCP_DBL* Ziin) override;
	OCP_DBL xiPhase_OGW(const OCP_DBL Pin, const OCP_DBL* Ziin);

	// return rho
	OCP_DBL rhoPhase(OCP_DBL Pin, OCP_DBL T, OCP_DBL* Ziin) override;
	OCP_DBL rhoPhase_OGW(OCP_DBL Pin, OCP_DBL* Ziin);


	// return gamma
	OCP_DBL gammaPhaseO(OCP_DBL Pin, OCP_DBL Pbbin) override;
	OCP_DBL gammaPhaseG(OCP_DBL Pin) override;
	OCP_DBL gammaPhaseW(OCP_DBL Pin) override;
	OCP_DBL gammaPhaseO_OW(OCP_DBL Pin);
	OCP_DBL gammaPhaseO_OGW(OCP_DBL Pin, OCP_DBL Pbbin);
	OCP_DBL gammaPhaseOG(OCP_DBL Pin, OCP_DBL Tin, OCP_DBL* Ziin) override { ERRORcheck("should not be used in BLKOIL"); exit(0); };

private:

	int									Mode;
	ReservoirTable<OCP_DBL>				PVCO;
	ReservoirTable<OCP_DBL>				PVDG;
	ReservoirTable<OCP_DBL>				PVTW;
    ReservoirTable<OCP_DBL>              PVDO;

	// Auxiliary parameters for Table interpolation
	int									len{ 0 };
	vector<OCP_DBL>						data;
	vector<OCP_DBL>						cdata;


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR
	OCP_DBL								Std_RhoO, Std_GammaO;
	OCP_DBL								Std_RhoG, Std_GammaG;
	OCP_DBL								Std_RhoW, Std_GammaW;

};
