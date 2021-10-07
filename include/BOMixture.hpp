/*! \file    BoMixture.hpp
 *  \brief   BOMixture class declaration
 *  \author  Shizhe Li
 *  \date    Oct/07/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOMIXTURE_HEADER__
#define __BOMIXTURE_HEADER__



#include "Mixture.hpp"
#include "ReservoirTable.hxx"
#include "ParamReservoir.hpp"

/// BOMixture is Inherited class of Mixture, it's used for black oil model.
class BOMixture : public Mixture
{
public:
	BOMixture() = default;
	BOMixture(const ParamReservoir& rs_param, const USI& PVTmode, const USI& i);

	bool empty_PVDG() const override{ return PVDG.isempty(); }

	// Flash
	void Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin, const OCP_DBL* Sjin, const OCP_DBL& Vpore, const OCP_DBL* Ziin) override;
	void BOFlash_Sj_W(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
	void BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
	void BOFlash_Sj_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);

	void Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
	void BOFlash_Ni_W(const OCP_DBL& Pin, const OCP_DBL* Niin);
	void BOFlash_Ni_OW(const OCP_DBL& Pin, const OCP_DBL* Niin);
	void BOFlash_Ni_OGW(const OCP_DBL& Pin, const OCP_DBL* Niin);

	// return Xi  molar density
	OCP_DBL xiPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin) override;
	OCP_DBL xiPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin);

	// return rho
	OCP_DBL rhoPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin) override;
	OCP_DBL rhoPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin);


	// return gamma
	OCP_DBL gammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override;
	OCP_DBL gammaPhaseG(const OCP_DBL& Pin) override;
	OCP_DBL gammaPhaseW(const OCP_DBL& Pin) override;
	OCP_DBL gammaPhaseO_OW(const OCP_DBL& Pin);
	OCP_DBL gammaPhaseO_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin);
	OCP_DBL gammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override { ERRORcheck("should not be used in BLKOIL"); exit(0); };

private:
	/// indicates the case of black oil, it's decided by user input.
	/// for example, PHASE_OW implies that only water phase and oil phase could be existing,
	/// which will determine which PVT tables will be used.
	USI									Mode;	
	ReservoirTable<OCP_DBL>				PVCO;	///< 
	ReservoirTable<OCP_DBL>				PVDG;
	ReservoirTable<OCP_DBL>				PVTW;
    ReservoirTable<OCP_DBL>             PVDO;

	// Auxiliary parameters for Table interpolation
	USI									len{ 0 };
	vector<OCP_DBL>						data;
	vector<OCP_DBL>						cdata;


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR
	OCP_DBL								Std_RhoO, Std_GammaO;
	OCP_DBL								Std_RhoG, Std_GammaG;
	OCP_DBL								Std_RhoW, Std_GammaW;

};

#endif