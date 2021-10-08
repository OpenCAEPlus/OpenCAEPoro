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


// OpenCAEPoro header files
#include "Mixture.hpp"
#include "ReservoirTable.hxx"
#include "ParamReservoir.hpp"

/// BOMixture is inherited class of Mixture, it's used for black oil model.
class BOMixture : public Mixture
{
public:
	BOMixture() = default;
	BOMixture(const ParamReservoir& rs_param, const USI& PVTmode, const USI& i);

	/// judge if table PVDG is empty.
	bool empty_PVDG() const override{ return PVDG.isempty(); }

	// Flash
	/// flash calculation with saturations of phase, pressure and buble point pressure in bulks.
	/// temperature is unnecessary now, Ziin, which represents the proportion of components is useless.
	void Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin, const OCP_DBL* Sjin, const OCP_DBL& Vpore, const OCP_DBL* Ziin) override;
	/// flash calculation with saturations while PVTmode is PHASE_W, where only water phase exists.
	void BOFlash_Sj_W(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
	/// flash calculation with saturations while PVTmode is PHASE_OW, where only water phase and oil phase could exist.
	void BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
	/// flash calculation with saturations while PVTmode is PHASE_OGW, where water phase, oil phase and gas phase could exist.
	/// (to do)in fact, the phasecase where if dissolved gas exists should be distinguished.
	void BOFlash_Sj_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
	/// flash calculation with saturations of phase, pressure in bulks. temperature is unnecessary now.
	void Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
	/// flash calculation with moles of components while PVTmode is PHASE_W, where only water phase exists.
	void BOFlash_Ni_W(const OCP_DBL& Pin, const OCP_DBL* Niin);
	/// flash calculation with moles of components while PVTmode is PHASE_OW, where only water phase and oil phase could exist.
	void BOFlash_Ni_OW(const OCP_DBL& Pin, const OCP_DBL* Niin);
	/// flash calculation with moles of components while PVTmode is PHASE_OGW, where water phase, oil phase and gas phase could exist.
	/// (to do)in fact, the phasecase where if dissolved gas exists should be distinguished.
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
	ReservoirTable<OCP_DBL>				PVCO;	///< PVT table for live oil (with dissolved gas).
	ReservoirTable<OCP_DBL>				PVDG;	///< PVT table for dry gas.
	ReservoirTable<OCP_DBL>				PVTW;	///< PVT table for water.
    ReservoirTable<OCP_DBL>             PVDO;	///< PVT table for dead oil (without dissolved gas).

	// Auxiliary parameters for Table interpolation
	USI									len{ 0 };	///< maximum number of columns of tables among all above.
	vector<OCP_DBL>						data;		///< container used to store the results of values of interpolation of PVT tables.
	vector<OCP_DBL>						cdata;		///< container used to store the results of slopes of interpolation of PVT tables.


	// Std_Gamma* = Std_Rho* * GRAVITY_FACTOR.
	// only one of rho and gamma is needed, the other will be calculated from it.
	OCP_DBL								Std_RhoO;	///< mass density of oil phase in standard condition.
	OCP_DBL								Std_GammaO; ///< Std_RhoO * gravity factor.
	OCP_DBL								Std_RhoG;	///< mass density of gas phase in standard condition.
	OCP_DBL								Std_GammaG; ///< Std_RhoG * gravity factor.
	OCP_DBL								Std_RhoW;	///< mass density of water phase in standard condition.
	OCP_DBL								Std_GammaW;	///< Std_RhoW * gravity factor.

};

#endif