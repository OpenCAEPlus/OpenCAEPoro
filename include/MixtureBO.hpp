/*! \file    MixtureBO.hpp
 *  \brief   MixtureBO class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BOMIXTURE_HEADER__
#define __BOMIXTURE_HEADER__

#include <cmath>

// OpenCAEPoro header files
#include "Mixture.hpp"
#include "OCPTable.hpp"

/// BOMixture is inherited class of Mixture, it's used for black oil model.
class BOMixture : public Mixture
{
public:
    BOMixture() = default;
    BOMixture(const ParamReservoir& rs_param, const USI& PVTmode, const USI& i);

    /// judge if table PVDG is empty.
    bool IsEmpty_PVDG() const override { return PVDG.IsEmpty(); }

    // Flash
    /// flash calculation with saturations of phase, pressure and buble point pressure
    /// in bulks. temperature is unnecessary now, Ziin, which represents the proportion
    /// of components is useless.
    void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
                  const OCP_DBL* Sjin, const OCP_DBL& Vpore,
                  const OCP_DBL* Ziin) override;
    /// flash calculation with saturations while PVTmode is PHASE_W, where only water
    /// phase exists.
    void BOFlash_Sj_W(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
    /// flash calculation with saturations while PVTmode is PHASE_OW, where only water
    /// phase and oil phase could exist.
    void BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore);
    /// flash calculation with saturations while PVTmode is PHASE_OGW, where water
    /// phase, oil phase and gas phase could exist. (to do)in fact, the phasecase where
    /// if dissolved gas exists should be distinguished.
    void BOFlash_Sj_ODGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL* Sjin,
                         const OCP_DBL& Vpore);
    /// flash calculation with saturations of phase, pressure in bulks. temperature is
    /// unnecessary now.
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    /// flash calculation with moles of components while PVTmode is PHASE_W, where only
    /// water phase exists.
    void BOFlash_Ni_W(const OCP_DBL& Pin, const OCP_DBL* Niin);
    /// flash calculation with moles of components while PVTmode is PHASE_OW, where only
    /// water phase and oil phase could exist.
    void BOFlash_Ni_OW(const OCP_DBL& Pin, const OCP_DBL* Niin);
    /// flash calculation with moles of components while PVTmode is PHASE_OGW, where
    /// water phase, oil phase and gas phase could exist. (to do)in fact, the phasecase
    /// where if dissolved gas exists should be distinguished.
    void BOFlash_Ni_ODGW(const OCP_DBL& Pin, const OCP_DBL* Niin);

    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
                        const OCP_DBL* Niin) override;
    void BOFlash_Ni_W_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin);
    void BOFlash_Ni_OW_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin);
    void BOFlash_Ni_ODGW_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin);

    // return xi  molar density
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
    OCP_DBL XiPhase_OW(const OCP_DBL& Pin, const OCP_DBL* Ziin);
    OCP_DBL XiPhase_ODGW(const OCP_DBL& Pin, const OCP_DBL* Ziin);

    // return rho
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
                     const OCP_DBL* Ziin) override;
    OCP_DBL RhoPhase_OW(const OCP_DBL& Pin, const OCP_DBL* Ziin);
    OCP_DBL RhoPhase_ODGW(const OCP_DBL& Pin, const OCP_DBL* Ziin);

    // return gamma
    OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override;
    OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override;
    OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override;
    OCP_DBL GammaPhaseO_OW(const OCP_DBL& Pin);
    OCP_DBL GammaPhaseO_ODGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin);
    OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
                         const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Should not be used in Black Oil mode!");
    };

private:
    /// indicates the case of black oil, it's decided by user input.
    /// for example, PHASE_OW implies that only water phase and oil phase could be
    /// existing, which will determine which PVT tables will be used.
    USI      mode;
    OCPTable PVCO; ///< PVT table for live oil (with dissolved gas).
    OCPTable PVDG; ///< PVT table for dry gas.
    OCPTable PVTW; ///< PVT table for water.
    OCPTable PVDO; ///< PVT table for dead oil (without dissolved gas).

    // Auxiliary parameters for Table interpolation
    USI             len{0}; ///< maximum number of columns of tables among all above.
    vector<OCP_DBL> data;   ///< container used to store the results of values of
                            ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata;  ///< container used to store the results of slopes of
                            ///< interpolation of PVT tables.

    // std_Gamma* = std_Rho* * GRAVITY_FACTOR.
    // only one of rho and gamma is needed, the other will be calculated from it.
    OCP_DBL std_RhoO;   ///< mass density of oil phase in standard condition.
    OCP_DBL std_GammaO; ///< std_RhoO * gravity factor.
    OCP_DBL std_RhoG;   ///< mass density of gas phase in standard condition.
    OCP_DBL std_GammaG; ///< std_RhoG * gravity factor.
    OCP_DBL std_RhoW;   ///< mass density of water phase in standard condition.
    OCP_DBL std_GammaW; ///< std_RhoW * gravity factor.
};

#endif /* end if __BOMIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/