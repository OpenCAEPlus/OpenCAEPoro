/*! \file    MixtureThermal.hpp
 *  \brief   MixtureThermal class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURETHERMAL_HEADER__
#define __MIXTURETHERMAL_HEADER__

#include <cmath>

 // OpenCAEPoro header files
#include "Mixture.hpp"
#include "OCPTable.hpp"

/// MixtureThermal is inherited class of Mixture, it's used for thermal model.
/// K-value Model
class MixtureThermal : public Mixture
{
public:
    MixtureThermal() = default;

    // usless in Thermal model
    USI GetFtype() override { OCP_ABORT("Should not be used in Thermal mode!"); return 100; }
    OCP_SIN GetMinEigenSkip() override { OCP_ABORT("Should not be used in Thermal mode!"); return 0; }
    OCP_BOOL GetFlagSkip() override { OCP_ABORT("Should not be used in Thermal mode!");  return OCP_FALSE; }
    OCP_DBL GetSurTen() override { OCP_ABORT("Should not be used in Thermal mode!"); return 0; }
    OCP_DBL GetErrorPEC() override { OCP_ABORT("Should not be used in Thermal mode!"); return 0; }
    OCP_ULL GetSSMSTAiters() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetNRSTAiters() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetSSMSPiters() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetNRSPiters() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetRRiters() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetSSMSTAcounts() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetNRSTAcounts() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetSSMSPcounts() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetNRSPcounts() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }
    OCP_ULL GetRRcounts() override {
        OCP_ABORT("Should not be used in Thermal mode!"); return 0;
    }

protected:


};


class MixtureThermal_K01 : public MixtureThermal
{
public:

    MixtureThermal_K01() = default;
    MixtureThermal_K01(const ParamReservoir& param, const USI& tarId);

    /// flash calculation with saturation of phases.
    void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
        const OCP_DBL* Sjin, const OCP_DBL& Vpore,
        const OCP_DBL* Ziin) override {};
    void InitFlashDer(const OCP_DBL& Pin, const OCP_DBL& Pbbin,
        const OCP_DBL& Tin, const OCP_DBL* Sjin,
        const OCP_DBL& Vpore, const OCP_DBL* Ziin) override {};
    void InitFlashDer_n(const OCP_DBL& Pin, const OCP_DBL& Pbbin,
        const OCP_DBL& Tin, const OCP_DBL* Sjin,
        const OCP_DBL& Vpore, const OCP_DBL* Ziin) override {};
    /// Flash calculation with moles of components.
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override {};
    /// Flash calculation with moles of components and Calculate the derivative
    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override {};
    void FlashDeriv_n(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const OCP_DBL* Sjin, const OCP_DBL* xijin,
        const OCP_DBL* njin, const USI& ftype, const USI* phaseExistin,
        const USI& lastNP, const OCP_DBL* lastKs) override {};
    /// Return molar density of phase, it's used to calculate the molar density of
    /// injection fluids in injection wells.
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Ziin, const USI& tarPhase) override {};

    /// return mass density of phase.
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Pbb, const OCP_DBL& Tin,
        const OCP_DBL* Ziin, const USI& tarPhase) override {};

    void CalProdWeight(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin,
        const vector<OCP_BOOL>& prodPhase, vector<OCP_DBL>& prodWeight) override {};

protected:

    OCP_DBL Pref{ PRESSURE_STD };   ///< reference pressure
    OCP_DBL Tref{ TEMPERATURE_STD };   ///< reference temperature

    vector<OCP_DBL>   xi_ref; ///< component molar density at reference temperature and reference pressure, lb/ft3
    vector<OCP_DBL>   cp;     ///< component compressibility, 1/psi
    vector<OCP_DBL>   ct1;    ///< the first thermal expansion coefficient, 1/F
    vector<OCP_DBL>   ct2;    ///< the second thermal expansion coefficient, 1/F
    vector<OCP_DBL>   cpt;    ///< the coefficient of density dependence on temperature and pressure, 1/psi-F
    vector<OCP_DBL>   cpl1;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
    vector<OCP_DBL>   cpl2;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
    vector<OCP_DBL>   cpl3;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
    vector<OCP_DBL>   cpl4;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
    vector<OCP_DBL>   cpg1;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F
    vector<OCP_DBL>   cpg2;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^2
    vector<OCP_DBL>   cpg3;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^3
    vector<OCP_DBL>   cpg4;   ///< coefficients in the component liquid enthalpy calculations, Btu/lbmol/F^4
    vector<OCP_DBL>   hvapr;  ///< coefficients in the component gas enthalpy calculations, Btu/lbmol
    vector<OCP_DBL>   hvr;    ///< coefficients in the vaporization enthalpy calculations
    vector<OCP_DBL>   ev;     ///< coefficients in the vaporization enthalpy calculations
    vector<OCP_DBL>   avisc;  ///< coefficients in water and oil viscosity correlation formulae
    vector<OCP_DBL>   bvisc;  ///< coefficients in water and oil viscosity correlation formulae
    vector<OCP_DBL>   avg;    ///< coefficients Ak in gas viscosity correlation formulae
    vector<OCP_DBL>   bvg;    ///< coefficients Bk in gas viscosity correlation formulae
    OCPTable         visc;    ///< viscosity-versus-temperature dependence
};



#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/