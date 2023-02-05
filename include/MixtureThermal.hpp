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

/// MixtureThermal is inherited class of Mixture, it's used for ifThermal model.
/// K-value Model
class MixtureThermal : public Mixture
{
public:
    MixtureThermal() = default;
    void Allocate()
    {
        Ni.resize(numCom);
        phaseExist.resize(numPhase);
        S.resize(numPhase);
        vj.resize(numPhase);
        nj.resize(numPhase);
        xij.resize(numPhase * numCom);
        rho.resize(numPhase);
        xi.resize(numPhase);
        mu.resize(numPhase);
        vfi.resize(numCom);
        rhoP.resize(numPhase);
        rhoT.resize(numPhase);
        rhox.resize(numPhase * numCom);
        xiP.resize(numPhase);
        xiT.resize(numPhase);
        xix.resize(numPhase * numCom);
        muP.resize(numPhase);
        muT.resize(numPhase);
        mux.resize(numPhase * numCom);
        dXsdXp.resize((numCom + 1) * numPhase * (numCom + 2));
        Ufi.resize(numCom);
        H.resize(numPhase);
        HT.resize(numPhase);
        Hx.resize(numPhase * numCom);
    }
    void SetupOptionalFeatures(OptionalFeatures& optFeatures,
                               const OCP_USI&    numBulk) override{};
    // usless in Thermal model
    OCP_DBL GetErrorPEC() override
    {
        OCP_ABORT("Should not be used in Thermal mode!");
        return 0;
    }
    void OutMixtureIters() const override{};
};

class MixtureThermal_K01 : public MixtureThermal
{
public:
    MixtureThermal_K01() = default;
    MixtureThermal_K01(const ParamReservoir& param, const USI& tarId);
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
    /// flash calculation with saturation of phases.
    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override;
    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override;
    void InitFlashFIMn(const OCP_DBL& Pin,
                       const OCP_DBL& Pbbin,
                       const OCP_DBL& Tin,
                       const OCP_DBL* Sjin,
                       const OCP_DBL& Vpore,
                       const OCP_DBL* Ziin,
                       const OCP_USI& bId) override{};
    /// Flash calculation with moles of components.
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
    /// Flash calculation with moles of components and Calculate the derivative
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override;
    void FlashFIMn(const OCP_DBL& Pin,
                   const OCP_DBL& Tin,
                   const OCP_DBL* Niin,
                   const OCP_DBL* Sjin,
                   const OCP_DBL* xijin,
                   const OCP_DBL* njin,
                   const USI*     phaseExistin,
                   const USI&     lastNP,
                   const OCP_USI& bId) override{};
    /// Return molar density of phase, it's used to calculate the molar density of
    /// injection fluids in injection wells.
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Ziin,
                    const USI&     tarPhase) override;

    /// return mass density of phase.
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const OCP_DBL* Ziin,
                     const USI&     tarPhase) override;

    // for Well
    void CalProdWeight(const OCP_DBL&         Pin,
                       const OCP_DBL&         Tin,
                       const OCP_DBL*         Niin,
                       const vector<OCP_DBL>& prodPhase,
                       vector<OCP_DBL>&       prodWeight) override;

    void CalProdRate(const OCP_DBL&   Pin,
                     const OCP_DBL&   Tin,
                     const OCP_DBL*   Niin,
                     vector<OCP_DBL>& prodRate) override;

    void    SetupWellOpt(WellOpt&                  wellopt,
                         const vector<SolventINJ>& sols,
                         const OCP_DBL&            Psurf,
                         const OCP_DBL&            Tsurf) override;
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override;

protected:
    void CalEnthalpy();

protected:
    OCP_DBL Pref{PRESSURE_STD};    ///< Reference pressure
    OCP_DBL Tref{TEMPERATURE_STD}; ///< Reference temperature

    vector<OCP_DBL> MWc;    ///< Molecular Weight of components
    vector<OCP_DBL> MWp;    ///< Molecular Weight of phase
    vector<OCP_DBL> Tcrit;  ///< Critical temperature of hydrocarbon components
    vector<OCP_DBL> xi_ref; ///< Component molar density at reference temperature and
                            ///< reference pressure, lb/ft3
    vector<OCP_DBL> cp;     ///< Component compressibility, 1/psi
    vector<OCP_DBL> ct1;    ///< The first thermal expansion coefficient, 1/F
    vector<OCP_DBL> ct2;    ///< The second thermal expansion coefficient, 1/F
    vector<OCP_DBL> cpt; ///< The coefficient of density dependence on temperature and
                         ///< pressure, 1/psi-F
    OCP_BOOL liquid_based{OCP_TRUE}; ///< Calculated enthalpy of fluid based liquid
    OCP_BOOL gas_based{OCP_FALSE};   ///< Calculated enthalpy of fluid based gas
    OCP_BOOL simple_hvap{OCP_FALSE}; ///< Calculated enthalpy of fluid with simple
                                     ///< enthalpy of evaporation
    vector<OCP_DBL> cpl1;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F
    vector<OCP_DBL> cpl2;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^2
    vector<OCP_DBL> cpl3;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^3
    vector<OCP_DBL> cpl4;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^4
    vector<OCP_DBL> cpg1;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F
    vector<OCP_DBL> cpg2;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^2
    vector<OCP_DBL> cpg3;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^3
    vector<OCP_DBL> cpg4;            ///< Coefficients in the component liquid enthalpy
                                     ///< calculations, Btu/lbmol/F^4
    vector<OCP_DBL>
        hvapr; ///< Coefficients in the component gas enthalpy calculations, Btu/lbmol
    vector<OCP_DBL> hvr; ///< Coefficients in the vaporization enthalpy calculations
    vector<OCP_DBL> ev;  ///< Coefficients in the vaporization enthalpy calculations
    vector<OCP_DBL>
        avisc; ///< Coefficients in water and oil viscosity correlation formulae
    vector<OCP_DBL>
        bvisc; ///< Coefficients in water and oil viscosity correlation formulae
    vector<OCP_DBL> avg; ///< Coefficients Ak in gas viscosity correlation formulae
    vector<OCP_DBL> bvg; ///< Coefficients Bk in gas viscosity correlation formulae
    OCP_BOOL        useViscTab; ///< If use viscosity to calculate viscosity
    OCPTable        visc;       ///< Viscosity-versus-temperature dependence
    vector<OCP_DBL> data;  ///< Auxiliary variables used to calculate viscosity with
                           ///< interpolating tab
    vector<OCP_DBL> cdata; ///< Auxiliary variables used to calculate viscosity with
                           ///< interpolating tab
};

#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
