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

    // return gamma
    virtual OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override { OCP_ABORT("Should not be used here!"); return 0; };
    virtual OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override { OCP_ABORT("Should not be used here!"); return 0; };
    virtual OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override { OCP_ABORT("Should not be used here!"); return 0; };
    OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Should not be used in Thermal mode!");  return 0;
    };

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
        const OCP_DBL* Ziin) override {};

    /// return mass density of phase.
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Ziin) override {};

    /// return gamma of oil phase, gamma equals to mass density times gravity factor.
    OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override {};
    /// return gamma of water phase, gamma equals to mass density times gravity factor.
    OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override {};
    /// return gamma of gas phase, gamma equals to mass density times gravity factor.
    OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override {};
    /// return gamma of hydrocarbon mixture, gamma equals to mass density times gravity
    /// factor.
    OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Ziin) override {};

protected:

    vector<OCP_DBL>  xi_ref;
    vector<OCP_DBL>  cp;
    vector<OCP_DBL>  ct1;
    vector<OCP_DBL>  ct2;
    vector<OCP_DBL>  cpt;
    vector<OCP_DBL>  cpl1;
    vector<OCP_DBL>  cpl2;
    vector<OCP_DBL>  cpl3;
    vector<OCP_DBL>  cpl4;
    vector<OCP_DBL>  cpg1;
    vector<OCP_DBL>  cpg2;
    vector<OCP_DBL>  cpg3;
    vector<OCP_DBL>  cpg4;
    vector<OCP_DBL>  hvapr;
    vector<OCP_DBL>  hvr;
    vector<OCP_DBL>  ev;
    vector<OCP_DBL>  avisc;
    vector<OCP_DBL>  bvisc;
    vector<OCP_DBL>  avg;
    vector<OCP_DBL>  bvg;
    OCPTable         visc;
};



#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
