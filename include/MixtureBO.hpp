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
    void BOMixtureInit(const ParamReservoir& rs_param);

    // return gamma
    virtual OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override { OCP_ABORT("Should not be used here!"); return 0; };
    virtual OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override { OCP_ABORT("Should not be used here!"); return 0; };
    virtual OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override { OCP_ABORT("Should not be used here!"); return 0; };
    OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
                         const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Should not be used in Black Oil mode!");  return 0;
    };

    // usless in BLKOIL
    USI GetFtype() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 100; }
    OCP_SIN GetMinEigenSkip() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0; }
    bool GetFlagSkip() override { OCP_ABORT("Should not be used in Black Oil mode!");  return false; }
    OCP_DBL GetSurTen() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0; }
    OCP_ULL GetSSMSTAiters() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0;
    }
    OCP_ULL GetNRSTAiters() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0;
    }
    OCP_ULL GetSSMSPiters() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0;
    }
    OCP_ULL GetNRSPiters() override { OCP_ABORT("Should not be used in Black Oil mode!"); return 0;
    }

protected:
    // USI mixtureType; ///< indicates the type of mixture, black oil or compositional or
                     ///< others.
    // std_Gamma* = std_Rho* * GRAVITY_FACTOR.
    // only one of rho and gamma is needed, the other will be calculated from it.
    OCP_DBL std_RhoO;   ///< mass density of oil phase in standard condition.
    OCP_DBL std_GammaO; ///< std_RhoO * gravity factor.
    OCP_DBL std_RhoG;   ///< mass density of gas phase in standard condition.
    OCP_DBL std_GammaG; ///< std_RhoG * gravity factor.
    OCP_DBL std_RhoW;   ///< mass density of water phase in standard condition.
    OCP_DBL std_GammaW; ///< std_RhoW * gravity factor.
};

///////////////////////////////////////////////
// BOMixture_W
///////////////////////////////////////////////

class BOMixture_W : public BOMixture
{
public:

    BOMixture_W() = default;
    BOMixture_W(const ParamReservoir& rs_param, const USI& i) { OCP_ABORT("Not Completed!"); };

    void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
        const OCP_DBL* Sjin, const OCP_DBL& Vpore,
        const OCP_DBL* Ziin) override {
        OCP_ABORT("Not Completed!");
    };
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override {
        OCP_ABORT("Not Completed!");
    };
    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override {
        OCP_ABORT("Not Completed!");
    };
    void FlashDeriv_n(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const OCP_DBL* Sjin, const OCP_DBL* xijin,
        const OCP_DBL* njin, const USI& ftype, const USI* phaseExistin, 
        const USI& lastNP, const OCP_DBL* lastKs) override {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override { OCP_ABORT("Not Completed!"); return 0; };
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override { OCP_ABORT("Not Completed!"); return 0; };
    OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override { OCP_ABORT("Not Completed!"); return 0; };

private:
    OCPTable PVTW;
};

///////////////////////////////////////////////
// BOMixture_OW
///////////////////////////////////////////////

class BOMixture_OW : public BOMixture
{
public:
    BOMixture_OW() = default;
    BOMixture_OW(const ParamReservoir& rs_param, const USI& i);

    void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
        const OCP_DBL* Sjin, const OCP_DBL& Vpore,
        const OCP_DBL* Ziin) override;
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
    void FlashDeriv_n(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const OCP_DBL* Sjin, const OCP_DBL* xijin,
        const OCP_DBL* njin, const USI& ftype, const USI* phaseExistin, 
        const USI& lastNP, const OCP_DBL* lastKs) override {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
    OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override;
    OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override;

private:
    OCPTable PVDO; ///< PVT table for dead oil
    OCPTable PVTW; ///< PVT table for water.
    vector<OCP_DBL> data;   ///< container used to store the results of values of
                            ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata;  ///< container used to store the results of slopes of
                            ///< interpolation of PVT tables.

};

///////////////////////////////////////////////
// BOMixture_ODGW
///////////////////////////////////////////////

class BOMixture_ODGW : public BOMixture
{
public:
    BOMixture_ODGW() = default;
    BOMixture_ODGW(const ParamReservoir& rs_param, const USI& i);

    void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
        const OCP_DBL* Sjin, const OCP_DBL& Vpore,
        const OCP_DBL* Ziin) override;
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
    void FlashDeriv_n(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const OCP_DBL* Sjin, const OCP_DBL* xijin,
        const OCP_DBL* njin, const USI& ftype, const USI* phaseExistin, 
        const USI& lastNP, const OCP_DBL* lastKs) override {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
    OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override;
    OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override;
    OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override;

private:
    OCPTable PVCO; ///< PVT table for live oil (with dissolved gas).
    OCPTable PVDG; ///< PVT table for dry gas.
    OCPTable PVTW; ///< PVT table for water.
    vector<OCP_DBL> data;   ///< container used to store the results of values of
                            ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata;  ///< container used to store the results of slopes of
                            ///< interpolation of PVT tables.
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