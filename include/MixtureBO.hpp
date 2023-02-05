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

#ifndef __MIXTUREBO_HEADER__
#define __MIXTUREBO_HEADER__

#include <cmath>

// OpenCAEPoro header files
#include "Mixture.hpp"
#include "OCPTable.hpp"

/// BOMixture is inherited class of Mixture, it's used for black oil model.
class BOMixture : public Mixture
{
public:
    BOMixture() = default;
    void SetupOptionalFeatures(OptionalFeatures& optFeatures,
                               const OCP_USI&    numBulk) override{};
    void BOMixtureInit(const ParamReservoir& rs_param);
    void InitFlashFIMn(const OCP_DBL& Pin,
                       const OCP_DBL& Pbbin,
                       const OCP_DBL& Tin,
                       const OCP_DBL* Sjin,
                       const OCP_DBL& Vpore,
                       const OCP_DBL* Ziin,
                       const OCP_USI& bId) override
    {
        OCP_ABORT("Not Used!");
    };

    // For Well
    void CalProdWeight(const OCP_DBL&         Pin,
                       const OCP_DBL&         Tin,
                       const OCP_DBL*         Niin,
                       const vector<OCP_DBL>& prodPhase,
                       vector<OCP_DBL>&       prodWeight) override
    {
        prodWeight = prodPhase;
    }

    void CalProdRate(const OCP_DBL&   Pin,
                     const OCP_DBL&   Tin,
                     const OCP_DBL*   Niin,
                     vector<OCP_DBL>& prodRate) override
    {
        prodRate.assign(Niin, Niin + numCom);
    };
    OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) override
    {
        OCP_ABORT("Can not be used in Black Oil Model!");
    }

    OCP_DBL GetErrorPEC() override
    {
        OCP_ABORT("Should not be used in Black Oil mode!");
        return 0;
    }
    void OutMixtureIters() const override{};

protected:
    // USI mixtureType; ///< indicates the type of mixture, black oil or compositional
    // or
    ///< others.
    // std_Gamma* = std_Rho* * GRAVITY_FACTOR.
    // only one of rho and gamma is needed, the other will be calculated from it.
    OCP_DBL std_RhoO; ///< The density of oil at surface conditions : lb/ft3
    OCP_DBL std_RhoG; ///< The density of gas at surface conditions : lb/ft3
    OCP_DBL std_RhoW; ///< The density of water at surface conditions : lb/ft3
};

///////////////////////////////////////////////
// BOMixture_W
///////////////////////////////////////////////

class BOMixture_W : public BOMixture
{
public:
    BOMixture_W() = default;
    BOMixture_W(const ParamReservoir& rs_param, const USI& i)
    {
        OCP_ABORT("Not Completed!");
    };
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override
    {
        OCP_ABORT("Not Completed!");
    }
    void InitFlashIMPEC(const OCP_DBL& Pin,
                        const OCP_DBL& Pbbin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Sjin,
                        const OCP_DBL& Vpore,
                        const OCP_DBL* Ziin,
                        const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void InitFlashFIM(const OCP_DBL& Pin,
                      const OCP_DBL& Pbbin,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Sjin,
                      const OCP_DBL& Vpore,
                      const OCP_DBL* Ziin,
                      const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void FlashFIM(const OCP_DBL& Pin,
                  const OCP_DBL& Tin,
                  const OCP_DBL* Niin,
                  const OCP_DBL* Sjin,
                  const USI&     lastNP,
                  const OCP_DBL* xijin,
                  const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    };
    void FlashFIMn(const OCP_DBL& Pin,
                   const OCP_DBL& Tin,
                   const OCP_DBL* Niin,
                   const OCP_DBL* Sjin,
                   const OCP_DBL* xijin,
                   const OCP_DBL* njin,
                   const USI*     phaseExistin,
                   const USI&     lastNP,
                   const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Ziin,
                    const USI&     tarPhase) override
    {
        OCP_ABORT("Not Completed!");
        return 0;
    };
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const OCP_DBL* Ziin,
                     const USI&     tarPhase) override
    {
        OCP_ABORT("Not Completed!");
        return 0;
    };

    // for Well
    void SetupWellOpt(WellOpt&                  opt,
                      const vector<SolventINJ>& sols,
                      const OCP_DBL&            Psurf,
                      const OCP_DBL&            Tsurf) override
    {
        OCP_ABORT("Not Completed!");
    };

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
    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
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
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
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
                   const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Ziin,
                    const USI&     tarPhase) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const OCP_DBL* Ziin,
                     const USI&     tarPhase) override;

    // for Well
    void SetupWellOpt(WellOpt&                  opt,
                      const vector<SolventINJ>& sols,
                      const OCP_DBL&            Psurf,
                      const OCP_DBL&            Tsurf) override;

private:
    OCPTable        PVDO;  ///< PVT table for dead oil
    OCPTable        PVTW;  ///< PVT table for water.
    vector<OCP_DBL> data;  ///< container used to store the results of values of
                           ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata; ///< container used to store the results of slopes of
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

    void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override;
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
    void FlashIMPEC(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Niin,
                    const USI&     lastNP,
                    const OCP_DBL* xijin,
                    const OCP_USI& bId) override;
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
                   const OCP_USI& bId) override
    {
        OCP_ABORT("Not Completed!");
    }
    OCP_DBL XiPhase(const OCP_DBL& Pin,
                    const OCP_DBL& Tin,
                    const OCP_DBL* Ziin,
                    const USI&     tarPhase) override;
    OCP_DBL RhoPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Pbb,
                     const OCP_DBL& Tin,
                     const OCP_DBL* Ziin,
                     const USI&     tarPhase) override;

    // for Well
    void SetupWellOpt(WellOpt&                  opt,
                      const vector<SolventINJ>& sols,
                      const OCP_DBL&            Psurf,
                      const OCP_DBL&            Tsurf) override;

private:
    OCPTable        PVCO;  ///< PVT table for live oil (with dissolved gas).
    OCPTable        PVDG;  ///< PVT table for dry gas.
    OCPTable        PVTW;  ///< PVT table for water.
    vector<OCP_DBL> data;  ///< container used to store the results of values of
                           ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata; ///< container used to store the results of slopes of
                           ///< interpolation of PVT tables.
};

#endif /* end if __MIXTUREBO_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/