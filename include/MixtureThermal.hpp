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
    MixtureThermal(const ComponentsParam& param, const USI& i);

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



#endif /* end if __MIXTURETHERMAL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
