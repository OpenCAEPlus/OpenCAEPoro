/*! \file    WellOpt.cpp
 *  \brief   WellOpt class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoro header files
#include "WellOpt.hpp"

WellOpt::WellOpt(const WellOptParam& optParam)
{
    if (optParam.type == "INJ") {
        type = INJ;
    } else if (optParam.type == "PROD") {
        type = PROD;
    } else {
        OCP_ABORT("Wrong well type!");
    }

    if (type == INJ) {
        fluidType = optParam.fluidType;
        if (fluidType == "WAT" || fluidType == "WATER") {
            fluidType = "WAT";
        }
    }

    if (optParam.state == "OPEN") {
        state = OPEN;
    } else if (optParam.state == "CLOSE") {
        state = CLOSE;
    } else {
        OCP_ABORT("Wrong state type!");
    }

    if (optParam.optMode == "RATE") {
        optMode = RATE_MODE;
    } else if (optParam.optMode == "ORAT") {
        optMode = ORATE_MODE;
    } else if (optParam.optMode == "GRAT") {
        optMode = GRATE_MODE;
    } else if (optParam.optMode == "WRAT") {
        optMode = WRATE_MODE;
    } else if (optParam.optMode == "LRAT") {
        optMode = LRATE_MODE;
    } else if (optParam.optMode == "BHP") {
        optMode = BHP_MODE;
    } else {
        OCP_ABORT("Wrong well option mode!");
    }

    initOptMode = optMode;
    maxRate     = optParam.maxRate;
    maxBHP      = optParam.maxBHP;
    minBHP      = optParam.minBHP;
    injTemp     = optParam.injTemp;
}

OCP_BOOL WellOpt::operator!=(const WellOpt& opt) const
{
    if (this->type != opt.type) return OCP_TRUE;
    if (this->state != opt.state) return OCP_TRUE;
    if (this->optMode != opt.optMode) return OCP_TRUE;
    if (this->initOptMode != opt.initOptMode) return OCP_TRUE;
    if (fabs(this->maxRate - opt.maxRate) > TINY) return OCP_TRUE;
    if (fabs(this->maxBHP - opt.maxBHP) > TINY) return OCP_TRUE;
    if (fabs(this->minBHP - opt.minBHP) > TINY) return OCP_TRUE;
    for (USI i = 0; i < injZi.size(); i++) {
        if (fabs(injZi[i] - opt.injZi[i]) > TINY) return OCP_TRUE;
    }
    for (USI i = 0; i < this->prodPhaseWeight.size(); i++) {
        if (fabs(this->prodPhaseWeight[i] - opt.prodPhaseWeight[i]) > TINY)
            return OCP_TRUE;
    }
    if (this->injProdPhase != opt.injProdPhase) return OCP_TRUE;
    if (fabs(this->injTemp - opt.injTemp) > TINY) return OCP_TRUE;
    return OCP_FALSE;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
