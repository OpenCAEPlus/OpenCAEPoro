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


WellOpt::WellOpt(const WellOptParam& Optparam)
{
    if (Optparam.type == "INJ") {
        type = INJ;
    }
    else if (Optparam.type == "PROD") {
        type = PROD;
    }
    else {
        OCP_ABORT("Wrong well type!");
    }

    if (type == INJ) {
        fluidType = Optparam.fluidType;
        if (fluidType == "WAT" || fluidType == "WATER") {
            fluidType = "WAT";
        }
    }

    if (Optparam.state == "OPEN") {
        state = OPEN;
    }
    else if (Optparam.state == "CLOSE") {
        state = CLOSE;
    }
    else {
        OCP_ABORT("Wrong state type!");
    }

    if (Optparam.optMode == "RATE") {
        optMode = RATE_MODE;
    }
    else if (Optparam.optMode == "ORAT") {
        optMode = ORATE_MODE;
    }
    else if (Optparam.optMode == "GRAT") {
        optMode = GRATE_MODE;
    }
    else if (Optparam.optMode == "WRAT") {
        optMode = WRATE_MODE;
    }
    else if (Optparam.optMode == "LRAT") {
        optMode = LRATE_MODE;
    }
    else if (Optparam.optMode == "BHP") {
        optMode = BHP_MODE;
    }
    else {
        OCP_ABORT("Wrong well option mode!");
    }

    initOptMode = optMode;
    maxRate = Optparam.maxRate;
    maxBHP = Optparam.maxBHP;
    minBHP = Optparam.minBHP;
}

OCP_BOOL WellOpt::operator!=(const WellOpt& Opt) const
{
    if (this->type != Opt.type) return OCP_TRUE;
    if (this->state != Opt.state) return OCP_TRUE;
    if (this->optMode != Opt.optMode) return OCP_TRUE;
    if (this->initOptMode != Opt.initOptMode) return OCP_TRUE;
    if (fabs(this->maxRate - Opt.maxRate) > TINY) return OCP_TRUE;
    if (fabs(this->maxBHP - Opt.maxBHP) > TINY) return OCP_TRUE;
    if (fabs(this->minBHP - Opt.minBHP) > TINY) return OCP_TRUE;
    for (USI i = 0; i < injZi.size(); i++) {
        if (fabs(injZi[i] - Opt.injZi[i]) > TINY) return OCP_TRUE;
    }
    for (USI i = 0; i < this->prodPhaseWeight.size(); i++) {
        if (fabs(this->prodPhaseWeight[i] - Opt.prodPhaseWeight[i]) > TINY) return OCP_TRUE;
    }
    if (this->injProdPhase != Opt.injProdPhase) return OCP_TRUE;
    if (fabs(this->Tinj - Opt.Tinj) > TINY) return OCP_TRUE;
    return OCP_FALSE;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
