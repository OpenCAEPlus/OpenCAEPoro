/*! \file    Rock.cpp
 *  \brief   Rock class declaration
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoro header files
#include "Rock.hpp"

///////////////////////////////////////////////
// Rock_Linear
///////////////////////////////////////////////

void Rock_Linear::CalPoro(const OCP_DBL& P,
                          const OCP_DBL& poroInit,
                          OCP_DBL&       poro,
                          OCP_DBL&       dPorodP) const
{
    OCP_DBL dP = (P - Pref);
    poro       = poroInit * (1 + (cp1 + cp2 / 2 * dP) * dP);
    dPorodP    = poroInit * (cp1 + cp2 * dP);
}

///////////////////////////////////////////////
// RockT
///////////////////////////////////////////////

void RockT::CalRockHT(const OCP_DBL& T, OCP_DBL& Hr, OCP_DBL& dHrdT) const
{
    const OCP_DBL Ta  = T + CONV5;
    const OCP_DBL Tra = Tref + CONV5;
    Hr                = hcp1 * (Ta - Tra) + 0.5 * hcp2 * (Ta * Ta - Tra * Tra);
    dHrdT             = hcp1 + hcp2 * Ta;
}

///////////////////////////////////////////////
// RockT_Linear
///////////////////////////////////////////////

void RockT_Linear ::CalPoroT(const OCP_DBL& P,
                             const OCP_DBL& T,
                             const OCP_DBL& poroInit,
                             OCP_DBL&       poro,
                             OCP_DBL&       dPorodP,
                             OCP_DBL&       dPorodT,
                             OCP_DBL&       RockV,
                             OCP_DBL&       dRockVdP,
                             OCP_DBL&       dRockVdT) const
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    poro             = poroInit * (1 + (cp * dP - ct * dT + cpt * dP * dT));
    dPorodP          = poroInit * (cp + cpt * dT);
    dPorodT          = poroInit * (-ct + cpt * dP);

    if (ConstRock) {
        RockV    = 1 - poroInit;
        dRockVdP = 0.0;
        dRockVdT = 0.0;
    } else {
        RockV    = 1 - poro;
        dRockVdP = -dPorodP;
        dRockVdT = -dPorodT;
    }
}

///////////////////////////////////////////////
// RockT_Exp
///////////////////////////////////////////////

void RockT_Exp::CalPoroT(const OCP_DBL& P,
                         const OCP_DBL& T,
                         const OCP_DBL& poroInit,
                         OCP_DBL&       poro,
                         OCP_DBL&       dPorodP,
                         OCP_DBL&       dPorodT,
                         OCP_DBL&       RockV,
                         OCP_DBL&       dRockVdP,
                         OCP_DBL&       dRockVdT) const
{
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    poro             = poroInit * exp(cp * dP - ct * dT + cpt * dP * dT);

    dPorodP = poro * (cp + cpt * dT);
    dPorodT = poro * (-ct + cpt * dP);
    if (ConstRock) {
        RockV    = 1 - poroInit;
        dRockVdP = 0.0;
        dRockVdT = 0.0;
    } else {
        RockV    = 1 - poro;
        dRockVdP = -dPorodP;
        dRockVdT = -dPorodT;
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
