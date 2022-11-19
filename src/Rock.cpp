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


void Rock_Linear::CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit, OCP_DBL& poro, OCP_DBL& dPorodP) const
{
	OCP_DBL dP = (P - Pref);
	poro = poroInit * (1 + (cp1  + cp2 / 2 * dP) * dP);
	dPorodP = poroInit * (cp1 + cp2 * dP);
}


///////////////////////////////////////////////
// RockT_Linear
///////////////////////////////////////////////

void RockT_Linear :: CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro,
							OCP_DBL& dPorodP, OCP_DBL& dPorodT, OCP_DBL& dRockVdP, OCP_DBL& dRockVdT) const
{
	poro = poroInit * (1 + (cp * (P - Pref) - ct * (T - Tref) + cpt * (P - Pref) * (T - Tref)));
}


///////////////////////////////////////////////
// RockT_Exp
///////////////////////////////////////////////

void RockT_Exp::CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro,
						OCP_DBL& dPorodP, OCP_DBL& dPorodT, OCP_DBL& dRockVdP, OCP_DBL& dRockVdT) const
{
	OCP_DBL dP = P - Pref;
	OCP_DBL dT = T - Tref;
	poro = poroInit * exp(cp * dP - ct * dT + cpt * dP * dT);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
