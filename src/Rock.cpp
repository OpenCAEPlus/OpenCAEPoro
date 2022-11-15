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
// Rock_Linear01
///////////////////////////////////////////////


void Rock_Linear01::CalPoro(const OCP_DBL& poroInit, const OCP_DBL& P, OCP_DBL& poro, OCP_DBL& dPorodP) const
{
	poro = poroInit * (1 + Rc1 * (P - Pref));
	dPorodP = poroInit * Rc1;
}



///////////////////////////////////////////////
// Rock_Linear02
///////////////////////////////////////////////


void Rock_Linear02::CalPoro(const OCP_DBL& poroInit, const OCP_DBL& P, OCP_DBL& poro, OCP_DBL& dPorodP) const
{
	poro = poroInit * (1 + Rc1 * (P - Pref) + Rc2 / 2 * (P - Pref) * (P - Pref));
	dPorodP = poroInit * (Rc1 + Rc2 * (P - Pref));
}




/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
