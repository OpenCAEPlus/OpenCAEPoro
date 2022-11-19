/*! \file    MixtureThermal_K.cpp
 *  \brief   MixtureThermal_K class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */



// OpenCAEPoro header files
#include "MixtureThermal.hpp"
#include "OCPTable.hpp"



MixtureThermal_K01 :: MixtureThermal_K01(const ParamReservoir& param, const USI& tarId)
{
	numCom = param.numCom; // water is included
	Allocate();

}




/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
