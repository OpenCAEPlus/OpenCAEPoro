/*! \file    MixtureBO.cpp
 *  \brief   MixtureBO class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureBO.hpp"

void BOMixture::BOMixtureInit(const ParamReservoir& rs_param)
{
    Allocate();

    if (rs_param.density.activity) {
        std_RhoO = rs_param.density.data[0];
        std_RhoW = rs_param.density.data[1];
        std_RhoG = rs_param.density.data[2];
    } else {
        std_RhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
        std_RhoW = RHOW_STD * rs_param.gravity.data[1];
        std_RhoG = RHOAIR_STD * rs_param.gravity.data[2];
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/