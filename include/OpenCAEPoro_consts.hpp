/*! \file    OpenCAEPoro_consts.hpp
 *  \brief   datatype, consts
 *  \author  Shizhe Li
 *  \date    Oct/05/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPENCAEPORO_CONSTS_HEADER__
#define __OPENCAEPORO_CONSTS_HEADER__


#include <iostream>

#define ERRORcheck(exp)                                                                \
    std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in "           \
              << __FILE__ << std::endl;

// data Type
typedef unsigned int USI;
typedef unsigned int OCP_USI;
typedef int          OCP_INT;
typedef double       OCP_DBL;


// Method
const USI IMPES = 0;
const USI FIM   = 1;

// general consts
const OCP_DBL TINY = 1E-8;
const OCP_DBL PI   = 3.141592653;

// pysical consts
const OCP_DBL GRAVITY_FACTOR = 0.00694444; // 0.00694444 ft2 psi / lb
const OCP_DBL RHOW_STD       = 62.3664;    // lb / ft3
const OCP_DBL RHOAIR_STD     = 0.076362;   // lb / ft3
const OCP_DBL PRESSURE_STD   = 14.6959;    // psia   =   1 atm

// Units consts
const OCP_DBL CONV1 = 5.61458;    // 1 bbl = 5.61458 ft3
const OCP_DBL CONV2 = 1.12712E-3; // Darcy constant in Field

// Mixture Type
const USI BLKOIL   = 1;
const USI EoS_PVTW = 2;

// Phase
const USI PHASE_W   = 1;
const USI PHASE_GW  = 2;
const USI PHASE_OW  = 3;
const USI PHASE_OGW = 4;
const USI PHASE_OG  = 5;

// Well params
const USI  INJ        = 0;
const USI  PROD       = 1;
const bool CLOSE      = 0;
const bool OPEN       = 1;
const USI  HORIZONTAL = 0;
const USI  VERTICAL   = 1;
// Well opt param
const USI RATE_MODE  = 0;
const USI ORATE_MODE = 1;
const USI GRATE_MODE = 2;
const USI WRATE_MODE = 3;
const USI LRATE_MODE = 4;
const USI BHP_MODE   = 5;
// Fluid type
const USI OIL     = 0;
const USI GAS     = 1;
const USI WATER   = 2;
const USI SOLVENT = 3;


#endif