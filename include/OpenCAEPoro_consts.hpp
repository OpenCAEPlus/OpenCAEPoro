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

 // Standard header files
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
const USI IMPES = 1;
const USI FIM   = 2;

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
const USI  INJ        = 1;
const USI  PROD       = 2;
const bool CLOSE      = false;
const bool OPEN       = true;
const USI  HORIZONTAL = 1;
const USI  VERTICAL   = 2;
// Well opt param
const USI RATE_MODE  = 1;
const USI ORATE_MODE = 2;
const USI GRATE_MODE = 3;
const USI WRATE_MODE = 4;
const USI LRATE_MODE = 5;
const USI BHP_MODE   = 6;
// Fluid type
const USI OIL     = 1;
const USI GAS     = 2;
const USI WATER   = 3;
const USI SOLVENT = 4;


#endif