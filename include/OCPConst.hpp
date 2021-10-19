/*! \file    OCPConst.hpp
 *  \brief   Definition of build-in datatypes and consts
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
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

// OpenCAEPoro header files
#include "UtilError.hpp"

// Check and reprot
#define ERRORcheck(exp)                                                                \
    std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in "           \
              << __FILE__ << std::endl;

// Build-in data type
typedef unsigned int USI;     ///< Unsigned integer
typedef unsigned int OCP_USI; ///< Long unsigned integer
typedef int          OCP_INT; ///< Long integer
typedef double       OCP_DBL; ///< Double precision
typedef float        OCP_SIN; ///< Single precision

// Solution methods
const USI IMPES = 1;
const USI FIM   = 2;

// General consts
const OCP_DBL MAX_TIME_STEP = 365;         ///< Maximal time stepsize
const OCP_DBL MIN_TIME_STEP = 1E-6;        ///< Minimal time stepsize
const OCP_DBL TINY          = 1E-8;        ///< Small constant
const OCP_DBL PI            = 3.141592653; ///< Pi

// Physical consts
const OCP_DBL GRAVITY_FACTOR = 0.00694444; ///< 0.00694444 ft2 psi / lb
const OCP_DBL RHOW_STD       = 62.3664;    // lb / ft3
const OCP_DBL RHOAIR_STD     = 0.076362;   // lb / ft3
const OCP_DBL PRESSURE_STD   = 14.6959;    ///< 14.6959 psia = 1 atm

// Unit conversion consts
const OCP_DBL CONV1 = 5.61458;    ///< 1 bbl = 5.61458 ft3
const OCP_DBL CONV2 = 1.12712E-3; ///< Darcy constant in Field

// Mixture types
const USI BLKOIL   = 1;
const USI EoS_PVTW = 2;

// Phase types
const USI PHASE_W   = 1;
const USI PHASE_GW  = 2;
const USI PHASE_OW  = 3;
const USI PHASE_OGW = 4;
const USI PHASE_OG  = 5;

// Well params
// perforation direction
const USI X_DIRECTION = 1;
const USI Y_DIRECTION = 2;
const USI Z_DIRECTION = 3;

const USI  INJ        = 1;
const USI  PROD       = 2;
const bool CLOSE      = false;
const bool OPEN       = true;
const USI  HORIZONTAL = 1;
const USI  VERTICAL   = 2;

// Well option params
const USI RATE_MODE  = 1;
const USI ORATE_MODE = 2;
const USI GRATE_MODE = 3;
const USI WRATE_MODE = 4;
const USI LRATE_MODE = 5;
const USI BHP_MODE   = 6;

// Fluid types
const USI OIL     = 0;
const USI GAS     = 1;
const USI WATER   = 2;
const USI SOLVENT = 3;

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/