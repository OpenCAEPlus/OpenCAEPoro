/*! \file    ThermalSolver.hpp
 *  \brief   ThermalSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __THERMALSOLVER_HEADER__
#define __THERMALSOLVER_HEADER__

 // OpenCAEPoro header files
#include "ThermalMethod.hpp"

/// ThermalSolver class for fluid solution method.
class ThermalSolver
{
public:


private:
    USI           method = FIM;
    LinearSystem  LSolver;
    LinearSystem  auxLSolver;
    OCP_IMPEC     impec;
    OCP_FIM       fim;
};

#endif /* end if __THERMALSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/