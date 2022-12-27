/*! \file    PhasePermeability.cpp
 *  \brief   PhasePermeability class declaration
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */


#include "PhasePermeability.hpp"


/////////////////////////////////////////////////////////////////////
// Miscible For Compositional Model
/////////////////////////////////////////////////////////////////////


void Miscible::InputParam(const Miscstr& misterm) {
    const USI len = misterm.surTenRef.size();
    if (len > 0) {
        surTenRef = misterm.surTenRef[0];
        surTenPc = 1;
        if (len > 2) {
            surTenPc = misterm.surTenRef[2] / surTenRef;
        }
        Fkexp = 0.25;
    }
    else {
        OCP_ABORT("No data in MISCSTR!");
    }
}


void Miscible::CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp)
{
    Fk[n] = min(1.0, pow(surTen[n] / surTenRef, Fkexp));
    Fp[n] = min(surTenPc, surTen[n] / surTenRef);
    fk = Fk[n];
    fp = Fp[n];
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/