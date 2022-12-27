/*! \file    PhasePermeability.hpp
 *  \brief   PhasePermeability class declaration
 *  \author  Shizhe Li
 *  \date    Dec/26/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PHASEPERMEABILITY_HEADER__
#define __PHASEPERMEABILITY_HEADER__



#include "OCPConst.hpp"
#include "ParamReservoir.hpp"


#include<vector>


using namespace std;


/////////////////////////////////////////////////////////////////////
// Miscible For Compositional Model
/////////////////////////////////////////////////////////////////////


class Miscible
{

public:

    void InputParam(const Miscstr& misterm);
    void CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp);
    void ResetTolastTimeStep() { surTen = lsurTen; }
    void UpdateLastTimeStep() { lsurTen = surTen; }
    OCP_DBL GetSurTen(const OCP_USI& n)const { return surTen[n]; }
    OCP_DBL GetFk(const OCP_USI& n)const { return Fk[n]; }
    OCP_DBL GetFp(const OCP_USI& n)const { return Fp[n]; }

protected:
    /// Miscible treatment of hydrocarbons, only used in compositional Model.
    OCP_BOOL    ifMiscible{ OCP_FALSE };  
    /// The reference surface tension - flow is immiscible when the surface tension 
    //  is greater than or equal to this value
    OCP_DBL     surTenRef;
    OCP_DBL     surTenPc;       ///< Maximum surface tension for capillary pressure / surTenRef
    OCP_DBL     Fkexp;          ///< Exponent set used to calculate Fk
    vector<OCP_DBL> surTen;     ///< Surface tensions between hydrocarbon phases.
    vector<OCP_DBL> Fk;         ///< The relative permeability interpolation parameter
    vector<OCP_DBL> Fp;         ///< The capillary pressure interpolation parameter.

    // Last time step
    vector<OCP_DBL> lsurTen;           ///< last surTen.
};





#endif /* end if __PHASEPERMEABILITY_HEADER__ */


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/