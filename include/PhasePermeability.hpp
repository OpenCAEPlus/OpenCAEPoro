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

#include "Grid.hpp"
#include "OCPConst.hpp"
#include "ParamReservoir.hpp"

#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////////////
// Miscible For Compositional Model
/////////////////////////////////////////////////////////////////////

class Miscible
{

public:
    /// Input param from input file
    void InputParam(const Miscstr& misterm);
    /// Allocate memory for Miscible term
    void Setup(const OCP_USI& numBulk);
    /// Return ifUseMiscible
    OCP_BOOL IfUseMiscible() const { return ifUseMiscible; }
    /// Assign value to surTen
    void AssignValue(const OCP_USI& n, const OCP_DBL& v) { surTen[n] = v; }
    /// Calculate Fk, Fp and return if miscible
    OCP_BOOL CalFkFp(const OCP_USI& n, OCP_DBL& fk, OCP_DBL& fp);
    /// Reset Miscible term to last time step
    void ResetTolastTimeStep() { surTen = lsurTen; }
    /// Update Miscible term at last time step
    void UpdateLastTimeStep() { lsurTen = surTen; }
    /// Return surTen
    OCP_DBL GetSurTen(const OCP_USI& n) const { return surTen[n]; }
    /// Return Fk
    OCP_DBL GetFk(const OCP_USI& n) const { return Fk[n]; }
    /// Return Fp
    OCP_DBL GetFp(const OCP_USI& n) const { return Fp[n]; }

protected:
    OCP_BOOL ifSetup{OCP_FALSE}; ///< Only one setup is needed.
    /// Miscible treatment of hydrocarbons, only used in compositional Model.
    OCP_BOOL ifUseMiscible{OCP_FALSE};
    /// The reference surface tension - flow is immiscible when the surface tension
    //  is greater than or equal to this value
    OCP_DBL surTenRef;
    OCP_DBL surTenPc; ///< Maximum surface tension for capillary pressure / surTenRef
    OCP_DBL Fkexp;    ///< Exponent set used to calculate Fk
    vector<OCP_DBL> surTen; ///< Surface tensions between hydrocarbon phases.
    vector<OCP_DBL> Fk;     ///< The relative permeability interpolation parameter
    vector<OCP_DBL> Fp;     ///< The capillary pressure interpolation parameter.

    // Last time step
    vector<OCP_DBL> lsurTen; ///< last surTen.
};

/////////////////////////////////////////////////////////////////////
// Scale The Water-Oil Capillary Pressure Curves (From SWATINIT)
/////////////////////////////////////////////////////////////////////

class ScalePcow
{
    // For Output
    friend class Out4RPT;

public:
    /// Setup ScalePcow term
    void Setup(const Grid& myGrid);
    /// Return ifScale
    OCP_BOOL IfScale() const { return ifScale; }
    /// Assign value to scaleVal
    void    AssignScaleValue(const OCP_USI& n, const OCP_DBL& v) { scaleVal[n] = v; }
    OCP_DBL GetSwInit(const OCP_USI& n) const { return swatInit[n]; }
    /// Return scaleVal
    OCP_DBL GetScaleVal(const OCP_USI& n) const { return scaleVal[n]; }

protected:
    OCP_BOOL ifSetup{OCP_FALSE}; ///< Only one setup is needed.
    OCP_BOOL ifScale{OCP_FALSE}; ///< If true, then Scale will be used
    vector<OCP_DBL>
        scaleVal; ///< Scale values for Pcow, it will be calculated from swatInit
    vector<OCP_DBL> swatInit; ///< Initial water distribution
};

#endif /* end if __PHASEPERMEABILITY_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/26/2022      Create file                          */
/*----------------------------------------------------------------------------*/