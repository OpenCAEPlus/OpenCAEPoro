/*! \file    Rock.hpp
 *  \brief   Rock class declaration
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCK_HEADER__
#define __ROCK_HEADER__

#include <math.h>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "OCPTable.hpp"
#include "ParamReservoir.hpp"

class Rock
{
public:
    Rock() = default;

    // Calculate porosity and d porosity / d P for isothermal model
    virtual void CalPoro(const OCP_DBL& P,
                         const OCP_DBL& poroInit,
                         OCP_DBL&       poro,
                         OCP_DBL&       dPorodP) const = 0;

    // Calculate porosity and d porosity / d P for ifThermal model
    virtual void CalPoroT(const OCP_DBL& P,
                          const OCP_DBL& T,
                          const OCP_DBL& poroInit,
                          OCP_DBL&       poro,
                          OCP_DBL&       dPorodP,
                          OCP_DBL&       dPorodT,
                          OCP_DBL&       RockV,
                          OCP_DBL&       dRockVdP,
                          OCP_DBL&       dRockVdT) const = 0;

    virtual void CalRockHT(const OCP_DBL& T, OCP_DBL& Hr, OCP_DBL& dHrdT) const = 0;
};

class RockIso : public Rock
{
public:
    RockIso() = default;
    void CalPoroT(const OCP_DBL& P,
                  const OCP_DBL& T,
                  const OCP_DBL& poroInit,
                  OCP_DBL&       poro,
                  OCP_DBL&       dPorodP,
                  OCP_DBL&       dPorodT,
                  OCP_DBL&       RockV,
                  OCP_DBL&       dRockVdP,
                  OCP_DBL&       dRockVdT) const override
    {
        OCP_ABORT("Not Used!");
    };
    void CalRockHT(const OCP_DBL& T, OCP_DBL& Hr, OCP_DBL& dHrdT) const override
    {
        OCP_ABORT("Not Used!");
    };
};

class Rock_Linear : public RockIso
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp1 * (P - Pref) + cp2 / 2 * (P - Pref) * (P - Pref))
public:
    Rock_Linear() = default;
    Rock_Linear(const RockParam& param)
    {
        Pref = param.Pref;
        cp1  = param.cp1;
        cp2  = param.cp2;
    };
    void CalPoro(const OCP_DBL& P,
                 const OCP_DBL& poroInit,
                 OCP_DBL&       poro,
                 OCP_DBL&       dPorodP) const override;

protected:
    OCP_DBL Pref;
    OCP_DBL cp1;
    OCP_DBL cp2;
};

class RockT : public Rock
{
public:
    RockT() = default;
    void Assign(const RockParam& param)
    {
        Pref      = param.Pref;
        Tref      = param.Tref;
        cp        = param.cp1;
        ct        = param.ct;
        cpt       = param.cpt;
        ConstRock = param.ConstRock;
        hcp1      = param.HCP1;
        hcp2      = param.HCP2;
    }
    void CalPoro(const OCP_DBL& P,
                 const OCP_DBL& poroInit,
                 OCP_DBL&       poro,
                 OCP_DBL&       dPorodP) const override
    {
        OCP_ABORT("Not Used!");
    };
    void CalRockHT(const OCP_DBL& T, OCP_DBL& Hr, OCP_DBL& dHrdT) const override;

protected:
    OCP_DBL  Pref;      ///< Reference pressure at initial porosity.
    OCP_DBL  Tref;      ///< Reference temperature at initial porosity.
    OCP_DBL  cp;        ///< Compressibility factor of rock in reservoir.
    OCP_DBL  ct;        ///< Expansion factor of rock in reservoir, ifThermal only
    OCP_DBL  cpt;       ///< cross items, ifThermal only
    OCP_BOOL ConstRock; ///< if true, rock volume remains const, else, bulk volume
                        ///< remains const
    OCP_DBL hcp1;       ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
    OCP_DBL hcp2;       ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
};

class RockT_Linear : public RockT
{
    // poro = poroInit * (1 + phi)
    // poro = poroInit * (1 + cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    RockT_Linear() = default;
    RockT_Linear(const RockParam& param) { Assign(param); };
    void CalPoroT(const OCP_DBL& P,
                  const OCP_DBL& T,
                  const OCP_DBL& poroInit,
                  OCP_DBL&       poro,
                  OCP_DBL&       dPorodP,
                  OCP_DBL&       dPorodT,
                  OCP_DBL&       RockV,
                  OCP_DBL&       dRockVdP,
                  OCP_DBL&       dRockVdT) const override;
};

class RockT_Exp : public RockT
{
    // poro = poroInit * (1 + exp(phi))
    // poro = poroInit * exp(cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
    RockT_Exp() = default;
    RockT_Exp(const RockParam& param) { Assign(param); };
    void CalPoroT(const OCP_DBL& P,
                  const OCP_DBL& T,
                  const OCP_DBL& poroInit,
                  OCP_DBL&       poro,
                  OCP_DBL&       dPorodP,
                  OCP_DBL&       dPorodT,
                  OCP_DBL&       RockV,
                  OCP_DBL&       dRockVdP,
                  OCP_DBL&       dRockVdT) const override;
};

#endif /* end if __ROCK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
