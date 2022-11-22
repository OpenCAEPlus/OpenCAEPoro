/*! \file    WellOpt.hpp
 *  \brief   WellOpt class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLOPT_HEADER__
#define __WELLOPT_HEADER__

// OpenCAEPoro header files
#include "ParamWell.hpp"

/// WellOpt describes the operation mode of a well.
/// usually it changes over time, specifically, each attributes could be changed
/// including the well type.
class WellOpt
{
    friend class Well;
    friend class AllWells;
    friend class Out4RPT;

public:
    /// Default constructor.
    WellOpt() = default;

    /// Constructor well operation mode using params.
    WellOpt(const WellOptParam& Optparam);

    /// overload inequality
    OCP_BOOL operator !=(const WellOpt& Opt) const;

private:
    USI type{ 0 }; ///< type of well, Inj or Prod.
    /// indicate which type of fluids will be injected, water, gas, or other solvent.
    /// it's decided by users and only useful for injection well.
    string fluidType;
    OCP_BOOL   state{ OCP_FALSE }; ///< state of well, close or open.
    USI optMode; ///< control mode of well: constant pressure, or constant flow rate of
                 ///< specified fluids.
    USI initOptMode; ///< Init opt mode during current step
    /// it gives the upper limit of flow rate of specified fluids if the well is under
    /// the control of constant pressure. it gives the flow rate of specified fluids if
    /// the well is under the control of constant flow rate.
    OCP_DBL maxRate;
    /// used for injection well.
    /// it gives the upper limit of well pressure if the well is under the control of
    /// constant flow rate. it gives the pressure of well if the well is under the
    /// control of constant pressure.
    OCP_DBL maxBHP;
    /// used for production well.
    /// it gives the lower limit of well pressure if the well is under the control of
    /// constant flow rate. it gives the pressure of well if the well is under the
    /// control of constant pressure.
    OCP_DBL minBHP;
    /// it's decided by users input.
    /// for injection well, it describes the components of injected fluids.
    /// for production well, it gives the the components of fluids which we are
    /// interested in.
    vector<OCP_DBL> injZi;
    USI injProdPhase;     ///< label the phase of injecting fluid if possible
    vector<OCP_DBL> prodPhaseWeight;
    OCP_DBL factorINJ;            ///< unit factor: Mscf -> lbmol for comps, Mscf -> Mscf in blackoil
    OCP_DBL Tinj;             ///< temperature of inj fluid F
                              // for Reinjection
    OCP_BOOL reInj{ OCP_FALSE }; ///< if OCP_TRUE, reinjection happens
    USI reinjPhase; ///< phase of Reinjection fluid
    OCP_DBL factor; ///< one moles Group production fluid has factor mole reinjection fluid
    vector<USI> connWell; ///< Well which connects to current Well
};


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/22/2022      Create file                          */
/*----------------------------------------------------------------------------*/
