/*! \file    Reservoir.hpp
 *  \brief   Reservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__

// OpenCAEPoro header files
#include "AllWells.hpp"
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "Grid.hpp"
#include "OptionalFeatures.hpp"
#include "ParamRead.hpp"

/// Reservoir is the core component in our simulator, it contains the all reservoir
/// information, and all operations on it.
///
/// Reservoir has four Core components.
/// Grids contains the basic informations of all grids as a database of reservoir.
/// Bulk only stores active grids, which defines the area used for calculation.
/// AllWells contains the well information, it's used to manage operations related to
/// wells. BulkConn contains connections between bulks(active grids).
class Reservoir
{
    friend class OCPControl;
    friend class Summary;
    friend class CriticalInfo;
    friend class Out4RPT;
    friend class Out4VTK;

    // temp
    friend class IsoT_IMPEC;
    friend class IsoT_FIM;
    friend class IsoT_FIMn;
    friend class IsoT_AIMc;
    friend class T_FIM;
    friend class Solver;

    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:
    /// Input param from internal param data structure, which stores the params from
    /// input files.
    void InputParam(ParamRead& param);
    /// Setup static information for reservoir with input params for Isothermal model
    void SetupIsoT();
    /// Setup static information for reservoir with input params for Thermal model
    void SetupT();
    /// Apply the control of ith critical time point.
    void ApplyControl(const USI& i);
    /// Calculate num of Injection, Production
    void CalIPRT(const OCP_DBL& dt);
    /// Calculate Maximum Change of some reference variables for IMPEC
    void CalMaxChange();
    /// Return the num of Bulk
    OCP_USI GetBulkNum() const { return bulk.GetBulkNum(); }
    /// Return the num of Well
    USI GetWellNum() const { return allWells.GetWellNum(); }
    /// Return the num of Components
    USI GetComNum() const { return bulk.GetComNum(); }

protected:
    Grid             grid;        ///< Init Grid info.
    Bulk             bulk;        ///< Active grid info.
    AllWells         allWells;    ///< Wells class info.
    BulkConn         conn;        ///< Bulk's connection info.
    OptionalFeatures optFeatures; ///< optional features.

public:
    /// Calculate the CFL number, including bulks and wells for IMPEC
    OCP_DBL CalCFL(const OCP_DBL& dt) const;
    /// Return NRdPmax
    OCP_DBL GetNRdPmax() { return bulk.GetNRdPmax(); }
    /// Return NRdSmax
    OCP_DBL GetNRdSmax(OCP_USI& index) { return bulk.CalNRdSmax(index); }
    /// Return NRdNmax
    OCP_DBL GetNRdNmax() { return bulk.GetNRdNmax(); }
    void    PrintSolFIM(const string& outfile) const;
    void    OutInfoFinal() const { bulk.OutMixtureIters(); }
};

#endif /* end if __RESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/