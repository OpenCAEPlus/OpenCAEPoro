/*! \file    WellPerf.hpp
 *  \brief   WellPerf class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PERFORATION_HEADER__
#define __PERFORATION_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

/// Perforation describe the connections between wells and bulks.
class Perforation
{
    friend class Well;
    friend class Out4RPT;

public:
    /// Default constructor.
    Perforation() = default;

    /// Set state of perf
    void SetState(const OCP_BOOL& flag) { state = flag; };
    /// Return the location of perf: index of bulk
    OCP_USI Location() const { return location; }

private:
    USI      I;        ///< I-index of Perforation in grid.
    USI      J;        ///< J-index of Perforation in grid.
    USI      K;        ///< K-index of Perforation in grid.
    OCP_BOOL state;    ///< True: perforation is open. False: perforation is close.
    OCP_USI  location; ///< Index of bulks which connects to current perforation.
    OCP_DBL  depth;    ///< Depth of bulks which connects to current perforation.
    OCP_DBL  P;        ///< Pressure in current perforation.

    OCP_DBL WI;     ///< Connection transmissibility factor, it can be provided directly
                    ///< from the users.
    OCP_DBL radius; ///< Well radius.
    OCP_DBL kh;     ///< Effective permeability times net thickness of the connection.
    OCP_DBL skinFactor; ///< Skin factor.
    USI     direction;  ///< Direction of the well penetrating the grid block

    /// Multiplier factor for transmissibility of current perforation.
    /// It equals to 0 (close) or 1 (open) now.
    OCP_DBL multiplier;
    /// Molar density of fluid in current perforation. It's used in injection well,
    /// where the fluid consists only single phase.
    mutable OCP_DBL xi;
    vector<OCP_DBL> qi_lbmol; ///< Flow rate of moles of components from into/out
                              ///< current perforation.
    vector<OCP_DBL> transj;   ///< Transmissibility of phase in current perforation.
    OCP_DBL         transINJ;
    vector<OCP_DBL> qj_ft3; ///< Flow rate of volume of phase from into/out current
                            ///< perforation.
    OCP_DBL qt_ft3;         ///< Flow rate of volume of fluids from into/out current
                            ///< perforation.
};

#endif /* end if __PERFORATION_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/