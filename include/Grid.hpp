/*! \file    Grid.hpp
 *  \brief   Grid class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __GRID_HEADER__
#define __GRID_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "GridInitInfo.hpp"
#include "Bulk.hpp"
#include "UtilOutput.hpp"

using namespace std;


/// Basic information of computational grid, including the rock properties.
//  Note: All grid cells are stored here, you can regard it as a database of the
//  reservoir. Considering there exist inactive cells (whose porosity or cell volume is
//  too small) and activeness status might change, the Grid class is necessary. The Grid
//  class is static while simulating, active grids will be stored in bulks, which is
//  "area" for calculating
class Grid
{
    friend class Reservoir;
    friend class BulkConn;
    friend class AllWells;
    friend class Well;
    friend class Out4RPT;
    friend class Out4VTK;
    friend class Summary;
    friend class CriticalInfo;
    friend class OCPControl;

    friend class Solver;
    friend class OCP_IMPEC;
    friend class OCP_FIM;
    friend class OCP_FIMn;

public:
    /// Default constructor.
    Grid() = default;
    /// Input parameters from the internal param structure.
    void InputParam(const ParamReservoir& rs_param);
    /// Setup the grid information and calculate the properties for isothermal model
    void SetupIsoT();
    /// Setup the grid information and calculate the properties for ifThermal model
    void SetupT();

protected:

    GridInitInfo    initInfo;       ///< Initial grid info
    Bulk            bulk;           ///< Fluid grid info

};

#endif /* end if __GRID_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Shizhe Li           Nov/18/2021      Add Connections between Grids        */
/*  Chensong Zhang      Jan/16/2022      Finish Doxygen                       */
/*----------------------------------------------------------------------------*/