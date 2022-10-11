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
#include "CornerGrid.hpp"
#include "OCPConst.hpp"
#include "ParamReservoir.hpp"
#include "UtilOutput.hpp"

using namespace std;




/// Effective area of intersection surfaces with neighboring cells.
class GPair
{
public:
    GPair() = default;
    GPair(const OCP_USI& Id, const OCP_DBL& Area)
        : id(Id)
        , area(Area){};
    static bool lessG(const GPair& G1, const GPair& G2) { return G1.id < G2.id; }

    OCP_USI id;   ///< Id of a neighboring cell
    OCP_DBL area; ///< Effective area between this cell and the neighbor indicated by id
};

/// Active cell indicator and its index among active cells.
//  Note: GB_Pair contains two variables, which indicates if a grid cell is active or
//  not and its index among active cells if it is active.
class GB_Pair
{
public:
    /// Default constructor.
    GB_Pair() = default;

    /// Constructor with given information. TODO: needed???
    GB_Pair(bool act, OCP_USI i)
        : activity(act)
        , index(i){};

    /// Return whether this cell is active or not.
    bool IsAct() const { return activity; }

    /// Return the index of this cell among active cells.
    OCP_USI GetId() const { return index; }

private:
    bool    activity; ///< Activeness of a grid cell.
    OCP_USI index;    ///< Active index of grid if active
};

/// Basic information of computational grid, including the rock properties.
//  Note: All grid cells are stored here, you can regard it as a database of the
//  reservoir. Considering there exist inactive cells (whose porosity or cell volume is
//  too small) and activeness status might change, the Grid class is necessary. The Grid
//  class is static while simulating, active grids will be stored in bulks, which is
//  "area" for calculating???
class Grid
{
    friend class Bulk;
    friend class BulkConn;
    friend class Well;

public:
    /// Default constructor.
    Grid() = default;
    /// Input parameters from the internal param structure.
    void InputParam(const ParamReservoir& rs_param);
    /// Setup the grid information and calculate the properties.
    void Setup();

    /// Setup an orthogonal grid.
    void SetupOrthogonalGrid();
    /// Setup the neighboring info for an orthogonal grid.
    void SetupNeighborOrthogonalGrid();
    /// Calculate Akd for an orthogonal grid.
    OCP_DBL CalAkdOrthogonalGrid(const OCP_USI& bId, const OCP_USI& eId,
                                const USI& direction);
    /// Calculate the depth and volume for an orthogonal grid.
    void CalDepthVOrthogonalGrid();

    /// Setup a corner-point grid.
    void SetupCornerGrid();
    /// Setup the neighboring info for a corner-point grid.
    void SetupNeighborCornerGrid(const COORD& CoTmp);
    /// Calculate Akd for a corner-point grid.
    OCP_DBL CalAkdCornerGrid(const GeneralConnect& conn);

    /// Calculate the activeness of grid cells.
    void CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Mapping from grid cells to bulks (active cells).
    const GB_Pair& MapG2B(const OCP_USI& i) const { return activeMap_G2B[i]; }
    /// Return nx of grid cell.
    OCP_USI GetGridNx() const { return nx; }
    /// Return ny of grid cell.
    OCP_USI GetGridNy() const { return ny; }
    /// Return nz of grid cell.
    OCP_USI GetGridNz() const { return nz; }
    /// Return the num of grid cells.
    OCP_USI GetGridNum() const { return numGrid; }
    /// Return the num of connections.
    OCP_USI GetConnNum() const { return numConn; }
    /// Return the num of bulks (active cells).
    OCP_USI GetActiveGridNum() const { return activeGridNum; }
    /// Return the index of active cell (i, j, k).
    OCP_USI GetActIndex(const USI& i, const USI& j, const USI& k) const;
    /// Return the 3D coordinate for object grid with Grid index
    void GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const;
    /// Return the 3D coordinate for object grid with bulk(active grids) index
    void GetIJKBulk(USI& i, USI& j, USI& k, const OCP_USI& n) const;
    void CalSomeInfo()const;
private:
    USI     nx;      ///< Number of cells in x-direction
    USI     ny;      ///< Number of cells in y-direction
    USI     nz;      ///< Number of cells in z-direction
    OCP_USI numGrid; ///< Number of all cells = nx * ny * nz
    OCP_USI numConn; ///< Number of connections

    USI                   gridType;  ///< Type of grid.
    vector<vector<GPair>> gNeighbor; ///< Neighboring information of grid.

    // Orthogonal grid
    vector<OCP_DBL> tops;  ///< Depth of top surface of the reservoir: nx*ny
    vector<OCP_DBL> depth; ///< Depth of center of grid cells: numGrid.
    vector<OCP_DBL> dx;    ///< Size of cell in x-direction: numGrid.
    vector<OCP_DBL> dy;    ///< Size of cell in y-direction: numGrid.
    vector<OCP_DBL> dz;    ///< Size of cell in z-direction: numGrid.

    // Corner-point grid
    vector<OCP_DBL> coord; ///< Lines of a corner-point grid.
    vector<OCP_DBL> zcorn; ///< ZValues of a corner-point grid.

    // Rock properties
    vector<OCP_DBL> v;    ///< Volume of cells: numGrid.
    vector<OCP_DBL> ntg;  ///< Net to gross ratio of cells: numGrid
    vector<OCP_DBL> poro; ///< Initial porosity of rock cells: numGrid
    vector<OCP_DBL> kx;   ///< Absolute permeability in x-direction: numGrid
    vector<OCP_DBL> ky;   ///< Absolute permeability in y-direction: numGrid
    vector<OCP_DBL> kz;   ///< Absolute permeability in z-direction: numGrid

    // Initial Properties
    vector<OCP_DBL>   SwatInit; ///< Initial water saturation

    // Region
    vector<USI> SATNUM; ///< Identify SAT region: numGrid.
    vector<USI> PVTNUM; ///< Identify PVT region for the blackoil model: numGrid.
    vector<USI> ACTNUM; ///< Indicate activity of grid from input file: numGrid. 0 = inactive, 1 = active.

    // Active grid cells
    OCP_USI         activeGridNum; ///< Num of active grid.
    vector<OCP_USI> activeMap_B2G; ///< Mapping from active grid to grid: activeGridNum = numBulk
    vector<GB_Pair> activeMap_G2B; ///< Mapping from grid to active grid: numGrid.

private:
    // Auxiliary variable
    USI             numDigutIJK;  ///< number of digits of maximum nx,ny,nz

public:
    void CalNumDigutIJK(); ///< only used in Structural grid
    USI GetNumDigitIJK() const { return numDigutIJK; }
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