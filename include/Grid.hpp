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

using namespace std;

class GPair
{
public:
    GPair() = default;
    GPair(const OCP_USI& Id, const OCP_DBL& Area)
        : id(Id)
        , area(Area){};
    static bool lessG(const GPair& G1, const GPair& G2) { return G1.id < G2.id; }
    OCP_USI     id;   ///< Id of some neighbor
    OCP_DBL     area; ///< Effective area between current grid and the id th neighbor
};

/// GB_Pair contains two variables, which indicates if a grid is active and what its
/// active index is if active.
class GB_Pair
{
public:
    /// Default constructor.
    GB_Pair() = default;

    /// TODO: Add Doxygen
    GB_Pair(bool act, OCP_USI i)
        : activity(act)
        , index(i){};

    /// return activity of some grid.
    bool IsAct() const { return activity; }

    /// return active index of some grid if active.
    OCP_USI GetId() const { return index; }

private:
    bool    activity; ///< activity of grid
    OCP_USI index;    ///< active index of grid if active
};

/// Grid contains basic information of grids of reservoir, the rock properties in grids
/// are also included. all grids are stored here, you can regard it as a database of
/// reservoir. considering the existence of inactive grids(whose volume of pores is too
/// small or effective resource is too little) or switch of activity of grid, Grid class
/// is necessary. Grid class is static while simulating, active grids will be stored in
/// bulks, which is "area" for calculating.
class Grid
{
    friend class Bulk;
    friend class BulkConn;
    friend class Well;

public:
    Grid() = default;
    /// Input the param from internal param structure to Grid.
    void InputParam(const ParamReservoir& rs_param);
    /// calculate the properties of Grid.
    void Setup();
    /// Setup Othogonal grid
    void SetupOthogonalGrid();
    /// Setup Neighbors for Othogonal grid
    void SetupNeighborOthogonalGrid();
    /// Calculate Akd for Othogonal Grid
    OCP_DBL CalAkdOthogonalGrid(const OCP_USI& bId, const OCP_USI& eId,
                                const USI& direction);
    /// Calculate the depth and volume of grids.
    void CalDepthVOthogonalGrid();
    /// Setup Corner Point Grid
    void SetupCornerGrid();
    /// Setup Neighbors for Othogonal grid
    void SetupNeighborCornerGrid(const COORD& CoTmp);
    /// Calculate Akd for Corner Point Grid
    OCP_DBL CalAkdCornerGrid(const GeneralConnect& conn);

    /// Calculate the active grid. If the volume, or the proportion of the effective
    /// parts is too small, then the grid is inactive, which means this grid dosen't
    /// paeticipate in the simumlation. Other rule can be given.
    void CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Mapping from grids to bulks.
    const GB_Pair& MapG2B(const OCP_USI& i) const { return activeMap_G2B[i]; }
    /// Return nx of grid.
    OCP_USI GetGridNx() const { return nx; }
    /// Return ny of grid.
    OCP_USI GetGridNy() const { return ny; }
    /// Return nz of grid.
    OCP_USI GetGridNz() const { return nz; }
    /// Return the num of all grids.
    OCP_USI GetGridNum() const { return numGrid; }
    /// Return the num of all connections.
    OCP_USI GetConnNum() const { return numConn; }
    /// Return thr num of bulks(active grids)
    OCP_USI GetActiveGridNum() const { return activeGridNum; }
    /// Return the index of active grid with (i, j, k).
    OCP_USI GetActIndex(const USI& i, const USI& j, const USI& k) const;

private:
    USI     nx;      ///< num of bulks along x-direction
    USI     ny;      ///< num of bulks along y-direction
    USI     nz;      ///< num of bulks along z-direction
    OCP_USI numGrid; ///< num of grids, nx * ny * nz
    OCP_USI numConn; ///< num of connection

    USI                   gridType;
    vector<vector<GPair>> gNeighbor; ///< Neighbors of grids.
    // Orthogonal Grid
    vector<OCP_DBL> tops;  ///< depth of top face of topest gird: nx*ny
    vector<OCP_DBL> depth; ///< depth of center of grid: numBulk.
    vector<OCP_DBL> dx;    ///< size of gird along the x direction: numBulk.
    vector<OCP_DBL> dy;    ///< size of grid along the y direction: numBulk.
    vector<OCP_DBL> dz;    ///< size of grid along the z direction: numBulk.
    // CornerPoint Grid
    vector<OCP_DBL> coord; ///< Lines of Corner Point Grid.
    vector<OCP_DBL> zcorn; ///< ZValues of Corner Point Grid.

    vector<OCP_DBL> v;    ///< volume of grid: numBulk.
    vector<OCP_DBL> ntg;  ///< net to gross of grid: numBulk
    vector<OCP_DBL> poro; ///< initial porosity of rock: numBulk
    vector<OCP_DBL> kx;   ///< Absolute permeability of rock along x direction: numBulk
    vector<OCP_DBL> ky;   ///< Absolute permeability of rock along y direction: numBulk
    vector<OCP_DBL> kz;   ///< Absolute permeability of rock along z direction: numBulk

    // Region
    vector<USI> SATNUM; ///< used to identify SAT region: numBulk.
    vector<USI> PVTNUM; ///< used to identify PVT region in blackoil model: numBulk.

    OCP_USI activeGridNum; ///< num of active grid.
    vector<OCP_USI>
        activeMap_B2G; ///< a index map form active grid to grid: activeGridNum.
    vector<GB_Pair> activeMap_G2B; ///< a index map form grid to active grid: numBulk.
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
/*----------------------------------------------------------------------------*/