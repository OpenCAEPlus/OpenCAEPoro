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
#include "ParamOutput.hpp"
#include "ParamReservoir.hpp"
#include "UtilOutput.hpp"

using namespace std;

/// Active cell indicator and its index among active cells.
//  Note: GB_Pair contains two variables, which indicates if a grid cell is active or
//  not and its index among active cells if it is active.
class GB_Pair
{
public:
    /// Default constructor.
    GB_Pair() = default;
    /// Constructor with given information.
    GB_Pair(const OCP_BOOL& act, const OCP_USI& i)
        : activity(act)
        , index(i){};

    /// Return whether this cell is active or not.
    OCP_BOOL IsAct() const { return activity; }
    /// Return the index of this cell among active cells.
    OCP_USI GetId() const { return index; }

private:
    OCP_BOOL activity; ///< Activeness of a grid cell.
    OCP_USI  index;    ///< Active index of grid if active
};

/// Effective area of intersection surfaces with neighboring cells.
class GPair
{
public:
    GPair() = default;
    GPair(const OCP_USI& Id,
          const USI&     Direct,
          const OCP_DBL& AreaB,
          const OCP_DBL& AreaE)
        : id(Id)
        , direction(Direct)
        , areaB(AreaB)
        , areaE(AreaE){};
    static OCP_BOOL lessG(const GPair& G1, const GPair& G2) { return G1.id < G2.id; }

    OCP_USI id;        ///< Id of a neighboring cell
    USI     direction; ///< direction: 1-x, 2-y, 3-z
    OCP_DBL
        areaB; ///< Effective intersection area between this cell and the neighbor, self
    OCP_DBL areaE; ///< Effective intersection area between this cell and the neighbor,
                   ///< neighbor
};

class OCPpolyhedron
{
public:
    OCPpolyhedron() = default;
    OCPpolyhedron(const USI& n)
        : numPoints(n)
    {
        Points.reserve(numPoints);
    };
    vector<Point3D> Points;
    USI             numPoints;
};

///< Record the initial grid information, all of grids are contained
class Grid
{
    friend class Bulk;
    friend class BulkConn;
    friend class Well;
    friend class ScalePcow;
    friend class Out4RPT;
    friend class Out4VTK;

    /////////////////////////////////////////////////////////////////////
    // Input Param and Setup
    /////////////////////////////////////////////////////////////////////

public:
    /// Input parameters from the internal param structure.
    void InputParam(const ParamReservoir& rs_param, const ParamOutput& output_param);

    /// Setup for Isothermal model
    void SetupIsoT();
    /// Setup for thermal model
    void SetupT();
    /// Setup the grid information and calculate the properties.
    void Setup();

protected:
    /// Setup orthogonal grid.
    void SetupOrthogonalGrid();
    /// Calculate the depth and volume for orthogonal grid.
    void CalDepthVOrthogonalGrid();
    /// Setup the neighboring info for an orthogonal grid.
    void SetupNeighborOrthogonalGrid();

    /// Setup corner-point grid.
    void SetupCornerGrid();
    /// Setup dx,dy,dz,depth, v for a corner-point grid.
    void SetupBasicCornerGrid(const OCP_COORD& CoTmp);
    /// Setup the neighboring info for a corner-point grid.
    void SetupNeighborCornerGrid(const OCP_COORD& CoTmp);

    /// Calculate the activity of grid cells.
    void CalActiveGridIsoT(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Calculate the activity of grid cells for ifThermal model
    void CalActiveGridT(const OCP_DBL& e1, const OCP_DBL& e2);

    /// Setup Grid location for Structured grid
    void SetupGridLocation();

public:
    OCP_USI GetGridNum() const { return numGrid; }
    OCP_INT GetActIndex(const USI& I, const USI& J, const USI& K) const;

protected:
    // Grid type
    USI     gridType; ///< Orthogonal or Corner grid
    OCP_USI numGrid;  ///< Number of all cells

    // structured grid
    USI nx; ///< Number of cells in x-direction
    USI ny; ///< Number of cells in y-direction
    USI nz; ///< Number of cells in z-direction

    // Corner-point grid
    vector<OCP_DBL> coord; ///< Lines of a corner-point grid.
    vector<OCP_DBL> zcorn; ///< ZValues of a corner-point grid.

    // Grid location
    // Isothermal model: useless now.
    // Thermal model: only top face and bottom face will be recorded.
    vector<USI> gLocation; ///< Top face, bottom face, side face, numGrid.

    // Orthogonal grid
    vector<OCP_DBL> tops; ///< Depth of center of grid cells: numGrid.

    // General informations
    vector<OCP_DBL> dx;    ///< Size of cell in x-direction: numGrid.
    vector<OCP_DBL> dy;    ///< Size of cell in y-direction: numGrid.
    vector<OCP_DBL> dz;    ///< Size of cell in z-direction: numGrid.
    vector<OCP_DBL> v;     ///< Volume of cells: numGrid.
    vector<OCP_DBL> depth; ///< Depth of center of grid cells: numGrid.

    // Rock properties
    vector<OCP_DBL> ntg;    ///< Net to gross ratio of cells: numGrid
    vector<OCP_DBL> poro;   ///< Initial porosity of rock cells: numGrid
    vector<OCP_DBL> kx;     ///< Absolute permeability in x-direction: numGrid
    vector<OCP_DBL> ky;     ///< Absolute permeability in y-direction: numGrid
    vector<OCP_DBL> kz;     ///< Absolute permeability in z-direction: numGrid
    vector<OCP_DBL> thconr; ///< Rock if Thermal conductivity: numGrid

    // Region
    vector<USI> SATNUM;  ///< Identify SAT region: numGrid.
    vector<USI> PVTNUM;  ///< Identify PVT region for the blackoil model: numGrid.
    vector<USI> ACTNUM;  ///< Indicate activity of grid from input file: numGrid. 0 =
                         ///< inactive, 1 = active.
    vector<USI> ROCKNUM; ///< index of rock table for each grid: numGrid

    // Initial Properties
    vector<OCP_DBL> SwatInit; ///< Initial water saturation

    // Connections
    vector<vector<GPair>> gNeighbor; ///< Neighboring information of grid.

    // Active grid cells
    OCP_USI activeGridNum; ///< Num of active grid.
    vector<OCP_USI>
        map_Act2All; ///< Mapping from active grid to all grid: activeGridNum.
    vector<GB_Pair> map_All2Act; ///< Mapping from grid to active all grid: numGrid.
    // Fluid grid cells
    OCP_USI         fluidGridNum; ///< Num of fluid grids.
    vector<GB_Pair> map_All2Flu;  ///< Mapping from all grid to fluid grid: numGrid.

    /////////////////////////////////////////////////////////////////////
    // Output
    /////////////////////////////////////////////////////////////////////

public:
    void     GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const;
    void     GetIJKBulk(USI& i, USI& j, USI& k, const OCP_USI& n) const;
    OCP_BOOL IfUseVtk() const
    {
        return useVTK;
    }                                    ///< return if use vtk format for outputing
    void SetHexaherdronGridOrthogonal(); ///< setup polyhedronGrid for orthogonal grid
    void SetHexaherdronGridCorner(
        const OCP_COORD& mycord); ///< setup polyhedronGrid for corner grid
    /// Setup grid tag
    void SetupGridTag();
    void OutputBaiscInfo() const; ///< Calculate and return basic informations for grid
    void CalNumDigutIJK();        ///< only used in structured grid
    USI  GetNumDigitIJK() const { return numDigutIJK; } ///< Return numDigutIJK
protected:
    OCP_BOOL              useVTK{OCP_FALSE}; ///< If output in vtk format
    vector<OCPpolyhedron> polyhedronGrid;    ///< Coordinates of grid points
    vector<USI>           gridTag;     ///< Tag of grid: dead, live(fluid), live(rock)
    USI                   numDigutIJK; ///< number of digits of maximum nx,ny,nz
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