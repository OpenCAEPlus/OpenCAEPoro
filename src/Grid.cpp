/*! \file    Grid.cpp
 *  \brief   Grid class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Grid.hpp"

/////////////////////////////////////////////////////////////////////
// Input Param and Setup
/////////////////////////////////////////////////////////////////////

void Grid::InputParam(const ParamReservoir& rs_param, const ParamOutput& output_param)
{
    // dimensions
    nx      = rs_param.dimens.nx;
    ny      = rs_param.dimens.ny;
    nz      = rs_param.dimens.nz;
    numGrid = rs_param.numGrid;

    // CornerPoint Grid or Orthogonal Grid
    if (!rs_param.coord.empty()) {
        // CornerPoint Grid
        gridType = CORNER_GRID;
        coord    = rs_param.coord;
        zcorn    = rs_param.zcorn;
    } else {
        // Orthogonal Grid
        gridType = ORTHOGONAL_GRID;
        tops     = rs_param.tops;
        dx       = rs_param.dx;
        dy       = rs_param.dy;
        dz       = rs_param.dz;
    }

    // General Information
    ntg    = rs_param.ntg;
    poro   = rs_param.poro;
    kx     = rs_param.permX;
    ky     = rs_param.permY;
    kz     = rs_param.permZ;
    thconr = rs_param.thconr;

    // Regions
    SATNUM.resize(numGrid, 0);
    if (rs_param.SATNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            SATNUM[i] = round(rs_param.SATNUM.data[i]) - 1;
        }
    }
    PVTNUM.resize(numGrid, 0);
    if (rs_param.PVTNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            PVTNUM[i] = round(rs_param.PVTNUM.data[i]) - 1;
        }
    }
    ACTNUM.resize(numGrid, 1);
    if (rs_param.ACTNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            ACTNUM[i] = round(rs_param.ACTNUM.data[i]);
        }
    }
    ROCKNUM.resize(numGrid, 0);
    if (rs_param.ROCKNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            ROCKNUM[i] = round(rs_param.ROCKNUM.data[i]);
        }
    }

    // Initial Properties
    SwatInit = rs_param.Swat;

    // Output
    useVTK = output_param.outVTKParam.useVTK;
}

void Grid::SetupIsoT()
{
    Setup();
    CalActiveGridIsoT(1E-6, 1E-6);
    SetupGridTag();
}

void Grid::SetupT()
{
    Setup();
    CalActiveGridT(1E-6, 1E-6);
    SetupGridLocation();
    SetupGridTag();
}

void Grid::Setup()
{
    switch (gridType) {
        case ORTHOGONAL_GRID:
            SetupOrthogonalGrid();
            break;
        case CORNER_GRID:
            SetupCornerGrid();
            break;
        default:
            OCP_ABORT("WRONG Grid Type!");
    }

    CalNumDigutIJK();
    OutputBaiscInfo();
}

void Grid::SetupOrthogonalGrid()
{
    // x -> y -> z
    CalDepthVOrthogonalGrid();
    SetupNeighborOrthogonalGrid();
    // for output
    SetHexaherdronGridOrthogonal();
}

void Grid::CalDepthVOrthogonalGrid()
{
    depth.resize(numGrid, 0);
    const OCP_USI nxny = nx * ny;
    // 0th layer
    for (USI j = 0; j < ny; j++) {
        for (USI i = 0; i < nx; i++) {
            OCP_USI id = j * nx + i;
            depth[id]  = tops[id] + dz[id] / 2;
        }
    }
    // 1th - (nz-1)th layer
    for (USI k = 1; k < nz; k++) {
        OCP_USI knxny = k * nxny;
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                OCP_USI id = knxny + j * nx + i;
                depth[id]  = depth[id - nxny] + dz[id - nxny] / 2 + dz[id] / 2;
            }
        }
    }

    v.resize(numGrid);
    for (OCP_USI i = 0; i < numGrid; i++) v[i] = dx[i] * dy[i] * dz[i];
}

void Grid::SetupNeighborOrthogonalGrid()
{
    gNeighbor.resize(numGrid);
    // PreAllocate
    for (OCP_USI n = 0; n < numGrid; n++) {
        gNeighbor[n].reserve(6);
    }

    // Begin Id and End Id in Grid, bIdg < eIdg
    OCP_USI       bIdg, eIdg;
    OCP_DBL       areaB, areaE;
    USI           direction;
    const OCP_USI nxny = nx * ny;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {

                bIdg = k * nxny + j * nx + i;
                // right  --  x-direction
                if (i < nx - 1) {
                    direction = 1;
                    eIdg      = bIdg + 1;
                    areaB     = 2 * dy[bIdg] * dz[bIdg] / dx[bIdg];
                    areaE     = 2 * dy[eIdg] * dz[eIdg] / dx[eIdg];
                    gNeighbor[bIdg].push_back(GPair(eIdg, direction, areaB, areaE));
                    gNeighbor[eIdg].push_back(GPair(bIdg, direction, areaE, areaB));
                }
                // front  --  y-direction
                if (j < ny - 1) {
                    direction = 2;
                    eIdg      = bIdg + nx;
                    areaB     = 2 * dz[bIdg] * dx[bIdg] / dy[bIdg];
                    areaE     = 2 * dz[eIdg] * dx[eIdg] / dy[eIdg];
                    gNeighbor[bIdg].push_back(GPair(eIdg, direction, areaB, areaE));
                    gNeighbor[eIdg].push_back(GPair(bIdg, direction, areaE, areaB));
                }
                // down --   z-direction
                if (k < nz - 1) {
                    direction = 3;
                    eIdg      = bIdg + nxny;
                    areaB     = 2 * dx[bIdg] * dy[bIdg] / dz[bIdg];
                    areaE     = 2 * dx[eIdg] * dy[eIdg] / dz[eIdg];
                    gNeighbor[bIdg].push_back(GPair(eIdg, direction, areaB, areaE));
                    gNeighbor[eIdg].push_back(GPair(bIdg, direction, areaE, areaB));
                }
            }
        }
    }
}

void Grid::SetupCornerGrid()
{
    OCP_COORD coordTmp;
    coordTmp.Allocate(nx, ny, nz);
    coordTmp.InputData(coord, zcorn);
    coordTmp.SetupCornerPoints();
    SetupBasicCornerGrid(coordTmp);
    SetupNeighborCornerGrid(coordTmp);
    // for output
    SetHexaherdronGridCorner(coordTmp);
}

void Grid::SetupBasicCornerGrid(const OCP_COORD& CoTmp)
{
    dx    = CoTmp.dx;
    dy    = CoTmp.dy;
    dz    = CoTmp.dz;
    v     = CoTmp.v;
    depth = CoTmp.depth;
}

void Grid::SetupNeighborCornerGrid(const OCP_COORD& CoTmp)
{
    gNeighbor.resize(numGrid);
    // PreAllocate
    for (OCP_USI n = 0; n < numGrid; n++) {
        gNeighbor[n].reserve(10);
    }

    OCP_USI bIdg, eIdg;
    OCP_DBL areaB, areaE;
    USI     direction;
    for (OCP_USI n = 0; n < CoTmp.numConn; n++) {
        const GeneralConnect& ConnTmp = CoTmp.connect[n];
        direction                     = ConnTmp.directionType;
        bIdg                          = ConnTmp.begin;
        eIdg                          = ConnTmp.end;
        areaB                         = ConnTmp.Ad_dd_begin;
        areaE                         = ConnTmp.Ad_dd_end;
        gNeighbor[bIdg].push_back(GPair(eIdg, direction, areaB, areaE));
        gNeighbor[eIdg].push_back(GPair(bIdg, direction, areaE, areaB));
    }
}

/////////////////////////////////////////////////////////////////////
// Basic Grid and Rock Information
/////////////////////////////////////////////////////////////////////

/// If porosity or volume of the grid cell is too small, then the cell is inactive.
//  Note: Inactive cells do NOT participate simumlation; other rules can be given.
void Grid::CalActiveGridIsoT(const OCP_DBL& e1, const OCP_DBL& e2)
{
    map_Act2All.reserve(numGrid);
    map_All2Act.resize(numGrid);
    OCP_USI count = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (ACTNUM[n] == 0 || poro[n] * ntg[n] < e1 || v[n] < e2) {
            map_All2Act[n] = GB_Pair(OCP_FALSE, 0);
            ACTNUM[n]      = 0;
            continue;
        }
        map_Act2All.push_back(n);
        map_All2Act[n] = GB_Pair(OCP_TRUE, count);
        count++;
    }
    activeGridNum = count;
    if (numGrid > activeGridNum) {
        cout << "  Number of inactive cells is " << (numGrid - activeGridNum) << " ("
             << (numGrid - activeGridNum) * 100.0 / numGrid << "%)" << endl;
    }

    // fluid grid = active grid
    fluidGridNum = activeGridNum;
    map_All2Flu  = map_All2Act;
}

void Grid::CalActiveGridT(const OCP_DBL& e1, const OCP_DBL& e2)
{
    map_Act2All.reserve(numGrid);
    map_All2Act.resize(numGrid);
    map_All2Flu.resize(numGrid);
    OCP_USI activeCount = 0;
    OCP_USI fluidCount  = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (ACTNUM[n] == 0 || v[n] < e1) {
            map_All2Act[n] = GB_Pair(OCP_FALSE, 0);
            map_All2Flu[n] = GB_Pair(OCP_FALSE, 0);
            ACTNUM[n]      = 0;
        } else {
            if (poro[n] * ntg[n] < e2) {
                map_All2Flu[n] = GB_Pair(OCP_FALSE, 0);
            } else {
                map_All2Flu[n] = GB_Pair(OCP_TRUE, fluidCount);
                fluidCount++;
            }
            map_Act2All.push_back(n);
            map_All2Act[n] = GB_Pair(OCP_TRUE, activeCount);
            activeCount++;
        }
    }
    activeGridNum = activeCount;
    if (numGrid > activeGridNum) {
        cout << "  Number of inactive cells is " << (numGrid - activeGridNum) << " ("
             << (numGrid - activeGridNum) * 100.0 / numGrid << "%)" << endl;
    }

    fluidGridNum = fluidCount;
    // temp
    // output num of fluid grid
}

void Grid::SetupGridLocation()
{
    gLocation.resize(numGrid);
    // Be careful if there's only one floor, then top face is also bottom face.
    // then one index is not enough
    // top face
    for (OCP_USI n = 0; n < nx * ny; n++) {
        gLocation[n] = 1;
    }
    // bottom face
    for (OCP_USI n = nx * ny * (nz - 1); n < nx * ny * nz; n++) {
        gLocation[n] = 2;
    }
}

OCP_INT Grid::GetActIndex(const USI& I, const USI& J, const USI& K) const
{
    const OCP_USI& n = K * nx * ny + J * nx + I;
    if (map_All2Flu[n].IsAct()) {
        return map_All2Act[n].GetId();
    } else {
        return -1;
    }
}

/////////////////////////////////////////////////////////////////////
// Output
/////////////////////////////////////////////////////////////////////

void Grid::GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const
{
    k = n / (nx * ny) + 1;
    j = (n - (k - 1) * nx * ny) / nx + 1;
    i = n - (k - 1) * nx * ny - (j - 1) * nx + 1;
}

void Grid::GetIJKBulk(USI& i, USI& j, USI& k, const OCP_USI& n) const
{
    GetIJKGrid(i, j, k, map_Act2All[n]);
}

void Grid::SetHexaherdronGridOrthogonal()
{
    // x,y-coordinate begins from 0
    if (!useVTK) return;

    polyhedronGrid.reserve(numGrid);
    OCPpolyhedron tmpP(8);
    OCP_DBL       tmpX, tmpY;
    OCP_USI       id;

    for (USI k = 0; k < nz; k++) {
        tmpY = 0;
        for (USI j = 0; j < ny; j++) {
            tmpX = 0;
            for (USI i = 0; i < nx; i++) {
                id = k * nx * ny + j * nx + i;
                tmpP.Points.push_back(Point3D(tmpX, tmpY, depth[id] + dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX + dx[id], tmpY, depth[id] + dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX + dx[id], tmpY + dy[id], depth[id] + dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX, tmpY + dy[id], depth[id] + dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX, tmpY, depth[id] - dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX + dx[id], tmpY, depth[id] - dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX + dx[id], tmpY + dy[id], depth[id] - dz[id] / 2));
                tmpP.Points.push_back(
                    Point3D(tmpX, tmpY + dy[id], depth[id] - dz[id] / 2));

                polyhedronGrid.push_back(tmpP);
                tmpP.Points.clear();
                tmpX += dx[id];
            }
            tmpY += dy[id];
        }
    }
}

void Grid::SetHexaherdronGridCorner(const OCP_COORD& mycord)
{
    if (!useVTK) return;

    polyhedronGrid.reserve(numGrid);
    OCPpolyhedron tmpP(8);

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p4);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p5);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p6);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p7);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p0);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p1);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p2);
                tmpP.Points.push_back(mycord.cornerPoints[i][j][k].p3);
                polyhedronGrid.push_back(tmpP);
                tmpP.Points.clear();
            }
        }
    }
}

void Grid::SetupGridTag()
{
    if (!useVTK) return;

    gridTag.resize(numGrid);
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (map_All2Act[n].IsAct()) {
            if (map_All2Flu[n].IsAct()) {
                gridTag[n] = 2;
            } else {
                gridTag[n] = 1;
            }
        } else {
            gridTag[n] = 0;
        }
    }
}

void Grid::OutputBaiscInfo() const
{
    OCP_DBL depthMax = 0;
    OCP_DBL depthMin = 1E8;
    OCP_DBL dxMax    = 0;
    OCP_DBL dxMin    = 1E8;
    OCP_DBL dyMax    = 0;
    OCP_DBL dyMin    = 1E8;
    OCP_DBL dzMax    = 0;
    OCP_DBL dzMin    = 1E8;

    for (OCP_USI n = 0; n < numGrid; n++) {
        if (depthMax < depth[n]) {
            depthMax = depth[n];
        }
        if (depthMin > depth[n]) {
            depthMin = depth[n];
        }
        if (dxMax < dx[n]) {
            dxMax = dx[n];
        }
        if (dxMin > dx[n]) {
            dxMin = dx[n];
        }
        if (dyMax < dy[n]) {
            dyMax = dy[n];
        }
        if (dyMin > dy[n]) {
            dyMin = dy[n];
        }
        if (dzMax < dz[n]) {
            dzMax = dz[n];
        }
        if (dzMin > dz[n]) {
            dzMin = dz[n];
        }
    }

    cout << "\n---------------------" << endl
         << "GRID"
         << "\n---------------------" << endl;
    cout << "  depthMax = " << depthMax << endl
         << "  depthMin = " << depthMin << endl
         << "  dxMax    = " << dxMax << endl
         << "  dxMin    = " << dxMin << endl
         << "  dyMax    = " << dyMax << endl
         << "  dyMin    = " << dyMin << endl
         << "  dzMax    = " << dzMax << endl
         << "  dzMin    = " << dzMin << endl;
}

void Grid::CalNumDigutIJK()
{
    OCP_ASSERT((nx > 0 && ny > 0 && nz > 0), "Wrong Dimension!");
    numDigutIJK = 1;
    if (log10(nx) >= numDigutIJK) numDigutIJK = ceil(log10(nx) + 1);
    if (log10(ny) >= numDigutIJK) numDigutIJK = ceil(log10(ny) + 1);
    if (log10(nz) >= numDigutIJK) numDigutIJK = ceil(log10(nz) + 1);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/16/2022      Fix Doxygen                          */
/*----------------------------------------------------------------------------*/