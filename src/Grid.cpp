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

void Grid::InputParam(const ParamReservoir& rs_param)
{
    nx      = rs_param.dimens.nx;
    ny      = rs_param.dimens.ny;
    nz      = rs_param.dimens.nz;
    numGrid = rs_param.numGrid;

    if (!rs_param.coord.empty()) {
        // CornerPoint Grid
        gridType = CORNER_GRID;

        coord = rs_param.coord;
        zcorn = rs_param.zcorn;
    } else {
        // Orthogonal Grid
        gridType = ORTHOGONAL_GRID;

        numConn = 3 * nx * ny * nz - nx * ny - ny * nz - nz * nx;
        tops    = rs_param.tops;
        dx      = rs_param.dx;
        dy      = rs_param.dy;
        dz      = rs_param.dz;
    }

    ntg      = rs_param.ntg;
    poro     = rs_param.poro;
    kx       = rs_param.permX;
    ky       = rs_param.permY;
    kz       = rs_param.permZ;
    SwatInit = rs_param.Swat;

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
}

void Grid::Setup(const OCP_BOOL& myVTK)
{
    useVTK = myVTK;
    CalNumDigutIJK();
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

    // test
    CalSomeInfo();
}

void Grid::SetupOrthogonalGrid()
{
    // x -> y -> z
    SetupNeighborOrthogonalGrid();
    CalDepthVOrthogonalGrid();
    CalActiveGrid(1E-6, 1E-6);
    
    // for output
    SetHexaherdronGridOrthogonal();
}

void Grid::SetupNeighborOrthogonalGrid()
{

    gNeighbor.resize(numGrid);
    // PreAllocate
    for (OCP_USI n = 0; n < numGrid; n++) {
        gNeighbor[n].reserve(6);
    }

    // Begin Id and End Id in Grid, bIdg < eIdg
    OCP_USI bIdg, eIdg;
    OCP_DBL area;
    OCP_USI nxny = nx * ny;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {

                bIdg = k * nxny + j * nx + i;
                // right  --  x-direction
                if (i < nx - 1) {
                    eIdg = bIdg + 1;
                    area = CalAkdOrthogonalGrid(bIdg, eIdg, 1);
                    gNeighbor[bIdg].push_back(GPair(eIdg, area));
                    gNeighbor[eIdg].push_back(GPair(bIdg, area));
                }

                // front  --  y-direction
                if (j < ny - 1) {
                    eIdg = bIdg + nx;
                    area = CalAkdOrthogonalGrid(bIdg, eIdg, 2);
                    gNeighbor[bIdg].push_back(GPair(eIdg, area));
                    gNeighbor[eIdg].push_back(GPair(bIdg, area));
                }

                // down --   z-direction
                if (k < nz - 1) {
                    eIdg = bIdg + nxny;
                    area = CalAkdOrthogonalGrid(bIdg, eIdg, 3);
                    gNeighbor[bIdg].push_back(GPair(eIdg, area));
                    gNeighbor[eIdg].push_back(GPair(bIdg, area));
                }
            }
        }
    }

    OCP_FUNCNAME;
}

OCP_DBL
Grid::CalAkdOrthogonalGrid(const OCP_USI& bId, const OCP_USI& eId, const USI& direction)
{
    OCP_DBL T1;
    OCP_DBL T2;
    switch (direction) {
        case 1:
            // x-direction
            T1 = kx[bId] * ntg[bId] * dy[bId] * dz[bId] / dx[bId];
            T2 = kx[eId] * ntg[eId] * dy[eId] * dz[eId] / dx[eId];
            return (2 / (1 / T1 + 1 / T2));
            break;
        case 2:
            // y-direction
            T1 = ky[bId] * ntg[bId] * dz[bId] * dx[bId] / dy[bId];
            T2 = ky[eId] * ntg[eId] * dz[eId] * dx[eId] / dy[eId];
            return (2 / (1 / T1 + 1 / T2));
            break;
        case 3:
            // z-direction -- no ntg
            T1 = kz[bId] * dx[bId] * dy[bId] / dz[bId];
            T2 = kz[eId] * dx[eId] * dy[eId] / dz[eId];
            return (2 / (1 / T1 + 1 / T2));
            break;
        default:
            OCP_ABORT("Wrong Direction!");
    }
}

void Grid::CalDepthVOrthogonalGrid()
{
    depth.resize(numGrid, 0);
    OCP_USI nxny = nx * ny;
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

void Grid::SetupCornerGrid()
{
    OCP_COORD coordTmp;
    coordTmp.Allocate(nx, ny, nz);
    coordTmp.InputData(coord, zcorn);
    // coordTmp.CalConn();
    coordTmp.SetupCornerPoints();
    SetupNeighborCornerGrid(coordTmp);
    CalActiveGrid(1E-6, 1E-6);

    // for output
    SetHexaherdronGridCorner(coordTmp);
}

void Grid::SetupNeighborCornerGrid(const OCP_COORD& CoTmp)
{

    dx    = CoTmp.dx;
    dy    = CoTmp.dy;
    dz    = CoTmp.dz;
    v     = CoTmp.v;
    depth = CoTmp.depth;

    gNeighbor.resize(numGrid);
    // PreAllocate
    for (OCP_USI n = 0; n < numGrid; n++) {
        gNeighbor[n].reserve(10);
    }

    // test
    // cout << "Grid : " << CoTmp.numConn << endl;

    OCP_USI bIdg, eIdg;
    OCP_DBL area;
    for (OCP_USI n = 0; n < CoTmp.numConn; n++) {
        const GeneralConnect& ConnTmp = CoTmp.connect[n];
        bIdg                          = ConnTmp.begin;
        eIdg                          = ConnTmp.end;
        area                          = CalAkdCornerGrid(ConnTmp);
        gNeighbor[bIdg].push_back(GPair(eIdg, area));
        gNeighbor[eIdg].push_back(GPair(bIdg, area));

        // USI I, J, K;
        // GetIJKGrid(I, J, K, bIdg);
        // cout << "(" << setw(3) << I << "," << setw(3) << J << "," << setw(3) << K <<
        // ")    "; cout << setw(6) << bIdg; cout << "       "; GetIJKGrid(I, J, K,
        // eIdg); cout << "(" << setw(3) << I << "," << setw(3) << J << "," << setw(3)
        // << K << ")    "; cout << setw(6) << eIdg; cout << setw(20) << fixed <<
        // setprecision(4) << area;

        // cout << endl;
    }
}

OCP_DBL Grid::CalAkdCornerGrid(const GeneralConnect& conn)
{
    OCP_USI bId   = conn.begin;
    OCP_USI eId   = conn.end;
    OCP_DBL bArea = conn.Ad_dd_begin;
    OCP_DBL eArea = conn.Ad_dd_end;
    OCP_DBL T1, T2;

    switch (conn.directionType) {
        case 1:
            T1 = ntg[bId] * kx[bId] * bArea;
            T2 = ntg[eId] * kx[eId] * eArea;
            return 1 / (1 / T1 + 1 / T2);
            break;
        case 2:
            T1 = ntg[bId] * ky[bId] * bArea;
            T2 = ntg[eId] * ky[eId] * eArea;
            return 1 / (1 / T1 + 1 / T2);
            break;
        case 3:
            T1 = kz[bId] * bArea;
            T2 = kz[eId] * eArea;
            return 1 / (1 / T1 + 1 / T2);
            break;
        default:
            OCP_ABORT("Wrong Direction Type!");
    }
}

/// If porosity or volume of the grid cell is too small, then the cell is inactive.
//  Note: Inactive cells do NOT participate simumlation; other rules can be given.
void Grid::CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2)
{
    activeMap_B2G.reserve(numGrid);
    activeMap_G2B.resize(numGrid);
    OCP_USI count = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (ACTNUM[n] == 0 || poro[n] * ntg[n] < e1 || v[n] < e2) {
            activeMap_G2B[n] = GB_Pair(OCP_FALSE, 0);
            continue;
        }
        activeMap_B2G.push_back(n);
        activeMap_G2B[n] = GB_Pair(OCP_TRUE, count);
        count++;
    }
    activeGridNum = count;
    cout << (numGrid - activeGridNum) * 100.0 / numGrid << "% ("
         << (numGrid - activeGridNum) << ") of grid cell is inactive" << endl;
}

/// Return id of the active cell and abort if the cell is inactive!
OCP_USI Grid::GetActIndex(const USI& i, const USI& j, const USI& k) const
{
    OCP_USI id = k * nx * ny + j * nx + i;
    if (id > numGrid) {
        OCP_ABORT("Id is out of Range!");
    }
    OCP_BOOL activity = activeMap_G2B[id].IsAct();
    if (!activity) {
        OCP_ABORT("(" + to_string(i) + "," + to_string(j) + "," + to_string(k) +
                  ") is inactive");
    }
    id = activeMap_G2B[id].GetId();
    return id;
}

// temp
void Grid::GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const
{
    // i,j,k begin from 1
    // n must be the index of grids instead bulks
    k = n / (nx * ny) + 1;
    j = (n - (k - 1) * nx * ny) / nx + 1;
    i = n - (k - 1) * nx * ny - (j - 1) * nx + 1;
}

void Grid::GetIJKBulk(USI& i, USI& j, USI& k, const OCP_USI& n) const
{
    GetIJKGrid(i, j, k, activeMap_B2G[n]);
}

void Grid::CalSomeInfo() const
{
    // test
    OCP_DBL depthMax = 0;
    OCP_DBL depthMin = 1E8;
    OCP_DBL dxMax    = 0;
    OCP_DBL dxMin    = 1E8;
    OCP_DBL dyMax    = 0;
    OCP_DBL dyMin    = 1E8;
    OCP_DBL dzMax    = 0;
    OCP_DBL dzMin    = 1E8;

    for (OCP_USI n = 0; n < numGrid; n++) {
        // if (!activeMap_G2B[nn].IsAct())
        //     continue;
        // OCP_USI n = activeMap_G2B[nn].GetId();
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

    cout << "---------------------" << endl
         << "GRID" << endl
         << "---------------------" << endl;
    cout << "  depthMax = " << depthMax << endl
         << "  depthMin = " << depthMin << endl
         << "  dxMax    = " << dxMax << endl
         << "  dxMin    = " << dxMin << endl
         << "  dyMax    = " << dyMax << endl
         << "  dyMin    = " << dyMin << endl
         << "  dzMax    = " << dzMax << endl
         << "  dzMin    = " << dzMin << endl;
}


void Grid::SetHexaherdronGridOrthogonal()
{
    // x,y-coordinate begins from 0

    if (!useVTK) return;

    polyhedronGrid.reserve(numGrid);
    OCPpolyhedron tmpP(8);
    OCP_DBL tmpX, tmpY;
    OCP_USI id;

    for (USI k = 0; k < nz; k++) {
        tmpY = 0;
        for (USI j = 0; j < ny; j++) {
            tmpX = 0;
            for (USI i = 0; i < nx; i++) {
                id = k * nx * ny + j * nx + i;
                tmpP.Points.push_back(Point3D(tmpX, tmpY, depth[id] + dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX + dx[id], tmpY, depth[id] + dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX + dx[id], tmpY + dy[id], depth[id] + dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX, tmpY + dy[id], depth[id] + dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX, tmpY, depth[id] - dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX + dx[id], tmpY, depth[id] - dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX + dx[id], tmpY + dy[id], depth[id] - dz[id] / 2));
                tmpP.Points.push_back(Point3D(tmpX, tmpY + dy[id], depth[id] - dz[id] / 2));

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



void Grid::CalNumDigutIJK()
{
    OCP_ASSERT((nx > 0 && ny > 0 && nz > 0), "Wrong Dimension!");
    numDigutIJK = 1;
    if (log10(nx) > numDigutIJK) numDigutIJK = ceil(log10(nx));
    if (log10(ny) > numDigutIJK) numDigutIJK = ceil(log10(ny));
    if (log10(nz) > numDigutIJK) numDigutIJK = ceil(log10(nz));
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