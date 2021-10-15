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

void Grid::Setup()
{
    CalDepthV();
    CalActiveGrid(1E-6, 1E-6);
}

void Grid::CalDepthV()
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
    cout << "Grid::calDepthV" << endl;
}

void Grid::CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2)
{
    activeMap_B2G.reserve(numGrid);
    activeMap_G2B.resize(numGrid);
    OCP_USI count = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        if (poro[n] * ntg[n] < e1 || dx[n] * dy[n] * dz[n] < e2) {
            activeMap_G2B[n] = GB_Pair(false, 0);
            continue;
        }
        activeMap_B2G.push_back(n);
        activeMap_G2B[n] = GB_Pair(true, count);
        count++;
    }
    activeGridNum = count;
    cout << (numGrid - activeGridNum) / numGrid << "%  grids is inactive" << endl;
    cout << "Grid::calActiveBulk" << endl;
}

void Grid::InputParam(const ParamReservoir& rs_param)
{
    nx      = rs_param.dimens.nx;
    ny      = rs_param.dimens.ny;
    nz      = rs_param.dimens.nz;
    numGrid = rs_param.numGrid;
    numConn = 3 * nx * ny * nz - nx * ny - ny * nz - nz * nx;

    tops = rs_param.tops;
    dx   = rs_param.dx;
    dy   = rs_param.dy;
    dz   = rs_param.dz;

    ntg  = rs_param.ntg;
    poro = rs_param.poro;
    kx   = rs_param.permX;
    ky   = rs_param.permY;
    kz   = rs_param.permZ;

    SATNUM.resize(numGrid, 0);
    if (rs_param.SATNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            SATNUM[i] = (USI)(rs_param.SATNUM.data[i]) - 1;
        }
    }
    PVTNUM.resize(numGrid, 0);
    if (rs_param.PVTNUM.activity) {
        for (OCP_USI i = 0; i < numGrid; i++) {
            PVTNUM[i] = (USI)(rs_param.PVTNUM.data[i]) - 1;
        }
    }
    cout << "Grid::InputParam" << endl;
}

OCP_USI Grid::GetActIndex(const USI& i, const USI& j, const USI& k) const
{
    OCP_USI id       = k * nx * ny + j * nx + i;
    bool    activity = activeMap_G2B[id].GetAct();
    if (!activity) {
        ERRORcheck("(" + to_string(i) + "," + to_string(j) + "," + to_string(k) +
                   ") is inactive");
        exit(0);
    }
    id = activeMap_G2B[id].GetId();
    return id;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/