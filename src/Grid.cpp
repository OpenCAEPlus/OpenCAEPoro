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
    initInfo.InputParam(rs_param);
    bulk.InputParam(rs_param);
}


void Grid::SetupIsoT()
{
    initInfo.Setup();
    initInfo.CalActiveGridIsoT(1E-6, 1E-6);
    bulk.SetupIsoT(initInfo);
}


void Grid::SetupT()
{
    initInfo.Setup();
    initInfo.CalActiveGridT(1E-6, 1E-6);
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