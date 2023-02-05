/*! \file    UtilOutput.cpp
 *  \brief   Utilities for reading output file
 *  \author  Shizhe Li
 *  \date    Oct/11/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "UtilOutput.hpp"

string GetIJKformat(const string& i, const string& j, const string& k, const USI& s)
{
    std::ostringstream IJKinfo;
    IJKinfo << "(" << setw(s) << i << ", " << setw(s) << j << ", " << setw(s) << k
            << ")";
    return IJKinfo.str();
}

string GetIJKformat(const USI& i, const USI& j, const USI& k, const USI& s)
{
    return GetIJKformat(to_string(i), to_string(j), to_string(k), s);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/11/2022      Create file                          */
/*----------------------------------------------------------------------------*/