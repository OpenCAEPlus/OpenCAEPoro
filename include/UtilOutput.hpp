/*! \file    UtilOutput.hpp
 *  \brief   Supply basic tools used to output files.
 *  \author  Shizhe Li
 *  \date    Oct/11/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILOUTPUT_HEADER__
#define __UTILOUTPUT_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

template <typename T>
constexpr auto OCP_SEP01(T n)
{
    return string(n, '-');
}

template <typename T>
constexpr auto OCP_SEP02(T n)
{
    return string(n, '=');
}

string GetIJKformat(const string& i, const string& j, const string& k, const USI& s);

string GetIJKformat(const USI& i, const USI& j, const USI& k, const USI& s);

#endif /* end if __UTILOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/11/2022      Create file                          */
/*----------------------------------------------------------------------------*/
