/*! \file    UtilInput.hpp
 *  \brief   Supply basic tools used to input files.
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILINPUT_HEADER__
#define __UTILINPUT_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

/// TODO: Replace it with error log
#define ParamCheck1(exp)                                                               \
    std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in "           \
              << __FILE__;

/// Map_str2int is used to map string to integer, which used to match the keyword
/// efficiently in input file in the switch structure.
constexpr inline long long Map_Str2Int(const char* mystr, const USI& len)
{
    long long res = 0;
    long long t   = 100;
    for (USI i = 0; i < len; i++) {
        res += (int)mystr[len - 1 - i] * t;
        t *= 100;
    }
    return res;
}

/// ReadLine is the core function while inputting the file. It will capture the next
/// line which is meanningful, for example, not blank line or comments, and then gets
/// rid of some useless characters such as space, commas. Finally, the segments of rest
/// string will be stored in result. And if return OCP_FALSE, it indicates we have reach
/// the end of file.
OCP_BOOL ReadLine(ifstream& ifs, vector<string>& result);

/// DealDefault is used to deal with the expression with asterisk, for example
/// m*n  -> <n,...,n> size m ,  m* -> <DEFAULT,..., DEFAULT> size m.
void DealDefault(vector<string>& result);

/// DealData change a series of product of integers into two arrays.
/// For example, 8*1  16*2  8*3  16*4  -> obj <8, 16, 8, 16> & val <1, 2, 3, 4>.
template <typename T>
void DealData(const vector<string>& vbuf, vector<OCP_USI>& obj, vector<T>& region)
{
    obj.resize(0);
    region.resize(0);
    for (auto& str : vbuf) {
        auto pos = str.find('*');
        if (pos != string::npos) {
            USI     len = str.size();
            OCP_USI num = stoi(str.substr(0, pos));
            USI     val = stoi(str.substr(pos + 1, len - (pos + 1)));
            obj.push_back(num);
            region.push_back(val);
        }
    }
}

#endif /* end if __UTILINPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/