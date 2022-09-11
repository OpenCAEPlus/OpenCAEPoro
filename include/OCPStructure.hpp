/*! \file    OCPStructure.hpp
 *  \brief   Some Structure in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPSTRUCTURE_HEADER__
#define __OCPSTRUCTURE_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

/// FIM
class ResFIM
{
public:
    void SetZero()
    {
        fill(res.begin(), res.end(), 0);
        maxRelRes_v   = 0;
        maxRelRes_mol = 0;
    }

    vector<double> res;
    OCP_DBL        maxRelRes0_v;
    OCP_DBL        maxRelRes_v;
    OCP_DBL        maxRelRes_mol;

    // use negative number to represent well number (ToDo)
    OCP_INT        maxId_v;
    OCP_INT        maxId_mol;
};

#endif