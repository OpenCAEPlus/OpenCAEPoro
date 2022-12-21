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


class OCPRes
{
public:
    void Setup_IsoT(const OCP_USI& nb, const OCP_USI& nw, const OCP_USI& nc) {
        OCP_USI reslen = (nb + nw) * (nc + 1);
        res.resize(reslen);
        resRelV.resize(nb);
        resRelN.resize(nb);
    }
    void SetZero()
    {
        fill(res.begin(), res.end(), 0);
        fill(resRelV.begin(), resRelV.end(), 0);
        fill(resRelN.begin(), resRelN.end(), 0);
        maxRelRes_v   = 0;
        maxRelRes_mol = 0;
        maxWellRelRes_mol = 0;
        maxId_v = 0;
        maxId_mol = 0;
    }
    void SetInitRes() {
        maxRelRes0_v = maxRelRes_v;
    }

    vector<OCP_DBL> res;        ///< residual for all equations for each bulk
    vector<OCP_DBL> resRelV;    ///< 2-norm of relative residual wrt. pore volume for all equations of each bulk
    vector<OCP_DBL> resRelN;    ///< 2-norm of relative residual wrt. total moles for mass conserve equations of each bulk
    vector<OCP_DBL> resRelT;    ///< 2-norm of relative residual wrt. total energy for energy conserve equations of each bulk
    OCP_DBL        maxRelRes0_v;    ///< (initial) maximum relative residual wrt. pore volume for each bulk
    OCP_DBL        maxRelRes_v;     ///< (iterations) maximum relative residual wrt. pore volume for each bulk
    OCP_DBL        maxRelRes_mol;   ///< maximum relative residual wrt. total moles for each bulk
    OCP_DBL        maxWellRelRes_mol;   ///< maximum relative residual wrt. total moles for each well

    // use negative number to represent well number (ToDo)
    OCP_INT        maxId_v;     ///< index of bulk who has maxRelRes_v
    OCP_INT        maxId_mol;   ///< index of bulk who has maxRelRes_mol
};

#endif