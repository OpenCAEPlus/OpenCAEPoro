/*! \file    LinearSystem.hpp
 *  \brief   Class declaration for internal linear systems
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#pragma once

#include <vector>

#include "OCPConst.hpp"

/// CSRx matrix type based on FASP++
template <typename T> class MAT_Faspxx
{
public:
    void ClearData();

    OCP_USI              NRow;
    OCP_USI              NCol;
    OCP_USI              Nnz;
    std::vector<T>       val;
    std::vector<OCP_USI> colId;
    std::vector<OCP_USI> RowPtr;
    std::vector<USI>     diagPtr;
};

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/