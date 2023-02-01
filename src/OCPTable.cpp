/*! \file    OCPTable.cpp
 *  \brief   OCPTable class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPTable.hpp"

OCPTable::OCPTable(const USI& row, const USI& col)
{
    nRow = row;
    nCol = col;
    bId  = 0;
    data.resize(nCol);
    for (USI j = 0; j < nCol; j++) {
        data[j].resize(nRow);
    }
}

OCPTable::OCPTable(const vector<vector<OCP_DBL>>& src) { this->Setup(src); }

void OCPTable::Setup(const std::vector<std::vector<OCP_DBL>>& src)
{
    data = src;
    nCol = data.size();
    nRow = data[0].size();
    bId  = nRow / 2;
}

OCP_INT OCPTable::GetRowZero(const USI& mycol) const
{
    if (mycol > nCol) OCP_ABORT("wrong specified column!");
    for (USI i = 0; i < nRow; i++) {
        if (data[mycol][i] >= TINY) {
            return i - 1;
        }
    }
    return nRow - 1;
}

/// Careful: the memory outdata and slope have not be allocated before
USI OCPTable::Eval_All(const USI&       j,
                       const OCP_DBL&   val,
                       vector<OCP_DBL>& outdata,
                       vector<OCP_DBL>& slope)
{
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 0; k < nCol; k++) {
                    slope[k] = (data[k][bId + 1] - data[k][bId]) /
                               (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope[k] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        for (USI k = 0; k < nCol; k++) {
            slope[k]   = 0;
            outdata[k] = data[k].back();
        }
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 0; k < nCol; k++) {
                    slope[k] = (data[k][bId + 1] - data[k][bId]) /
                               (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope[k] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        for (USI k = 0; k < nCol; k++) {
            slope[k]   = 0;
            outdata[k] = data[k].front();
        }
    }
    return bId;
}

USI OCPTable::Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata)
{
    const USI j    = 0;
    OCP_DBL   tmpk = 0;
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 1; k < nCol; k++) {
                    tmpk = (data[k][bId + 1] - data[k][bId]) /
                           (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + tmpk * (val - data[j][bId]);
                }
                return bId;
            }
        }
        for (USI k = 1; k < nCol; k++) {
            outdata[k - 1] = data[k].back();
        }
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 1; k < nCol; k++) {
                    tmpk = (data[k][bId + 1] - data[k][bId]) /
                           (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + tmpk * (val - data[j][bId]);
                }
                return bId;
            }
        }
        for (USI k = 1; k < nCol; k++) {
            outdata[k - 1] = data[k].front();
        }
    }
    return bId;
}

OCP_DBL OCPTable::Eval(const USI& j, const OCP_DBL& val, const USI& destj)
{
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId       = i - 1;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        return data[destj].back();
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        return data[destj].front();
    }
}

OCP_DBL OCPTable::Eval(const USI& j, const OCP_DBL& val, const USI& destj, OCP_DBL& myK)
{
    // becareful when the memory outdata and slope have not be allocated before

    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                myK = (data[destj][bId + 1] - data[destj][bId]) /
                      (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + myK * (val - data[j][bId]));
            }
        }
        return data[destj].back();
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                myK = (data[destj][bId + 1] - data[destj][bId]) /
                      (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + myK * (val - data[j][bId]));
            }
        }
        return data[destj].front();
    }
}

OCP_DBL OCPTable::Eval_Inv(const USI& j, const OCP_DBL& val, const USI& destj)
{
    // becareful when the memory outdata and slope have not be allocated before

    if (val > data[j][bId]) {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val <= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        return data[destj].front();
    } else {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val >= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId] - data[destj][bId - 1]) /
                            (data[j][bId] - data[j][bId - 1]);
                return (data[destj][bId - 1] + k * (val - data[j][bId - 1]));
            }
        }
        return data[destj].back();
    }
}

void OCPTable::Display() const
{
    cout << "\n---------------------" << endl
         << "Pressure Distribution"
         << "\n---------------------" << endl;
    cout << "   Depth         "
         << "   Poil         "
         << "   Pgas         "
         << "   Pwat         " << endl;

    cout << fixed << setprecision(3);
    for (USI i = 0; i < nRow; i++) {
        for (USI j = 0; j < nCol; j++) {
            cout << data[j][i] << "\t";
        }
        cout << "\n";
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/