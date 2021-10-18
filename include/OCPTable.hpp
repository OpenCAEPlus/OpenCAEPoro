/*! \file    OCPTable.hpp
 *  \brief   OCPTable class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_Table_HEADER__
#define __OCP_Table_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

/// OCP_Table is a template Table class, which used to deal with everything about table
/// in OpenCAEPoro such as PVT table, saturation table.
template <typename T> class OCP_Table
{
public:
    /// Default constructor.
    OCP_Table() = default;

    /// Constructor a table of fixed size.
    OCP_Table(const USI& row, const USI& col);

    /// Setup tables from existing data of table.
    void Setup(const vector<vector<T>>& src);

    /// judge if table is empty.
    bool IsEmpty() const { return data.empty(); }

    /// return the column num of table.
    USI GetCol() const { return nCol; }

    /// push v into the last column of table.
    void PushCol(const vector<T>& v) { data.push_back(v); }

    /// return the jth column in table to modify or use.
    vector<T>& GetCol(const USI& j) { return data[j]; }

    /// Setup row nums and col nums of tables, initialize the bId.
    void SetRowCol()
    {
        nRow = data[0].size();
        nCol = data.size();
        bId  = nRow / 2;
    }

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// all columns.
    USI Eval_All(const USI& j, const T& val, vector<T>& outdata, vector<T>& slope);

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// the target column.
    T Eval(const USI& j, const T& val, const USI& destj);

    /// interpolate the specified monotonically decreasing column in table to evaluate
    /// the target column.

    T Eval_Inv(const USI& j, const T& val, const USI& destj);

    /// Display the data of table on screen.
    void Display() const;

private:
    USI               nRow; ///< number of rows of the table
    USI               nCol; ///< number of columns of the table
    USI               bId;  ///< the starting point of rows when interpolating
    vector<vector<T>> data; ///< data of the table, data[i] is the ith column.
};

template <typename T> OCP_Table<T>::OCP_Table(const USI& row, const USI& col)
{
    nRow = row;
    nCol = col;
    bId  = 0;
    data.resize(nCol);
    for (USI j = 0; j < nCol; j++) {
        data[j].resize(nRow);
    }
}

template <typename T> void OCP_Table<T>::Setup(const std::vector<std::vector<T>>& src)
{
    data = src;
    nCol = data.size();
    nRow = data[0].size();
    bId  = nRow / 2;
}

template <typename T>
inline USI OCP_Table<T>::Eval_All(const USI& j, const T& val, vector<T>& outdata,
                                  vector<T>& slope)
{
    // becareful when the memory outdata and slope have not be allocated before

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

template <typename T>
inline T OCP_Table<T>::Eval(const USI& j, const T& val, const USI& destj)
{
    // becareful when the memory outdata and slope have not be allocated before

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

template <typename T>
inline T OCP_Table<T>::Eval_Inv(const USI& j, const T& val, const USI& destj)
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

template <typename T> void OCP_Table<T>::Display() const
{
    for (USI i = 0; i < nRow; i++) {
        for (USI j = 0; j < nCol; j++) {
            cout << data[j][i] << "\t";
        }
        cout << "\n";
    }
}

#endif /* end if __OCP_TABLE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/