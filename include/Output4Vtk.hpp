/*! \file    Output4Vtk.hpp
 *  \brief   Output reservoir information in vtk format
 *  \author  Shizhe Li
 *  \date    Oct/19/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OUTPUT4VTK_HEADER__
#define __OUTPUT4VTK_HEADER__

#include <fstream>
#include <string>
#include <vector>

#include "Grid.hpp"

using namespace std;

// Basic datatype
typedef OCP_DBL VTK_DBL;
typedef OCP_SIN VTK_SIN;
typedef OCP_USI VTK_USI;
typedef OCP_ULL VTK_ULL;

// Basic Keyword
const string  VTK_HEADER            = "# vtk DataFile Version 3.0";
const string  VTK_ASCII             = "ASCII";
const string  VTK_DATASET           = "DATASET";
const string  VTK_UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";
const string  VTK_POINTS            = "POINTS";
const string  VTK_CELLS             = "CELLS";
const string  VTK_CELL_TYPES        = "CELL_TYPES";
const string  VTK_CELL_DATA         = "CELL_DATA";
const string  VTK_POINT_DATA        = "POINT_DATA";
const VTK_USI VTK_MAX_TITLE_LENGTH  = 256;
const string  VTK_LOOKUP_TABLE      = "LOOKUP_TABLE";
const string  VTK_DEFAULT           = "default";
const string  VTK_SCALARS           = "SCALARS";

// Basic Cell Type
const VTK_USI VTK_HEXAHEDRON = 12;
const VTK_USI VTK_POLY_LINE  = 4;

const string VTK_FLOAT        = "float";
const string VTK_UNSIGNED_INT = "unsigned_int";

class Output4Vtk
{
    friend class Out4VTK;

public:
    void Init(const string&  myFile,
              const string&  shortInfo,
              const string&  myCodeWay,
              const string&  girdType,
              const VTK_USI& nG,
              const VTK_USI& nW); ///< create a new file and write basic information
    void OutputPOINTS(const string&                myFile,
                      const vector<OCPpolyhedron>& myHexGrid,
                      const vector<OCPpolyhedron>& myHexWell,
                      const string&                dataType) const;
    void OutputCELLS(const string&                myFile,
                     const vector<OCPpolyhedron>& myHexGrid,
                     const vector<OCPpolyhedron>& myHexWell) const;
    void OutputCELL_TYPES(const string&                myFile,
                          const vector<OCPpolyhedron>& myHex,
                          const vector<OCPpolyhedron>& myHexWell) const;
    template <typename T>
    void OutputCELL_DATA_SCALARS(const string&          myFile,
                                 const string&          dataName,
                                 const string&          dataType,
                                 const T*               gridVal,
                                 const USI&             gap,
                                 const vector<GB_Pair>& gbPair,
                                 const bool&            useActive,
                                 const T*               wellVal) const;
    void BeginCellData() const { cellData = true; };

private:
    mutable bool cellData{false};
    VTK_USI      numGrid;
    VTK_USI      numWell;
    VTK_USI      numCell;
};

template <typename T>
void Output4Vtk::OutputCELL_DATA_SCALARS(const string&          myFile,
                                         const string&          dataName,
                                         const string&          dataType,
                                         const T*               gridVal,
                                         const USI&             gap,
                                         const vector<GB_Pair>& gbPair,
                                         const bool&            useActive,
                                         const T*               wellVal) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open file: " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    if (cellData) {
        myVtk << VTK_CELL_DATA << " " << numCell << "\n";
        cellData = false;
    }

    myVtk << VTK_SCALARS << " " << dataName << " " << dataType << " " << 1 << "\n";
    myVtk << VTK_LOOKUP_TABLE << " " << VTK_DEFAULT << "\n";

    // Grid
    if (useActive) {
        for (VTK_USI n = 0; n < numGrid; n++) {
            if (gbPair[n].IsAct()) {
                myVtk << gridVal[gbPair[n].GetId() * gap] << "\n";
            } else {
                myVtk << 0 << "\n";
            }
        }
    } else {
        for (VTK_USI n = 0; n < numGrid; n++) {
            myVtk << gridVal[n * gap] << "\n";
        }
    }

    // Well
    for (USI w = 0; w < numWell; w++) {
        myVtk << wellVal[w] << "\n";
    }

    myVtk << "\n";
    myVtk.close();
}

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
