/*! \file    Output4Vtk.cpp
 *  \brief   Output reservoir information in vtk format
 *  \author  Shizhe Li
 *  \date    Oct/19/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Output4Vtk.hpp"

void Output4Vtk::Init(const string&  myFile,
                      const string&  shortInfo,
                      const string&  myCodeWay,
                      const string&  girdType,
                      const VTK_USI& nG,
                      const VTK_USI& nW)
{
    ofstream myVtk(myFile);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open file: " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    // Header
    myVtk << VTK_HEADER << "\n";

    // Title
    if (shortInfo.size() > VTK_MAX_TITLE_LENGTH) {
        OCP_WARNING("Length of title is beyond the limit: 256");
        myVtk << "Invalid short info, Too many characters\n";
    } else {
        myVtk << shortInfo << "\n";
    }

    // Code
    myVtk << myCodeWay << "\n";

    // Grid Type
    myVtk << VTK_DATASET << " " << girdType << "\n";

    myVtk << "\n";
    myVtk.close();

    // Init numGrid and numGrid
    numGrid = nG;
    numWell = nW;
    numCell = numGrid + numWell;
}

void Output4Vtk::OutputPOINTS(const string&                myFile,
                              const vector<OCPpolyhedron>& myHexGrid,
                              const vector<OCPpolyhedron>& myHexWell,
                              const string&                dataType) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open file: " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    // Grid Points
    VTK_USI numWellPoints = numGrid * 8;
    for (auto& w : myHexWell) {
        numWellPoints += w.numPoints;
    }
    myVtk << VTK_POINTS << " " << numWellPoints << " " << VTK_FLOAT << "\n";

    for (VTK_USI i = 0; i < numGrid; i++) {
        for (USI j = 0; j < myHexGrid[i].numPoints; j++) {
            myVtk << setw(6) << myHexGrid[i].Points[j].x << "   " << setw(6)
                  << myHexGrid[i].Points[j].y << "   " << setw(6)
                  << myHexGrid[i].Points[j].z << "\n";
        }
    }

    // Well Points
    for (VTK_USI w = 0; w < numWell; w++) {
        for (USI j = 0; j < myHexWell[w].numPoints; j++) {
            myVtk << setw(6) << myHexWell[w].Points[j].x << "   " << setw(6)
                  << myHexWell[w].Points[j].y << "   " << setw(6)
                  << myHexWell[w].Points[j].z << "\n";
        }
    }

    myVtk << "\n";
    myVtk.close();
}

void Output4Vtk::OutputCELLS(const string&                myFile,
                             const vector<OCPpolyhedron>& myHexGrid,
                             const vector<OCPpolyhedron>& myHexWell) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open file: " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    USI numSize = numCell;
    for (VTK_USI n = 0; n < numGrid; n++) {
        numSize += myHexGrid[n].numPoints;
    }
    for (USI w = 0; w < numWell; w++) {
        numSize += myHexWell[w].numPoints;
    }

    myVtk << VTK_CELLS << " " << numCell << " " << numSize << "\n";

    // EASY output!
    // Grid Cell
    VTK_USI tmp = 0;
    for (VTK_USI n = 0; n < numGrid; n++) {
        myVtk << 8 << "  ";
        for (VTK_USI j = 0; j < 8; j++) {
            myVtk << j + tmp << "   ";
        }
        myVtk << "\n";
        tmp += 8;
    }

    // Well Cell
    for (USI w = 0; w < numWell; w++) {
        myVtk << myHexWell[w].numPoints << "  ";
        for (VTK_USI j = 0; j < myHexWell[w].numPoints; j++) {
            myVtk << j + tmp << "   ";
        }
        myVtk << "\n";
        tmp += myHexWell[w].numPoints;
    }

    myVtk << "\n";
    myVtk.close();
}

void Output4Vtk::OutputCELL_TYPES(const string&                myFile,
                                  const vector<OCPpolyhedron>& myHexGrid,
                                  const vector<OCPpolyhedron>& myHexWell) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open file: " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    myVtk << VTK_CELL_TYPES << " " << numCell << "\n";

    // EASY output!
    // Grid Cell
    for (VTK_USI n = 0; n < numGrid; n++) {
        myVtk << VTK_HEXAHEDRON << "\n";
    }

    // Well Cell
    for (USI w = 0; w < numWell; w++) {
        myVtk << VTK_POLY_LINE << "\n";
    }

    myVtk << "\n";
    myVtk.close();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*  Chensong Zhang      Feb/05/2023      Update output in vtk files           */
/*----------------------------------------------------------------------------*/
