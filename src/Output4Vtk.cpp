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




void Output4Vtk::Init(const string& myFile, const string& shortInfo, const string& myCodeWay, const string& girdType) const
{
    ofstream myVtk(myFile);
    if (!myVtk.is_open()) {
        OCP_ABORT("Fail to Open " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    // Header
    myVtk << VTK_HEADER << "\n";

    // Title
    if (shortInfo.size() > VTK_MAX_TITLE_LENGTH) {
        OCP_WARNING("length of title is beyond the limit: 256");
        myVtk << "Invalid short info, Too many characters" << "\n";
    }
    else {
        myVtk << shortInfo << "\n";
    }
    
    // Code
    myVtk << myCodeWay << "\n";

    // Grid Type
    myVtk << VTK_DATASET << " " << girdType << "\n";

    myVtk << "\n";
    myVtk.close();
}


void Output4Vtk::OutputPOINTS(const string& myFile, const vector<OCPpolyhedron>& myHex, const string& dataType) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    const VTK_USI numGrid = myHex.size();
    myVtk << VTK_POINTS << " " << numGrid * 8 << " " << VTK_FLOAT << "\n";

    for (VTK_USI i = 0; i < numGrid; i++) {
        for (USI j = 0; j < myHex[i].numPoints; j++) {
            myVtk << setw(6) << myHex[i].Points[j].x << "   "
                << setw(6) << myHex[i].Points[j].y << "   "
                << setw(6) << myHex[i].Points[j].z << "\n";
        }
    }

    myVtk << "\n";
    myVtk.close();
}


void Output4Vtk::OutputCELLS(const string& myFile, const vector<OCPpolyhedron>& myHex) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    const VTK_USI numGrid = myHex.size();
    myVtk << VTK_CELLS << " " << numGrid << " " << numGrid + numGrid * 8 << "\n";
    
    // EASY output!
    VTK_USI tmp = 0;
    for (VTK_USI i = 0; i < numGrid; i++) {
        myVtk << 8 << "  ";
        for (VTK_USI j = 0; j < 8; j++) {
            myVtk << j + tmp << "   ";
        }
        myVtk << "\n";
        tmp += 8;
    }

    myVtk << "\n";
    myVtk.close();
}


void Output4Vtk::OutputCELL_TYPES(const string& myFile, const vector<OCPpolyhedron>& myHex) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);

    const VTK_USI numGrid = myHex.size();
    myVtk << VTK_CELL_TYPES << " " << numGrid << "\n";

    // EASY output!
    for (VTK_USI i = 0; i < numGrid; i++) {
        myVtk << VTK_HEXAHEDRON << "\n";
    }

    myVtk << "\n";
    myVtk.close();
}


void Output4Vtk::OutputCELL_DATA_SCALARS(const string& myFile, const string& dataName, const string& dataType,
    const VTK_DBL* val, const USI& gap, const vector<GB_Pair>& gbPair) const
{
    ofstream myVtk;
    myVtk.open(myFile, ios::app);
    if (!myVtk.is_open()) {
        OCP_ABORT("Can not open " + myFile);
    }

    ios::sync_with_stdio(false);
    myVtk.tie(0);


    const VTK_USI numGrid = gbPair.size();
    if (cellData) {
        myVtk << VTK_CELL_DATA << " " << numGrid << "\n";
        cellData = false;
    }  
    
    myVtk << VTK_SCALARS << " " << dataName << " " << dataType << " " << 1 << "\n";
    myVtk << VTK_LOOKUP_TABLE << " " << VTK_DEFAULT << "\n";

    
    for (VTK_USI n = 0; n < numGrid; n++) {
        if (gbPair[n].IsAct()) {
            myVtk << val[gbPair[n].GetId() * gap] << "\n";
        }
        else {
            myVtk << val[0] << "\n"; //tmp
        }
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
/*----------------------------------------------------------------------------*/