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

    // Header
    myVtk << VTK_HEADER << endl;

    // Title
    if (shortInfo.size() > VTK_MAX_TITLE_LENGTH) {
        OCP_WARNING("length of title is beyond the limit: 256");
        myVtk << "Invalid short info, Too many characters" << endl;
    }
    else {
        myVtk << shortInfo << endl;
    }
    
    // Code
    myVtk << myCodeWay << endl;

    // Grid Type
    myVtk << VTK_DATASET << " " << girdType << endl;

    myVtk.close();
}








/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/19/2022      Create file                          */
/*----------------------------------------------------------------------------*/