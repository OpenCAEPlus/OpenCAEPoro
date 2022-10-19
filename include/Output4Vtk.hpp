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

#include <vector>
#include <string>
#include <fstream>
#include "UtilError.hpp"

using namespace std;


// Basic datatype
typedef double VTK_DOUBLE;
typedef float VTK_FLOAT;
typedef unsigned int VTK_USI;
typedef unsigned long long VTK_ULL;


// Basic Keyword
const string VTK_HEADER = "# DataFile Version 3.0";
const string VTK_ASCII = "ASCII";
const string VTK_DATASET = "DATASET";
const string VTK_UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";
const string VTK_POINTS = "POINTS";
const string VTK_CELLS = "CELLS";
const string VTK_CELL_TYPES = "CELL_TYPES";
const string VTK_CELL_DATA = "CELL_DATA";
const string VTK_POINT_DATA = "POINT_DATA";
const VTK_USI    VTK_MAX_TITLE_LENGTH = 256;

class Output4Vtk
{
	friend class Out4VTK;

public:
	void Init(const string& myFile, const string& shortInfo, const string& myCodeWay, const string& girdType) const; ///< create a new file and write basic information
	void OutputPOINTS();
	void OutputCELLS();
	void OutputCELL_TYPES();
	void OutputPOINT_DATA();
	void OutputCELL_DATA();

};



#endif


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Oct/19/2022      Create file                          */
 /*----------------------------------------------------------------------------*/