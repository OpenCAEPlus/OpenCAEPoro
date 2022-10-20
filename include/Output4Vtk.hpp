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
#include "Grid.hpp"
#include "UtilError.hpp"

using namespace std;


// Basic datatype
typedef OCP_DBL VTK_DBL;
typedef OCP_SIN VTK_SIN;
typedef OCP_USI VTK_USI;
typedef OCP_ULL VTK_ULL;


// Basic Keyword
const string VTK_HEADER = "# vtk DataFile Version 3.0";
const string VTK_ASCII = "ASCII";
const string VTK_DATASET = "DATASET";
const string VTK_UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";
const string VTK_POINTS = "POINTS";
const string VTK_CELLS = "CELLS";
const string VTK_CELL_TYPES = "CELL_TYPES";
const string VTK_CELL_DATA = "CELL_DATA";
const string VTK_POINT_DATA = "POINT_DATA";
const VTK_USI    VTK_MAX_TITLE_LENGTH = 256;
const string VTK_LOOKUP_TABLE = "LOOKUP_TABLE";
const string VTK_DEFAULT = "default";
const string VTK_SCALARS = "SCALARS";

// Basic Cell Type
const VTK_USI VTK_HEXAHEDRON = 12;

const string VTK_FLOAT = "float";

class Output4Vtk
{
	friend class Out4VTK;

public:
	void Init(const string& myFile, const string& shortInfo, const string& myCodeWay, const string& girdType) const; ///< create a new file and write basic information
	void OutputPOINTS(const string& myFile, const vector<OCPpolyhedron>& myHex, const string& dataType) const;
	void OutputCELLS(const string& myFile, const vector<OCPpolyhedron>& myHex) const;
	void OutputCELL_TYPES(const string& myFile, const vector<OCPpolyhedron>& myHex) const;
	void OutputPOINT_DATA();
	void OutputCELL_DATA_SCALARS(const string& myFile, const string& dataName, const string& dataType,
		const vector<VTK_DBL>& val, const vector<GB_Pair>& gbPair) const;

};



#endif


 /*----------------------------------------------------------------------------*/
 /*  Brief Change History of This File                                         */
 /*----------------------------------------------------------------------------*/
 /*  Author              Date             Actions                              */
 /*----------------------------------------------------------------------------*/
 /*  Shizhe Li           Oct/19/2022      Create file                          */
 /*----------------------------------------------------------------------------*/