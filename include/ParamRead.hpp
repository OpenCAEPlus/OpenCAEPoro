/*! \file    ParamRead.hpp
 *  \brief   ParamRead class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMREAD_HEADER__
#define __PARAMREAD_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoro header files
#include "ParamControl.hpp"
#include "ParamOutput.hpp"
#include "ParamReservoir.hpp"
#include "ParamWell.hpp"
#include "ReadTool.hpp"

using namespace std;

/// ParamRead is pre-processing part in our simulator, which is responsible for
/// inputting params from files supplied by users, it's almost compatible with Eclipse
/// but has own rules for easy to use. it is extensible friendly. ParamRead contains
/// four main parts: read reservoir params, read well params, read control params and
/// read output param. These work is assigned to four class respectively.
class ParamRead
{
public:
    /// Input file, it must give the full path: absolute path or relative path.
    string inputFile;
    string workDir;  ///< Current work directory.
    string fileName; ///< File name of input file.

    ParamReservoir param_Rs;      ///< Read the reservoir params.
    ParamWell      param_Well;    ///< Read the well params.
    ParamControl   param_Control; ///< Read the control params.
    ParamOutput    param_Output;  ///< Read the output params.

    /// The general interface for reading input file.
    void ReadInputFile(const string& file);
    /// Get the current work dir and the input file name from the full file path.
    void GetDirAndName();
    /// Read the input file.
    void ReadFile(const string& file);

    /// Initialize the param inputting.
    void Init();
    /// Input the keyword INCLUDE. INCLUDE contains other input files needed to be read.
    /// All these files has identical format.
    void InputINCLUDE(ifstream& ifs);

    /// Check if the param input is right.
    void CheckParam();
};

#endif /* end if __PARAMREAD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/