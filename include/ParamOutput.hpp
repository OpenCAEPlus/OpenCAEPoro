/*! \file    ParamOutput.hpp
 *  \brief   ParamOutput class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMOUTPUT_HEADER__
#define __PARAMOUTPUT_HEADER__

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "ReadTool.hpp"

/// A structure of three-dimensional coordinates.
class COOIJK
{
public:
    COOIJK() = default;
    COOIJK(const USI& i, const USI& j, const USI& k)
        : I(i)
        , J(j)
        , K(k){};
    USI I;
    USI J;
    USI K;
};

/// Used to stores the contents of keyword whose contents are in form of coordinates.
class Type_B_o
{
public:
    bool           activity{false};
    vector<COOIJK> obj;
};

/// Used to stores the contents of keyword whose contents are in form of string.
class Type_A_o
{
public:
    bool           activity{false};
    vector<string> obj;
};

/// OutputSummary contains all information about SUMMARY, these information tells
/// which results will be print to the summary output file.
class OutputSummary
{

public:
    bool FPR{false};  ///< Field average Pressure.
    bool FOPR{false}; ///< Field oil production rate.
    bool FOPT{false}; ///< Field total oil production.
    bool FGPR{false}; ///< Field gas production rate.
    bool FGPt{false}; ///< Field total gas production.
    bool FWPR{false}; ///< Field water production rate.
    bool FWPT{false}; ///< Field total water production.
    bool FGIR{false}; ///< Field gas injection rate.
    bool FGIT{false}; ///< Field total gas injection.
    bool FWIR{false}; ///< Field water injection rate.
    bool FWIT{false}; ///< Field total water injection.

    Type_A_o WOPR; ///< Well oil production rate.
    Type_A_o WOPT; ///< Well total oil production rate.
    Type_A_o WGPR; ///< Well gas production rate.
    Type_A_o WGPT; ///< Well total gas production.
    Type_A_o WWPR; ///< Well water production rate.
    Type_A_o WWPT; ///< Well total water production.
    Type_A_o WGIR; ///< Well gas injection rate.
    Type_A_o WGIT; ///< Well total gas injection.
    Type_A_o WWIR; ///< Well water injection rate.
    Type_A_o WWIT; ///< Well total water injection.
    Type_A_o WBHP; ///< Well pressure.

    Type_B_o BPR; ///< Bulk pressure.
};

/// OutputDetail is a part of ParamOutput, it's used to control the output of detailed
/// information of reservoir. For example, pressure in every grid, due to excessive data
/// volume, this option is usually used in Debug Model.
class OutputDetail
{
public:
    bool    PRE{ false };   ///< Pressure of grids.
    bool    PGAS{ false };  ///< Gas pressure of grids.
    bool    PWAT{ false };  ///< Water pressure of grids.
    bool    SOIL{ false };  ///< Oil saturation of grids.
    bool    SGAS{ false };  ///< Gas saturation of grids.
    bool    SWAT{ false };  ///< Water saturation of grids.
    bool    DENO{ false };  ///< Oil density saturation of grids.
    bool    DENG{ false };  ///< Gas density saturation of grids.
    bool    DENW{ false };  ///< Water density saturation of grids.
};


/// ParamOutput is an internal structure used to stores the information of outputting
/// from input files. It is an intermediate interface and independent of the main
/// simulator. After all file inputting finishs, the params in it will pass to
/// corresponding modules.
class ParamOutput
{
public:
    OutputSummary summary; ///< See OutputSummary.
    OutputDetail detailInfo;    ///< See OutputDetail.

    /// Input the keyword SUMMARY, which contains many subkeyword, indicating which
    /// results are interested by user. After the simulation, these results will be
    /// output into a summary file.
    void InputSUMMARY(ifstream& ifs);
    /// Input the subkeyword in SUMMARY, the contents in these keyword is in the form of
    /// string.
    void InputType_A(ifstream& ifs, Type_A_o& obj);
    /// Input the subkeyword in SUMMARY, the contents in these keyword is in the form of
    /// coordinates.
    void InputType_B(ifstream& ifs, Type_B_o& obj);

    /// Input the keyword RPTSCHED, which tells which detailed information will be output
    /// to the RPTfile.
    void InputRPTSCHED(ifstream& ifs);
};

#endif /* end if __PARAMOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/