/*! \file    ParamRead.cpp
 *  \brief   ParamRead class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamRead.hpp"
#include "UtilError.hpp"

/// Initialize param_Rs, param_Well, and param_Control.
void ParamRead::Init()
{
    param_Rs.Init();
    param_Well.Init();
    param_Control.Init(workDir);
}

/// Get workDir and fileName from inputFile.
void ParamRead::GetDirAndName()
{
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)
    // for Window file system
    OCP_INT pos = inputFile.find_last_of('\\') + 1;
    workDir     = inputFile.substr(0, pos);
    fileName    = inputFile.substr(pos, inputFile.size() - pos);
#else
    // for Linux and Mac OSX file system
    OCP_INT pos = inputFile.find_last_of('/') + 1;
    workDir     = inputFile.substr(0, pos);
    fileName    = inputFile.substr(pos, inputFile.size() - pos);
#endif
}

/// The general interface for reading input file.
void ParamRead::ReadInputFile(const string& file)
{
    inputFile = file;
    GetDirAndName();
    Init();
    ReadFile(inputFile);
    CheckParam();
}

// TODO: What's this?
void ParamRead::ReadFile(const string& file)
{
    ifstream ifs(file, ios::in);
    if (!ifs) {
        OCP_MESSAGE("Trying to open input file: " << (file));
        OCP_ABORT("Cannot open the input file!");
    }

    while (!ifs.eof()) {
        vector<string> vbuf;
        if (!ReadLine(ifs, vbuf)) break;
        string keyword = vbuf[0];

        switch (Map_Str2Int(&keyword[0], keyword.size())) {
            case Map_Str2Int("BLACKOIL", 8):
                param_Rs.blackOil = true;
                break;

            case Map_Str2Int("COMPS", 5):
                param_Rs.InputCOMPS(ifs);
                break;

            case Map_Str2Int("OIL", 3):
                param_Rs.oil = true;
                break;

            case Map_Str2Int("GAS", 3):
                param_Rs.gas = true;
                break;

            case Map_Str2Int("WATER", 5):
                param_Rs.water = true;
                break;

            case Map_Str2Int("DISGAS", 6):
                param_Rs.disGas = true;
                break;

            case Map_Str2Int("DIMENS", 6):
                param_Rs.InputDIMENS(ifs);
                param_Rs.DisplayDIMENS();
                break;

            case Map_Str2Int("RTEMP", 5):
                param_Rs.InputRTEMP(ifs);
                break;

            case Map_Str2Int("EQUALS", 6):
                param_Rs.InputEQUALS(ifs);
                break;

            case Map_Str2Int("DX", 2):
            case Map_Str2Int("DY", 2):
            case Map_Str2Int("DZ", 2):
            case Map_Str2Int("NTG", 3):
            case Map_Str2Int("PORO", 4):
            case Map_Str2Int("TOPS", 4):
            case Map_Str2Int("PERMX", 5):
            case Map_Str2Int("PERMY", 5):
            case Map_Str2Int("PERMZ", 5):
            case Map_Str2Int("PRESSURE", 8):
            case Map_Str2Int("Ni", 2):
                param_Rs.InputGRID(ifs, keyword);
                break;

            case Map_Str2Int("COPY", 4):
                param_Rs.InputCOPY(ifs);
                break;

            case Map_Str2Int("MULTIPLY", 8):
                param_Rs.InputMULTIPLY(ifs);
                break;

            case Map_Str2Int("SWOF", 4):
            case Map_Str2Int("SGOF", 4):
            case Map_Str2Int("PVCO", 4):
            case Map_Str2Int("PVDG", 4):
            case Map_Str2Int("PVTW", 4):
            case Map_Str2Int("PBVD", 4):
                param_Rs.InputTABLE(ifs, keyword);
                break;

            case Map_Str2Int("ROCK", 4):
                param_Rs.InputROCK(ifs);
                break;

            case Map_Str2Int("GRAVITY", 7):
                param_Rs.InputGRAVITY(ifs);
                break;

            case Map_Str2Int("DENSITY", 7):
                param_Rs.InputDENSITY(ifs);
                break;

            case Map_Str2Int("EQUIL", 5):
                param_Rs.InputEQUIL(ifs);
                break;

            case Map_Str2Int("TABDIMS", 7):
                param_Rs.InputTABDIMS(ifs);
                break;

            case Map_Str2Int("SATNUM", 6):
            case Map_Str2Int("PVTNUM", 6):
                param_Rs.InputRegion(ifs, keyword);
                break;

            case Map_Str2Int("INCLUDE", 7):
                ReadINCLUDE(ifs);
                break;

            case Map_Str2Int("METHOD", 6):
                param_Control.InputMETHOD(ifs);
                break;

            case Map_Str2Int("TUNING", 6):
                param_Control.InputTUNING(ifs);
                break;

            case Map_Str2Int("WELSPECS", 8):
                param_Well.InputWELSPECS(ifs);
                break;

            case Map_Str2Int("COMPDAT", 7):
                param_Well.InputCOMPDAT(ifs);
                break;

            case Map_Str2Int("WCONINJE", 8):
                param_Well.InputWCONINJE(ifs);
                break;

            case Map_Str2Int("WCONPROD", 8):
                param_Well.InputWCONPROD(ifs);
                break;

            case Map_Str2Int("TSTEP", 5):
                param_Well.InputTSTEP(ifs);
                param_Control.criticalTime = param_Well.criticalTime;
                break;

            case Map_Str2Int("WELTARG", 7):
            case Map_Str2Int("WELLTARG", 8):
                param_Well.InputWELTARG(ifs);
                break;

            case Map_Str2Int("SUMMARY", 7):
                param_Output.InputSUMMARY(ifs);
                break;

            case Map_Str2Int("RPTSCHED", 8):
                param_Output.InputRPTSCHED(ifs);
                break;

            default:
                OCP_WARNING("Skipping unknonw keywords!");
                break;
        }
    }

    ifs.close();
}

/// Read INCLUDE files; these files should have identical format.
void ParamRead::ReadINCLUDE(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    ReadFile(workDir + vbuf[0]);
}

/// Check param_Rs and param_Well parameters.
void ParamRead::CheckParam()
{
    param_Rs.CheckParam();
    param_Well.CheckParam();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/