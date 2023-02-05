/*! \file    ParamOutput.cpp
 *  \brief   ParamOutput class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamOutput.hpp"

void ParamOutput::InputSUMMARY(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;
        string keyword = vbuf[0];

        switch (Map_Str2Int(&keyword[0], keyword.size())) {
            case Map_Str2Int("FPR", 3):
                summary.FPR = OCP_TRUE;
                break;

            case Map_Str2Int("FTR", 3):
                summary.FTR = OCP_TRUE;
                break;

            // Field
            case Map_Str2Int("FOPR", 4):
                summary.FOPR = OCP_TRUE;
                break;

            case Map_Str2Int("FOPT", 4):
                summary.FOPT = OCP_TRUE;
                break;

            case Map_Str2Int("FGPR", 4):
                summary.FGPR = OCP_TRUE;
                break;

            case Map_Str2Int("FGPT", 4):
                summary.FGPt = OCP_TRUE;
                break;

            case Map_Str2Int("FWPR", 4):
                summary.FWPR = OCP_TRUE;
                break;

            case Map_Str2Int("FWPT", 4):
                summary.FWPT = OCP_TRUE;
                break;

            case Map_Str2Int("FGIR", 4):
                summary.FGIR = OCP_TRUE;
                break;

            case Map_Str2Int("FGIT", 4):
                summary.FGIT = OCP_TRUE;
                break;

            case Map_Str2Int("FWIR", 4):
                summary.FWIR = OCP_TRUE;
                break;

            case Map_Str2Int("FWIT", 4):
                summary.FWIT = OCP_TRUE;
                break;

            // Well
            case Map_Str2Int("WOPR", 4):
                InputType_A(ifs, summary.WOPR);
                break;

            case Map_Str2Int("WOPT", 4):
                InputType_A(ifs, summary.WOPT);
                break;

            case Map_Str2Int("WGPR", 4):
                InputType_A(ifs, summary.WGPR);
                break;

            case Map_Str2Int("WGPT", 4):
                InputType_A(ifs, summary.WGPT);
                break;

            case Map_Str2Int("WWPR", 4):
                InputType_A(ifs, summary.WWPR);
                break;

            case Map_Str2Int("WWPT", 4):
                InputType_A(ifs, summary.WWPT);
                break;

            case Map_Str2Int("WGIR", 4):
                InputType_A(ifs, summary.WGIR);
                break;

            case Map_Str2Int("WGIT", 4):
                InputType_A(ifs, summary.WGIT);
                break;

            case Map_Str2Int("WWIR", 4):
                InputType_A(ifs, summary.WWIR);
                break;

            case Map_Str2Int("WWIT", 4):
                InputType_A(ifs, summary.WWIT);
                break;

            case Map_Str2Int("WBHP", 4):
                InputType_A(ifs, summary.WBHP);
                break;

            case Map_Str2Int("DG", 2):
                InputType_A(ifs, summary.DG);
                break;

            case Map_Str2Int("BPR", 3):
                InputType_B(ifs, summary.BPR);
                break;

            case Map_Str2Int("SOIL", 4):
                InputType_B(ifs, summary.SOIL);
                break;

            case Map_Str2Int("SGAS", 4):
                InputType_B(ifs, summary.SGAS);
                break;

            case Map_Str2Int("SWAT", 4):
                InputType_B(ifs, summary.SWAT);
                break;
            default:
                break;
        }
    }
    // cout << "SUMMARY" << endl;
}

void ParamOutput::InputType_A(ifstream& ifs, Type_A_o& obj)
{
    obj.activity = OCP_TRUE;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") {
        obj.obj.push_back("All");
    } else {
        OCP_INT len = vbuf.size();
        for (OCP_INT i = 0; i < len - 1; i++) {
            obj.obj.push_back(vbuf[i]);
        }
        if (vbuf.back() != "/") obj.obj.push_back(vbuf.back());

        while (ReadLine(ifs, vbuf)) {
            if (vbuf[0] == "/") break;

            OCP_INT len = vbuf.size();
            for (OCP_INT i = 0; i < len - 1; i++) {
                obj.obj.push_back(vbuf[i]);
            }
            if (vbuf.back() != "/") obj.obj.push_back(vbuf.back());
        }
    }
    // cout << "Type_A" << endl;
}

void ParamOutput::InputType_B(ifstream& ifs, Type_B_o& obj)
{

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        obj.activity = OCP_TRUE;
        DealDefault(vbuf);
        USI i = stoi(vbuf[0]);
        USI j = stoi(vbuf[1]);
        USI k = stoi(vbuf[2]);

        obj.obj.push_back(COOIJK(i, j, k));
    }
    // cout << "Type_B" << endl;
}

void ParamOutput::InputRPTSCHED(ifstream& ifs, const string& keyword)
{
    BasicGridPropertyParam* tmpBgpp;
    if (keyword == "RPTSCHED") {
        outRPTParam.useRPT = OCP_TRUE;
        tmpBgpp            = &outRPTParam.bgp;
    } else if (keyword == "VTKSCHED") {
        outVTKParam.useVTK = OCP_TRUE;
        tmpBgpp            = &outVTKParam.bgp;
    } else {
        return;
    }

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        USI len = vbuf.size();

        for (USI i = 0; i < len; i++) {

            string keyword = vbuf[i];

            switch (Map_Str2Int(&keyword[0], keyword.size())) {
                case Map_Str2Int("PRES", 4):
                case Map_Str2Int("PRESSURE", 8):
                    tmpBgpp->PRE = OCP_TRUE;
                    break;
                case Map_Str2Int("PGAS", 4):
                    tmpBgpp->PGAS = OCP_TRUE;
                    break;
                case Map_Str2Int("PWAT", 4):
                    tmpBgpp->PWAT = OCP_TRUE;
                    break;
                case Map_Str2Int("SOIL", 4):
                    tmpBgpp->SOIL = OCP_TRUE;
                    break;
                case Map_Str2Int("SGAS", 4):
                    tmpBgpp->SGAS = OCP_TRUE;
                    break;
                case Map_Str2Int("SWAT", 4):
                    tmpBgpp->SWAT = OCP_TRUE;
                    break;
                case Map_Str2Int("DENO", 4):
                    tmpBgpp->DENO = OCP_TRUE;
                    break;
                case Map_Str2Int("DENG", 4):
                    tmpBgpp->DENG = OCP_TRUE;
                    break;
                case Map_Str2Int("DENW", 4):
                    tmpBgpp->DENW = OCP_TRUE;
                    break;
                case Map_Str2Int("KRO", 3):
                    tmpBgpp->KRO = OCP_TRUE;
                    break;
                case Map_Str2Int("KRG", 3):
                    tmpBgpp->KRG = OCP_TRUE;
                    break;
                case Map_Str2Int("KRW", 3):
                    tmpBgpp->KRW = OCP_TRUE;
                    break;
                case Map_Str2Int("BOIL", 4):
                    tmpBgpp->BOIL = OCP_TRUE;
                    break;
                case Map_Str2Int("BGAS", 4):
                    tmpBgpp->BGAS = OCP_TRUE;
                    break;
                case Map_Str2Int("BWAT", 4):
                    tmpBgpp->BWAT = OCP_TRUE;
                    break;
                case Map_Str2Int("VOIL", 4):
                    tmpBgpp->VOIL = OCP_TRUE;
                    break;
                case Map_Str2Int("VGAS", 4):
                    tmpBgpp->VGAS = OCP_TRUE;
                    break;
                case Map_Str2Int("VWAT", 4):
                    tmpBgpp->VWAT = OCP_TRUE;
                    break;
                case Map_Str2Int("XMF", 3):
                    tmpBgpp->XMF = OCP_TRUE;
                    break;
                case Map_Str2Int("YMF", 3):
                    tmpBgpp->YMF = OCP_TRUE;
                    break;
                case Map_Str2Int("PCW", 3):
                    tmpBgpp->PCW = OCP_TRUE;
                    break;
                default:
                    break;
            }
        }
    }
    // cout << keyword << endl;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/