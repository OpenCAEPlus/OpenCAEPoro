/*! \file    ParamWell.cpp
 *  \brief   ParamWell class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamWell.hpp"

WellOptParam::WellOptParam(string intype, vector<string>& vbuf)
{
    type = intype;
    if (type == "INJ") {
        fluidType = vbuf[1];
        state     = vbuf[2];
        optMode   = vbuf[3];
        if (vbuf[4] == "1*") {
            if (optMode == "BHP")
                maxRate = 1E10;
            else {
                cout << "### ERROR: Inj Rate is missing in WCONINJE!" << endl;
            }
        } else {
            maxRate = stod(vbuf[4]);
        }
        maxBHP = stod(vbuf[5]);
    } else if (type == "PROD") {
        state   = vbuf[1];
        optMode = vbuf[2];
        if (vbuf[3] == "1*") {
            if (optMode == "BHP")
                maxRate = 1E10;
            else {
                cout << "### ERROR: Prod Rate is missing in WCONINJE!" << endl;
            }
        } else {
            maxRate = stod(vbuf[3]);
        }
        minBHP = stod(vbuf[4]);
    } else {
        OCP_ABORT("Wrong Well Type!");
    }
}

WellParam::WellParam(vector<string>& info)
{
    name = info[0];
    if (info[1] != "DEFAULT") group = info[1];
    I = stoi(info[2]);
    J = stoi(info[3]);
    if (info[4] != "DEFAULT") depth = stod(info[4]);
}

Solvent::Solvent(const vector<string>& vbuf)
{
    name    = vbuf[0];
    USI len = vbuf.size();
    for (USI i = 0; i < len - 1; i++) {
        if (vbuf[i + 1] == "/") break;
        comRatio.push_back(stod(vbuf[i + 1]));
    }
}

void ParamWell::InputWELSPECS(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        well.push_back(WellParam(vbuf));
    }
    // cout << "WELSPECS" << endl;
}

void ParamWell::InputCOMPDAT(ifstream& ifs)
{
    USI            num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src   = vbuf[0];
        string::size_type pos   = src.find("*");
        OCP_BOOL          match = (pos != string::npos);
        if (match) {
            src.erase(pos);
        }
        OCP_BOOL tmp = OCP_FALSE;

        for (USI w = 0; w < num; w++) {
            if (match)
                tmp = (well[w].name.substr(0, pos) == src);
            else
                tmp = (well[w].name == src);

            if (tmp) {

                USI k1 = stoi(vbuf[3]);
                USI k2 = stoi(vbuf[4]);

                for (USI k = k1; k <= k2; k++) {
                    if (vbuf[1] == "DEFAULT" || vbuf[2] == "DEFAULT") {
                        well[w].I_perf.push_back(well[w].I);
                        well[w].J_perf.push_back(well[w].J);
                    } else {
                        well[w].I_perf.push_back(stoi(vbuf[1]));
                        well[w].J_perf.push_back(stoi(vbuf[2]));
                    }
                    well[w].K_perf.push_back(k);

                    if (vbuf[5] != "DEFAULT")
                        well[w].WI.push_back(stod(vbuf[5]));
                    else
                        well[w].WI.push_back(-1.0);

                    if (vbuf[6] != "DEFAULT")
                        well[w].diameter.push_back(stod(vbuf[6]));
                    else
                        well[w].diameter.push_back(1.0);

                    if (vbuf[7] != "DEFAULT")
                        well[w].kh.push_back(stod(vbuf[7]));
                    else
                        well[w].kh.push_back(-1.0);

                    if (vbuf[8] != "DEFAULT")
                        well[w].skinFactor.push_back(stod(vbuf[8]));
                    else
                        well[w].skinFactor.push_back(0.0);

                    if (vbuf[9] != "DEFAULT")
                        well[w].direction.push_back(vbuf[9]);
                    else
                        well[w].direction.push_back("z");
                }
            }
        }
    }
    // cout << "COMPDAT" << endl;
}

void ParamWell::InputWCONINJE(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++) {
                if (well[w].name.find(src) != string::npos) {
                    well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
                }
            }
        } else {
            for (USI w = 0; w < num; w++) {
                if (well[w].name == src) {
                    well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
                }
            }
        }
    }
    // cout << "WCONINJE" << endl;
}

void ParamWell::InputWCONPROD(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++)
                if (well[w].name.find(src) != string::npos)
                    well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
        } else {
            for (USI w = 0; w < num; w++)
                if (well[w].name == src)
                    well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
        }
    }
    // cout << "WCONPROD" << endl;
}

void ParamWell::InputTSTEP(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealDefault(vbuf);
        OCP_INT len = vbuf.size();
        for (OCP_INT i = 0; i < len - 1; i++) {
            OCP_DBL t = criticalTime.back() + stod(vbuf[i]);
            criticalTime.push_back(t);
        }
        if (vbuf.back() != "/") {
            OCP_DBL t = criticalTime.back() + stod(vbuf.back());
            criticalTime.push_back(t);
        }
    }
}

void ParamWell::InputWELTARG(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    cout << "\n---------------------" << endl
         << "WELTARG"
         << "\n---------------------" << endl;

    const USI      d   = criticalTime.size() - 1;
    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (const auto& s : vbuf) {
            cout << s << "   ";
        }
        cout << endl;

        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }
        OCP_BOOL succMatch = OCP_FALSE;

        for (USI w = 0; w < num; w++) {
            if (fuzzyMatch)
                succMatch = (well[w].name.find(src) != string::npos);
            else
                succMatch = (well[w].name == src);

            if (succMatch) {
                if (well[w].optParam.size() == 0) {
                    OCP_ABORT("No Well Control Defined in Advance!");
                }
                WellOptPair tar = well[w].optParam.back();
                tar.d           = d;
                tar.opt.optMode = vbuf[1];
                OCP_DBL val     = stod(vbuf[2]);
                if (vbuf[1] == "BHP") {
                    if (tar.opt.type == "INJ")
                        tar.opt.maxBHP = val;
                    else
                        tar.opt.minBHP = val;
                } else {
                    tar.opt.maxRate = val;
                }
                well[w].optParam.push_back(tar);
            }
        }
    }
}

void ParamWell::InputWTEMP(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    cout << "\n---------------------" << endl
         << "WTEMP"
         << "\n---------------------" << endl;

    const USI      d       = criticalTime.size() - 1;
    const USI      numWell = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (const auto& s : vbuf) {
            cout << s << "   ";
        }
        cout << endl;

        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }
        OCP_BOOL succMatch = OCP_FALSE;

        for (USI w = 0; w < numWell; w++) {
            if (fuzzyMatch)
                succMatch = (well[w].name.find(src) != string::npos);
            else
                succMatch = (well[w].name == src);

            if (succMatch) {
                if (well[w].optParam.size() == 0) {
                    OCP_ABORT("No Well Control Defined in Advance!");
                }
                WellOptPair tar = well[w].optParam.back();
                tar.d           = d;
                tar.opt.injTemp = stod(vbuf[1]);
                well[w].optParam.push_back(tar);
            }
        }
    }
}

void ParamWell::InputUNWEIGHT(ifstream& ifs)
{
    assert(criticalTime.size() > 0);

    cout << "\n---------------------" << endl
         << "UNWEIGHT"
         << "\n---------------------" << endl;

    const USI      num = well.size();
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (const auto& s : vbuf) {
            cout << s << "   ";
        }
        cout << endl;

        DealDefault(vbuf);
        string            src        = vbuf[0];
        string::size_type pos        = src.find("*");
        const OCP_BOOL    fuzzyMatch = (pos != string::npos);
        if (fuzzyMatch) {
            src.erase(pos);
        }

        if (fuzzyMatch) {
            for (USI w = 0; w < num; w++)
                if (well[w].name.find(src) != string::npos)
                    well[w].ifUseUnweight = OCP_TRUE;
        } else {
            for (USI w = 0; w < num; w++)
                if (well[w].name == src) well[w].ifUseUnweight = OCP_TRUE;
        }
    }
}

void ParamWell::InputWELLSTRE(ifstream& ifs)
{
    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;
        solSet.push_back(Solvent(vbuf));
    }
    // cout << "WELLSTRE" << endl;
}

void ParamWell::InputPSURF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    Psurf = stod(vbuf[0]);

    cout << "\n---------------------" << endl
         << "PSURF"
         << "\n---------------------" << endl;
    cout << Psurf << endl;
}

void ParamWell::InputTSURF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    Tsurf = stod(vbuf[0]);

    cout << "\n---------------------" << endl
         << "TSURF"
         << "\n---------------------" << endl;
    cout << Tsurf << endl;
}

// check
void ParamWell::CheckParam(const OCP_BOOL& boModel) const
{
    if (boModel) {
        CheckINJFluid();
    }
    CheckPerf();
}

void ParamWell::CheckINJFluid() const
{
    USI wellnum = well.size();
    USI len;
    for (USI w = 0; w < wellnum; w++) {
        len = well[w].optParam.size();
        for (USI i = 0; i < len; i++) {
            const WellOptParam& tmpOpt = well[w].optParam[i].opt;
            if (tmpOpt.type == "INJ") {
                if (tmpOpt.fluidType == "WAT" || tmpOpt.fluidType == "WATER") {
                } else if (tmpOpt.fluidType == "GAS") {
                } else {
                    OCP_ABORT("Wrong FluidType!");
                }
            }
        }
    }
}

void ParamWell::CheckPerf() const
{
    USI wellnum = well.size();
    USI perfnum;
    for (USI w = 0; w < wellnum; w++) {
        perfnum = well[w].I_perf.size();
        if (well[w].J_perf.size() != perfnum) {
            OCP_ABORT("Wrong perforation size J_perf!");
        }
        if (well[w].K_perf.size() != perfnum) {
            OCP_ABORT("Wrong perforation size K_perf!");
        }
        if (well[w].diameter.size() != perfnum) {
            OCP_ABORT("Wrong perforation diameter size!");
        }
        if (well[w].WI.size() != perfnum) {
            OCP_ABORT("Wrong perforation WI size!");
        }
        if (well[w].kh.size() != perfnum) {
            OCP_ABORT("Wrong perforation kh size!");
        }
        if (well[w].skinFactor.size() != perfnum) {
            OCP_ABORT("Wrong perforation skinFactor size!");
        }
        if (well[w].direction.size() != perfnum) {
            OCP_ABORT("Wrong perforation direction size!");
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Oct/27/2021      Unify error messages                 */
/*----------------------------------------------------------------------------*/