/*! \file    ParamReservoir.cpp
 *  \brief   ParamReservoir class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamReservoir.hpp"

/// Find pointer to the specified variable.
vector<OCP_DBL>* ParamReservoir::FindPtr(const string& varName)
{
    vector<OCP_DBL>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("DX", 2):
            dx.reserve(numGrid);
            myPtr = &dx;
            break;

        case Map_Str2Int("DY", 2):
            dy.reserve(numGrid);
            myPtr = &dy;
            break;

        case Map_Str2Int("DZ", 2):
            dz.reserve(numGrid);
            myPtr = &dz;
            break;

        case Map_Str2Int("COORD", 5):
            coord.reserve((dimens.nx + 1) * (dimens.ny + 1) * 6);
            myPtr = &coord;
            break;

        case Map_Str2Int("ZCORN", 5):
            zcorn.reserve(numGrid * 8);
            myPtr = &zcorn;
            break;

        case Map_Str2Int("PORO", 4):
            poro.reserve(numGrid);
            myPtr = &poro;
            break;

        case Map_Str2Int("NTG", 3):
            ntg.reserve(numGrid);
            myPtr = &ntg;
            break;

        case Map_Str2Int("PERMX", 5):
            permX.reserve(numGrid);
            myPtr = &permX;
            break;

        case Map_Str2Int("PERMY", 5):
            permY.reserve(numGrid);
            myPtr = &permY;
            break;

        case Map_Str2Int("PERMZ", 5):
            permZ.reserve(numGrid);
            myPtr = &permZ;
            break;

        case Map_Str2Int("THCONR", 6):
            thconr.reserve(numGrid);
            myPtr = &thconr;
            break;

        case Map_Str2Int("TOPS", 4):
            tops.reserve(dimens.nx * dimens.ny);
            myPtr = &tops;
            break;

        case Map_Str2Int("PRESSURE", 8):
            P.reserve(numGrid);
            myPtr = &P;
            break;

        case Map_Str2Int("Ni", 2):
            Ni.reserve(numGrid);
            myPtr = &Ni;
            break;

        case Map_Str2Int("SWATINIT", 8):
            Swat.reserve(numGrid);
            myPtr     = &Swat;
            ScalePcow = OCP_TRUE;
            break;

        case Map_Str2Int("SATNUM", 6):
            SATNUM.activity = OCP_TRUE;
            SATNUM.data.reserve(numGrid);
            myPtr = &SATNUM.data;
            break;

        case Map_Str2Int("PVTNUM", 6):
            PVTNUM.activity = OCP_TRUE;
            PVTNUM.data.reserve(numGrid);
            myPtr = &PVTNUM.data;
            break;

        case Map_Str2Int("ROCKNUM", 7):
            ROCKNUM.activity = OCP_TRUE;
            ROCKNUM.data.reserve(numGrid);
            myPtr = &ROCKNUM.data;
            break;
    }

    return myPtr;
}

/// Find pointer to the specified table.
TableSet* ParamReservoir::FindPtr_T(const string& varName)
{
    TableSet* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("SWFN", 4):
            myPtr = &SWFN_T;
            break;

        case Map_Str2Int("SWOF", 4):
            myPtr = &SWOF_T;
            break;

        case Map_Str2Int("SGFN", 4):
            myPtr = &SGFN_T;
            break;

        case Map_Str2Int("SGOF", 4):
            myPtr = &SGOF_T;
            break;

        case Map_Str2Int("SOF3", 4):
            myPtr = &SOF3_T;
            break;

        case Map_Str2Int("PBVD", 4):
            myPtr = &PBVD_T;
            break;

        case Map_Str2Int("PVCO", 4):
            myPtr = &PVCO_T;
            break;

        case Map_Str2Int("PVDO", 4):
            myPtr = &PVDO_T;
            break;

        case Map_Str2Int("PVDG", 4):
            myPtr = &PVDG_T;
            break;

        case Map_Str2Int("PVTW", 4):
            myPtr = &PVTW_T;
            break;

        case Map_Str2Int("ZMFVD", 5):
            myPtr = &ZMFVD_T;
            break;

        case Map_Str2Int("TEMPVD", 6):
            myPtr = &TEMPVD_T;
            break;
    }

    return myPtr;
}

/// Initialize tables and other reservoir parameters.
void ParamReservoir::Init()
{
    InitTable();

    gravity.data.resize(3);
    gravity.data[0] = 45.5;   // oil
    gravity.data[1] = 1.0;    // pure water
    gravity.data[2] = 0.7773; // air

    density.data.resize(3);
    density.data[0] = 37.457;    // The density of oil at surface conditions: lb/ft3
    density.data[1] = 62.366416; // The density of water at surface conditions: lb/ft3
    density.data[2] = 0.062428;  // The density of gas at surface conditions: lb/ft3

    rsTemp = 60.0;
}

/// Initialize tables.
void ParamReservoir::InitTable()
{
    SWFN_T.name     = "SWFN";
    SWFN_T.colNum   = 3;
    SWOF_T.name     = "SWOF";
    SWOF_T.colNum   = 4;
    SGFN_T.name     = "SGFN";
    SGFN_T.colNum   = 3;
    SGOF_T.name     = "SGOF";
    SGOF_T.colNum   = 4;
    SOF3_T.name     = "SOF3";
    SOF3_T.colNum   = 3;
    PBVD_T.name     = "PBVD";
    PBVD_T.colNum   = 2;
    PVCO_T.name     = "PVCO";
    PVCO_T.colNum   = 6;
    PVDO_T.name     = "PVDO";
    PVDO_T.colNum   = 3;
    PVDG_T.name     = "PVDG";
    PVDG_T.colNum   = 3;
    PVTW_T.name     = "PVTW";
    PVTW_T.colNum   = 5;
    ZMFVD_T.name    = "ZMFVD";  // colnum equals numCom(hydrocarbon) + 1
    TEMPVD_T.name   = "TEMPVD"; // colnum equals 2
    TEMPVD_T.colNum = 2;
}

/// TODO: Add Doxygen
template <typename T>
void ParamReservoir::setVal(vector<T>& obj, const T& val, const vector<USI>& index)
{
    USI     Nx   = dimens.nx;
    USI     Ny   = dimens.ny;
    OCP_USI NxNy = Nx * Ny;
    OCP_USI id   = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id      = k * NxNy + j * Nx + i;
                obj[id] = val;
            }
        }
    }
}

/// TODO: Add Doxygen
template <typename T>
void ParamReservoir::CopyVal(vector<T>&         obj,
                             const vector<T>&   src,
                             const vector<USI>& index)
{
    USI     Nx   = dimens.nx;
    USI     Ny   = dimens.ny;
    OCP_USI NxNy = Nx * Ny;
    OCP_USI id   = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id      = k * NxNy + j * Nx + i;
                obj[id] = src[id];
            }
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::MultiplyVal(vector<OCP_DBL>&   obj,
                                 const OCP_DBL&     val,
                                 const vector<USI>& index)
{
    USI     Nx   = dimens.nx;
    USI     Ny   = dimens.ny;
    OCP_USI NxNy = Nx * Ny;
    OCP_USI id   = 0;

    for (USI k = index[4]; k <= index[5]; k++) {
        for (USI j = index[2]; j <= index[3]; j++) {
            for (USI i = index[0]; i <= index[1]; i++) {
                id = k * NxNy + j * Nx + i;
                obj[id] *= val;
            }
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::InputCOMPS(ifstream& ifs)
{
    comps = OCP_TRUE;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    numCom           = stoi(vbuf[0]);
    comsParam.numCom = numCom;
    comsParam.Init();

    cout << endl << "COMPS" << endl;
    cout << numCom << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputDIMENS(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    dimens.nx = stoi(vbuf[0]);
    dimens.ny = stoi(vbuf[1]);
    dimens.nz = stoi(vbuf[2]);
    numGrid   = dimens.nx * dimens.ny * dimens.nz;

    DisplayDIMENS();
}

/// TODO: Add Doxygen
void ParamReservoir::DisplayDIMENS()
{
    cout << "\n---------------------" << endl
         << "DIMENS"
         << "\n---------------------" << endl;
    cout << "   " << dimens.nx << "  " << dimens.ny << "  " << dimens.nz << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputRTEMP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    rsTemp = stod(vbuf[0]);
    cout << "RTEMP\n" << rsTemp << endl << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputEQUALS(ifstream& ifs)
{
    cout << "\n---------------------" << endl
         << "EQUALS"
         << "\n---------------------" << endl;

    vector<USI>    index(6, 0);
    vector<string> vbuf;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto v : vbuf) {
            if (v != "/") cout << setw(10) << v;
        }
        cout << "\n";

        index[0] = 0, index[1] = dimens.nx - 1;
        index[2] = 0, index[3] = dimens.ny - 1;
        index[4] = 0, index[5] = dimens.nz - 1;

        string  objName = vbuf[0];
        OCP_DBL val     = stod(vbuf[1]);

        DealDefault(vbuf);

        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }
        if (index[0] < 0 || index[2] < 0 || index[4] < 0 || index[1] > dimens.nx - 1 ||
            index[3] > dimens.ny - 1 || index[5] > dimens.nz - 1) {
            OCP_ABORT("WRONG Range in " + objName + " in EQUALS!");
        }

        vector<OCP_DBL>* objPtr = FindPtr(objName);

        if (objPtr != nullptr) {
            if (objName == "TOPS") {
                objPtr->resize(dimens.nx * dimens.ny);
                index[4] = index[5] = 0;
            } else {
                objPtr->resize(numGrid);
            }
            setVal(*objPtr, val, index);
        } else {
            OCP_ABORT("Wrong object name: " + objName);
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::InputGRID(ifstream& ifs, string& keyword)
{
    vector<OCP_DBL>* objPtr = nullptr;

    objPtr = FindPtr(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto& str : vbuf) {
            // if m*n occurs, then push back n  m times
            auto pos = str.find('*');
            if (pos == string::npos) {
                objPtr->push_back(stod(str));
            } else {
                USI     len = str.size();
                OCP_USI num = stoi(str.substr(0, pos));
                OCP_DBL val = stod(str.substr(pos + 1, len - (pos + 1)));
                for (USI i = 0; i < num; i++) objPtr->push_back(val);
            }
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::InputCOPY(ifstream& ifs)
{
    cout << "\n---------------------" << endl
         << "COPY"
         << "\n---------------------" << endl;

    vector<string> vbuf;
    vector<USI>    index(6, 0);

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto v : vbuf) {
            if (v != "/") cout << setw(10) << v;
        }
        cout << "\n";

        index[0] = 0, index[1] = dimens.nx - 1;
        index[2] = 0, index[3] = dimens.ny - 1;
        index[4] = 0, index[5] = dimens.nz - 1;

        string srcName = vbuf[0];
        string objName = vbuf[1];
        DealDefault(vbuf);
        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }

        vector<OCP_DBL>* srcPtr = FindPtr(srcName);
        vector<OCP_DBL>* objPtr = FindPtr(objName);
        if (srcPtr != nullptr && objPtr != nullptr) {
            objPtr->resize(srcPtr->size());
            CopyVal(*objPtr, *srcPtr, index);
        } else {
            OCP_ABORT("Wrong object names: " + srcName + ", " + objName);
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::InputMULTIPLY(ifstream& ifs)
{
    vector<string> vbuf;
    vector<USI>    index(6, 0);

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        index[0] = 0, index[1] = dimens.nx - 1;
        index[2] = 0, index[3] = dimens.ny - 1;
        index[4] = 0, index[5] = dimens.nz - 1;

        string  objName = vbuf[0];
        OCP_DBL val     = stod(vbuf[1]);

        DealDefault(vbuf);
        for (USI n = 2; n < 8; n++) {
            if (vbuf[n] != "DEFAULT") index[n - 2] = stoi(vbuf[n]) - 1;
        }

        vector<OCP_DBL>* objPtr = FindPtr(objName);
        if (objPtr != nullptr) {
            if (objName == "TOPS") {
                index[4] = index[5] = 0;
            }
            MultiplyVal(*objPtr, val, index);
        } else {
            OCP_ABORT("Wrong object name: " + objName);
        }
    }
}

/// TODO: Add Doxygen
void ParamReservoir::InputTABLE(ifstream& ifs, const string& tabName)
{
    TableSet* obj;
    obj = FindPtr_T(tabName);
    if (obj == nullptr) {
        OCP_ABORT("Wrong table name :" + tabName);
    }

    USI col = obj->colNum;
    if (tabName == "ZMFVD") {
        if (!comps) OCP_ABORT("COMPS isn't set correctly!");
        obj->colNum = numCom + 1;
        col         = obj->colNum;
    }
    vector<vector<OCP_DBL>> tmpTab(col);

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (USI i = 0; i < col; i++) {
            tmpTab[i].push_back(stod(vbuf[i]));
        }

        if (vbuf.back() == "/") {
            obj->data.push_back(tmpTab);
            for (USI j = 0; j < col; j++) {
                tmpTab[j].clear();
            }
        }
    }
    if (!tmpTab[0].empty()) obj->data.push_back(tmpTab);

    obj->DisplayTable();
}

/// Read data from the ROCK keyword.
void ParamReservoir::InputROCK(ifstream& ifs)
{
    cout << "\n---------------------" << endl
         << "ROCK"
         << "\n---------------------" << endl;

    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        RockParam rock;
        rock.type = vbuf[0];
        rock.Pref = stod(vbuf[1]);
        rock.cp1  = stod(vbuf[2]);

        if (rock.type == "LINEAR02") {
            if (vbuf.size() > 3 && vbuf[3] != "/") {
                rock.cp2 = stod(vbuf[3]);
            } else {
                rock.cp2 = rock.cp1;
            }
        }
        rockSet.push_back(rock);

        cout << "   " << rock.type << "   " << rock.Pref << "   " << rock.cp1 << "   "
             << rock.cp2 << endl;
    }
}

/// Read data from the ROCK keyword.
void ParamReservoir::InputROCKT(ifstream& ifs)
{
    cout << "\n---------------------" << endl
         << "ROCKT"
         << "\n---------------------" << endl;

    RockParam      rock;
    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        USI index = 0;
        USI len   = vbuf.size();
        while (index < len) {
            if (vbuf[index] == "*PORFORM") {
                rock.type = vbuf[index + 1];
            } else if (vbuf[index] == "*PRPOR") {
                rock.Pref = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*TRPOR") {
                rock.Tref = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CPOR") {
                rock.cp1 = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CTPOR") {
                rock.ct = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CPTPOR") {
                rock.cpt = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*VOLCONST") {
                if (vbuf[index + 1] == "BULK") rock.ConstRock = OCP_FALSE;
            } else if (vbuf[index] == "*CP1") {
                rock.HCP1 = stod(vbuf[index + 1]);
            } else if (vbuf[index] == "*CP2") {
                rock.HCP2 = stod(vbuf[index + 1]);
            }
            index += 2;
        }
    }
    rockSet.push_back(rock);

    cout << "*PORFORM   " << rock.type << endl;
    cout << "*PRPOR     " << rock.Pref << endl;
    cout << "*TRPOR     " << rock.Tref << endl;
    cout << "*CPOR      " << rock.cp1 << endl;
    cout << "*CTPOR     " << rock.ct << endl;
    cout << "*CPTPOR    " << rock.cpt << endl;
    cout << "*VOLCONST  " << (rock.ConstRock ? "ROCK" : "BULK") << endl;
    cout << "*CP1       " << rock.HCP1 << endl;
    cout << "*CP2       " << rock.HCP2 << endl;
}

void ParamReservoir::InputHLOSS(ifstream& ifs)
{
    cout << "\n---------------------" << endl
         << "HLOSSPROR"
         << "\n---------------------" << endl;

    hLoss.ifHLoss = OCP_TRUE;

    vector<string> vbuf;
    while (true) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") break;

        USI index = 0;
        USI len   = vbuf.size();
        while (index < len) {
            if (vbuf[index] == "*OVERBUR") {
                hLoss.obC = stod(vbuf[index + 1]);
                hLoss.obK = stod(vbuf[index + 2]);
            } else if (vbuf[index] == "*UNDERBUR") {
                hLoss.ubC = stod(vbuf[index + 1]);
                hLoss.ubK = stod(vbuf[index + 2]);
            }
            index += 3;
        }
    }
    cout << "*OVERBUR   " << hLoss.obC << "   " << hLoss.obK << endl;
    cout << "*UNDERBUR  " << hLoss.ubC << "   " << hLoss.ubK << endl;
}

/// Read data from the MISCSTR keyword.
void ParamReservoir::InputMISCSTR(ifstream& ifs)
{
    if (!comsParam.miscible) {
        OCP_WARNING("MISCIBLE has not been declared. Keyword ignored!");
    } else {
        vector<string> vbuf;
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") return;
        if (vbuf.back() == "/") vbuf.pop_back();

        USI len = vbuf.size();
        for (USI i = 0; i < len; i++) {
            miscstr.surTenRef.push_back(stod(vbuf[i]));
        }
    }
    cout << "\n---------------------" << endl
         << "MISCSTR"
         << "\n---------------------" << endl;
    for (auto& v : miscstr.surTenRef) cout << v << "   ";
    cout << endl;
}

/// Read data from the GRAVITY keyword.
void ParamReservoir::InputGRAVITY(ifstream& ifs)
{
    gravity.activity = OCP_TRUE;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;
    DealDefault(vbuf);
    OCP_ASSERT(vbuf.size() == 4, "Wrong Keyword GRAVITY!");
    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            gravity.data[i] = stod(vbuf[i]);
        }
    }

    cout << "\n---------------------" << endl
         << "GRAVITY"
         << "\n---------------------" << endl;
    cout << "   " << gravity.data[0] << "  " << gravity.data[1] << "  "
         << gravity.data[2] << endl;
}

/// Read data from the DENSITY keyword.
void ParamReservoir::InputDENSITY(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    DealDefault(vbuf);
    OCP_ASSERT(vbuf.size() == 3, "Wrong Keyword DENSITY!");
    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            density.activity = OCP_TRUE;
            density.data[i]  = stod(vbuf[i]);
        }
    }

    cout << "\n---------------------" << endl
         << "DENSITY"
         << "\n---------------------" << endl;
    cout << density.data[0] << "  " << density.data[1] << "  " << density.data[2]
         << endl;
}

/// Read data from the THCONO, THCONG, THCONW
void ParamReservoir::InputTHCON(ifstream& ifs, const string& keyword)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (keyword == "THCONO") {
        thcono = stod(vbuf[0]);
    } else if (keyword == "THCONG") {
        thcong = stod(vbuf[0]);
    } else if (keyword == "THCONW") {
        thconw = stod(vbuf[0]);
    }

    cout << "THCONO\n" << thcono << endl << endl;
    cout << "THCONG\n" << thcong << endl << endl;
    cout << "THCONW\n" << thconw << endl << endl;
}

/// Read data from the EQUIL keyword.
void ParamReservoir::InputEQUIL(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    EQUIL.resize(6, 0);
    DealDefault(vbuf);
    for (USI i = 0; i < 6; i++) {
        if (vbuf[i] != "DEFAULT") EQUIL[i] = stod(vbuf[i]);
    }

    cout << "\n---------------------" << endl
         << "EQUIL"
         << "\n---------------------" << endl;
    cout << "   ";
    for (USI i = 0; i < 6; i++) cout << EQUIL[i] << "  ";
    cout << endl;
}

/// Read data from the TABDIMS keyword.
void ParamReservoir::InputTABDIMS(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);

    if (vbuf.size() < 3) {
        OCP_ABORT("Input the number of Saturation tables, PVT tables, and Rock tables "
                  "in turn!");
    }

    NTSFUN = stoi(vbuf[0]);
    NTPVT  = stoi(vbuf[1]);
    NTROOC = stoi(vbuf[2]);

    cout << "\n---------------------" << endl
         << "TABDIMS"
         << "\n---------------------" << endl;
    cout << "   " << NTSFUN << "   " << NTPVT << "   " << NTROOC << endl;
}

/// Region information like SATNUM to decide which grid belongs to which saturation
/// region, so corresponding saturation table will be used.
void ParamReservoir::InputRegion(ifstream& ifs, const string& keyword)
{
    Type_A_r<OCP_DBL>* ptr = &PVTNUM;
    USI                lim = NTPVT;

    if (keyword == "SATNUM") {
        ptr = &SATNUM;
        lim = NTSFUN;
    } else if (keyword == "ACTNUM") {
        ptr = &ACTNUM;
    } else if (keyword == "ROCKNUM") {
        ptr = &ROCKNUM;
    }

    ptr->activity = OCP_TRUE;
    ptr->data.reserve(numGrid);
    vector<string>  vbuf;
    vector<OCP_USI> obj;
    vector<USI>     region;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto& str : vbuf) {
            // if m*n occurs, then push back n m times
            auto pos = str.find('*');
            if (pos == string::npos) {
                ptr->data.push_back(stod(str));
            } else {
                USI     len = str.size();
                OCP_USI num = stoi(str.substr(0, pos));
                OCP_DBL val = stod(str.substr(pos + 1, len - (pos + 1)));
                for (USI i = 0; i < num; i++) ptr->data.push_back(val);
            }
        }
    }

    cout << "Number of Tables = " << lim << endl;
}

/// Check consistency of input parameters.
void ParamReservoir::CheckParam()
{
    CheckGrid();
    CheckEQUIL();
    CheckDenGra();
    CheckPhase();
    CheckRegion();
    CheckRock();
}

/// Check data dimension for potential problems.
void ParamReservoir::CheckGrid()
{
    if (coord.size() == 0) {
        if (tops.size() != dimens.nx * dimens.ny) OCP_ABORT("Wrong TOPS size!");
        if (dx.size() != numGrid) OCP_ABORT("Wrong DX size!");
        if (dy.size() != numGrid) OCP_ABORT("Wrong DY size!");
        if (dz.size() != numGrid) OCP_ABORT("Wrong DZ size!");
    } else {
        if (coord.size() != (dimens.nx + 1) * (dimens.ny + 1) * 6)
            OCP_ABORT("Wrong COORD size!");
        if (zcorn.size() != numGrid * 8) OCP_ABORT("Wrong ZCORN size!");
    }

    if (poro.size() != numGrid) OCP_ABORT("Wrong PORO size!");
    if (permX.size() != numGrid) OCP_ABORT("Wrong PERMX size!");
    if (permY.size() != numGrid) OCP_ABORT("Wrong PERMY size!");
    if (permZ.size() != numGrid) OCP_ABORT("Wrong PERMZ size!");

    if (ntg.size() != numGrid) {
        ntg.resize(numGrid, 1);
        cout << "Set net-to-gross ratio to 1.0!" << endl;
    }
}

/// Check rock keyword.
void ParamReservoir::CheckRock()
{
    if (rockSet.size() != NTROOC) {
        OCP_ABORT("Wrong ROCK or ROCKT!");
    }
}

/// Check EQUIL keyword.
void ParamReservoir::CheckEQUIL() const
{
    if (EQUIL.empty()) OCP_ABORT("EQUIL is missing!");
}

/// TODO: Add Doxygen
void ParamReservoir::CheckDenGra() const
{
    if (density.activity && gravity.activity) {
        OCP_ABORT("Both DENSITY and GRAVITY have been given, just one can be used!");
    }
}

/// TODO: Add Doxygen
void ParamReservoir::CheckPhase() const
{
    if (blackOil && disGas && (!gas && !oil)) {
        OCP_ABORT("DISGAS can only be used only if OIL and GAS are both present!");
    }
}

/// Check tables: Different tables will be used under different conditions.
void ParamReservoir::CheckPhaseTab() const
{
    if (!blackOil && !comps) OCP_ABORT("Unknown model: Use BLACKOIL or COMPS!");

    if (water && oil && SWOF_T.data.empty()) OCP_ABORT("SWOF is missing!");
    if (gas && oil && SGOF_T.data.empty()) OCP_ABORT("SGOF is missing!");
    if (water && PVTW_T.data.empty()) OCP_ABORT("PVTW is missing!");

    if (blackOil) {
        if (oil && disGas && PVCO_T.data.empty()) OCP_ABORT("PVCO is missing!");
        if (oil && (!disGas) && PVDO_T.data.empty()) OCP_ABORT("PVDO is missing!");
        if (gas && PVDG_T.data.empty()) OCP_ABORT("PVDG is missing!");
    }
}

/// TODO: Add Doxygen
void ParamReservoir::CheckRegion() const
{
    if (SATNUM.activity && SATNUM.data.size() != numGrid) {
        OCP_ABORT("Missing data in SATNUM!");
    }
    if (PVTNUM.activity && PVTNUM.data.size() != numGrid) {
        OCP_ABORT("Missing data in PVTNUM!");
    }
    if (ACTNUM.activity && ACTNUM.data.size() != numGrid) {
        OCP_ABORT("Missing data in ACTNUM!");
    }
    if (ROCKNUM.activity && ROCKNUM.data.size() != numGrid) {
        OCP_ABORT("Missing data in ROCKNUM!");
    }
}

/// TODO: Add Doxygen
void ParamReservoir::CheckEqlRegion() const
{
    if (PBVD_T.data.size() > 1) {
        OCP_ABORT("More than one equilibrium region is not supported!");
    }
}

/// TODO: Add Doxygen
void TableSet::DisplayTable() const
{
    cout << "\n---------------------" << endl
         << name << "\n---------------------" << endl;

    for (USI n = 0; n < data.size(); n++) {
        if (refName.size() > n) {
            cout << refName[n] << "   ";
            cout << refData[n] << endl;
        }

        const USI len = data[n][0].size();
        for (USI i = 0; i < len; i++) {
            for (USI j = 0; j < colNum; j++) {
                cout << setw(10) << data[n][j][i];
            }
            cout << "\n";
        }
    }
}

void ComponentParam::Init()
{
    // Init LBC coefficient
    LBCcoef.resize(5);
    LBCcoef[0] = 0.1023;
    LBCcoef[1] = 0.023364;
    LBCcoef[2] = 0.058533;
    LBCcoef[3] = -0.040758;
    LBCcoef[4] = 0.0093324;
}

Type_A_r<vector<OCP_DBL>>* ComponentParam::FindPtr01(const string& varName)
{
    Type_A_r<vector<OCP_DBL>>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("TCRIT", 5):
            myPtr = &Tc;
            break;

        case Map_Str2Int("PCRIT", 5):
            myPtr = &Pc;
            break;

        case Map_Str2Int("VCRIT", 5):
            myPtr = &Vc;
            break;

        case Map_Str2Int("ZCRIT", 5):
            myPtr = &Zc;
            break;

        case Map_Str2Int("MW", 2):
            myPtr = &MW;
            break;

        case Map_Str2Int("ACF", 3):
            myPtr = &Acf;
            break;

        case Map_Str2Int("OMEGAA", 6):
            myPtr = &OmegaA;
            break;

        case Map_Str2Int("OMEGAB", 6):
            myPtr = &OmegaB;
            break;

        case Map_Str2Int("SSHIFT", 6):
            myPtr = &Vshift;
            break;

        case Map_Str2Int("PARACHOR", 8):
            myPtr = &Parachor;
            break;

        case Map_Str2Int("VCRITVIS", 8):
            myPtr = &Vcvis;
            break;

        case Map_Str2Int("ZCRITVIS", 8):
            myPtr = &Zcvis;
            break;

        case Map_Str2Int("MOLDEN", 6):
            myPtr = &molden;
            break;

        case Map_Str2Int("CP", 2):
            myPtr = &cp;
            break;

        case Map_Str2Int("CT1", 3):
            myPtr = &ct1;
            break;

        case Map_Str2Int("CT2", 3):
            myPtr = &ct2;
            break;

        case Map_Str2Int("CPT", 3):
            myPtr = &cpt;
            break;

        case Map_Str2Int("CPL1", 4):
            myPtr = &cpl1;
            break;

        case Map_Str2Int("CPL2", 4):
            myPtr = &cpl2;
            break;

        case Map_Str2Int("CPL3", 4):
            myPtr = &cpl3;
            break;

        case Map_Str2Int("CPL4", 4):
            myPtr = &cpl4;
            break;

        case Map_Str2Int("CPG1", 4):
            myPtr = &cpg1;
            break;

        case Map_Str2Int("CPG2", 4):
            myPtr = &cpg2;
            break;

        case Map_Str2Int("CPG3", 4):
            myPtr = &cpg3;
            break;

        case Map_Str2Int("CPG4", 4):
            myPtr = &cpg4;
            break;

        case Map_Str2Int("HVAPR", 5):
            myPtr = &hvapr;
            break;

        case Map_Str2Int("HVR", 3):
            myPtr = &hvr;
            break;

        case Map_Str2Int("EV", 2):
            myPtr = &ev;
            break;

        case Map_Str2Int("AVSIC", 5):
            myPtr = &avisc;
            break;

        case Map_Str2Int("BVSIC", 5):
            myPtr = &bvisc;
            break;

        case Map_Str2Int("AVG", 3):
            myPtr = &avg;
            break;

        case Map_Str2Int("BVG", 3):
            myPtr = &bvg;
            break;
    }

    return myPtr;
}

void ComponentParam::InputRefPR(ifstream& ifs, const string& keyword)
{
    OCP_ASSERT(NTPVT > 0, "NTPVT has not been set!");

    vector<OCP_DBL>* objPtr = nullptr;
    objPtr                  = FindPtr02(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }

    vector<string> vbuf;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);

        if (vbuf[0] == "/") break;

        for (auto& v : vbuf) {
            if (v != "/") objPtr->push_back(stod(v));
            if (objPtr->size() >= NTPVT) break;
        }
        if (objPtr->size() >= NTPVT) break;
    }
    cout << keyword << endl;
    for (USI i = 0; i < NTPVT; i++) {
        cout << objPtr->at(i) << "   ";
    }
    cout << "\n/" << endl << endl;
}

vector<OCP_DBL>* ComponentParam::FindPtr02(const string& varName)
{
    vector<OCP_DBL>* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("PRSR", 4):
            myPtr = &Pref;
            break;

        case Map_Str2Int("TEMR", 4):
            myPtr = &Tref;
            break;
    }

    return myPtr;
}

void ComponentParam::InputCOMPONENTS(ifstream& ifs, const string& keyword)
{
    OCP_ASSERT((numCom > 0) && (NTPVT > 0), "NPNC has not been set!");

    Type_A_r<vector<OCP_DBL>>* objPtr = nullptr;
    objPtr                            = FindPtr01(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unknown keyword!");
    }
    objPtr->activity = OCP_TRUE;

    vector<string>  vbuf;
    vector<OCP_DBL> tmp;
    USI             nReg = 0;

    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            nReg++;
            objPtr->data.push_back(tmp);
            if (nReg >= NTPVT) break;
            tmp.clear();
            continue;
        }
        for (auto& v : vbuf) {
            if (v != "/") {
                tmp.push_back(stod(v));
            }
        }
        if (vbuf.back() == "/") {
            nReg++;
            objPtr->data.push_back(tmp);
            tmp.clear();
            if (nReg >= NTPVT) break;
        }
    }

    cout << keyword << endl;
    for (USI i = 0; i < NTPVT; i++) {
        for (auto& v : objPtr->data[i]) {
            cout << v << endl;
        }
        cout << "/" << endl;
    }
    cout << endl;
}

void ComponentParam::InputCNAMES(ifstream& ifs)
{
    OCP_ASSERT(numCom > 0, "numCom has not been set!");

    vector<string> vbuf;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            break;
        }
        for (auto& v : vbuf) {
            if (v != "/") Cname.push_back(v);
        }
        if (vbuf.back() == "/") break;
    }

    OCP_FUNCNAME;
    cout << "CNAMES" << endl;
    for (USI i = 0; i < numCom; i++) {
        cout << Cname[i] << "   ";
    }
    cout << endl << endl;
}

void ComponentParam::InputLBCCOEF(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    DealDefault(vbuf);
    for (USI i = 0; i < 5; i++) {
        if (vbuf[i] != "DEFAULT") LBCcoef[i] = stod(vbuf[i]);
    }

    OCP_FUNCNAME;
    cout << "LBCCOEF" << endl;
    for (USI i = 0; i < 5; i++) {
        cout << LBCcoef[i] << "   ";
    }
    cout << endl << endl;
}

/// Input Binary Interaction Coefficients Matrix
void ComponentParam::InputBIC(ifstream& ifs)
{
    OCP_ASSERT((numCom > 0) && (NTPVT > 0), "numCom or NTPVT has not been set!");

    BIC.resize(NTPVT);

    vector<string> vbuf;
    USI            nReg = 0;
    while (OCP_TRUE) {
        ReadLine(ifs, vbuf);
        if (vbuf[0] == "/") {
            nReg++;
            if (nReg >= NTPVT) break;
            continue;
        }
        for (auto& v : vbuf) {
            if (v != "/") {
                BIC[nReg].push_back(stod(v));
                cout << setw(10) << BIC[nReg].back();
            }
        }
        cout << endl;
        if (vbuf.back() == "/") {
            nReg++;
            if (nReg >= NTPVT) break;
        }
    }
}

void ComponentParam::InputVISCTAB(ifstream& ifs)
{
    vector<string>          vbuf;
    vector<vector<OCP_DBL>> tmp;
    USI                     ncol = numCom + 1; // temp + comps
    tmp.resize(ncol);
    OCP_BOOL flag = OCP_TRUE;

    while (OCP_TRUE) {
        if (flag) {
            ReadLine(ifs, vbuf);
            flag = OCP_FALSE;
        }

        if (vbuf[0] == "/") break;

        if (vbuf[0] == "ATPRES") {
            viscTab.refName.push_back("ATPRES");
            viscTab.refData.push_back(stod(vbuf[1]));
            ReadLine(ifs, vbuf);
        }
        // Read Table
        while (OCP_TRUE) {
            for (USI i = 0; i < ncol; i++) {
                tmp[i].push_back(stod(vbuf[i]));
            }
            ReadLine(ifs, vbuf);
            if (vbuf[0] == "ATPRES" || vbuf[0] == "/") {
                viscTab.data.push_back(tmp);
                tmp.clear();
                tmp.resize(ncol);
                break;
            }
        }
        if (vbuf[0] == "/") break;
    }
    viscTab.name   = "VISCTAB";
    viscTab.colNum = ncol;
    // output
    viscTab.DisplayTable();

    cout << "/" << endl;
}

/// TODO: Add Doxygen
void ComponentParam::InputSSMSTA(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    int len = vbuf.size();
    for (int i = 0; i < len; i++) {
        SSMparamSTA.push_back(vbuf[i]);
    }
    OCP_FUNCNAME;
    for (int i = 0; i < len; i++) {
        cout << SSMparamSTA[i] << "   ";
    }
    cout << endl << endl;
}

/// TODO: Add Doxygen
void ComponentParam::InputNRSTA(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (int i = 0; i < 2; i++) {
        NRparamSTA.push_back(vbuf[i]);
    }
    OCP_FUNCNAME;
    for (int i = 0; i < 2; i++) {
        cout << NRparamSTA[i] << "   ";
    }
    cout << endl << endl;
}

/// TODO: Add Doxygen
void ComponentParam::InputSSMSP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        SSMparamSP.push_back(vbuf[i]);
    }
    OCP_FUNCNAME;
    for (USI i = 0; i < 2; i++) {
        cout << SSMparamSP[i] << "   ";
    }
    cout << endl << endl;
}

/// TODO: Add Doxygen
void ComponentParam::InputNRSP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        NRparamSP.push_back(vbuf[i]);
    }
    OCP_FUNCNAME;
    for (USI i = 0; i < 2; i++) {
        cout << NRparamSP[i] << "   ";
    }
    cout << endl << endl;
}

/// TODO: Add Doxygen
void ComponentParam::InputRR(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    for (USI i = 0; i < 2; i++) {
        RRparam.push_back(vbuf[i]);
    }
    OCP_FUNCNAME;
    for (USI i = 0; i < 2; i++) {
        cout << RRparam[i] << "   ";
    }
    cout << endl << endl;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update output and Doxygen            */
/*----------------------------------------------------------------------------*/