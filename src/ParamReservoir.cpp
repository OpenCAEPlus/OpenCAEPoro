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

/// TODO: Add Doxygen
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

        case Map_Str2Int("SATNUM", 6):
            SATNUM.activity = true;
            SATNUM.data.reserve(numGrid);
            myPtr           = &SATNUM.data;
            break;

        case Map_Str2Int("PVTNUM", 6):
            PVTNUM.activity = true;
            PVTNUM.data.reserve(numGrid);
            myPtr           = &PVTNUM.data;
            break;
    }

    return myPtr;
}

/// TODO: Add Doxygen
TableSet* ParamReservoir::FindPtr_T(const string& varName)
{
    TableSet* myPtr = nullptr;

    switch (Map_Str2Int(&varName[0], varName.size())) {
        case Map_Str2Int("SWOF", 4):
            myPtr = &SWOF_T;
            break;

        case Map_Str2Int("SGOF", 4):
            myPtr = &SGOF_T;
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
    }

    return myPtr;
}

/// TODO: Add Doxygen
void ParamReservoir::Init()
{
    InitTable();

    rsTemp = 60;

    gravity.data.resize(3);
    gravity.data[0] = 45.5;   // oil
    gravity.data[1] = 1.0;    // pure water
    gravity.data[2] = 0.7773; // air

    density.data.resize(3);
    density.data[0] = 37.457;    // oil
    density.data[1] = 62.366416; // pure water
    density.data[2] = 0.062428;  // air

    rock.Pref = 14.7;
    rock.Cr   = 3.406E-6;
}

/// TODO: Add Doxygen
void ParamReservoir::InitTable()
{
    SWOF_T.name   = "SWOF";
    SWOF_T.colNum = 4;
    SGOF_T.name   = "SGOF";
    SGOF_T.colNum = 4;
    PBVD_T.name   = "PBVD";
    PBVD_T.colNum = 2;
    PVCO_T.name   = "PVCO";
    PVCO_T.colNum = 6;
    PVDO_T.name   = "PVDO";
    PVDO_T.colNum = 3;
    PVDG_T.name   = "PVDG";
    PVDG_T.colNum = 3;
    PVTW_T.name   = "PVTW";
    PVTW_T.colNum = 5;
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
void ParamReservoir::CopyVal(vector<T>& obj, const vector<T>& src,
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
void ParamReservoir::MultiplyVal(vector<OCP_DBL>& obj, const OCP_DBL& val,
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
    comps = true;
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    numCom = stoi(vbuf[0]);
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
}

/// TODO: Add Doxygen
void ParamReservoir::DisplayDIMENS()
{
    std::cout << "DIMENS" << endl;
    std::cout << dimens.nx << "  " << dimens.ny << "  " << dimens.nz << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputRTEMP(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    rsTemp = stod(vbuf[0]);
    cout << "RTEMP" << endl;
    cout << rsTemp << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputEQUALS(ifstream& ifs)
{
    vector<USI>    index(6, 0);
    vector<string> vbuf;

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

    std::cout << SATNUM.activity << endl;
    std::cout << PVTNUM.activity << endl;
    cout << "EQUALS" << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputGRID(ifstream& ifs, string& keyword)
{
    vector<OCP_DBL>* objPtr = nullptr;

    objPtr = FindPtr(keyword);
    if (objPtr == nullptr) {
        OCP_ABORT("Unmatched keyword!");
    }

    vector<string> vbuf;
    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        for (auto str : vbuf) {
            objPtr->push_back(stod(str));
        }
    }
    std::cout << &permX << endl;
    std::cout << &permY << endl;
    std::cout << &permZ << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputCOPY(ifstream& ifs)
{
    vector<string> vbuf;
    vector<USI>    index(6, 0);

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

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
    std::cout << permX[0] << endl;
    std::cout << permY[0] << endl;
    std::cout << permZ[0] << endl;
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
    std::cout << permX[0] << endl;
    std::cout << permY[0] << endl;
    std::cout << permZ[0] << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputTABLE(ifstream& ifs, const string& tabName)
{
    TableSet* obj;
    obj = FindPtr_T(tabName);
    if (obj == nullptr) {
        OCP_ABORT("Wrong table name :" + tabName);
    }

    USI                     col = obj->colNum;
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

/// TODO: Add Doxygen
void ParamReservoir::InputROCK(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    rock.Pref = stod(vbuf[0]);
    rock.Cr   = stod(vbuf[1]);

    std::cout << "ROCK" << endl;
    std::cout << rock.Pref << "  " << rock.Cr << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputGRAVITY(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            gravity.activity = true;
            gravity.data[i]  = stod(vbuf[i]);
        }
    }

    std::cout << "GRAVITY" << endl;
    std::cout << gravity.data[0] << "  " << gravity.data[1] << "  " << gravity.data[2]
              << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputDENSITY(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    if (vbuf[0] == "/") return;

    DealDefault(vbuf);
    for (USI i = 0; i < 3; i++) {
        if (vbuf[i] != "DEFAULT") {
            density.activity = true;
            density.data[i]  = stod(vbuf[i]);
        }
    }

    std::cout << "DENSITY" << endl;
    std::cout << density.data[0] << "  " << density.data[1] << "  " << density.data[2]
              << endl;
}

/// TODO: Add Doxygen
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

    std::cout << "EQUIL" << endl;
    for (USI i = 0; i < 6; i++) std::cout << EQUIL[i] << "  ";
    std::cout << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputTABDIMS(ifstream& ifs)
{
    vector<string> vbuf;
    ReadLine(ifs, vbuf);
    NTSFUN = stoi(vbuf[0]);
    NTPVT  = stoi(vbuf[1]);
    cout << "TABDIMS" << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::InputRegion(ifstream& ifs, const string& keyword)
{
    Type_A_r<OCP_DBL>* ptr = &PVTNUM;
    USI                lim = NTPVT;

    if (keyword == "SATNUM") {
        ptr = &SATNUM;
        lim = NTSFUN;
    }

    ptr->activity = true;
    ptr->data.reserve(numGrid);
    vector<string>  vbuf;
    vector<OCP_USI> obj;
    vector<USI>     region;

    while (ReadLine(ifs, vbuf)) {
        if (vbuf[0] == "/") break;

        DealData(vbuf, obj, region);

#ifdef DEBUG
        // check region
        for (auto r : region) {
            if (r > lim) {
                OCP_ABORT("Region is out of range!");
            }
        }
#endif // DEBUG

        USI len = obj.size();
        for (USI i = 0; i < len; i++) {
            for (OCP_USI j = 0; j < obj[i]; j++) {
                ptr->data.push_back(region[i]);
            }
        }
    }

    cout << "Number of tables = " << lim << endl;
    cout << &SATNUM << endl << &PVTNUM << endl;
}

/// TODO: Add Doxygen
void ParamReservoir::CheckParam()
{
    CheckGrid();
    CheckEQUIL();
    CheckDenGra();
    CheckPhase();
    CheckPhaseTab();
    CheckRegion();
}

/// Check data dimension for potential problems.
void ParamReservoir::CheckGrid()
{
    if (coord.size() == 0) {
        if (tops.size() != dimens.nx * dimens.ny) OCP_ABORT("Wrong TOPS size!");
        if (dx.size() != numGrid) OCP_ABORT("Wrong DX size!");
        if (dy.size() != numGrid) OCP_ABORT("Wrong DY size!");
        if (dz.size() != numGrid) OCP_ABORT("Wrong DZ size!");
    }
    else {
        if (coord.size() != (dimens.nx+1)*(dimens.ny+1)*6) OCP_ABORT("Wrong COORD size!");
        if (zcorn.size() != numGrid*8) OCP_ABORT("Wrong ZCORN size!");
    }
    
    if (poro.size() != numGrid) OCP_ABORT("Wrong PORO size!");
    if (permX.size() != numGrid) OCP_ABORT("Wrong PERMX size!");
    if (permY.size() != numGrid) OCP_ABORT("Wrong PERMY size!");
    if (permZ.size() != numGrid) OCP_ABORT("Wrong PERMZ size!");
    if (ntg.size() != numGrid) {
        OCP_WARNING("Wrong Ntg size; set it to 1."); // TODO: Should not appear every time!
        ntg.resize(numGrid, 1);
    }
}

/// TODO: Add Doxygen
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
    if (disGas && (!gas && !oil)) {
        OCP_ABORT("DISGAS can only be used only if OIL and GAS are both present!");
    }
}

/// TODO: Add Doxygen
void ParamReservoir::CheckPhaseTab() const
{
    if (!blackOil && !comps) OCP_ABORT("Unknown model: choose BLACKOIL or COMPS!");

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
}

/// TODO: Add Doxygen
void ParamReservoir::CheckEqlRegion() const
{
    if (PBVD_T.data.size() > 1) {
        OCP_ABORT("Only one equilibration region is supported!");
    }
}

/// TODO: Add Doxygen
void TableSet::DisplayTable() const
{
    cout << "---------------------\n";
    cout << "TABLE: " << name << "\n";
    cout << "---------------------\n";
    for (auto v : data) {
        const USI len = v[0].size();
        for (USI i = 0; i < len; i++) {
            for (USI j = 0; j < colNum; j++) {
                cout << v[j][i] << "\t";
            }
            cout << "\n";
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
/*----------------------------------------------------------------------------*/