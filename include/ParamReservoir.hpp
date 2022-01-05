/*! \file    ParamReservoir.hpp
 *  \brief   ParamReservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMRESERVOIR_HEADER__
#define __PARAMRESERVOIR_HEADER__

// Standard header files
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"

using namespace std;

/// TableSet is used to store a series of tables which have the same type. For example,
/// a series of PVTW Table, each PVTW table is two-dimension.
class TableSet
{
public:
    void DisplayTable() const; ///< Print table

public:
    string                          name;   ///< Name of table.
    USI                             colNum; ///< Number of columns of table.
    vector<vector<vector<OCP_DBL>>> data;   ///< All table with the same name.
};

/// Dimens contains the dimensions of grids.
class Dimens
{
public:
    USI nx; ///< Num of bulks along x-direction.
    USI ny; ///< Num of bulks along y-direction.
    USI nz; ///< Num of bulks along z-direction.
};

/// Rock class contains information about the keyword ROCK.
class Rock
{
public:
    OCP_DBL Pref; ///< Reference pressure at initial porosity.
    OCP_DBL Cr;   ///< Compressibility factor of rock in reservoir.
};

/// A internal structure used to store some params for reservoir, it can tell
/// if these params are given by users.
template <typename T> class Type_A_r
{
public:
    bool      activity{false}; ///< If false, this param is not given.
    vector<T> data;            ///< Data of param.
};

/// EoSParam contains the params for Compositional Model and functions to read them
class EoSparam
{
public:
    // Basic params
    void InputNCNP(ifstream& ifs);
    void InputZI(ifstream& ifs);
    void InputCOM(ifstream& ifs);
    void InputBIP(ifstream& ifs);
    // Method params
    void InputSSMSTA(ifstream& ifs); 
    void InputNRSTA(ifstream& ifs);
    void InputSSMSP(ifstream& ifs); 
    void InputNRSP(ifstream& ifs); 
    void InputRR(ifstream& ifs);  

public:
    USI numComp{ 0 };   ///< num of components, water is excluded.
    USI numPhase{ 0 };  ///< num of phase, water is excluded.
    vector<OCP_DBL> zi;          ///< molar fraction of components in initial reservoir, water is excluded.
    vector<vector<string>>  COM; ///< Coponents information
    vector<vector<OCP_DBL>> BIP; ///< Binary interaction

    vector<string>  SSMparamSTA; ///< Params for Solving Phase Spliting with SSM
    vector<string>	NRparamSTA; ///< Params for Solving Phase Spliting with NR
    vector<string>  SSMparamSP; ///< Params for Solving Phase Spliting with SSM
    vector<string>	NRparamSP; ///< Params for Solving Phase Spliting with NR
    vector<string>	RRparam;  ///< Params for Solving Rachford-Rice equations
};


/// ParamReservoir is an internal structure used to stores the information of
/// reservoir(except wells) from input files. It is an intermediate interface and
/// independent of the main simulator. After all file inputting finishs, the params in
/// it will pass to corresponding modules.
class ParamReservoir
{

public:
    // Grid
    // Cartesian
    Dimens  dimens;  ///< Dimension of grid: the number of grids along x,y,z direction.
    OCP_USI numGrid; ///< Num of grids.
    vector<OCP_DBL> tops; ///< Depth of the top surface of the uppermost grids.
    vector<OCP_DBL> dx;   ///< Size along the x - direction for each grid.
    vector<OCP_DBL> dy;   ///< Size along the y - direction for each grid.
    vector<OCP_DBL> dz;   ///< Size along the z - direction for each grid.

    // Corner point geometry
    vector<OCP_DBL> coord;
    vector<OCP_DBL> zcorn;

    OCP_DBL rsTemp; ///< Temperature for reservoir.

    // Rock
    vector<OCP_DBL> ntg;   ///< Net to gross for each grid.
    vector<OCP_DBL> poro;  ///< Porosity for each grid.
    vector<OCP_DBL> permX; ///< Permeability along the x - direction for each grid.
    vector<OCP_DBL> permY; ///< Permeability along the y-direction for each grid.
    vector<OCP_DBL> permZ; ///< Permeability along the z-direction for each grid.
    Rock rock; ///< Contains the compressibility factor and reference pressure at
               ///< initial porosity.

    // If P and Ni are given, then calculation of initial equilibration is unneeded.
    vector<OCP_DBL> P;  ///< Initial pressure of components in each grid.
    vector<OCP_DBL> Ni; ///< Initial moles of components in each grid.

    // phase property
    Type_A_r<OCP_DBL> density; ///< Density of oil, water, gas in standard conditions.
    Type_A_r<OCP_DBL> gravity; ///< Gravity of oil, water, gas in standard conditions.

    // Model and Phase
    bool blackOil{false}; ///< If ture, blackoil model will be used.
    bool comps{false};    ///< If true, compositional model will be used.
    bool oil{false};      ///< If true, oil phase could exist.
    bool gas{false};      ///< If true, gas phase could exist.
    bool water{false};    ///< If true, water phase could exist.
    bool disGas{false};   ///< If true, dissolve gas could exist in oil phase.

    // Compositional Model
    /// Contains the initial proportion of components. It's used in compositional model.
    EoSparam EoSp;

    // SAT Region & PVT Region
    USI               NTSFUN{1}; ///< Num of SAT regions.
    USI               NTPVT{1};  ///< Num of PVT regions.
    Type_A_r<OCP_DBL> SATNUM;    ///< Records the index of SAT region for each grid.
    Type_A_r<OCP_DBL> PVTNUM;    ///< Records the index of PVT region for each grid.

    // Saturation table & buble point pressure
    TableSet        SWOF_T; ///< Table set of SWOF.
    TableSet        SGOF_T; ///< Table set of SGOF.
    TableSet        PBVD_T; ///< Table set of PBVD.
    vector<OCP_DBL> EQUIL;  ///< See ParamEQUIL.

    // PVT property
    USI      numPhase; ///< Number of phases
    USI      numCom;   ///< Number of components, used in Compositional Model
    TableSet PVCO_T;   ///< Table set of PVCO.
    TableSet PVDO_T;   ///< Table set of PVDO.
    TableSet PVDG_T;   ///< Table set of PVDG.
    TableSet PVTW_T;   ///< Table set of PVTW.

    // internal method
    /// Find corresponding variable according to the name of variable.
    /// It is used for the basic properties of reservoir such as DX.
    vector<OCP_DBL>* FindPtr(const string& varName);
    /// Find corresponding variable according to the name of variable.
    /// It is used for the scope of the table.
    TableSet* FindPtr_T(const string& varName);

    /// Initialize the default value in reservoir, such as temperature, density, table.
    void Init();
    /// Initialize the tables' name and num of colum.
    void InitTable();

    /// It's used in InputEQUALS, assigning values in batches.
    template <typename T>
    void setVal(vector<T>& obj, const T& val, const vector<USI>& index);

    /// It's used in InputCOPY, copying the value of one variable to another.
    template <typename T>
    void CopyVal(vector<T>& obj, const vector<T>& src, const vector<USI>& index);

    /// It's used in InputMULTIPLY, multipling the value of a certain range of a
    /// variable by a coefficient.
    void MultiplyVal(vector<OCP_DBL>& obj, const OCP_DBL& val,
                     const vector<USI>& index);

    /// Input the keyword: COMPS. COMPS is used in compositional model, which gives the
    /// num of components.
    void InputCOMPS(ifstream& ifs);

    /// Input the keyword: DIMENS. DIMENS contain the dimension of grids of reservoir.
    /// It gives the num of grids along the x,y,z direction.
    void InputDIMENS(ifstream& ifs);
    /// Display the dimens, it's used to chech input.
    void DisplayDIMENS();

    /// Input the keyword: RTEMP. RTEMP gives the temperature of reservoir.
    void InputRTEMP(ifstream& ifs);

    /// Input the keyword: EQUALS. EQUALS contains many keywords about grids which has
    /// special input format. These keywords contains DX, TOPS, PORO and so on. You can
    /// assign values to them in batches
    void InputEQUALS(ifstream& ifs);

    /// Input the keyword about grids, actually, it's a supplement for EQUALS.
    /// It supplies another way to input the params in EQUALS.
    void InputGRID(ifstream& ifs, string& keyword);

    /// Input the keyword: COPY. COPY could copy the value of one variable to another.
    void InputCOPY(ifstream& ifs);

    /// Input the keyword: MULTIPLY. MULTIIPLY could multiply the value of a certain
    /// range of a variable by a coefficient.
    void InputMULTIPLY(ifstream& ifs);

    /// Input PVTtable and SATtable such as SWOF, PVCO.
    void InputTABLE(ifstream& ifs, const string& tabName);

    /// Input the keyword: ROCK. ROCK contains the compressibility factor and reference
    /// pressure at initial porosity.
    void InputROCK(ifstream& ifs);

    /// Input the reference gravity of oil, water, and air in standard condition.
    void InputGRAVITY(ifstream& ifs);

    /// Input the reference density of oil, water, and air in standard condition.
    void InputDENSITY(ifstream& ifs);

    /// Input the keyword: EQUIL. EQUIL contains initial information of reservoir.
    /// See ParamEQUIL.
    void InputEQUIL(ifstream& ifs);

    // SATNUM & PVTNUM  -- Region
    /// Input the keyword: TABDIMS. It contains the num of saturation region and PVT
    /// region.
    void InputTABDIMS(ifstream& ifs);
    /// input the keyword: SATNUM, PVTNUM. For example, SATNUM decides which grid
    /// belongs to which saturation region, so corresponding saturation table will
    /// be used.
    void InputRegion(ifstream& ifs, const string& keyword);

    // Input EoSparam
    // Basic params
    void InputNCNP(ifstream& ifs) { EoSp.InputNCNP(ifs); };
    void InputZI(ifstream& ifs) { EoSp.InputZI(ifs); };
    void InputCOM(ifstream& ifs) { EoSp.InputCOM(ifs); };
    void InputBIP(ifstream& ifs) { EoSp.InputBIP(ifs); };
    // Method params
    void InputSSMSTA(ifstream& ifs) { EoSp.InputSSMSTA(ifs); };
    void InputNRSTA(ifstream& ifs) { EoSp.InputNRSTA(ifs); };
    void InputSSMSP(ifstream& ifs) { EoSp.InputSSMSP(ifs); };
    void InputNRSP(ifstream& ifs) { EoSp.InputNRSP(ifs); };
    void InputRR(ifstream& ifs) { EoSp.InputRR(ifs); };

    // check
    /// Check the reservoir param from input file.
    void CheckParam();
    /// Check the size of properties of grids.
    void CheckGrid();
    /// Check if keyword EQUIL is given.
    void CheckEQUIL() const;
    /// Check if density and gravity are both input, only one of them is needed.
    void CheckDenGra() const;
    /// Check the existence of disgas, it could only exist when both oil and gas exist.
    void CheckPhase() const;
    /// Check the existence of PVTtable and SATtable. Different tables will be used in
    /// different conditions.
    void CheckPhaseTab() const;
    /// Check if each grid is assigned to an area or all defaulted.
    void CheckRegion() const;
    /// (To do) In current program, only initialization of equilibration of only one
    /// region is realized.
    void CheckEqlRegion() const;
};

#endif /* end if __PARAMRESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/