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
/// a series of Function Table, each Function table is two-dimension.
class TableSet
{
public:
    void DisplayTable() const; ///< Print table

public:
    string                          name;    ///< Name of table.
    USI                             colNum;  ///< Number of columns of table.
    vector<string>                  refName; ///< name of refData
    vector<OCP_DBL>                 refData; ///< refData for each sub data.
    vector<vector<vector<OCP_DBL>>> data;    ///< All table with the same name.
};

/// Dimens contains the dimensions of grids.
class Dimens
{
public:
    USI nx; ///< Num of bulks along x-direction.
    USI ny; ///< Num of bulks along y-direction.
    USI nz; ///< Num of bulks along z-direction.
};

class HLoss
{
public:
    OCP_BOOL ifHLoss{OCP_FALSE}; ///< If use Heat loss
    OCP_DBL  obC;                ///< Volumetric heat capacity of overburden rock
    OCP_DBL  obK;                ///< Thermal conductivity of overburden rock
    OCP_DBL  ubC;                ///< Volumetric heat capacity of underburden rock
    OCP_DBL  ubK;                ///< Thermal conductivity of underburden rock
};

/// RockParam class contains information about the keyword ROCK.
class RockParam
{
public:
    string   type{"LINEAR"}; ///< LINEAR or EXPONENT for porosity model
    OCP_DBL  Pref{14.7};     ///< Reference pressure at initial porosity.
    OCP_DBL  Tref{60};       ///< Reference temperature at initial porosity.
    OCP_DBL  cp1{3.406E-6};  ///< Compressibility factor of rock in reservoir.
    OCP_DBL  cp2{0};         ///< 2 order Compressibility factor of rock in reservoir.
    OCP_DBL  ct{0};          ///< Expansion factor of rock in reservoir, ifThermal only
    OCP_DBL  cpt{0};         ///< cross items, ifThermal only
    OCP_BOOL ConstRock{OCP_TRUE}; ///< if true, rock volume remains const, else, bulk
                                  ///< volume remains const
    OCP_DBL HCP1{35}; ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
    OCP_DBL HCP2{0};  ///< coefficients of the rock enthalpy formula, Btu/ft^3 - F
};

/// A internal structure used to store some params for reservoir, it can tell
/// if these params are given by users.
template <typename T>
class Type_A_r
{
public:
    OCP_BOOL  activity{OCP_FALSE}; ///< If OCP_FALSE, this param is not given.
    vector<T> data;                ///< Data of param.
};

/// ComponentParam contains information of components
class ComponentParam
{
public:
    // Basic params
    /// Init Params
    void Init();
    /// Input the information of hydrocarbon components
    void InputCOMPONENTS(ifstream& ifs, const string& keyword);
    /// Find corresponding variable according to the name of variable.
    /// It is used for the basic properties of hydrocarbon components such as TCRIT
    Type_A_r<vector<OCP_DBL>>* FindPtr01(const string& varName);
    /// input reference pressure, temperature
    void InputRefPR(ifstream& ifs, const string& keyword);
    /// Find corresponding variable according to the name of variable
    vector<OCP_DBL>* FindPtr02(const string& varName);
    /// Input the names of hydrocarbon components
    void InputCNAMES(ifstream& ifs);
    /// Input LBC coefficients for viscosity calculation
    void InputLBCCOEF(ifstream& ifs);
    /// Input the Binary interaction of components
    void InputBIC(ifstream& ifs);
    /// Input VISCTAB
    void InputVISCTAB(ifstream& ifs);
    // Method params
    void InputSSMSTA(ifstream& ifs);
    void InputNRSTA(ifstream& ifs);
    void InputSSMSP(ifstream& ifs);
    void InputNRSP(ifstream& ifs);
    void InputRR(ifstream& ifs);

public:
    USI            NTPVT{1};       ///< num of EoS region, constant now.
    USI            numCom{0};      ///< num of components, water is excluded.
    USI            numPhase{2};    ///< num of phase, water is excluded, constant now.
    vector<string> Cname;          ///< Name of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Tc;  ///< Critical temperature of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Pc;  ///< Critical pressure of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Vc;  ///< Critical volume of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Zc;  ///< Critical Z-factor of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> MW;  ///< Molecular Weight of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Acf; ///< Acentric factor of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> OmegaA;   ///< OMEGA_A of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> OmegaB;   ///< OMEGA_B of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Vshift;   ///< Volume shift of hydrocarbon components
    Type_A_r<vector<OCP_DBL>> Parachor; ///< PARACHOR of hydrocarbon components
    // for viscosity calculation
    Type_A_r<vector<OCP_DBL>>
        Vcvis; ///< Critical volume used for viscosity calculations only.
    Type_A_r<vector<OCP_DBL>>
                    Zcvis; ///< Critical Z-factor used for viscosity calculations only.
    vector<OCP_DBL> LBCcoef;     ///< LBC coefficients for viscosity calculation
    vector<vector<OCP_DBL>> BIC; ///< Binary interaction
    // ifThermal only
    Type_A_r<vector<OCP_DBL>> molden; ///< component molar density at reference
                                      ///< temperature and reference pressure, lb/ft3
    Type_A_r<vector<OCP_DBL>> cp;     ///< component compressibility, 1/psi
    Type_A_r<vector<OCP_DBL>> ct1;  ///< the first ifThermal expansion coefficient, 1/F
    Type_A_r<vector<OCP_DBL>> ct2;  ///< the second ifThermal expansion coefficient, 1/F
    Type_A_r<vector<OCP_DBL>> cpt;  ///< the coefficient of density dependence on
                                    ///< temperature and pressure, 1/psi-F
    Type_A_r<vector<OCP_DBL>> cpl1; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F
    Type_A_r<vector<OCP_DBL>> cpl2; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^2
    Type_A_r<vector<OCP_DBL>> cpl3; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^3
    Type_A_r<vector<OCP_DBL>> cpl4; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^4
    Type_A_r<vector<OCP_DBL>> cpg1; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F
    Type_A_r<vector<OCP_DBL>> cpg2; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^2
    Type_A_r<vector<OCP_DBL>> cpg3; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^3
    Type_A_r<vector<OCP_DBL>> cpg4; ///< coefficients in the component liquid enthalpy
                                    ///< calculations, Btu/lbmol/F^4
    Type_A_r<vector<OCP_DBL>>
        hvapr; ///< coefficients in the component gas enthalpy calculations, Btu/lbmol
    Type_A_r<vector<OCP_DBL>>
        hvr; ///< coefficients in the vaporization enthalpy calculations
    Type_A_r<vector<OCP_DBL>>
        ev; ///< coefficients in the vaporization enthalpy calculations
    Type_A_r<vector<OCP_DBL>>
        avisc; ///< coefficients in water and oil viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>>
        bvisc; ///< coefficients in water and oil viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>>
        avg; ///< coefficients Ak in gas viscosity correlation formulae
    Type_A_r<vector<OCP_DBL>>
        bvg; ///< coefficients Bk in gas viscosity correlation formulae
    /// viscosity-versus-temperature dependence,  This table can specify the
    /// viscosity-versus-temperature-pressure dependence.
    TableSet viscTab;

    vector<OCP_DBL> Pref; ///< reference pressure
    vector<OCP_DBL> Tref; ///< reference temperature

    // Miscibility
    OCP_BOOL miscible{OCP_FALSE}; ///< Miscible treatment of hydrocarbons, used in
                                  ///< compositional Model.

    vector<string> SSMparamSTA; ///< Params for Solving Phase Spliting with SSM
    vector<string> NRparamSTA;  ///< Params for Solving Phase Spliting with NR
    vector<string> SSMparamSP;  ///< Params for Solving Phase Spliting with SSM
    vector<string> NRparamSP;   ///< Params for Solving Phase Spliting with NR
    vector<string> RRparam;     ///< Params for Solving Rachford-Rice equations
};

class Miscstr
{
public:
    vector<OCP_DBL> surTenRef;
    // 0th entry: reference surface tension - flow is immiscible when the surface
    // tension is greater than or equal to this value. 1th entry: maximum surface
    // tension expected, it should be greater than surTenRef. 2th entry: maximum surface
    // tension used to scale the input capillary pressure curves.
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
    vector<OCP_DBL> coord; ///< TODO: Add Doxygen.
    vector<OCP_DBL> zcorn; ///< TODO: Add Doxygen.

    // RockParam
    vector<OCP_DBL>   ntg;     ///< Net to gross for each grid.
    vector<OCP_DBL>   poro;    ///< Porosity for each grid.
    vector<OCP_DBL>   permX;   ///< Permeability along the x - direction for each grid.
    vector<OCP_DBL>   permY;   ///< Permeability along the y-direction for each grid.
    vector<OCP_DBL>   permZ;   ///< Permeability along the z-direction for each grid.
    vector<OCP_DBL>   thconr;  ///< Rock ifThermal conductivity.
    OCP_DBL           rsTemp;  ///< Temperature for reservoir.
    vector<RockParam> rockSet; ///< a set of rock
    HLoss             hLoss;   ///< Heat loss property
    Miscstr           miscstr; ///< reference Miscibility surface tension

    // If P and Ni are given, then calculation of initial equilibration is unneeded.
    vector<OCP_DBL> P;    ///< Initial pressure of components in each grid.
    vector<OCP_DBL> Ni;   ///< Initial moles of components in each grid.
    vector<OCP_DBL> Swat; ///< Initial water saturation in each grid.
    OCP_BOOL        ScalePcow{OCP_FALSE}; ///< whether Pcow should be scaled.

    // phase property
    Type_A_r<OCP_DBL> density; ///< Density of oil, water, gas in standard conditions.
    Type_A_r<OCP_DBL> gravity; ///< Gravity of oil, water, gas in standard conditions.
    OCP_DBL           thcono{24}; ///< oil ifThermal conductivity
    OCP_DBL           thcong{24}; ///< gas ifThermal conductivity
    OCP_DBL           thconw{24}; ///< water ifThermal conductivity

    // Models
    OCP_BOOL blackOil{OCP_FALSE}; ///< If ture, blackoil model will be used.
    OCP_BOOL comps{OCP_FALSE};    ///< If OCP_TRUE, compositional model will be used.
    OCP_BOOL thermal{OCP_FALSE};  ///< If OCP_TRUE, ifThermal model will be used.
    OCP_BOOL oil{OCP_FALSE};      ///< If OCP_TRUE, oil phase could exist.
    OCP_BOOL gas{OCP_FALSE};      ///< If OCP_TRUE, gas phase could exist.
    OCP_BOOL water{OCP_FALSE};    ///< If OCP_TRUE, water phase could exist.
    OCP_BOOL disGas{OCP_FALSE}; ///< If OCP_TRUE, dissolve gas could exist in oil phase.

    ComponentParam comsParam; ///< information for components

    // SAT Region & PVT Region
    USI               NTSFUN{1}; ///< Num of SAT regions.
    USI               NTPVT{1};  ///< Num of PVT regions.
    USI               NTROOC{1}; ///< Num of Rock regions.
    Type_A_r<OCP_DBL> SATNUM;    ///< Records the index of SAT region for each grid.
    Type_A_r<OCP_DBL> PVTNUM;    ///< Records the index of PVT region for each grid.
    Type_A_r<OCP_DBL> ACTNUM;    ///< Records the index of Active region for each grid.
    Type_A_r<OCP_DBL> ROCKNUM;   ///< Records the index of ROCK region for each grid.

    // Saturation tables & bubble point pressure
    TableSet SWFN_T; ///< Table set of SWFN.
    TableSet SWOF_T; ///< Table set of SWOF.
    TableSet SGFN_T; ///< Table set of SGFN.
    TableSet SGOF_T; ///< Table set of SGOF.
    TableSet SOF3_T; ///< Table set of SOF3.
    TableSet PBVD_T; ///< Table set of PBVD.
    // initial zi vs depth
    TableSet        ZMFVD_T;  ///< Table set of ZMFVD
    TableSet        TEMPVD_T; ///< Table set of TEMPVD
    vector<OCP_DBL> EQUIL;    ///< See ParamEQUIL.

    // PVT properties
    USI numPhase; ///< Number of phases
    USI numCom; ///< Number of components(hydrocarbon components), used in Compositional
                ///< Model when input
    TableSet PVCO_T; ///< Table set of PVCO.
    TableSet PVDO_T; ///< Table set of PVDO.
    TableSet PVDG_T; ///< Table set of PVDG.
    TableSet PVTW_T; ///< Table set of PVTW.

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
    void
    MultiplyVal(vector<OCP_DBL>& obj, const OCP_DBL& val, const vector<USI>& index);

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
    /// Input Rock information for ifThermal model
    void InputROCKT(ifstream& ifs);
    /// Input heat loss property for overburden rock and underburden rock
    void InputHLOSS(ifstream& ifs);

    /// Input the Miscibility information
    void InputMISCSTR(ifstream& ifs);

    /// Input the reference gravity of oil, water, and air in standard condition.
    void InputGRAVITY(ifstream& ifs);

    /// Input the reference density of oil, water, and air in standard condition.
    void InputDENSITY(ifstream& ifs);

    /// Input the phase ifThermal conductivity
    void InputTHCON(ifstream& ifs, const string& keyword);

    /// EQUIL contains initial information of reservoir; see ParamEQUIL.
    void InputEQUIL(ifstream& ifs);

    // SATNUM & PVTNUM  -- Region
    /// TABDIMS contains the num of saturation region and PVT region.
    void InputTABDIMS(ifstream& ifs);

    /// Input the keyword: SATNUM and PVTNUM.
    void InputRegion(ifstream& ifs, const string& keyword);

    // Input ComponentParam
    // Basic params
    void InputCNAMES(ifstream& ifs) { comsParam.InputCNAMES(ifs); };
    void InputCOMPONENTS(ifstream& ifs, const string& keyword)
    {
        comsParam.InputCOMPONENTS(ifs, keyword);
    }
    void InputLBCCOEF(ifstream& ifs) { comsParam.InputLBCCOEF(ifs); }
    void InputBIC(ifstream& ifs) { comsParam.InputBIC(ifs); };
    void InputVISCTAB(ifstream& ifs) { comsParam.InputVISCTAB(ifs); }
    void InputRefPR(ifstream& ifs, const string& keyword)
    {
        comsParam.InputRefPR(ifs, keyword);
    };

    // Method params
    void InputSSMSTA(ifstream& ifs) { comsParam.InputSSMSTA(ifs); };
    void InputNRSTA(ifstream& ifs) { comsParam.InputNRSTA(ifs); };
    void InputSSMSP(ifstream& ifs) { comsParam.InputSSMSP(ifs); };
    void InputNRSP(ifstream& ifs) { comsParam.InputNRSP(ifs); };
    void InputRR(ifstream& ifs) { comsParam.InputRR(ifs); };

    // check
    /// Check the reservoir param from input file.
    void CheckParam();

    /// Check the size of properties of grids.
    void CheckGrid();

    /// Check Rock
    void CheckRock();

    /// Check if keyword EQUIL is given.
    void CheckEQUIL() const;

    /// Check if density and gravity are both input, only one of them is needed.
    void CheckDenGra() const;

    /// Check existence of disgas, it could only exist when both oil and gas exist.
    void CheckPhase() const;

    /// Check existence of PVTtable and SATtable.
    void CheckPhaseTab() const;

    /// Check if each grid is assigned to an area or all defaulted.
    void CheckRegion() const;

    /// (Todo) Initialization of equilibration of only one region is realized.
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
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/