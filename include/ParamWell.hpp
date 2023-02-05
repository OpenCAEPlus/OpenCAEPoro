/*! \file    ParamWell.hpp
 *  \brief   ParamWell class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMWELL_HEADER__
#define __PARAMWELL_HEADER__

// Standard header files
#include <cassert>
#include <fstream>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"

using namespace std;

/// WellOptParam contains all the param used for well operation. for a well, all the
/// params in it can be changed over the time.
class WellOptParam
{
public:
    WellOptParam(string intype, vector<string>& vbuf);
    // WCONINJE & WCONPROD
    string type;      ///< Type of well, injection or production?
    string fluidType; ///< Type of fluid into the injection well. (injection well only)
    string state;     ///< State of well, open or close?
    string optMode;   ///< Mode of well, Rate or BHP?

    OCP_DBL maxRate; ///< Maximum allowable flow rate into/out the well.
    OCP_DBL maxBHP;  ///< Maximum allowable pressure in the injection well.
    OCP_DBL minBHP;  ///< Minimum allowable pressure in the production well.
    OCP_DBL injTemp; ///< Temperature of injected fluid.
};

/// WellOptPair contains two parts, one is the operation mode of well, the other is the
/// beginning time at which the operation is applied. The beginning time is represented
/// by an index in critical time.
class WellOptPair
{
public:
    WellOptPair(USI i, string type, vector<string>& vbuf)
        : d(i)
        , opt(type, vbuf){};
    USI          d;
    WellOptParam opt;
};

/// TODO: Add Doxygen
class WellParam
{
public:
    WellParam(vector<string>& info);
    // static infomation
    // WELSPECS
    string  name;           ///< Name of Well
    string  group{"FEILD"}; ///< Group the well belongs to.
    USI     I;              ///< I index of well.
    USI     J;              ///< J index of well.
    OCP_DBL depth{-1.0};    ///< Depth of well.

    // COMPDAT ---- for all perforation.
    vector<USI>     I_perf;   ///< I-index of perforation in grid.
    vector<USI>     J_perf;   ///< J-index of perforation in grid.
    vector<USI>     K_perf;   ///< K-index of perforation in grid.
    vector<OCP_DBL> WI;       ///< Transmissibility connection factor.
    vector<OCP_DBL> diameter; ///< Diameter of perforations.
    vector<OCP_DBL> kh;
    vector<OCP_DBL> skinFactor;               ///< Skin factor.
    vector<string>  direction;                ///< Direction of perforations.
    OCP_BOOL        ifUseUnweight{OCP_FALSE}; ///< If use unweighted injector.

    // dynamic infomation
    vector<WellOptPair> optParam;
};

/// Describe the molar fraction of components of fluid injected to reservoir from INJ.
class Solvent
{
public:
    Solvent() = default;
    Solvent(const vector<string>& vbuf);
    string          name;
    vector<OCP_DBL> comRatio;
};

/// ParamWell is an internal structure used to stores the information of wells from
/// input files. It is an intermediate interface and independent of the main simulator.
/// After all file inputting finishes, the params in it will pass to corresponding
/// modules.
class ParamWell
{
public:
    vector<WellParam> well;         ///< Contains all the information of wells.
    vector<OCP_DBL>   criticalTime; ///< Records the critical time given by users.
    vector<Solvent>   solSet;       ///< Sets of Solvent.
    OCP_DBL           Psurf{PRESSURE_STD};    ///< Pressure in surface condition.
    OCP_DBL           Tsurf{TEMPERATURE_STD}; ///< Temperature in surface condition.

    /// Initialize the inputting the params of wells.
    void Init() { InitTime(); };
    /// Initialize the critical time.
    void InitTime() { criticalTime.push_back(0); };
    /// Input the well keyword WELSPECS. WELSPECS defines wells including well name,
    /// well location, well depth and son on.
    void InputWELSPECS(ifstream& ifs);
    /// Input the well keyword COMPDAT. COMPDAT contains the information of perforations
    /// of wells, for example the location, the trans or directions.
    void InputCOMPDAT(ifstream& ifs);
    /// Input the well keyword WCONINJE. WCONINJE describes the initial operation mode
    /// of injection well.
    void InputWCONINJE(ifstream& ifs);
    /// Input the well keyword WCONPROD. WCONPROD describes the initial operation mode
    /// of production well.
    void InputWCONPROD(ifstream& ifs);
    /// Input the keyword: TSTEP. TSTEP is used to divide the simulation time and
    /// ususally the time point is critical, at which for example, the operation mode of
    /// well will change. So the params of solving equations could be adjusted
    /// correspondingly.
    void InputTSTEP(ifstream& ifs);
    /// Input the well keyword: WELTARG. WELTARG is used to change the operation mode of
    /// well anytime. For example, the oil production rate is changed from 1000 stb/day
    /// to 1500 stb/day at the 100th day.
    void InputWELTARG(ifstream& ifs);
    /// Input the temperature of injected fluid
    void InputWTEMP(ifstream& ifs);
    /// Input injector type -- MOBWEIGHT(defaulted) or UNWEIGHT
    void InputUNWEIGHT(ifstream& ifs);
    /// Input well keyword: Solvent. It describes the molar fraction of components of
    /// fluid injected to reservoir from INJ.
    void InputWELLSTRE(ifstream& ifs);
    /// Input surface pressure
    void InputPSURF(ifstream& ifs);
    /// Input surface temperature
    void InputTSURF(ifstream& ifs);

    // check
    /// Check if wrong params are input.
    void CheckParam(const OCP_BOOL& boModel) const;
    /// Check if params of Perforation is wrong.
    void CheckPerf() const;
    void CheckINJFluid() const;
};

#endif /* end if __PARAMWELL_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/