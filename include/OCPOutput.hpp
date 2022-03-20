/*! \file    OCPOutput.hpp
 *  \brief   OCPOutput class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_OUTPUT_HEADER__
#define __OCP_OUTPUT_HEADER__

// Standard header files
#include <iomanip>
#include <iostream>

// OpenCAEPoro header files
#include "OCPControl.hpp"
#include "ParamOutput.hpp"
#include "Reservoir.hpp"

using namespace std;

/// 3D coordinate representation in OpenCAEPoro
class OCPIJK
{
public:
    OCPIJK() = default;
    OCPIJK(const USI& i, const USI& j, const USI& k)
        : I(i)
        , J(j)
        , K(k) {};
    OCPIJK(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
    };
    OCPIJK& operator=(const COOIJK& src)
    {
        I = src.I;
        J = src.J;
        K = src.K;
        return *this;
    }
    USI I, J, K;
};

/// TODO: Add Doxygen
template <typename T> class OCPType_Sum
{
public:
    OCPType_Sum& operator=(const Type_A_o& src)
    {
        activity = src.activity;
        obj      = src.obj;
        return *this;
    }
    OCPType_Sum& operator=(const Type_B_o& src)
    {
        activity = src.activity;
        obj.assign(src.obj.begin(), src.obj.end());
        return *this;
    }
    bool      activity{false};
    vector<T> obj;
    vector<USI>
        index; ///< Records the index of bulk or well, whose properties will be printed
};

/// The SumPair class is an auxiliary structure storing summary data to output.
class SumPair
{
public:
    SumPair(const string& item, const string& obj, const string& unit)
        : Item(item)
        , Obj(obj)
        , Unit(unit){};
    string          Item;
    string          Obj;
    string          Unit;
    vector<OCP_DBL> val;
};

/// The Summary class manages the output in the summary file.
//  Note: It contains the most interested information in each time step, which usually
//  will be convert to figures for later analysis.
class Summary
{
public:
    /// TODO: Add Doxygen
    void InputParam(const OutputSummary& summary_param);

    /// TODO: Add Doxygen
    void Setup(const Reservoir& reservoir, const OCP_DBL& totalTime);

    /// TODO: Add Doxygen
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl);

    /// Write output information to a file.
    void PrintInfo(const string& dir) const;

private:
    vector<SumPair> Sumdata; ///< Contains all information to be printed.

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

    OCPType_Sum<string> WOPR; ///< Well oil production rate.
    OCPType_Sum<string> WOPT; ///< Well total oil production.
    OCPType_Sum<string> WGPR; ///< Well gas production rate.
    OCPType_Sum<string> WGPT; ///< Well total gas production.
    OCPType_Sum<string> WWPR; ///< Well water production rate.
    OCPType_Sum<string> WWPT; ///< Well total water production.
    OCPType_Sum<string> WGIR; ///< Well gas injection rate.
    OCPType_Sum<string> WGIT; ///< Well total gas injection.
    OCPType_Sum<string> WWIR; ///< Well water injection rate.
    OCPType_Sum<string> WWIT; ///< Well total water injection.
    OCPType_Sum<string> WBHP; ///< Well pressure.
    OCPType_Sum<string> DG;   ///< Pressure difference between wells and perforations.

    OCPType_Sum<OCPIJK> BPR;  ///< Bulk pressure.
    OCPType_Sum<OCPIJK> SOIL; ///< Oil saturation of bulk.
    OCPType_Sum<OCPIJK> SGAS; ///< Gas saturation of bulk.
    OCPType_Sum<OCPIJK> SWAT; ///< Water saturation of bulk.
};

/// Collect important information of each time step for fast review.
class CriticalInfo
{
public:
    /// TODO: Add Doxygen
    void Setup(const OCP_DBL& totalTime);

    /// TODO: Add Doxygen
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl);

    /// TODO: Add Doxygen
    void PrintInfo(const string& dir) const;

private:
    vector<OCP_DBL> time;  ///< TODO: Add Doxygen
    vector<OCP_DBL> dt;    ///< TODO: Add Doxygen
    vector<OCP_DBL> dPmax; ///< TODO: Add Doxygen
    vector<OCP_DBL> dVmax; ///< TODO: Add Doxygen
    vector<OCP_DBL> dSmax; ///< TODO: Add Doxygen
    vector<OCP_DBL> dNmax; ///< TODO: Add Doxygen
    vector<OCP_DBL> cfl;   ///< TODO: Add Doxygen
};

/// Collect more detailed information of each time step.
class DetailInfo
{
public:
    void InputParam(const OutputDetail& detail_param);
    void Setup(const string& dir);
    void PrintInfo(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;

private:
    bool PRE{false};  ///< Pressure of grids.
    bool PGAS{false}; ///< Gas pressure of grids.
    bool PWAT{false}; ///< Water pressure of grids.
    bool SOIL{false}; ///< Oil saturation of grids.
    bool SGAS{false}; ///< Gas saturation of grids.
    bool SWAT{false}; ///< Water saturation of grids.
    bool DENO{false}; ///< Oil density of grids.
    bool DENG{false}; ///< Gas density of grids.
    bool DENW{false}; ///< Water density of grids.
    bool KRO{false}; ///< Oil relative permeability of grids.
    bool KRG{false}; ///< Gas relative permeability of grids.
    bool KRW{false}; ///< Water relative permeability of grids.
    bool BOIL{false}; ///< Oil reservoir molar densities of grids.
    bool BGAS{false}; ///< Gas reservoir molar densities of grids.
    bool BWAT{false}; ///< Water reservoir molar densities of grids.
    bool VOIL{false}; ///< Oil viscosity of grids.
    bool VGAS{false}; ///< Gas viscosity of grids.
    bool VWAT{false}; ///< Water viscosity of grids.
    bool XMF{false}; ///< liquid component mole fractions.
    bool YMF{false}; ///< gas component mole fractions.
};

/// The OCPOutput class manages different kinds of ways to output information.
//  Note: The most commonly used is the summary file, which usually gives the
//  information of bulks and wells in each time step, such as average pressure, oil
//  production rate of wells. If other information at critical dates is of interest, you
//  can chose the PRT file (TODO). Also, some infomation will be printed on the screen
//  at the critical dates to make sure the program is at the right way.
class OCPOutput
{
    friend class OpenCAEPoro;

public:
    void InputParam(const ParamOutput& paramOutput);
    void Setup(const Reservoir& reservoir, const OCPControl& ctrl);
    void SetVal(const Reservoir& reservoir, const OCPControl& ctrl);
    void PrintInfo() const;
    void PrintInfoSched(const Reservoir& rs, const OCPControl& ctrl,
                        const OCP_DBL& time) const;

private:
    string       wordDir;
    Summary      summary;
    CriticalInfo crtInfo;
    DetailInfo   dtlInfo;
};

#endif /* end if __OCPOUTPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/