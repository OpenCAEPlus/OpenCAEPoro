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
#include "UtilOutput.hpp"
#include "Output4Vtk.hpp"

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
    OCP_BOOL      activity{OCP_FALSE};
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

    OCP_BOOL FPR{OCP_FALSE};  ///< Field average Pressure.
    OCP_BOOL FOPR{OCP_FALSE}; ///< Field oil production rate.
    OCP_BOOL FOPT{OCP_FALSE}; ///< Field total oil production.
    OCP_BOOL FGPR{OCP_FALSE}; ///< Field gas production rate.
    OCP_BOOL FGPt{OCP_FALSE}; ///< Field total gas production.
    OCP_BOOL FWPR{OCP_FALSE}; ///< Field water production rate.
    OCP_BOOL FWPT{OCP_FALSE}; ///< Field total water production.
    OCP_BOOL FGIR{OCP_FALSE}; ///< Field gas injection rate.
    OCP_BOOL FGIT{OCP_FALSE}; ///< Field total gas injection.
    OCP_BOOL FWIR{OCP_FALSE}; ///< Field water injection rate.
    OCP_BOOL FWIT{OCP_FALSE}; ///< Field total water injection.

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

/// Basic grid properties for output
class BasicGridProperty
{
    friend class DetailInfo;
    friend class Out4VTK;

private:
    OCP_BOOL PRE{ OCP_FALSE };  ///< Pressure of grids.
    OCP_BOOL PGAS{ OCP_FALSE }; ///< Gas pressure of grids.
    OCP_BOOL PWAT{ OCP_FALSE }; ///< Water pressure of grids.
    OCP_BOOL SOIL{ OCP_FALSE }; ///< Oil saturation of grids.
    OCP_BOOL SGAS{ OCP_FALSE }; ///< Gas saturation of grids.
    OCP_BOOL SWAT{ OCP_FALSE }; ///< Water saturation of grids.
    OCP_BOOL DENO{ OCP_FALSE }; ///< Oil density of grids.
    OCP_BOOL DENG{ OCP_FALSE }; ///< Gas density of grids.
    OCP_BOOL DENW{ OCP_FALSE }; ///< Water density of grids.
    OCP_BOOL KRO{ OCP_FALSE }; ///< Oil relative permeability of grids.
    OCP_BOOL KRG{ OCP_FALSE }; ///< Gas relative permeability of grids.
    OCP_BOOL KRW{ OCP_FALSE }; ///< Water relative permeability of grids.
    OCP_BOOL BOIL{ OCP_FALSE }; ///< Oil reservoir molar densities of grids.
    OCP_BOOL BGAS{ OCP_FALSE }; ///< Gas reservoir molar densities of grids.
    OCP_BOOL BWAT{ OCP_FALSE }; ///< Water reservoir molar densities of grids.
    OCP_BOOL VOIL{ OCP_FALSE }; ///< Oil viscosity of grids.
    OCP_BOOL VGAS{ OCP_FALSE }; ///< Gas viscosity of grids.
    OCP_BOOL VWAT{ OCP_FALSE }; ///< Water viscosity of grids.
    OCP_BOOL XMF{ OCP_FALSE }; ///< liquid component mole fractions.
    OCP_BOOL YMF{ OCP_FALSE }; ///< gas component mole fractions.
    OCP_BOOL PCW{ OCP_FALSE }; ///< capilary pressure: Po - Pw.
};

/// Collect more detailed information of each time step.
class DetailInfo
{
public:
    void InputParam(const OutputDetail& detail_param);
    void Setup(const string& dir);
    void PrintInfo(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;

private:
    OCP_BOOL    useRPT{ OCP_FALSE };
    BasicGridProperty bgp;
};

class Out4VTK
{
public:
    void PrintVTK(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
private:

    OCP_BOOL     useVtk{ OCP_TRUE };
    mutable USI  index{ 0 };   ///< index of output file
    BasicGridProperty bgp;
    Output4Vtk  out4vtk;
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
    string       workDir;
    Summary      summary;
    CriticalInfo crtInfo;
    DetailInfo   dtlInfo;

    Out4VTK      out4VTK;
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