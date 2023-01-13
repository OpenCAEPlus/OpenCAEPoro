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
#include "Output4Vtk.hpp"
#include "ParamOutput.hpp"
#include "Reservoir.hpp"
#include "UtilOutput.hpp"
#include "UtilTiming.hpp"

#ifdef USE_METIS
#include "metis.h"
#endif

using namespace std;

/// 3D coordinate representation in OpenCAEPoro
class OCPIJK
{
public:
    OCPIJK() = default;
    OCPIJK(const USI& i, const USI& j, const USI& k)
        : I(i)
        , J(j)
        , K(k){};
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
template <typename T>
class OCPType_Sum
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
    OCP_BOOL  activity{OCP_FALSE};
    vector<T> obj;
    vector<USI>
        index; ///< Records the index of bulk or well, whose properties will be printed
};

/// The SumItem class is an auxiliary structure storing summary data to output.
class SumItem
{
public:
    SumItem(const string& item,
            const string& obj,
            const string& unit,
            const string& type)
        : Item(item)
        , Obj(obj)
        , Unit(unit)
        , Type(type){};
    string          Item;
    string          Obj;
    string          Unit;
    string          Type;
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
    vector<SumItem> Sumdata; ///< Contains all information to be printed.

    OCP_BOOL FPR{OCP_FALSE};  ///< Field average Pressure.
    OCP_BOOL FTR{OCP_FALSE};  ///< Field average Temperature.
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
    void PrintFastReview(const string& dir) const;

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
    friend class Out4RPT;
    friend class Out4VTK;

public:
    void SetBasicGridProperty(const BasicGridPropertyParam& param);

private:
    OCP_BOOL PRE{OCP_FALSE};  ///< Pressure of grids.
    OCP_BOOL PGAS{OCP_FALSE}; ///< Gas pressure of grids.
    OCP_BOOL PWAT{OCP_FALSE}; ///< Water pressure of grids.
    OCP_BOOL SOIL{OCP_FALSE}; ///< Oil saturation of grids.
    OCP_BOOL SGAS{OCP_FALSE}; ///< Gas saturation of grids.
    OCP_BOOL SWAT{OCP_FALSE}; ///< Water saturation of grids.
    OCP_BOOL DENO{OCP_FALSE}; ///< Oil density of grids.
    OCP_BOOL DENG{OCP_FALSE}; ///< Gas density of grids.
    OCP_BOOL DENW{OCP_FALSE}; ///< Water density of grids.
    OCP_BOOL KRO{OCP_FALSE};  ///< Oil relative permeability of grids.
    OCP_BOOL KRG{OCP_FALSE};  ///< Gas relative permeability of grids.
    OCP_BOOL KRW{OCP_FALSE};  ///< Water relative permeability of grids.
    OCP_BOOL BOIL{OCP_FALSE}; ///< Oil reservoir molar densities of grids.
    OCP_BOOL BGAS{OCP_FALSE}; ///< Gas reservoir molar densities of grids.
    OCP_BOOL BWAT{OCP_FALSE}; ///< Water reservoir molar densities of grids.
    OCP_BOOL VOIL{OCP_FALSE}; ///< Oil viscosity of grids.
    OCP_BOOL VGAS{OCP_FALSE}; ///< Gas viscosity of grids.
    OCP_BOOL VWAT{OCP_FALSE}; ///< Water viscosity of grids.
    OCP_BOOL XMF{OCP_FALSE};  ///< liquid component mole fractions.
    OCP_BOOL YMF{OCP_FALSE};  ///< gas component mole fractions.
    OCP_BOOL PCW{OCP_FALSE};  ///< capillary pressure: Po - Pw.
};

/// Collect more detailed information of each time step.
class Out4RPT
{
public:
    void InputParam(const OutputRPTParam& RPTparam);
    void Setup(const string& dir, const Reservoir& reservoir);
    void PrintRPT(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
    template <typename T>
    void PrintRPT_Scalar(ofstream&              ifs,
                         const string&          dataName,
                         const OCP_DBL&         days,
                         const T*               gridVal,
                         const USI&             gap,
                         const vector<GB_Pair>& gbPair,
                         const bool&            useActive,
                         const OCP_DBL&         alpha = 1.0) const;
    void GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const;

private:
    OCP_BOOL          useRPT{OCP_FALSE};
    OCP_USI           numGrid;
    OCP_USI           nx;
    OCP_USI           ny;
    USI               IJKspace;
    BasicGridProperty bgp;
};

template <typename T>
void Out4RPT::PrintRPT_Scalar(ofstream&              myRPT,
                              const string&          dataName,
                              const OCP_DBL&         days,
                              const T*               gridVal,
                              const USI&             gap,
                              const vector<GB_Pair>& gbPair,
                              const bool&            useActive,
                              const OCP_DBL&         alpha) const
{
    USI     I, J, K;
    OCP_USI bId;

    myRPT << OCP_SEP01(50) << "\n";
    myRPT << dataName << "                   " << fixed << setprecision(3) << days
          << "  DAYS";

    if (useActive) {
        for (OCP_USI n = 0; n < numGrid; n++) {
            if (n % nx == 0) myRPT << "\n";
            if (n % (nx * ny) == 0) myRPT << "\n\n";

            if (n % nx == 0) {
                GetIJKGrid(I, J, K, n);
                myRPT << GetIJKformat("*", to_string(J), to_string(K), IJKspace);
            }

            if (gbPair[n].IsAct()) {
                bId = gbPair[n].GetId();
                myRPT << setw(10) << fixed << setprecision(3)
                      << gridVal[bId * gap] * alpha;
            } else {
                myRPT << setw(10) << " --- ";
            }
        }
    } else {
        for (OCP_USI n = 0; n < numGrid; n++) {
            if (n % nx == 0) myRPT << "\n";
            if (n % (nx * ny) == 0) myRPT << "\n\n";

            if (n % nx == 0) {
                GetIJKGrid(I, J, K, n);
                myRPT << GetIJKformat("*", to_string(J), to_string(K), IJKspace);
            }
            myRPT << setw(10) << fixed << setprecision(3) << gridVal[n * gap] * alpha;
        }
    }

    myRPT << "\n\n\n";
}

class Out4VTK
{
public:
    void InputParam(const OutputVTKParam& VTKParam);
    void Setup(const string& dir, const Reservoir& rs, const USI& ndates);
    void PrintVTK(const string& dir, const Reservoir& rs, const OCP_DBL& days) const;
    OCP_BOOL IfOutputVTK() const { return useVTK; }

private:
    OCP_BOOL          useVTK{OCP_FALSE}; ///< If use vtk
    mutable USI       index{0};          ///< Index of output file
    BasicGridProperty bgp;               ///< Basic grid information
    Output4Vtk        out4vtk;           ///< Output for vtk

    // test for Parallel version
#ifdef USE_METIS
    mutable class MyMetisTest
    {
        friend class Out4VTK;

    private:
        OCP_BOOL useMetis{OCP_FALSE};

        OCP_USI nb;
        OCP_USI nw;
        OCP_USI ng;

        mutable idx_t         nvtxs;
        mutable idx_t         nedges;
        mutable idx_t         ncon{1};
        mutable idx_t         nparts{8};
        mutable vector<idx_t> xadj;
        mutable vector<idx_t> adjncy;
        mutable vector<idx_t> adjwgt;
        mutable vector<idx_t> vwgt;
        mutable vector<idx_t> part;
        mutable vector<USI>   partitions;
        const idx_t           MAXWEIGHT = 100000000;

    public:
        void Setup(const Reservoir& rs)
        {
            if (!useMetis) return;
            nb = rs.GetBulkNum();
            nw = rs.GetWellNum();
            ng = rs.grid.GetGridNum();

            nvtxs                           = nb + nw;
            nedges                          = 0;
            vector<vector<OCP_USI>> tmpConn = rs.conn.neighbor;
            vector<OCP_USI>         tmp;
            for (USI w = 0; w < rs.allWells.numWell; w++) {
                tmp.clear();
                for (USI p = 0; p < rs.allWells.GetWellPerfNum(w); p++) {
                    OCP_USI n = rs.allWells.wells[w].PerfLocation(p);
                    tmpConn[n].push_back(nb + w);
                    tmp.push_back(n);
                }
                tmpConn.push_back(tmp);
                nedges += rs.allWells.GetWellPerfNum(w) * 2;
            }

            const vector<vector<OCP_USI>>& myConn = tmpConn;

            for (OCP_USI n = 0; n < nb; n++) {
                nedges += rs.conn.neighborNum[n];
            }
            nedges -= nb;

            // generate xadj  adjncy  adjwgt
            xadj.resize(nvtxs + 1);
            xadj[0] = 0;
            adjncy.reserve(nedges);
            adjwgt.reserve(nedges);
            // Bulk
            for (OCP_USI n = 0; n < nb; n++) {
                xadj[n + 1] = xadj[n] + myConn[n].size() - 1;
                USI i       = 0;
                // left
                for (i = 0; i < rs.conn.selfPtr[n]; i++) {
                    adjncy.push_back(myConn[n][i]);
                    adjwgt.push_back(1);
                }
                // right
                for (i = rs.conn.selfPtr[n] + 1; i < rs.conn.neighborNum[n]; i++) {
                    adjncy.push_back(myConn[n][i]);
                    adjwgt.push_back(1);
                }
                // well
                for (i = rs.conn.neighborNum[n]; i < myConn[n].size(); i++) {
                    adjncy.push_back(myConn[n][i]);
                    adjwgt.push_back(MAXWEIGHT);
                }
            }
            // Well
            for (USI w = 0; w < nw; w++) {
                xadj[nb + w + 1] += xadj[nb + w] + myConn[nb + w].size();
                for (auto v : myConn[nb + w]) {
                    adjncy.push_back(v);
                    adjwgt.push_back(MAXWEIGHT);
                }
            }
            // Setup active/inactive grid flags
            // inactive grids' flags equal nparts
            part.resize(nvtxs, 0);
            partitions.resize(ng + nw, nparts);
        }

        void MyPartitionFunc(decltype(METIS_PartGraphKway)* METIS_PartGraphFunc) const
        {

            if (!useMetis) return;
            idx_t objval;
            int   ret = METIS_PartGraphFunc(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                            vwgt.data(), NULL, adjwgt.data(), &nparts,
                                            NULL, NULL, NULL, &objval, part.data());

            if (ret != rstatus_et::METIS_OK) {
                OCP_ABORT("METIS ERROR");
            }
            cout << "METIS_OK" << endl;
            cout << "objval: " << objval << endl;
        }

        void SetPartitions(const vector<USI> b2g) const
        {
            if (!useMetis) return;
            // bulk
            OCP_USI n = 0;
            for (n = 0; n < nb; n++) {
                partitions[b2g[n]] = part[n];
            }
            // well
            for (USI w = 0; w < nw; w++) {
                partitions[ng + w] = part[n];
                n++;
            }
        }
    } metisTest;
#endif // USE_METIS
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
    void     InputParam(const ParamOutput& paramOutput);
    void     Setup(const Reservoir& reservoir, const OCPControl& ctrl);
    void     SetVal(const Reservoir& reservoir, const OCPControl& ctrl);
    void     PrintInfo() const;
    void     PrintInfoSched(const Reservoir&  rs,
                            const OCPControl& ctrl,
                            const OCP_DBL&    time) const;
    OCP_BOOL IfOutputVTK() const { return out4VTK.IfOutputVTK(); }

private:
    string       workDir;
    Summary      summary;
    CriticalInfo crtInfo;
    Out4RPT      out4RPT;
    Out4VTK      out4VTK;

    mutable OCP_DBL outputTime{0}; ///< Total time for main output
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