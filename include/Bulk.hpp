/*! \file    Bulk.hpp
 *  \brief   Bulk class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULK_HEADER__
#define __BULK_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoro header files
#include "DenseMat.hpp"
#include "FlowUnit.hpp"
#include "Grid.hpp"
#include "LinearSystem.hpp"
#include "Mixture.hpp"
#include "MixtureBO.hpp"
#include "MixtureComp.hpp"
#include "OCPConst.hpp"
#include "OCPStructure.hpp"
#include "ParamReservoir.hpp"

using namespace std;

/// Initial reservoir infomation for calculating initial equilibration.
//  Note: This class includes reference depth and pressure at it, depth of contacts
//  between phases, and capillary pressure at phase contact surfaces.
class ParamEQUIL
{
    friend class Bulk;

private:
    OCP_DBL Dref; ///< The reference depth
    OCP_DBL Pref; ///< Pressure at the reference depth
    OCP_DBL DOWC; ///< The depth of oil-water contact surface
    OCP_DBL DGOC; ///< The depth of gas-oil contact surface
    OCP_DBL PcOW; ///< Capillary pressure at oil-water contact Pcow = Po - Pw
    OCP_DBL PcGO; ///< capillary pressure at gas-oil contact Pcgo = Pg - Po

    OCPTable PBVD; ///< PBVD Table: bubble point pressere vs depth
};

/// Physical information of each active reservoir bulk.
//  Note: Bulk contains main physical infomation of active grids. It describes the
//  actural geometric domain for simulating. Variables are stored bulk by bulk, and then
//  phase by phase, then component by component. The bulks are ordered in the alphabetic
//  order, i.e. the X-axis indices first, followed by the Y-axis and Z-axis indices.
//  Operations on each bulk are also defined here.
class Bulk
{
    friend class BulkConn;
    friend class Well;
    friend class DetailInfo;

    // temp
    friend class OCP_FIM;
    friend class Reservoir;
    friend class OCP_AIMt;

    /////////////////////////////////////////////////////////////////////
    // For general usage
    /////////////////////////////////////////////////////////////////////

public:
    Bulk() = default;

    /// Input param from internal data structure ParamReservoir.
    void InputParam(ParamReservoir& rs_param);
    /// Allocate memory for bulk data of grid.
    void Setup(const Grid& myGrid);
    /// Calculate initial equilibrium for blkoil model according to EQUIL.
    void InitSjPcBo(const USI& tabrow);
    /// Calculate initial equilibrium for compositional model according to EQUIL.
    void InitSjPcComp(const USI& tabrow, const Grid& myGrid);
    /// Perform flash calculation with saturations.
    void InitFlash(const bool& flag = false);
    /// Perform flash calculation with Ni.
    void Flash();
    /// Perform flash calculation with Ni and calculate derivatives.
    void FlashDeriv();
    /// Pass values from Flash to Bulk after Flash calculation.
    void PassFlashValue(const OCP_USI& n);
    /// Pass derivative values from Flash to Bulk after Flash calculation.
    void PassFlashValueDeriv(const OCP_USI& n);
    /// Reset variables in flash calculations.
    void ResetFlash();

    /// Calculate relative permeability and capillary pressure with saturation.
    void CalKrPc();
    /// Calculate relative permeability and capillary pressure and their derivatives.
    void CalKrPcDeriv();
    /// Calculate volume of pore with pressure.
    void CalVpore();
    /// Calculate average pressure in reservoir.
    OCP_DBL CalFPR() const;
    /// Calculate max change of some variables.
    void CalMaxChange();

    /// Check if error occurs in Setup.
    void CheckSetup() const;
    /// Check initial pore volume.
    void CheckInitVpore() const;
    /// Check pore volume.
    void CheckVpore() const;
    /// Check if negative P occurs, return false if so.
    bool CheckP() const;
    /// Check if negative Ni occurs, return false if so.
    bool CheckNi() const;
    /// Check if relative volume error is out of range, return false if so.
    bool CheckVe(const OCP_DBL& Vlim) const;
    /// Check if the sum of saturations is one.
    void CheckSat() const;
    /// Check difference from last time step, for Debug and Test.
    void CheckDiff();

    /// Return the number of bulks.
    OCP_USI GetBulkNum() const { return numBulk; }
    /// Return the number of components.
    USI GetComNum() const { return numCom; }
    /// Return the number of phases.
    USI GetPhaseNum() const { return numPhase; }
    /// Return the mixture mode.
    USI GetMixMode() const;
    /// Return flash results (it has not been used by far).
    const vector<Mixture*>& GetMixture() const { return flashCal; }
    /// Return pressure of the n-th bulk.
    OCP_DBL GetP(const OCP_USI& n) const { return P[n]; }
    /// Return oil saturation of the n-th bulk.
    OCP_DBL GetSOIL(const OCP_USI& n) const
    {
        return S[n * numPhase + phase2Index[OIL]];
    }
    /// Return gas saturation of the n-th bulk.
    OCP_DBL GetSGAS(const OCP_USI& n) const
    {
        return S[n * numPhase + phase2Index[GAS]];
    }
    /// Return water saturation of the n-th bulk.
    OCP_DBL GetSWAT(const OCP_USI& n) const
    {
        return S[n * numPhase + phase2Index[WATER]];
    }
    /// Return dPmax.
    OCP_DBL GetdPmax() const { return dPmax; }
    /// Return dNmax.
    OCP_DBL GetdNmax() const { return dNmax; }
    /// Return dSmax.
    OCP_DBL GetdSmax() const { return dSmax; }
    /// Return dVmax.
    OCP_DBL GetdVmax() const { return dVmax; }
    
    // Reset phaseNum to the ones of the last time step.
    void ResetphaseNum() { phaseNum = lphaseNum; }
    // Reset minEigenSkip to the ones of the last time step.
    void ResetminEigenSkip() { minEigenSkip = lminEigenSkip; }
    // Reset flagSkip to the ones of the last time step.
    void ResetflagSkip() { flagSkip = lflagSkip; }
    // Reset flagSkip to the ones of the last time step.
    void ResetziSkip() { ziSkip = lziSkip; }
    // Reset flagSkip to the ones of the last time step.
    void ResetPSkip() { PSkip = lPSkip; }
    // Reser Ks to the ones of the last time step.
    void ResetKs() { Ks = lKs; }
    
    /// Reset Nt to the ones of the last time step.
    void ResetNt() { Nt = lNt; }

    /// Reset P to the ones of the last time step.
    void ResetP() { P = lP; }
    /// Reset Pj to the ones of the last time step.
    void ResetPj() { Pj = lPj; }
    /// Reset Ni to the ones of the last time step.
    void ResetNi() { Ni = lNi; }
    /// Reset Vp to the ones of the last time step.
    void ResetVp() { rockVp = rockLVp; }
    void CalSomeInfo(const Grid& myGrid) const;
    
    /// Allocate memory for WellbulkId
    void AllocateWellBulkId(const USI& n) { wellBulkId.reserve(n); }
    void ClearWellBulkId() { wellBulkId.clear(); }
    

private:
    /////////////////////////////////////////////////////////////////////
    // General variables
    /////////////////////////////////////////////////////////////////////
    OCP_USI numBulk;  ///< Number of bulks (active grids).
    USI     numPhase; ///< Number of phase.
    USI     numCom;   ///< Number of component.
    USI     numCom_1; ///< numCom - 1

    // Initial proportion of each component for EoS : numCom - 1, water is excluded.
    vector<OCP_DBL>   initZi;   ///< Initial proportion of each component.
    
    vector<OCP_DBL>   SwatInit; ///< Initial water saturation.
    bool SwatInitExist{false};    ///< If SwatInit has been given.
    vector<OCP_DBL>  ScaleValuePcow; ///< Scale valus for Pcow.
    bool ScalePcow{false};  ///< whether Pcow should be scaled.
    

    USI               PVTmode;  ///< Identify PVT mode in blackoil model.
    vector<USI>       PVTNUM;   ///< Identify PVT region in blackoil model: numBulk.
    vector<Mixture*>  flashCal; ///< Flash calculation class.
    USI               SATmode;  ///< Identify SAT mode.
    vector<USI>       SATNUM;   ///< Identify SAT region: numBulk.
    vector<FlowUnit*> flow;     ///< Vector for capillary pressure, relative perm.
 
    // Skip stability analysis
    vector<USI>       phaseNum;     ///< Num of hydrocarbon phase in each bulk
    vector<USI>       lphaseNum;    ///< last phaseNum 
    /// minimal eigenvalue used to determinined if skip the stability analysis
    vector<OCP_SIN>   minEigenSkip; 
    vector<bool>      flagSkip;
    vector<OCP_DBL>   ziSkip;
    vector<OCP_DBL>   PSkip;
    vector<OCP_SIN>   lminEigenSkip;
    vector<bool>      lflagSkip;
    vector<OCP_DBL>   lziSkip;
    vector<OCP_DBL>   lPSkip;

    // phase split calculation
    vector<OCP_DBL>   Ks;  ///< Equilibrium constant in phase split calculation
    vector<OCP_DBL>   lKs; ///< last Ks

    /////////////////////////////////////////////////////////////////////
    // Basic model information
    /////////////////////////////////////////////////////////////////////
    ParamEQUIL EQUIL;    ///< Initial Equilibration.
    bool       blackOil; ///< If true, blackoil model will be used.
    bool       comps;    ///< If true, compositional model will be used.
    bool       oil;      ///< If true, oil phase could exist.
    bool       gas;      ///< If true, gas phase could exist.
    bool       water;    ///< If true, water phase could exist.
    bool       disGas;   ///< If true, dissolve gas in live oil could exist.
    bool       miscible; ///< Miscible treatment of hydrocarbons, used in compositional Model.

    /////////////////////////////////////////////////////////////////////
    // Basic physical variables
    /////////////////////////////////////////////////////////////////////
    OCP_DBL         T;          ///< Reservoir temperature: numBulk.
    vector<OCP_DBL> Pb;         ///< Bubble point pressere: numBulk.
    vector<OCP_DBL> P;          ///< Pressure: numBulk.
    vector<OCP_DBL> Pj;         ///< Pressure of phase: numPhase*numBulk.
    vector<OCP_DBL> Pc;         ///< Capillary pressure of phase: numPhase*numBulk.
    vector<bool>    phaseExist; ///< Existence of phase: numPhase*numBulk.
    vector<OCP_DBL> S;          ///< Saturation of phase j: numPhase*numBulk.
    vector<OCP_DBL> rho;        ///< Mass density of phase: numPhase*numBulk.
    vector<OCP_DBL> xi;         ///< Moles density of phase: numPhase*numBulk.
    vector<OCP_DBL> xij;        ///< Nij / Nj: numPhase*numCom*numBulk.
    vector<OCP_DBL> Ni;         ///< Moles of component: numCom*numBulk.
    vector<OCP_DBL> mu;         ///< Viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL> kr;         ///< Relative permeability of phase: numPhase*numBulk.
    vector<OCP_DBL> vj;         ///< Volume of phase: numPhase*numBulk.
    vector<OCP_DBL> vf;         ///< Total fluid volume: numBulk.
    vector<OCP_DBL> Nt;         ///< Total moles of components in bulks: numBulk.
    vector<OCP_DBL> lNt;        ///< last Nt
    // Note: Nij is the moles of component i in phase j, Nj is the moles of phase j.

    /////////////////////////////////////////////////////////////////////
    // Derivatives
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL> vfi; ///< dVf / dNi: numCom*numBulk.
    vector<OCP_DBL> vfp; ///< dVf / dP: numBulk.

    /////////////////////////////////////////////////////////////////////
    // Properties at the last time step, determined by methods.
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL> lP;          ///< Pressure: numBulk.
    vector<OCP_DBL> lPj;         ///< Pressure of phase: numPhase*numBulk.
    vector<OCP_DBL> lPc;         ///< Capillary pressure: numPhase*numBulk.
    vector<bool>    lphaseExist; ///< Existence of phases: numPhase*numBulk.
    vector<OCP_DBL> lS;          ///< Saturation of phase: numPhase*numBulk.
    vector<OCP_DBL> lrho;        ///< Mass density of phase: numPhase*numBulk.
    vector<OCP_DBL> lxi;         ///< Moles density of phase: numPhase*numBulk.
    vector<OCP_DBL> lxij;        ///< Nij / Nj: numPhase*numCom*numBulk.
    vector<OCP_DBL> lNi;         ///< Moles of component: numCom*numBulk.
    vector<OCP_DBL> lmu;         ///< Viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL> lkr;         ///< Relative permeability of phase: numPhase*numBulk.
    vector<OCP_DBL> lvj;         ///< Volume of phase: numPhase*numBulk.
    vector<OCP_DBL> lvf;         ///< Total fluid volume: numBulk.
    vector<OCP_DBL> lvfi;        ///< dVf / dNi: numCom*numBulk.
    vector<OCP_DBL> lvfp;        ///< dVf / dP: numBulk.
    vector<OCP_DBL> rockLVp;     ///< Pore volume: numBulk.

    vector<OCP_DBL> surTen;      ///< surface tensions of hydroncarbon phase.
    vector<OCP_DBL> Fk;          ///< The relative permeability interpolation parameter
    vector<OCP_DBL> Fp;          ///< The capillary pressure interpolation parameter
    

    /////////////////////////////////////////////////////////////////////
    // Reservoir rock infomation of each bulk (size = numBulk)
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL> dx;         ///< size of bulk along the x direction.
    vector<OCP_DBL> dy;         ///< size of bulk along the y direction.
    vector<OCP_DBL> dz;         ///< size of bulk along the z direction.
    vector<OCP_DBL> depth;      ///< depth of center of bulk.
    vector<OCP_DBL> ntg;        ///< net to gross of bulk.
    vector<OCP_DBL> rockVpInit; ///< init pore volume = Vgrid * ntg * poro_init.
    vector<OCP_DBL> rockVp;     ///< pore volume = Vgrid * ntg * poro.
    OCP_DBL         rockPref;   ///< reference pressure for initial rock volume.
    OCP_DBL         rockC1;     ///< rock compressibility term 1.
    OCP_DBL         rockC2;     ///< rock compressibility term 2.
    vector<OCP_DBL> rockKxInit; ///< initial rock permeability along the x direction.
    vector<OCP_DBL> rockKx;     ///< current rock permeability along the x direction.
    vector<OCP_DBL> rockKyInit; ///< initial rock permeability along the y direction.
    vector<OCP_DBL> rockKy;     ///< current rock permeability along the y direction.
    vector<OCP_DBL> rockKzInit; ///< initial rock permeability along the z direction.
    vector<OCP_DBL> rockKz;     ///< current rock permeability along the z direction.

    /////////////////////////////////////////////////////////////////////
    // Auxiliary variables
    /////////////////////////////////////////////////////////////////////
    vector<USI> index2Phase; ///< Identify phase name according to index: numPhase.
    // Note: For example, phase 0 is Oil.
    vector<USI> phase2Index; ///< Location of phase according to its name: numPhase.
    // Note: For example, `Oil' is at the i-th location.

    /////////////////////////////////////////////////////////////////////
    // Max range for some variables at the current time step
    /////////////////////////////////////////////////////////////////////
    OCP_DBL dPmax; ///< Max change in pressure.
    OCP_DBL dSmax; ///< Max change in saturation.
    OCP_DBL dNmax; ///< Max change in moles of component.
    OCP_DBL dVmax; ///< Max relative diff between fluid and pore volume.

    /////////////////////////////////////////////////////////////////////
    // For IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for IMPEC.
    void AllocateAuxIMPEC();
    /// Get solution from solver class after linear system is solved.
    void GetSolIMPEC(const vector<OCP_DBL>& u);
    /// Initialize the CFL number.
    void InitCFLIMPEC() const { fill(cfl.begin(), cfl.end(), 0); }
    /// Calculate the CFL number.
    OCP_DBL CalCFL01IMPEC() const;
    /// Update value of last step for IMPEC.
    void UpdateLastStepIMPEC();

private:
    mutable vector<OCP_DBL> cfl; ///< CFL number for each bulk

    /////////////////////////////////////////////////////////////////////
    // For FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for FIM.
    void AllocateAuxFIM();
    /// Get the solution for FIM after a Newton iteration.
    void GetSolFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
                   const OCP_DBL& dSmaxlim);
    /// Get the solution for FIM after a Newton iteration???
    OCP_DBL GetSol01FIM(const vector<OCP_DBL>& u);

    /// Calculate relative resiual for FIM.
    void CalRelResFIM(ResFIM& resFIM) const;
    /// Reset FIM.
    void ResetFIM();
    /// Update values of last step for FIM.
    void UpdateLastStepFIM();
    /// Return NRdPmax.
    OCP_DBL GetNRdPmax() const { return NRdPmax; }
    /// Return NRdSmax.
    OCP_DBL GetNRdSmax() const { return NRdSmax; }

private:
    // Derivatives for FIM
    vector<OCP_DBL> muP;       ///< d Mu   / d P: numPhase*numBulk.
    vector<OCP_DBL> xiP;       ///< d Xi   / d P: numPhase*numBulk.
    vector<OCP_DBL> rhoP;      ///< d Rho  / d P: numPhase*numBulk.
    vector<OCP_DBL> mux;       ///< d Muj  / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> xix;       ///< d Xi_j / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> rhox;      ///< d Rhoj / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> dPcj_dS;   ///< d Pcj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dKr_dS;    ///< d Krj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dSec_dPri; ///< d Secondary variable / d Primary variable.
    // Size: (numPhase + numPhase * numCom) * (numCom + 1) * numBulk

    // vars at last step


    OCP_DBL NRdPmax; ///< Max pressure difference in an NR step
    OCP_DBL NRdSmax; ///< Max saturation difference in an NR step(Predict)

public:
    // for debug!
    void OutputInfo(const OCP_USI& n) const;
    OCP_ULL GetSSMSTAiters()const { return flashCal[0]->GetSSMSTAiters(); }
    OCP_ULL GetNRSTAiters()const { return flashCal[0]->GetNRSTAiters(); }
    OCP_ULL GetSSMSPiters()const { return flashCal[0]->GetSSMSPiters(); }
    OCP_ULL GetNRSPiters()const { return flashCal[0]->GetNRSPiters(); }


    /////////////////////////////////////////////////////////////////////
    // For AIMs, AIMt
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for AIMt.
    void AllocateAuxAIM(const OCP_DBL& ratio);
    OCP_USI GetMaxFIMBulk()const { return maxNumFIMBulk; }
    /// Perform flash calculation with Ni and calculate derivatives.
    void FlashDerivAIM(const bool& IfAIMs);
    void PassFlashValueDerivAIM(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure and their derivatives.
    void CalKrPcDerivAIM(const bool& IfAIMs);
    /// Calculate relative resiual for local FIM.
    void CalRelResAIMt(ResFIM& resFIM) const;
    /// Get the solution for local FIM after a Newton iteration.
    void GetSolAIMt(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
        const OCP_DBL& dSmaxlim);

    /// Calculate relative resiual for AIMs, parts related to FIM are considered.
    void CalRelResAIMs(ResFIM& resFIM) const;
    void GetSolAIMs(const vector<OCP_DBL>& u, const OCP_DBL& dPmaxlim,
        const OCP_DBL& dSmaxlim);


private:
    vector<OCP_USI>   wellBulkId;   ///< Index of bulks which are penetrated by wells ans their K-neighbor
    vector<OCP_INT>   map_Bulk2FIM; ///< Stores the index of FIM bulk in equations
    vector<OCP_USI>   FIMBulk;      ///< index of bulks which performs FIM
    OCP_USI           numFIMBulk;      ///< current num of bulks which are performed by FIM
    OCP_USI           maxNumFIMBulk;   ///< max num of bulks which are performed by FIM
};

#endif /* end if __BULK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/