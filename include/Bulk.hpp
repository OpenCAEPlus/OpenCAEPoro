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
    friend class Reservoir;

    /////////////////////////////////////////////////////////////////////
    // For general usage
    /////////////////////////////////////////////////////////////////////

public:
    Bulk() = default;

    /// Input param from internal data structure ParamReservoir.
    void InputParam(ParamReservoir& rs_param);
    /// Allocate memory for General data.
    void Setup(const Grid& myGrid);
    /// Calculate initial equilibrium for blkoil model according to EQUIL.
    void InitSjPcBo(const USI& tabrow);
    /// Calculate initial equilibration for compositional model according to EQUIL.
    void InitSjPcComp(const USI& tabrow);

    /// Perform flash calculation with saturations.
    void FlashSj();
    /// Perform flash calculation with Ni.
    void FlashNi();
    /// Perform flash calculation with Ni and calculate derivatives.
    void FlashNiDeriv();
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

    /// Reset P to the ones of the last time step.
    void ResetP() { P = lP; }
    /// Reset Pj to the ones of the last time step.
    void ResetPj() { Pj = lPj; }
    /// Reset Ni to the ones of the last time step.
    void ResetNi() { Ni = lNi; }
    /// Reset Vp to the ones of the last time step.
    void ResetVp() { rockVp = rockLVp; }

private:
    /////////////////////////////////////////////////////////////////////
    // General variables
    /////////////////////////////////////////////////////////////////////
    OCP_USI numBulk;  ///< Number of bulks (active grids).
    USI     numPhase; ///< Number of phase.
    USI     numCom;   ///< Number of component.

    // Initial proportion of each component for EoS : numCom - 1, water is excluded.
    vector<OCP_DBL>   initZi;   ///< Initial proportion of each component.
    USI               PVTmode;  ///< Identify PVT mode in blackoil model.
    vector<USI>       PVTNUM;   ///< Identify PVT region in blackoil model: numBulk.
    vector<Mixture*>  flashCal; ///< Flash calculation class.
    USI               SATmode;  ///< Identify SAT mode.
    vector<USI>       SATNUM;   ///< Identify SAT region: numBulk.
    vector<FlowUnit*> flow;     ///< Vector for capillary pressure, relative perm.

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
    void InitCFLIMPEC() const { cfl.assign(numBulk * numPhase, 0); }
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
    // Derivatives For FIM
    vector<OCP_DBL> muP;       ///< dMu / dP: numPhase*numBulk.
    vector<OCP_DBL> xiP;       ///< dXi / dP: numPhase*numBulk.
    vector<OCP_DBL> rhoP;      ///< dRho / dP: numPhase*numBulk.
    vector<OCP_DBL> mux;       ///< dMuj / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> xix;       ///< dXi_j / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> rhox;      ///< dRhoj / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> dPcj_dS;   ///< d Pcj / dSk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dKr_dS;    ///< d Krj / dSk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dSec_dPri; ///< d Secondary var / d Primary var
    // Size: (numPhase + numPhase * numCom) * (numCom + 1) * numBulk

    OCP_DBL NRdPmax; ///< Max pressure difference in NR???
    OCP_DBL NRdSmax; ///< Max saturation difference in NR???
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