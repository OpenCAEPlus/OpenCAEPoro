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
#include "Rock.hpp"
#include "Grid.hpp"
#include "LinearSystem.hpp"
#include "Mixture.hpp"
#include "MixtureBO.hpp"
#include "MixtureComp.hpp"
#include "MixtureThermal.hpp"
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
    OCP_DBL  Dref; ///< The reference depth
    OCP_DBL  Pref; ///< Pressure at the reference depth
    OCP_DBL  DOWC; ///< The depth of oil-water contact surface
    OCP_DBL  DGOC; ///< The depth of gas-oil contact surface
    OCP_DBL  PcOW; ///< Capillary pressure at oil-water contact Pcow = Po - Pw
    OCP_DBL  PcGO; ///< capillary pressure at gas-oil contact Pcgo = Pg - Po
    OCPTable PBVD; ///< PBVD Table: bubble point pressure vs depth
};


class SkipStaAnaly
{
    friend class Bulk;
public:
    SkipStaAnaly() = default;
    void SetUseSkip(const OCP_BOOL& flag) { useSkip = flag; }
    OCP_BOOL IfUseSkip() const { return useSkip; }
    void Setup(const OCP_USI& numBulk, const USI& numCom_1) {
        flagSkip.resize(numBulk);
        PSkip.resize(numBulk);
        TSkip.resize(numBulk);
        minEigenSkip.resize(numBulk);       
        ziSkip.resize(numBulk * numCom_1);
    }
    OCP_BOOL IfSkip(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL& Ntin, 
        const OCP_DBL* Niin, const OCP_USI& n, const USI& numCom_1) const {
        if (flagSkip[n]) {
            if (fabs(1 - PSkip[n] / Pin) >= minEigenSkip[n] / 10) {
                return OCP_FALSE;
            }
            OCP_DBL Nt_w = Ntin - Niin[numCom_1];
            for (USI i = 0; i < numCom_1; i++) {
                if (fabs(Niin[i] / Nt_w - ziSkip[n*numCom_1 + i]) >= minEigenSkip[n] / 10) {
                    return OCP_FALSE;
                }
            }
            if (fabs(TSkip[n] - Tin) >= minEigenSkip[n] * 10) {
                return OCP_FALSE;
            }
            return OCP_TRUE;
        }
        else {
            return OCP_FALSE;
        }
    }

protected:
    OCP_BOOL        useSkip{ OCP_FALSE };
    vector<OCP_DBL> flagSkip;
    vector<OCP_DBL> minEigenSkip;
    vector<OCP_DBL> PSkip;
    vector<OCP_DBL> TSkip;
    vector<OCP_DBL> ziSkip;
};


/// Physical information of each active reservoir bulk.
//  Note: Bulk contains main physical infomation of active grids. It describes the
//  actual geometric domain for simulating. Variables are stored bulk by bulk, and then
//  phase by phase, then component by component. The bulks are ordered in the alphabetic
//  order, i.e. the X-axis indices first, followed by the Y-axis and Z-axis indices.
//  Operations on each bulk are also defined here.
class Bulk
{
    friend class BulkConn;
    friend class Well;
    friend class Out4RPT;
    friend class Out4VTK;

    // temp
    friend class Reservoir;
    friend class OCP_IMPEC;
    friend class OCP_FIM;
    friend class OCP_FIMn;
    friend class OCP_AIMc;

    /////////////////////////////////////////////////////////////////////
    // For general usage
    /////////////////////////////////////////////////////////////////////

public:
    Bulk() = default;

    /// Input param from internal data structure ParamReservoir.
    void InputParam(ParamReservoir& rs_param);
    void InputParamBLKOIL(ParamReservoir& rs_param);
    void InputParamCOMPS(const ParamReservoir& rs_param);
    void InputParamTHERMAL(const ParamReservoir& rs_param);
    void InputSatFunc(const ParamReservoir& rs_param);
    void InputRockFunc(const ParamReservoir& rs_param);
    void InputRockFuncT(const ParamReservoir& rs_param);
    /// Allocate memory for bulk data of grid.
    void Setup(const Grid& myGrid);
    /// Calculate initial equilibrium
    void InitSjPc(const Grid& myGrid, const USI& tabrow);
    /// Perform flash calculation with saturations.
    void InitFlash(const OCP_BOOL& flag = OCP_FALSE);
    /// Perform flash calculation with saturations and calculate derivatives.
    void InitFlashDer();
    void InitFlashDer_n();
    /// Perform flash calculation with Ni.
    void Flash();
    /// Perform flash calculation with Ni and calculate derivatives.
    void FlashDeriv();
    void FlashDeriv_n();
    /// Perform flash calculation with Ni in Black Oil Model
    void FlashDerivBLKOIL_n();
    /// Perform flash calculation with Ni in Compositional Model
    void FlashDerivCOMP_n();
    /// determine which flash type will be used
    USI CalFlashType(const OCP_USI& n, const OCP_BOOL& fimbulk) const;
    /// Pass values from Flash to Bulk after Flash calculation.
    void PassFlashValue(const OCP_USI& n);
    void PassFlashValueAIMc(const OCP_USI& n);
    /// Pass derivative values from Flash to Bulk after Flash calculation.
    void PassFlashValueDeriv(const OCP_USI& n);
    void PassFlashValueDeriv_n(const OCP_USI& n);
    void PassAdditionInfo(const OCP_USI& n, const USI& pvtnum);
    /// Reset variables in flash calculations.
    void ResetFlash();

    /// Calculate relative permeability and capillary pressure with saturation.
    void CalKrPc();
    /// Calculate relative permeability and capillary pressure and their derivatives.
    void CalKrPcDeriv();
    /// Calculate volume of pore with pressure.
    void CalRock();
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
    /// Check if negative P occurs, return OCP_FALSE if so.
    OCP_BOOL CheckP() const;
    /// Check if negative Ni occurs, return OCP_FALSE if so.
    OCP_BOOL CheckNi();
    /// Check if relative volume error is out of range, return OCP_FALSE if so.
    OCP_BOOL CheckVe(const OCP_DBL& Vlim) const;
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
    void ResetSkipStaAnalyTerm() { skipStaAnalyTerm = lskipStaAnalyTerm; }

    /// Reset Nt to the ones of the last time step.
    void ResetNt() { Nt = lNt; }

    /// Reset P to the ones of the last time step.
    void ResetP() { P = lP; }
    /// Reset Pj to the ones of the last time step.
    void ResetPj() { Pj = lPj; }
    /// Reset Ni to the ones of the last time step.
    void ResetNi() { Ni = lNi; }
    /// Reset Kr to the ones of the last time step.
    void ResetKr() { kr = lkr; }
    /// Reset Vp to the ones of the last time step.
    void ResetRock() { rockVp = lrockVp; poro = lporo; poroP = lporoP; }
    void CalSomeInfo(const Grid& myGrid) const;

    /// Allocate memory for wellBulkId
    void    AllocateWellBulkId(const USI& n) { wellBulkId.reserve(n); }
    void    ClearWellBulkId() { wellBulkId.clear(); }
    OCP_DBL CalNT()
    {
        NT = Dnorm1(numBulk, &Nt[0]);
        return NT;
    }

private:
    /////////////////////////////////////////////////////////////////////
    // General variables
    /////////////////////////////////////////////////////////////////////
    OCP_USI numBulk;  ///< Number of bulks (active grids).
    USI     numPhase; ///< Number of phase.
    USI     numCom;   ///< Number of component.
    USI     numCom_1; ///< numCom - 1

    // Initial proportion of each component for EoS : numCom - 1, water is excluded.
    vector<OCPTable> initZi_Tab; ///< InitZi tab set
    vector<OCPTable> initT_Tab;  ///< InitT  tab set

    vector<OCP_DBL> SwatInit;                 ///< Initial water saturation.
    OCP_BOOL        SwatInitExist{OCP_FALSE}; ///< If SwatInit has been given.
    vector<OCP_DBL> ScaleValuePcow;           ///< Scale values for Pcow.
    OCP_BOOL        ScalePcow{OCP_FALSE};     ///< whether Pcow should be scaled.

    USI               PVTmode;  ///< Identify PVT mode in black-oil model.
    vector<USI>       PVTNUM;   ///< Identify PVT region in black-oil model: numBulk.
    USI               NTPVT;    ///< num of PVT regions
    vector<Mixture*>  flashCal; ///< Flash calculation class.
    USI               SATmode;  ///< Identify SAT mode.
    vector<USI>       SATNUM;   ///< Identify SAT region: numBulk.
    USI               NTSFUN;   ///< num of SAT regions
    vector<FlowUnit*> flow;     ///< Vector for capillary pressure, relative perm.
    vector<vector<OCP_DBL>>
        satcm; ///< critical saturation when phase becomes mobile / immobile.
    USI               NTROCC;   ///< num of Rock regions
    vector<USI>       ROCKNUM;  ///< index of Rock table for each bulk
    vector<Rock*>     rock;     ///< rock model

    // Skip stability analysis
    vector<USI> phaseNum;   ///< Num of hydrocarbon phase in each bulk
    vector<USI> lphaseNum;  ///< last phaseNum
    vector<USI> NRphaseNum; ///< phaseNum in NR step

    SkipStaAnaly skipStaAnalyTerm;
    SkipStaAnaly lskipStaAnalyTerm;

    /////////////////////////////////////////////////////////////////////
    // Basic model information
    /////////////////////////////////////////////////////////////////////
    ParamEQUIL EQUIL;    ///< Initial Equilibration.
    OCP_BOOL   blackOil{ OCP_FALSE }; ///< If OCP_TRUE, black-oil model will be used.
    OCP_BOOL   comps{ OCP_FALSE };    ///< If OCP_TRUE, compositional model will be used.
    OCP_BOOL   thermal{ OCP_FALSE };
    OCP_BOOL   useEoS{ OCP_FALSE };   ///< If OCP_TRUE, then EoS model is used.
    OCP_BOOL   oil{ OCP_FALSE };      ///< If OCP_TRUE, oil phase could exist.
    OCP_BOOL   gas{ OCP_FALSE };      ///< If OCP_TRUE, gas phase could exist.
    OCP_BOOL   water{ OCP_FALSE };    ///< If OCP_TRUE, water phase could exist.
    OCP_BOOL   disGas{ OCP_FALSE };   ///< If OCP_TRUE, dissolve gas in live oil could exist.
    OCP_BOOL
    miscible{ OCP_FALSE }; ///< Miscible treatment of hydrocarbons, used in compositional Model.

    /////////////////////////////////////////////////////////////////////
    // Basic physical variables
    /////////////////////////////////////////////////////////////////////
    OCP_DBL          RTemp;      ///< Reservoir temperature.
    vector<OCP_DBL>  Pb;         ///< Bubble point pressure: numBulk.
    vector<OCP_DBL>  T;          ///< Temperature: numBulk. unit£ºF
    vector<OCP_DBL>  P;          ///< Pressure: numBulk.
    vector<OCP_DBL>  Pj;         ///< Pressure of phase: numPhase*numBulk.
    vector<OCP_DBL>  Pc;         ///< Capillary pressure of phase: numPhase*numBulk.
    vector<OCP_BOOL> phaseExist; ///< Existence of phase: numPhase*numBulk.
    vector<OCP_DBL>  S;          ///< Saturation of phase j: numPhase*numBulk.
    vector<OCP_DBL>  nj;         ///< moles number of phase j: numPhase*numBulk.
    vector<OCP_DBL>  rho;        ///< Mass density of phase: numPhase*numBulk.
    vector<OCP_DBL>  xi;         ///< Moles density of phase: numPhase*numBulk.
    vector<OCP_DBL>  xij;        ///< Nij / Nj: numPhase*numCom*numBulk.
    vector<OCP_DBL>  Ni;         ///< Moles of component: numCom*numBulk.
    vector<OCP_DBL>  mu;         ///< Viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL>  kr;         ///< Relative permeability of phase: numPhase*numBulk.
    vector<OCP_DBL>  vj;         ///< Volume of phase: numPhase*numBulk.
    vector<OCP_DBL>  vf;         ///< Total fluid volume: numBulk.
    vector<OCP_DBL>  Nt;         ///< Total moles of components in bulks: numBulk.
    OCP_DBL          NT;         ///< sum of Nt in all bulks
    // Note: Nij is the moles of component i in phase j, Nj is the moles of phase j.
    // for thermal model
    vector<OCP_DBL>  thconp;     ///< phase thermal conductivity: numPhase
    vector<OCP_DBL> vfT;  ///< d vf  / dT, numBulk
    vector<OCP_DBL> muT;  ///< d mu j  / dT, numPhase * numbulk
    vector<OCP_DBL> xiT;  ///< d xi j / dT, numPhase * numbulk
    vector<OCP_DBL> rhoT; ///< d rho j / dT, numPhase * numbulk
    vector<OCP_DBL> Uf;   ///< Internal energy of fluid, numBulk
    vector<OCP_DBL> Ufi;  ///< dUf / dNi, numCom * numBulk
    vector<OCP_DBL> Ufp;  ///< dUf / dP, numbulk
    vector<OCP_DBL> UfT;  ///< dUf / dT, numbulk
    vector<OCP_DBL> H;    ///< Enthalpy, numPhase * numbulk
    vector<OCP_DBL> HT;   ///< d Hj / d T, numPhase * numbulk
    vector<OCP_DBL> Hx;   ///< d Hj / d xij, numPhase * numCom * numbulk

    /////////////////////////////////////////////////////////////////////
    // Derivatives
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL> vfi; ///< dVf / dNi: numCom*numBulk.
    vector<OCP_DBL> vfp; ///< dVf / dP: numBulk.

    /////////////////////////////////////////////////////////////////////
    // Properties at the last time step, determined by methods.
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL>  lP;          ///< Pressure: numBulk.
    vector<OCP_DBL>  lPj;         ///< Pressure of phase: numPhase*numBulk.
    vector<OCP_DBL>  lPc;         ///< Capillary pressure: numPhase*numBulk.
    vector<OCP_BOOL> lphaseExist; ///< Existence of phases: numPhase*numBulk.
    vector<OCP_DBL>  lS;          ///< Saturation of phase: numPhase*numBulk.
    vector<OCP_DBL>  lnj;         ///< last nj: numPhase*numBulk.
    vector<OCP_DBL>  lrho;        ///< Mass density of phase: numPhase*numBulk.
    vector<OCP_DBL>  lxi;         ///< Moles density of phase: numPhase*numBulk.
    vector<OCP_DBL>  lxij;        ///< Nij / Nj: numPhase*numCom*numBulk.
    vector<OCP_DBL>  lNi;         ///< Moles of component: numCom*numBulk.
    vector<OCP_DBL>  lmu;         ///< Viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL>  lkr;         ///< Relative permeability of phase: numPhase*numBulk.
    vector<OCP_DBL>  lvj;         ///< Volume of phase: numPhase*numBulk.
    vector<OCP_DBL>  lvf;         ///< Total fluid volume: numBulk.
    vector<OCP_DBL>  lNt;         ///< last Nt
    vector<OCP_DBL>  lvfi;        ///< dVf / dNi: numCom*numBulk.
    vector<OCP_DBL>  lvfp;        ///< dVf / dP: numBulk.
    vector<OCP_DBL>  lrockVp;     ///< Pore volume: numBulk.

    vector<OCP_DBL> surTen; ///< surface tensions of hydrocarbon phase.
    vector<OCP_DBL> Fk;     ///< The relative permeability interpolation parameter
    vector<OCP_DBL> Fp;     ///< The capillary pressure interpolation parameter

    vector<OCP_DBL> lsurTen; ///< last surTen.

    /////////////////////////////////////////////////////////////////////
    // Reservoir rock infomation of each bulk (size = numBulk)
    /////////////////////////////////////////////////////////////////////
    vector<OCP_DBL> ntg;        ///< net to gross of bulk.
    vector<OCP_DBL> poroInit;   ///< initial rock porosity.
    vector<OCP_DBL> poro;       ///< rock porosity.
    vector<OCP_DBL> poroP;      ///< d poro / d P.
    vector<OCP_DBL> rockVntg; ///< init pore volume = Vgrid * ntg.
    vector<OCP_DBL> rockVp;     ///< pore volume = Vgrid * ntg * poro.
    vector<OCP_DBL> rockKx;     ///< current rock permeability along the x direction.
    vector<OCP_DBL> rockKy;     ///< current rock permeability along the y direction.
    vector<OCP_DBL> rockKz;     ///< current rock permeability along the z direction.

    // last step
    vector<OCP_DBL> lporo;      ///< last poro.
    vector<OCP_DBL> lporoP;     ///< last poroP.

    /////////////////////////////////////////////////////////////////////
    // Auxiliary variables
    /////////////////////////////////////////////////////////////////////
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
    // Error
    /////////////////////////////////////////////////////////////////////

    vector<OCP_DBL>         ePEC; ///< error for fugacity balance equations
    mutable vector<OCP_DBL> eN;   ///< error for mass conservation equations
    mutable vector<OCP_DBL> eV;   ///< error for volume conservation equations

    /////////////////////////////////////////////////////////////////////
    // For IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for IMPEC.
    void AllocateAuxIMPEC();
    /// Update P and Pj after linear system is solved.
    void GetSolIMPEC(const vector<OCP_DBL>& u);
    /// Initialize the CFL number.
    void SetCFL2Zero() const { fill(cfl.begin(), cfl.end(), 0); }
    /// Calculate the CFL number.
    OCP_DBL CalCFL() const;
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
    void AllocateAuxFIMn();
    /// Get the solution for FIM after a Newton iteration.
    void GetSolFIM(const vector<OCP_DBL>& u,
                   const OCP_DBL&         dPmaxlim,
                   const OCP_DBL&         dSmaxlim);
    void GetSolFIM_n(const vector<OCP_DBL>& u,
                     const OCP_DBL&         dPmaxlim,
                     const OCP_DBL&         dSmaxlim);
    /// Calculate relative residual for FIM.
    void CalRelResFIM(OCPRes& resFIM) const;
    // Show Res
    void ShowRes(const vector<OCP_DBL>& res) const;
    /// Reset FIM.
    void ResetFIM();
    void ResetFIMn();
    /// Update values of last step for FIM.
    void UpdateLastStepFIM();
    void UpdateLastStepFIMn();
    /// Calculate some auxiliary variable, for example, dSmax
    OCP_DBL CalNRdSmax(OCP_USI& index);

    /// Return NRdPmax.
    OCP_DBL GetNRdPmax();
    /// Return NRdSmaxP.
    OCP_DBL GetNRdSmaxP();
    OCP_DBL GetNRdNmax();
    void    CorrectNi(const vector<OCP_DBL>& res);

private:
    // Derivatives for FIM
    vector<OCP_DBL> muP;          ///< d Mu   / d P: numPhase*numBulk.
    vector<OCP_DBL> xiP;          ///< d Xi   / d P: numPhase*numBulk.
    vector<OCP_DBL> rhoP;         ///< d Rho  / d P: numPhase*numBulk.
    vector<OCP_DBL> mux;          ///< d Muj  / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> xix;          ///< d Xi_j / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> rhox;         ///< d Rhoj / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL> dPcj_dS;      ///< d Pcj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dKr_dS;       ///< d Krj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dSec_dPri;    ///< d Secondary variable / d Primary variable.
    vector<OCP_DBL> res_n;        ///< ...
    vector<OCP_DBL> resPc;        ///< a precalculated value
    USI             maxLendSdP;   ///< length of dSec_dPri.
    vector<USI>     bRowSizedSdP; ///< length of dSec_dPri in each bulk
    vector<OCP_USI> resIndex; ///< store the starting position of res_n of each bulk.

    // Auxiliary variable for dSec_dPr
    vector<OCP_BOOL> pSderExist; ///< Existence of  derivative of phase saturation
    vector<USI>      pVnumCom;   ///< num of variable components in the phase

    // vars at last step
    vector<OCP_DBL>  lmuP;          ///< last muP
    vector<OCP_DBL>  lxiP;          ///< last xiP
    vector<OCP_DBL>  lrhoP;         ///< last rhoP
    vector<OCP_DBL>  lmux;          ///< last mux
    vector<OCP_DBL>  lxix;          ///< last xix
    vector<OCP_DBL>  lrhox;         ///< last rhox
    vector<OCP_DBL>  ldPcj_dS;      ///< last Pcj_dS
    vector<OCP_DBL>  ldKr_dS;       ///< last dKr_dS
    vector<OCP_DBL>  ldSec_dPri;    ///< last dSec_dPri
    vector<OCP_DBL>  lres_n;        ///< last res_n
    vector<OCP_DBL>  lresPc;        ///< last lresPc;
    vector<USI>      lbRowSizedSdP; ///< last bRowSizedSdP
    vector<OCP_USI>  lresIndex;     ///< last res_n
    vector<OCP_BOOL> lpSderExist;   ///< last pSderExist
    vector<USI>      lpVnumCom;     ///< last pVnumCom

    vector<OCP_DBL> dSNR;  ///< saturation change between NR steps
    vector<OCP_DBL> dSNRP; ///< predicted saturation change between NR steps
    vector<OCP_DBL> dNNR;  ///< Ni change between NR steps
    vector<OCP_DBL> dPNR;  ///< dP change between NR steps

    OCP_DBL NRdSSP;    ///< difference between dSNR and dSNRP, 2-norm
    OCP_DBL maxNRdSSP; ///< max difference between dSNR and dSNRP
    OCP_USI index_maxNRdSSP;
    OCP_DBL NRdPmax;  ///< Max pressure difference in an NR step
    OCP_DBL NRdNmax;  ///< Max Ni difference in an NR step
    OCP_DBL NRdSmax;  ///< Max saturation difference in an NR step(Real)
    OCP_DBL NRdSmaxP; ///< Max saturation difference in an NR step(Predict)

    vector<OCP_DBL> NRstep; ///< NRstep for FIM

public:
    // for debug!
    void    OutputInfo(const OCP_USI& n) const;
    OCP_ULL GetSSMSTAiters() const { return flashCal[0]->GetSSMSTAiters(); }
    OCP_ULL GetNRSTAiters() const { return flashCal[0]->GetNRSTAiters(); }
    OCP_ULL GetSSMSPiters() const { return flashCal[0]->GetSSMSPiters(); }
    OCP_ULL GetNRSPiters() const { return flashCal[0]->GetNRSPiters(); }
    OCP_ULL GetRRiters() const { return flashCal[0]->GetRRiters(); }
    OCP_ULL GetSSMSTAcounts() const { return flashCal[0]->GetSSMSTAcounts(); }
    OCP_ULL GetNRSTAcounts() const { return flashCal[0]->GetNRSTAcounts(); }
    OCP_ULL GetSSMSPcounts() const { return flashCal[0]->GetSSMSPcounts(); }
    OCP_ULL GetNRSPcounts() const { return flashCal[0]->GetNRSPcounts(); }
    OCP_ULL GetRRcounts() const { return flashCal[0]->GetRRcounts(); }

    /////////////////////////////////////////////////////////////////////
    // For AIMc
    /////////////////////////////////////////////////////////////////////

public:
    // Print FIM Bulk
    void ShowFIMBulk(const OCP_BOOL& flag = OCP_FALSE) const;
    // Allocate auxiliary variable
    void AllocateAuxAIMc();
    /// Perform flash calculation with Ni for IMPEC bulk -- Update partial properties
    void FlashAIMc();
    /// Perform flash calculation with Ni for IMPEC bulk -- Update all properties
    void FlashAIMc01();
    /// Perform flash calculation with Ni and calculate derivatives.
    void FlashDerivAIMc();
    /// Calculate relative permeability and capillary pressure with saturation.
    void CalKrPcAIMc();
    /// Calculate relative permeability and capillary pressure and their derivatives.
    void CalKrPcDerivAIMc();
    void GetSolAIMc(const vector<OCP_DBL>& u,
                    const OCP_DBL&         dPmaxlim,
                    const OCP_DBL&         dSmaxlim);
    void UpdatePj();
    void ResetXijNR() { xijNR = lxij; }

private:
    vector<OCP_USI> wellBulkId;    ///< Index of bulks which are penetrated by wells and their K-neighbor
    class BulkTypeAIM
    {
    public:
        BulkTypeAIM() = default;
        void Init(const OCP_USI& nb) { indicator.resize(nb, -1); }
        void Set0() { fill(indicator.begin(), indicator.end(), -1); }
        void SetBulkType(const OCP_USI& n, const OCP_INT& flag) { indicator[n] = flag; }
        OCP_BOOL IfFIMbulk(const OCP_USI& n)const { return indicator[n] > 0; }
        OCP_BOOL IfIMPECbulk(const OCP_USI& n)const { return indicator[n] < 0; }
    private:
        vector<OCP_INT> indicator;  ///< Stores the index of FIM bulk in equations, FIM bulk: >=0; IMPEC bulk: <0;
    } bulkTypeAIM;
    vector<OCP_DBL> xijNR; ///< store the current NR step's xij in AIM
    
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