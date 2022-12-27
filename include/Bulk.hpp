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
#include "Grid.hpp"
#include "FlowUnit.hpp"
#include "Rock.hpp"
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
    // Input Param and Setup
    /////////////////////////////////////////////////////////////////////

public:
    /// Input param from internal data structure ParamReservoir.
    void InputParam(const ParamReservoir& rs_param);
    void InputParamBLKOIL(const ParamReservoir& rs_param);
    void InputParamCOMPS(const ParamReservoir& rs_param);
    void InputParamTHERMAL(const ParamReservoir& rs_param);
    void InputSatFunc(const ParamReservoir& rs_param);
    void InputRockFunc(const ParamReservoir& rs_param);
    void InputRockFuncT(const ParamReservoir& rs_param);

    /// Allocate memory for fluid grid for isothermal model
    void SetupIsoT(const Grid& myGrid);
    /// Allocate memory for fluid grid for ifThermal model
    void SetupT(const Grid& myGrid);
    /// Setup Optional Feature
    void SetupOptionalFeatures(const Grid& myGrid, OptionalFeatures& optFeatures);

    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Return the number of bulks.
    OCP_USI GetBulkNum() const { return numBulk; }
    /// Return the number of phases.
    USI GetPhaseNum() const { return numPhase; }
    /// Return the number of components.
    USI GetComNum() const { return numCom; }

protected:

    OCP_USI numBulk;  ///< Number of bulks (active grids).
    USI     numPhase; ///< Number of phase.
    USI     numCom;   ///< Number of component.
    USI     numComH;  ///< Number of HydronCarbon


    /////////////////////////////////////////////////////////////////////
    // Initial Properties
    /////////////////////////////////////////////////////////////////////
public:
    /// Calculate initial equilibrium
    void InitSjPc(const USI& tabrow);

protected:

    vector<OCPTable> initZi_Tab; ///< initial mole ratio of components vs. depth, table set
    vector<OCPTable> initT_Tab;  ///< initial temperature vs. depth, table set
    ParamEQUIL       EQUIL;      ///< Initial Equilibration.
    OCP_DBL          RTemp;      ///< Reservoir temperature.
    vector<OCP_DBL>  thconp;     ///< phase ifThermal conductivity: numPhase


    /////////////////////////////////////////////////////////////////////
    // Region
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for region num
    void AllocateRegion(const Grid& myGrid);
    /// Return flash.
    const vector<Mixture*>& GetMixture() const { return flashCal; }
    /// Output iterations in Mixture
    void OutMixtureIters()const { flashCal[0]->OutMixtureIters(); }

protected:

    USI               NTPVT;    ///< num of PVT regions
    USI               PVTmodeB; ///< Identify PVT mode in black-oil model.
    vector<USI>       PVTNUM;   ///< Identify PVT region in black-oil model: numBulk.
    vector<Mixture*>  flashCal; ///< Flash calculation class.

    USI               NTSFUN;   ///< num of SAT regions
    USI               SATmode;  ///< Identify SAT mode.
    vector<USI>       SATNUM;   ///< Identify SAT region: numBulk.
    vector<FlowUnit*> flow;     ///< Vector for capillary pressure, relative perm.
    vector<vector<OCP_DBL>> satcm; ///< critical saturation when phase becomes mobile / immobile.

    USI               NTROCC;   ///< num of Rock regions
    vector<USI>       ROCKNUM;  ///< index of Rock table for each bulk
    vector<Rock*>     rock;     ///< rock model
 

    /////////////////////////////////////////////////////////////////////
    // Basic PVT Model Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Return ifUseEoS.
    OCP_BOOL IfUseEoS() const { return ifUseEoS; };

protected:

    OCP_BOOL   ifBlackOil{ OCP_FALSE }; ///< If OCP_TRUE, black-oil model will be used.
    OCP_BOOL   ifComps{ OCP_FALSE };    ///< If OCP_TRUE, compositional model will be used.
    OCP_BOOL   ifThermal{ OCP_FALSE };  ///< Id OCP_TRUE, ifThermal model will be used.
    OCP_BOOL   ifUseEoS{ OCP_FALSE };   ///< If OCP_TRUE, then EoS model is used.
    OCP_BOOL   oil{ OCP_FALSE };      ///< If OCP_TRUE, oil phase could exist.
    OCP_BOOL   gas{ OCP_FALSE };      ///< If OCP_TRUE, gas phase could exist.
    OCP_BOOL   water{ OCP_FALSE };    ///< If OCP_TRUE, water phase could exist.
    OCP_BOOL   disGas{ OCP_FALSE };   ///< If OCP_TRUE, dissolve gas in live oil could exist.

    /////////////////////////////////////////////////////////////////////
    // Basic Grid and Basic Rock Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for Rock properties
    void AllocateGridRockIsoT(const Grid& myGrid);
    /// Initialize the information of rock
    void InitRock();
    /// Calculate volume of pore with pressure.
    void CalRock();

protected:

    vector<OCP_DBL> dx;     ///< Size of cell in x-direction: activeGridNum.
    vector<OCP_DBL> dy;     ///< Size of cell in y-direction: activeGridNum.
    vector<OCP_DBL> dz;     ///< Size of cell in z-direction: activeGridNum.
    vector<OCP_DBL> v;      ///< Volume of grids: activeGridNum.
    vector<OCP_DBL> depth;  ///< Depth of center of grid cells: activeGridNum.

    vector<OCP_DBL> ntg;        ///< net to gross of bulk.
    vector<OCP_DBL> poroInit;   ///< initial rock porosity.
    vector<OCP_DBL> rockVntg;   ///< init effective volume = Vgrid * ntg.
    vector<OCP_DBL> poro;       ///< rock porosity.
    vector<OCP_DBL> rockVp;     ///< pore volume = Vgrid * ntg * poro.
    vector<OCP_DBL> rockKx;     ///< current rock permeability along the x direction.
    vector<OCP_DBL> rockKy;     ///< current rock permeability along the y direction.
    vector<OCP_DBL> rockKz;     ///< current rock permeability along the z direction.
    vector<OCP_DBL> thconr;     ///< Rock ifThermal conductivity: activeGridNum.
    vector<OCP_DBL> vr;         ///< Volume of rock: activeGridNum.
    vector<OCP_DBL> Hr;         ///< Enthalpy of rock: activeGridNum.  

    // Last time step
    vector<OCP_DBL> lporo;      ///< last poro.
    vector<OCP_DBL> lrockVp;     ///< Pore volume: numBulk.
    vector<OCP_DBL> lvr;         ///< Last vr: activeGridNum.
    vector<OCP_DBL> lHr;         ///< Last Hr: activeGridNum.

    // Derivatives
    vector<OCP_DBL> poroP;      ///< d poro / d P.
    vector<OCP_DBL> poroT;      ///< d poro / d T.
    vector<OCP_DBL> vrP;        ///< d vr / d p, numbulk
    vector<OCP_DBL> vrT;        ///< dvr / dT: activeGridNum.
    vector<OCP_DBL> HrT;        ///< dHr / dT: activeGridNum.
    
    // Last time step   
    vector<OCP_DBL> lporoP;     ///< last poroP.    
    vector<OCP_DBL> lporoT;     ///< last poroT.
    vector<OCP_DBL> lvrP;       ///< last vrp.
    vector<OCP_DBL> lvrT;       ///< Last vrT.
    vector<OCP_DBL> lHrT;       ///< Last HrT.

 
    /////////////////////////////////////////////////////////////////////
    // Basic Fluid Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Calculate average pressure in reservoir.
    OCP_DBL CalFPR() const;
    /// Return pressure of the n-th bulk.
    OCP_DBL GetP(const OCP_USI& n) const { return P[n]; }
    /// Return oil saturation of the n-th bulk.
    OCP_DBL GetSOIL(const OCP_USI& n) const { return S[n * numPhase + phase2Index[OIL]];}
    /// Return gas saturation of the n-th bulk.
    OCP_DBL GetSGAS(const OCP_USI& n) const { return S[n * numPhase + phase2Index[GAS]];}
    /// Return water saturation of the n-th bulk.
    OCP_DBL GetSWAT(const OCP_USI& n) const { return S[n * numPhase + phase2Index[WATER]]; }

protected:
    vector<USI> phase2Index; ///< Location of phase according to its name: numPhase.
                             // Note: For example, `Oil' is at the i-th location.
    vector<USI>      phaseNum;   ///< Num of hydrocarbon phase in each bulk
    vector<OCP_DBL>  Nt;         ///< Total moles of components in bulks: numBulk.
    vector<OCP_DBL>  Ni;         ///< Moles of component: numCom*numBulk. 
    vector<OCP_DBL>  vf;         ///< Total fluid volume: numBulk.   
    vector<OCP_DBL>  T;          ///< Temperature: numBulk.
    vector<OCP_DBL>  P;          ///< Pressure: numBulk.
    vector<OCP_DBL>  Pb;         ///< Bubble point pressure: numBulk.
    vector<OCP_DBL>  Pj;         ///< Pressure of phase: numPhase*numBulk.
    vector<OCP_DBL>  Pc;         ///< Capillary pressure of phase: numPhase*numBulk.
    vector<OCP_BOOL> phaseExist; ///< Existence of phase: numPhase*numBulk.
    vector<OCP_DBL>  S;          ///< Saturation of phase: numPhase*numBulk.
    vector<OCP_DBL>  vj;         ///< Volume of phase: numPhase*numBulk.
    vector<OCP_DBL>  nj;         ///< moles number of phase: numPhase*numBulk.
    vector<OCP_DBL>  xij;        ///< Nij / Nj: numPhase*numCom*numBulk.
    vector<OCP_DBL>  rho;        ///< Mass density of phase: numPhase*numBulk.
    vector<OCP_DBL>  xi;         ///< Moles density of phase: numPhase*numBulk.
    vector<OCP_DBL>  mu;         ///< Viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL>  kr;         ///< Relative permeability of phase: numPhase*numBulk.
    vector<OCP_DBL>  Uf;         ///< Internal energy of fluid: numBulk
    vector<OCP_DBL>  H;          ///< Enthalpy of phase: numPhase*numBulk.
    vector<OCP_DBL>  kt;         ///< Coefficient of ifThermal diffusivity: activeGridNum.

    // Last time step
    vector<USI>      lphaseNum;   ///< last phaseNum
    vector<OCP_DBL>  lNt;         ///< last Nt
    vector<OCP_DBL>  lNi;         ///< last Ni
    vector<OCP_DBL>  lvf;         ///< last vf
    vector<OCP_DBL>  lT;          ///< last T
    vector<OCP_DBL>  lP;          ///< last P
    vector<OCP_DBL>  lPj;         ///< last Pj
    vector<OCP_DBL>  lPc;         ///< last Pc
    vector<OCP_BOOL> lphaseExist; ///< last phaseExist
    vector<OCP_DBL>  lS;          ///< last S
    vector<OCP_DBL>  lvj;         ///< last vj
    vector<OCP_DBL>  lnj;         ///< last nj
    vector<OCP_DBL>  lxij;        ///< last xij
    vector<OCP_DBL>  lrho;        ///< last rho
    vector<OCP_DBL>  lxi;         ///< last xi
    vector<OCP_DBL>  lmu;         ///< last mu
    vector<OCP_DBL>  lkr;         ///< last kr
    vector<OCP_DBL>  lUf;         ///< last Uf
    vector<OCP_DBL>  lH;          ///< last H
    vector<OCP_DBL>  lkt;         ///< last kt 

    // Derivatives  
    vector<OCP_DBL>  vfP;         ///< d vf   / d P: numBulk.
    vector<OCP_DBL>  vfT;         ///< d vf   / d T, numBulk
    vector<OCP_DBL>  vfi;         ///< d vf   / d Ni: numCom*numBulk.
    vector<OCP_DBL>  rhoP;        ///< d Rho  / d P: numPhase*numBulk.
    vector<OCP_DBL>  rhoT;        ///< d rhoj / d T: numPhase * numbulk 
    vector<OCP_DBL>  rhox;        ///< d Rhoj / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL>  xiP;         ///< d xi   / d P: numPhase*numBulk.
    vector<OCP_DBL>  xiT;         ///< d xij  / d T, numPhase * numbulk
    vector<OCP_DBL>  xix;         ///< d Xi_j / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL>  muP;         ///< d Mu   / d P: numPhase*numBulk.
    vector<OCP_DBL>  muT;         ///< d muj  / d T: numPhase * numbulk
    vector<OCP_DBL>  mux;         ///< d Muj  / d xij: numPhase*numCom*numBulk.
    vector<OCP_DBL>  dPcj_dS;     ///< d Pcj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL>  dKr_dS;      ///< d Krj  / d Sk: numPhase * numPhase * bulk.
    vector<OCP_DBL>  UfP;         ///< d Uf   / d P: numbulk
    vector<OCP_DBL>  UfT;         ///< d Uf   / d T: numbulk
    vector<OCP_DBL>  Ufi;         ///< d Uf   / d Ni: numCom * numBulk  
    vector<OCP_DBL>  HT;          ///< d Hj   / d T: numPhase * numbulk
    vector<OCP_DBL>  Hx;          ///< d Hj   / d xij: numPhase * numCom * numbulk
    vector<OCP_DBL>  ktp;         ///< d kt   / d P: numbulk
    vector<OCP_DBL>  ktT;         ///< d kt   / d T: activeGridNum.
    vector<OCP_DBL>  ktS;         ///< d kt   / d S: numPhase * numbulk 

    // Last time step 
    vector<OCP_DBL>  lvfP;        ///< last vfP
    vector<OCP_DBL>  lvfT;        ///< last vfT
    vector<OCP_DBL>  lvfi;        ///< last vfi
    vector<OCP_DBL>  lrhoP;       ///< last rhoP
    vector<OCP_DBL>  lrhoT;       ///< last rhoT
    vector<OCP_DBL>  lrhox;       ///< last rhox
    vector<OCP_DBL>  lxiP;        ///< last xiP
    vector<OCP_DBL>  lxiT;        ///< last xiT
    vector<OCP_DBL>  lxix;        ///< last xix
    vector<OCP_DBL>  lmuP;        ///< last muP
    vector<OCP_DBL>  lmuT;        ///< last muT
    vector<OCP_DBL>  lmux;        ///< last mux
    vector<OCP_DBL>  ldPcj_dS;    ///< last Pcj_dS
    vector<OCP_DBL>  ldKr_dS;     ///< last dKr_dS
    vector<OCP_DBL>  lUfP;        ///< last UfP
    vector<OCP_DBL>  lUfT;        ///< last UfT
    vector<OCP_DBL>  lUfi;        ///< last Ufi
    vector<OCP_DBL>  lHT;         ///< last HT
    vector<OCP_DBL>  lHx;         ///< last Hx
    vector<OCP_DBL>  lktP;        ///< last ktP
    vector<OCP_DBL>  lktT;        ///< last ktT
    vector<OCP_DBL>  lktS;        ///< last ktS 


    /////////////////////////////////////////////////////////////////////
    // Newton Iteration Information
    /////////////////////////////////////////////////////////////////////   

public:
    /// Calculate some auxiliary variable, for example, dSmax
    OCP_DBL CalNRdSmax(OCP_USI& index);
    /// Return NRdPmax.
    OCP_DBL GetNRdPmax() const { return NRdPmax; };
    /// Return NRdSmaxP.
    OCP_DBL GetNRdSmaxP() const { return NRdSmaxP; };
    /// Return NRdNmax.
    OCP_DBL GetNRdNmax() const { return NRdNmax; };
            
protected:

    vector<OCP_DBL> dSNR;  ///< saturation change between NR steps
    vector<OCP_DBL> dSNRP; ///< predicted saturation change between NR steps
    vector<OCP_DBL> dNNR;  ///< Ni change between NR steps
    vector<OCP_DBL> dPNR;  ///< P  change between NR steps

    OCP_DBL NRdSSP;        ///< difference between dSNR and dSNRP, 2-norm
    OCP_DBL maxNRdSSP;     ///< max difference between dSNR and dSNRP
    OCP_USI index_maxNRdSSP; ///< index of grid which has maxNRdSSP
    OCP_DBL NRdPmax;      ///< Max pressure difference in an NR step
    OCP_DBL NRdNmax;      ///< Max Ni difference in an NR step
    OCP_DBL NRdSmax;      ///< Max saturation difference in an NR step(Real)
    OCP_DBL NRdSmaxP;     ///< Max saturation difference in an NR step(Predict)

    vector<OCP_DBL> NRstep; ///< NRstep for FIM
    vector<USI> NRphaseNum;     ///< phaseNum in NR step


    /////////////////////////////////////////////////////////////////////
    // Important Indicator Variable and Check
    /////////////////////////////////////////////////////////////////////

public:
    /// Check if negative P occurs, return OCP_FALSE if so.
    OCP_BOOL CheckP() const;
    /// Check if negative Ni occurs, return OCP_FALSE if so.
    OCP_BOOL CheckNi();
    /// Check if relative volume error is out of range, return OCP_FALSE if so.
    OCP_BOOL CheckVe(const OCP_DBL& Vlim) const;

    /// Calculate the CFL number.
    OCP_DBL CalCFL() const;
    /// Return maxCFL
    OCP_DBL GetMaxCFL() const { return maxCFL; }

    /// Calculate max change of indicator variables.
    void CalMaxChange();
    /// Return dPmax.
    OCP_DBL GetdPmax() const { return dPmax; }
    /// Return dNmax.
    OCP_DBL GetdNmax() const { return dNmax; }
    /// Return dSmax.
    OCP_DBL GetdSmax() const { return dSmax; }
    /// Return dVmax.
    OCP_DBL GetdVmax() const { return dVmax; }

protected:

    OCP_DBL dPmax; ///< Max change in pressure during the current time step.
    OCP_DBL dSmax; ///< Max change in saturation during the current time step.
    OCP_DBL dNmax; ///< Max change in moles of component during the current time step.
    OCP_DBL dVmax; ///< Max relative diff between fluid and pore volume during the current time step.


    mutable vector<OCP_DBL> cfl; ///< CFL number for each bulk
    mutable OCP_DBL maxCFL;      ///< max CFL number


    /////////////////////////////////////////////////////////////////////
    // Error
    /////////////////////////////////////////////////////////////////////

public:
    void AllocateError();

protected:
    vector<OCP_DBL>         ePEC; ///< error for fugacity balance equations, EoS only now


protected:
    /////////////////////////////////////////////////////////////////////
    // Method-Specified Variable
    /////////////////////////////////////////////////////////////////////

    USI              maxLendSdP;    ///< length of dSec_dPri.
    vector<USI>      bRowSizedSdP;  ///< length of dSec_dPri in each bulk
    vector<OCP_DBL>  dSec_dPri;     ///< d Secondary variable / d Primary variable.
    vector<OCP_BOOL> pSderExist;    ///< Existence of derivative of phase saturation
    vector<USI>      pVnumCom;      ///< num of variable components in the phase
    vector<OCP_DBL>  res_n;         ///< residual for FIM_n
    vector<OCP_DBL>  resPc;         ///< a precalculated value
    vector<OCP_USI>  resIndex;      ///< store the starting position of res_n of each bulk.
    
    // Last time step
    vector<USI>      lbRowSizedSdP; ///< last bRowSizedSdP
    vector<OCP_DBL>  ldSec_dPri;    ///< last dSec_dPri
    vector<OCP_BOOL> lpSderExist;   ///< last pSderExist
    vector<USI>      lpVnumCom;     ///< last pVnumCom
    vector<OCP_DBL>  lres_n;        ///< last res_n
    vector<OCP_DBL>  lresPc;        ///< last lresPc;
    vector<OCP_USI>  lresIndex;     ///< last res_n
    
public:

    /////////////////////////////////////////////////////////////////////
    // For IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for IMPEC.
    void AllocateIMPEC_IsoT();
    /// Perform Flash with Sj and calculate values needed for IMPEC
    void InitFlashIMPEC();
    /// Perform Flash with Ni and calculate values needed for IMPEC
    void CalFlashIMPEC();
    /// Pass value needed for IMPEC from flash to bulk
    void PassFlashValueIMPEC(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure with saturation.
    void CalKrPcIMPEC();
    /// Update P and Pj after linear system is solved.
    void GetSolIMPEC(const vector<OCP_DBL>& u);
    /// Reset variables needed for IMPEC to last time step
    void ResetVal01IMPEC();
    /// Reset variables needed for IMPEC to last time step
    void ResetVal02IMPEC();
    /// Reset variables needed for IMPEC to last time step
    void ResetVal03IMPEC();
    /// Update value of last step for IMPEC.
    void UpdateLastStepIMPEC();


    /////////////////////////////////////////////////////////////////////
    // For FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for FIM.
    void AllocateFIM_IsoT();
    /// Perform Flash with Sj and calculate values needed for FIM
    void InitFlashFIM();
    /// Perform Flash with Ni and calculate values needed for FIM
    void CalFlashFIM();
    /// Pass value needed for FIM from flash to bulk
    void PassFlashValueFIM(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure needed for FIM
    void CalKrPcFIM();
    /// Get the solution for FIM after a Newton iteration.
    void GetSolFIM(const vector<OCP_DBL>& u,
                   const OCP_DBL&         dPmaxlim,
                   const OCP_DBL&         dSmaxlim);
    /// Calculate relative residual for FIM.
    void CalRelResFIM(OCPRes& resFIM) const;
    /// Reset FIM.
    void ResetFIM();
    /// Update values of last step for FIM.
    void UpdateLastStepFIM();

    /////////////////////////////////////////////////////////////////////
    // For FIMn
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for FIMn
    void AllocateFIMn_IsoT();
    /// Perform Flash with Sj and calculate values needed for FIMn
    void InitFlashFIMn();
    /// Perform Flash with Ni and calculate values needed for FIMn
    void CalFlashFIMn();
    /// Perform Flash with Ni and calculate values needed for FIMn in Black Oil Model
    void CalFlashFIMn_BLKOIL();
    /// Perform Flash with Ni and calculate values needed for FIMn in Compositional Model
    void CalFlashFIMn_COMP();
    /// Pass value needed for FIMn from flash to bulk
    void PassFlashValueFIMn(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure needed for FIMn
    void CalKrPcFIMn() { CalKrPcFIM(); }
    /// Get the solution for FIMn after a Newton iteration.
    void GetSolFIMn(const vector<OCP_DBL>& u,
        const OCP_DBL& dPmaxlim,
        const OCP_DBL& dSmaxlim);
    /// Calculate relative residual for FIM.
    void CalRelResFIMn(OCPRes& resFIM) const { CalRelResFIM(resFIM); }
    /// Reset FIMn
    void ResetFIMn();
    /// Update values of last step for FIMn
    void UpdateLastStepFIMn();


    /////////////////////////////////////////////////////////////////////
    // For AIMc
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for AIMc
    void AllocateAIMc_IsoT();
    /// Perform flash calculation with Ni for Explicit bulk -- Update partial properties
    void CalFlashAIMcEp();
    /// Perform flash calculation with Ni for Explicit bulk -- Update all properties
    void CalFlashAIMcEa();
    /// Perform flash calculation with Ni for Implicit bulk
    void CalFlashAIMcI();
    /// Pass flash value needed for Explicit bulk -- Update partial properties 
    void PassFlashValueAIMcEp(const OCP_USI& n);
    /// Pass flash value needed for Explicit bulk -- Update all properties 
    void PassFlashValueAIMcEa(const OCP_USI& n) { PassFlashValueIMPEC(n); }
    /// Pass value needed for Implicit bulk
    void PassFlashValueAIMcI(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure for Explicit bulk
    void CalKrPcAIMcE();
    /// Calculate relative permeability and capillary pressure for Implicit bulk
    void CalKrPcAIMcI();
    /// Get the solution for AIMc after a Newton iteration.
    void GetSolAIMc(const vector<OCP_DBL>& u,
                    const OCP_DBL&         dPmaxlim,
                    const OCP_DBL&         dSmaxlim);
    /// Reset AIMc
    void ResetAIMc();
    /// Update values of last step for AIMc
    void UpdateLastStepAIMc();
    /// Print Bulk which are implicit
    void ShowFIMBulk(const OCP_BOOL& flag = OCP_FALSE) const;
    /// clear wellBulkId
    void ClearWellBulkId() { wellBulkId.clear(); }
    /// push back an element for wellBulkId
    void AddWellBulkId(const OCP_USI& n) { wellBulkId.push_back(n); }

protected:
    vector<OCP_USI> wellBulkId;    ///< Index of bulks which are penetrated by wells and their K-neighbor
    class BulkTypeAIM
    {
    public:
        void Setup(const OCP_USI& nb) { indicator.resize(nb, -1); }
        void Init() { fill(indicator.begin(), indicator.end(), -1); }
        void SetBulkType(const OCP_USI& n, const OCP_INT& flag) { indicator[n] = flag; }
        OCP_BOOL IfFIMbulk(const OCP_USI& n)const { return indicator[n] > 0; }
        OCP_BOOL IfIMPECbulk(const OCP_USI& n)const { return indicator[n] < 0; }
    protected:
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