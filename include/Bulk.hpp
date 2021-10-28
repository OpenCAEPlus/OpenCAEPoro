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
#include "FlowUnit.hpp"
#include "Grid.hpp"
#include "LinearSolver.hpp"
#include "Mixture.hpp"
#include "MixtureBO.hpp"
#include "OCPConst.hpp"
#include "ParamReservoir.hpp"

using namespace std;

/// ParamEQUIL contains some initial reservoir infomation used to calculate initial
/// Equilibration, including reference depth and pressure at it, depth of contacts
/// between phases, and capillary at them.
class ParamEQUIL
{
    friend class Bulk;

private:
    OCP_DBL Dref;  ///< reference depth.
    OCP_DBL Pref;  ///< pressure in reference depth.
    OCP_DBL DOWC;  ///< depth of oil-water contact.
    OCP_DBL PcOWC; ///< capillary pressure in oil-water contact: Po - Pw.
    OCP_DBL DGOC;  ///< depth of gas-oil contact.
    OCP_DBL PcGOC; ///< capillary pressure in gas-oil contact: Pg - Po.

    OCPTable PBVD; ///< PBVD Table: buble point pressere vs depth.
};

/// Bulk contains main physical infomation of reservoir, but only active grids are
/// stored. acturally is's real geometric "area" for simulating. variables here are
/// stored bulks by bulks, and then phases by phases, components by components. the
/// bulks are ordered with the X axis index cycling fastest, followed by the Y and Z
/// axis indices defaulted. operations refered to single bulk are included here.
class Bulk
{
    friend class BulkConn;
    friend class Well;
    friend class DetailInfo;

public:
    Bulk() = default;

    /// Return the num of active bulk.
    OCP_USI GetBulkNum() const { return numBulk; }
    /// Input param from internal data structure of param: ParamReservoir.
    void InputParam(ParamReservoir& rs_param);
    /// Setup basic data from grid, Setup of grid must be finished before.
    void Setup(const Grid& myGrid);
    /// Calculate initial equilibration for blkoil model according to EQUIL.
    /// tabrow is maximum number of depth nodes in table of depth vs pressure.
    void InitSjPcBlk(const USI& tabrow);
    /// Calculate initial equilibration for compositional model according to EQUIL.
    /// tabrow is maximum number of depth nodes in table of depth vs pressure.
    void InitSjPcComp(const USI& tabrow);
    /// Assignment value for some variable, it's called when one time step finished.
    void UpdateLastStep();

    /// calculate max change of some physical variables to predict the size of next time
    /// step.
    void CalMaxChange();

    /// Flash calculation, initial saturation is needed in blackoil model, and initial
    /// zi is needed in compositional model.
    void FlashSj();
    /// Flash calculation, moles of component and pressure are needed both in blackoil
    /// model and compositional model.
    void FlashNi();
    void FlashNiDeriv();
    /// when flash calculation finished, values in flash class should be passed to ones
    /// in bulk class.
    void PassFlashValue(const OCP_USI& n);
    void PassFlashValueDeriv(const OCP_USI& n);
    /// Calculate relative permeability and capillary pressure with saturation.
    void CalKrPc();
    void CalKrPcDeriv();
    /// Calculate volume of pore with pressure.
    void CalVpore();
    
    /// pass flash to well, it hasn't been used now.
    const vector<Mixture*>& GetMixture() const { return flashCal; }

    /// return the mixture mode.
    USI GetMixMode() const;

    USI GetComNum()const { return numCom; }
    /// get solution from solver class after linear system is solved.
    void GetSolIMPEC(const vector<OCP_DBL>& u);
    void GetSolFIM(const vector<OCP_DBL>& u);

    /// calculate average pressure in reservoir.
    OCP_DBL CalFPR() const;
    /// return pressure in ith bulk.
    OCP_DBL GetP(const OCP_USI& n) const { return P[n]; }
    /// Return oil saturation in ith bulk.
    OCP_DBL GetSOIL(const OCP_USI& n) const { return S[n * numPhase + phase2Index[OIL]]; }
    /// Return gas saturation in ith bulk.
    OCP_DBL GetSGAS(const OCP_USI& n) const { return S[n * numPhase + phase2Index[GAS]]; }
    /// Return water saturation in ith bulk.
    OCP_DBL GetSWAT(const OCP_USI& n) const { return S[n * numPhase + phase2Index[WATER]]; }
    /// check if negative P occurs, return false if so.
    bool CheckP() const;
    /// Check if negative Ni occurs, return false if so.
    bool CheckNi() const;
    /// Check if relative volume error is outranged, return false if so.
    bool CheckVe(const OCP_DBL& Vlim) const;
    /// Reset P to the ones at last time step.
    void ResetP() { P = lP; }
    /// Reset Pj to the ones at last time step.
    void ResetPj() { Pj = lPj; }
    /// Reset Ni to the ones at last time step.
    void ResetNi() { Ni = lNi; }
    /// Reset Vp to the ones at last time step.
    void ResetVp() { rockVp = rockLVp; }
    /// Reset variables in flash calculations.
    void ResetFlash();
    /// Check difference from last time step.
    void CheckDiff();
    /// Check if the sum of saturations is one.
    void CheckSat() const;
    /// Return dPmax.
    OCP_DBL GetdPmax() const { return dPmax; }
    /// Return dNmax.
    OCP_DBL GetdNmax() const { return dNmax; }
    /// Return dSmax.
    OCP_DBL GetdSmax() const { return dSmax; }
    /// Return dVmax.
    OCP_DBL GetdVmax() const { return dVmax; }
    /// Initialize cfl number.
    void InitCFL() const { cfl.assign(numBulk * numPhase, 0); }
    /// Calculate the cfl number.
    OCP_DBL CalCFL(bool flag) const;

private:
    OCP_USI numBulk; ///< num of bulks (active grids).

    // Physical infomation
    USI numPhase; ///< num of phase.
    USI numCom;   ///< num of component.

    OCP_DBL         T;          ///< temperature: numBulk.
    vector<OCP_DBL> Pb;         ///< buble point pressere: numBulk.
    vector<OCP_DBL> P;          ///< pressure: numBulk.
    vector<OCP_DBL> Pj;         ///< pressure of phase: numPhase*numBulk.
    vector<OCP_DBL> Pc;         ///< capillary pressure of phase: numPhase*numBulk.
    vector<bool>    phaseExist; ///< existence of phase: numPhase*numBulk.
    vector<OCP_DBL> S;          ///< saturation of phase j: numPhase*numBulk.
    vector<OCP_DBL> rho;        ///< mass density of phase: numPhase*numBulk.
    vector<OCP_DBL> xi;         ///< moles density of phase: numPhase*numBulk.
    vector<OCP_DBL> cij; ///< Nij / Nj: numPhase*numCom*numBulk. Nij is the moles of
                         ///< component i in phase j, Nj is the moles of phase j.
    vector<OCP_DBL> Ni;  ///< moles of component: numCom*numBulk.
    vector<OCP_DBL> mu;  ///< viscosity of phase: numPhase*numBulk.
    vector<OCP_DBL> kr;  ///< relative permeability of phase: numPhase*numBulk.

    vector<OCP_DBL> vj;  ///< volume of phase: numPhase*numBulk.
    vector<OCP_DBL> vf;  ///< total fluid volume: numBulk.
    // Derivatives
    // For IMPEC
    vector<OCP_DBL> vfi; ///< dVf / dNi: numCom*numBulk.
    vector<OCP_DBL> vfp; ///< dVf / dP: numBulk.
    // For FIM
    vector<OCP_DBL> muP;    ///< dMu / dP: numPhase*numBulk.
    vector<OCP_DBL> xiP;    ///< dXi / dP: numPhase*numBulk.
    vector<OCP_DBL> rhoP;   ///< dRho / dP: numPhase*numBulk.
    vector<OCP_DBL> muC;    ///< dMuj / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> xiC;    ///< dXi_j / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> rhoC;   ///< dRhoj / dxij: numPhase*numCom*numBulk.
    vector<OCP_DBL> dSec_dPri; ///< d Second var / d Primary var: (numPhase + numPhase*numCom)*(numCom + 1)*numBulk
    vector<OCP_DBL> dPcj_dS; ///< d Pcj / dSk: numPhase * numPhase * bulk.
    vector<OCP_DBL> dKr_dS; ///< d Krj / dSk: numPhase * numPhase * bulk.

    /// used to identify phase name according to its index: numPhase.
    /// For example, 0th phase is Oil.
    vector<USI> index2Phase; 
    /// used to tell the location of phase according to its name: numPhase.
    /// For example, Oil is in ith location.
    vector<USI> phase2Index;


    vector<OCP_DBL> initZi; ///< initial proportion of each component for EoS : numCom -
                            ///< 1, water is excluded.
    USI         PVTmode;    ///< used to identify PVT mode in blackoil model.
    vector<USI> PVTNUM;     ///< used to identify PVT region in blackoil model: numBulk.
    vector<Mixture*>  flashCal; ///< a flash class used to flash calculation.
    USI               SATmode;  ///< used to identify SAT mode.
    vector<USI>       SATNUM;   ///< used to identify SAT region: numBulk.
    vector<FlowUnit*> flow;     ///< a flow class used to calculate capillary pressure,
                                ///< relative permeability.

    // Properties at last STEP
    vector<OCP_DBL> lP;  ///< Pressure at last time step: numBulk.
    vector<OCP_DBL> lPj; ///< Pressure of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lPc; ///< Capillary pressure of phase at last time step: numPhase*numBulk.
    vector<bool>    lphaseExist; ///< existence of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lS;  ///< Saturation of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lrho; ///< Mass density of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lxi; ///< Moles density of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lcij; ///< Nij / Nj at last time step: numPhase*numCom*numBulk.
    vector<OCP_DBL> lNi; ///< Moles of component at last time step: numCom*numBulk.
    vector<OCP_DBL> lmu; ///< Viscosity of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lkr; ///< Relative permeability of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lvj; ///< Volume of phase at last time step: numPhase*numBulk.
    vector<OCP_DBL> lvf; ///< Total fluid volume at last time step: numBulk.
    vector<OCP_DBL> lvfi; ///< dVf / dNi at last time step: numCom*numBulk.
    vector<OCP_DBL> lvfp; ///< dVf / dP at last time step: numBulk.
    vector<OCP_DBL> rockLVp; ///< volume of pore at last time step: numBulk.

    // CFL number
    mutable vector<OCP_DBL> cfl;

    // max change
    OCP_DBL dPmax; ///< max change in pressure during current time step.
    OCP_DBL dNmax; ///< max change in moles of component during current time step.
    OCP_DBL dVmax; ///< max relative error between fluid volume and pore volume during
                   ///< current time step.
    OCP_DBL dSmax; ///< max change in saturation during current time step.

    // Reservoir rock infomation
    vector<OCP_DBL> dx;    ///< size of bulk along the x direction: numBulk.
    vector<OCP_DBL> dy;    ///< size of bulk along the y direction: numBulk.
    vector<OCP_DBL> dz;    ///< size of bulk along the z direction: numBulk.
    vector<OCP_DBL> depth; ///< depth of center of bulk: numBulk.
    vector<OCP_DBL> ntg;   ///< net to gross of bulk: numBulk.
    vector<OCP_DBL>
        rockVpInit;         ///< initial pore volume: Vgrid * ntg * poro_init: numBulk.
    vector<OCP_DBL> rockVp; ///< pore volume: Vgrid * ntg * poro: numBulk.
    OCP_DBL         rockPref; ///< reference pressure for initial rock volume.
    OCP_DBL         rockC1;   ///< rock compressibility.
    OCP_DBL         rockC2;   ///< rock compressibility.
    vector<OCP_DBL>
        rockKxInit; ///< initial rock permeability along the x direction: numBulk.
    vector<OCP_DBL>
        rockKx; ///< current rock permeability along the x direction: numBulk.
    vector<OCP_DBL>
        rockKyInit; ///< initial rock permeability along the y direction: numBulk.
    vector<OCP_DBL>
        rockKy; ///< current rock permeability along the y direction: numBulk.
    vector<OCP_DBL>
        rockKzInit; ///< initial rock permeability along the z direction: numBulk.
    vector<OCP_DBL>
        rockKz; ///< current rock permeability along the z direction: numBulk.
    // vector<OCP_DBL>		Rock_Poro;			// current porosity
    // vector<OCP_DBL>		Rock_PoroInit;		// initial porosity

    // basic Reservoir model
    ParamEQUIL EQUIL;    ///< initial Equilibration.
    bool       blackOil; ///< if true, blackoil model will be used.
    bool       comps;    ///< if true, compositional model will be used.
    bool       oil;      ///< if true, oil phase could exist.
    bool       gas;      ///< if true, gas phase could exist.
    bool       water;    ///< if true, water phase could exist.
    bool       disGas;   ///< if true, dissolve gas in live oil could exist.
};

#endif /* end if __BULK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/