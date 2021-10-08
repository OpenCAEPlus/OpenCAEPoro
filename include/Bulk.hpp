/*! \file    Bulk.hpp
 *  \brief   Bulk class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2021
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
#include "BOMixture.hpp"
#include "FlowUnit.hpp"
#include "Grid.hpp"
#include "Mixture.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamReservoir.hpp"
#include "Solver.hxx"


using namespace std;


/// ParamEQUIL contains some initial reservoir infomation used to calculate initial Equilibration,
/// including reference depth and pressure at it, depth of contacts between phases, and capillary at them.
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
    
    ReservoirTable<OCP_DBL> PBVD; ///< PBVD Table: buble point pressere vs depth.
};


/// Bulk contains main physical infomation of reservoir, but only active grids are stored.
/// acturally is's real geometric "area" for simulating.
/// variables here are stored bulks by bulks, and then phases by phases, components by components.
/// the bulks are ordered with the X axis index cycling fastest, followed by the Y and Z axis indices defaulted.
/// operations refered to single bulk are included here.
class Bulk
{
    friend class Connection_BB;
    friend class Well;

public:
    Bulk() = default;

    /// return the num of active bulk.
    OCP_USI getBulkNum() const { return Num; }
    /// input param from internal data structure of param: ParamReservoir.
    void inputParam(ParamReservoir& rs_param);
    /// setup basic data from grid, setup of grid must be finished before. 
    void setup(const Grid& myGrid);
    /// calculate initial equilibration for blkoil model according to EQUIL.
    /// tabrow is maximum number of depth nodes in table of depth vs pressure.
    void initSjPc_blk(const USI& tabrow);
    /// calculate initial equilibration for compositional model according to EQUIL.
    /// tabrow is maximum number of depth nodes in table of depth vs pressure.
    void initSjPc_comp(const USI& tabrow);
    /// assignment value for some variable, it's called when one time step finished.
    void setLastStep()
    {
        lP       = P;
        lPj      = Pj;
        lNi      = Ni;
        lS       = S;
        Rock_lVp = Rock_Vp;
    }
    /// calculate max change of some physical variables to predict the size of next time step.
    void calMaxChange();

    /// Flash calculation, initial saturation is needed in blackoil model, and initial zi is needed in compositional model.
    void flash_Sj();
    /// Flash calculation, moles of component and pressure are needed both in blackoil model and compositional model.
    void flash_Ni();
    /// when flash calculation finished, values in flash class should be passed to ones in bulk class.
    void passFlashValue(const OCP_USI& n);

    /// calculate relative permeability and capillary pressure with saturation.
    void calKrPc();
    /// calculate volume of pore with pressure.
    void calVporo();

    /// pass flash to well, it hasn't been used now.
    const std::vector<Mixture*>& getMixture() const { return Flashcal; }

    /// return the mixture mode.
    USI mixMode() const;

    /// get solution from solver class after linear system is solved.
    void getSol_IMPES(const vector<OCP_DBL>& u);

    /// calculate average pressure in reservoir.
    OCP_DBL calFPR() const;
    /// return pressure in ith bulk.
    OCP_DBL getP(int n) const { return P[n]; }
    /// check if negative P occurs, return false if so.
    bool    checkP() const;
    /// check if negative Ni occurs, return false if so.
    bool    checkNi() const;
    /// check if relative volume error is outranged, return false if so.
    bool    checkVe(const OCP_DBL& Vlim) const;
    /// reset P to the ones at last time step.
    void    resetP() { P = lP; }
    /// reset Pj to the ones at last time step.
    void    resetPj() { Pj = lPj; }
    /// reset Ni to the ones at last time step.
    void    resetNi() { Ni = lNi; }
    /// reset Vp to the ones at last time step.
    void    resetVp() { Rock_Vp = Rock_lVp; }

    /// return dPmax.
    OCP_DBL getdPmax() const { return dPmax; }
    /// return dNmax.
    OCP_DBL getdNmax() const { return dNmax; }
    /// return dSmax.
    OCP_DBL getdSmax() const { return dSmax; }
    /// return dVmax.
    OCP_DBL getdVmax() const { return dVmax; }

private:

    OCP_USI Num; ///< num of active bulk.

    // Physical infomation
    USI Np; ///< num of phase.
    USI Nc; ///< num of component.

    OCP_DBL              T;          ///< temperature: Num.
    std::vector<OCP_DBL> Pbub;       ///< buble point pressere: Num.
    std::vector<OCP_DBL> P;          ///< pressure: Num.
    std::vector<OCP_DBL> Pj;         ///< pressure of phase: Np*Num.
    std::vector<OCP_DBL> Pc;         ///< capillary pressure of phase: Np*Num.
    std::vector<bool>    PhaseExist; ///< existence of phase: Np*Num.
    std::vector<OCP_DBL> S;          ///< saturation of phase j: Np*Num.
    std::vector<OCP_DBL> Rho;        ///< mass density of phase: Np*Num.
    std::vector<OCP_DBL> Xi;         ///< moles density of phase: Np*Num.
    std::vector<OCP_DBL> Cij;        ///< Nij / Nj: Np*Nc*Num. Nij is the moles of component i in phase j, Nj is the moles of phase j.
    std::vector<OCP_DBL> Ni;         ///< moles of component: Nc*Num.
    std::vector<OCP_DBL> Mu;         ///< viscosity of phase: Np*Num.
    std::vector<OCP_DBL> Kr;         ///< relative permeability of phase: Np*Num.

    std::vector<OCP_DBL> Vj;         ///< volume of phase: Np*Num.
    std::vector<OCP_DBL> Vf;         ///< total fluid volume: Num.
    std::vector<OCP_DBL> Vfi;        ///< dVf / dNi: Num.
    std::vector<OCP_DBL> Vfp;        ///< dVf / dP.

    vector<USI>            PhaseLabel;  ///< used to identify phase according to its index: Np.
    std::vector<OCP_DBL>   InitZi;      ///< initial proportion of each component for EoS : Nc - 1, water is excluded.
    USI                    PVTmode;     ///< used to identify PVT mode in blackoil model.
    std::vector<USI>       PVTNUM;      ///< used to identify PVT region in blackoil model: Num.
    std::vector<Mixture*>  Flashcal;    ///< a flash class used to flash calculation.
    USI                    SATmode;     ///< used to identify SAT mode.
    std::vector<USI>       SATNUM;      ///< used to identify SAT region: Num.
    std::vector<FlowUnit*> Flow;        ///< a flow class used to calculate capillary pressure, relative permeability.

    // infomation at last STEP
    std::vector<OCP_DBL> lP;            ///< pressure at last time step: Num.
    std::vector<OCP_DBL> lPj;           ///< capillary pressure at last time step: Np*Num.
    std::vector<OCP_DBL> lNi;           ///< moles of component at last time step: Nc*Num.
    std::vector<OCP_DBL> lS;            ///< saturation of phase at last time step: Np*Num.
    std::vector<OCP_DBL> Rock_lVp;      ///< volume of pore at last time step: Num.

    // max change
    OCP_DBL dPmax;                      ///< max change in pressure during current time step.
    OCP_DBL dNmax;                      ///< max change in moles of component during current time step.
    OCP_DBL dVmax;                      ///< max relative error between fluid volume and pore volume during current time step.
    OCP_DBL dSmax;                      ///< max change in saturation during current time step.

    // Reservoir rock infomation
    std::vector<OCP_DBL> Dx;            ///< size of bulk along the x direction: Num.
    std::vector<OCP_DBL> Dy;            ///< size of bulk along the y direction: Num.
    std::vector<OCP_DBL> Dz;            ///< size of bulk along the z direction: Num.
    std::vector<OCP_DBL> Depth;         ///< depth of center of bulk: Num.
    std::vector<OCP_DBL> Ntg;           ///< net to gross of bulk: Num.
    std::vector<OCP_DBL> Rock_VpInit;   ///< initial pore volume: Vgrid * ntg * poro_init: Num.
    std::vector<OCP_DBL> Rock_Vp;       ///< pore volume: Vgrid * ntg * poro: Num.
    OCP_DBL              Rock_Pref;     ///< reference pressure for initial rock volume.
    OCP_DBL              Rock_C1;       ///< rock compressibility.
    OCP_DBL              Rock_C2;       ///< rock compressibility.
    std::vector<OCP_DBL> Rock_KxInit;   ///< initial rock permeability along the x direction: Num.
    std::vector<OCP_DBL> Rock_Kx;       ///< current rock permeability along the x direction: Num.
    std::vector<OCP_DBL> Rock_KyInit;   ///< initial rock permeability along the y direction: Num.
    std::vector<OCP_DBL> Rock_Ky;       ///< current rock permeability along the y direction: Num.
    std::vector<OCP_DBL> Rock_KzInit;   ///< initial rock permeability along the z direction: Num.
    std::vector<OCP_DBL> Rock_Kz;       ///< current rock permeability along the z direction: Num.
    // std::vector<OCP_DBL>		Rock_Poro;			// current porosity
    // std::vector<OCP_DBL>		Rock_PoroInit;		// initial porosity

    // basic Reservoir model
    ParamEQUIL EQUIL;                   ///< initial Equilibration.
    bool       BLACKOIL;                ///< if true, blackoil model will be used.
    bool       COMPS;                   ///< if true, compositional model will be used.
    bool       Oil;                     ///< if true, indicates oil phase could exist.
    bool       Gas;                     ///< if true, indicates gas phase could exist.
    bool       Water;                   ///< if true, indicates water phase could exist. 
    bool       DisGas;                  ///< if true, indicates dissolve gas in live oil could exist. 
};

#endif


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/