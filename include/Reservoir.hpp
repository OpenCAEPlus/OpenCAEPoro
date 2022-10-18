/*! \file    Reservoir.hpp
 *  \brief   Reservoir class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIR_HEADER__
#define __RESERVOIR_HEADER__

// OpenCAEPoro header files
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "Grid.hpp"
#include "ParamRead.hpp"
#include "AllWells.hpp"

/// Reservoir is the core component in our simulator, it contains the all reservoir
/// information, and all operations on it.
///
/// Reservoir has four Core components.
/// Grids contains the basic informations of all grids as a database of reservoir.
/// Bulk only stores active grids, which defines the area used for calculation.
/// AllWells contains the well information, it's used to manage operations related to
/// wells. BulkConn contains connections between bulks(active grids).
class Reservoir
{
    friend class OCPControl;
    friend class Summary;
    friend class CriticalInfo;
    friend class DetailInfo;

    // temp
    friend class OCP_IMPEC;
    friend class OCP_FIM;
    friend class OCP_AIMs;
    friend class OCP_AIMt;
    friend class OCP_AIMc;
    friend class Solver;

    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:
    /// Input param from internal param data structure, which stores the params from
    /// input files.
    void InputParam(ParamRead& param);
    /// Setup static information for reservoir with input params.
    void Setup();
    /// Apply the control of ith critical time point.
    void ApplyControl(const USI& i);
    /// Calculate Well Properties at the beginning of each time step.
    void PrepareWell();
    /// Calculate Flux between Bulk and Wells.
    void CalWellFlux();
    /// Calculate Trans of Wells.
    void CalWellTrans();
    /// Calculate pore of Bulks.
    void CalVpore();
    /// Calculate Relative Permeability and Capillary for each Bulk
    void CalKrPc();
    /// Calculate Maximum Change of some reference variables for IMPEC
    void CalMaxChange();
    /// Calculate num of Injection, Production
    void CalIPRT(const OCP_DBL& dt);
    /// Check if abnormal Pressure occurs
    OCP_INT CheckP(const OCP_BOOL& bulkCheck = OCP_TRUE, const OCP_BOOL& wellCheck = OCP_TRUE);
    /// Check if abnormal Pressure occurs
    OCP_BOOL CheckNi();
    /// Check error between Fluids and Pores
    OCP_BOOL CheckVe(const OCP_DBL& Vlim) const;
    /// Return the num of Bulk
    OCP_USI GetBulkNum() const { return bulk.GetBulkNum(); }
    /// Return MaxNUMFIMBulk
    OCP_USI GetMaxFIMBulk() const { return bulk.GetMaxFIMBulk(); }
    /// Return the num of Well
    USI GetWellNum() const { return allWells.GetWellNum(); }
    /// Return the num of Components
    USI GetComNum() const { return bulk.GetComNum(); }
    void SetupWellBulk() { allWells.SetupWellBulk(bulk); }
    void GetNTQT(const OCP_DBL& dt);

private:
    Grid      grid;      ///< Grid class.
    Bulk      bulk;      ///< Bulk class.
    AllWells  allWells; ///< AllWells class.
    BulkConn  conn;      ///< BulkConn class.

    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for IMPEC
    void AllocateAuxIMPEC();
    /// Initialize the properties of Reservoir for IMPEC
    void InitIMPEC();
    /// Calcluate the CFL number, including bulks and wells for IMPEC
    OCP_DBL CalCFL(const OCP_DBL& dt);
    /// Calculate flux between bulks, bulks and wells
    void CalFLuxIMPEC();
    /// Calculate flux between bulks
    void CalConnFluxIMPEC();
    /// Calculate Ni according to Flux
    void MassConseveIMPEC(const OCP_DBL& dt);
    /// Calculate Flash For IMPEC
    void CalFlashIMPEC();
    /// Update value of last step for IMPEC
    void UpdateLastStepIMPEC();
    /// Allocate Maxmimum memory for internal Matirx for IMPEC
    void AllocateMatIMPEC(LinearSystem& myLS) const;
    /// Assemble Matrix for IMPEC
    void AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Return the Solution to Reservoir Pressure for IMPEC
    void GetSolutionIMPEC(const vector<OCP_DBL>& u);
    /// Reset Well for IMPEC
    void ResetWellIMPEC();
    /// Reset Capillary Pressure, Flux for IMPEC
    void ResetVal01IMPEC();
    /// Reset Capillary Pressure, Moles of Componnets, Flux for IMPEC
    void ResetVal02IMPEC();
    /// Reset Pressure, Capillary Pressure, Moles of components, Flux, Volume of Pores
    /// for IMPEC
    void ResetVal03IMPEC();

private:

    OCP_DBL cfl; ///< CFL number.

public:
    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

    /// Allocate memory for auxiliary variables used for FIM
    void AllocateAuxFIM();
    /// Initialize the properties of Reservoir for FIM
    void InitFIM();
    void InitFIM_n();
    /// Calculate Flash for FIM, some derivatives are needed
    void CalFlashDerivFIM();
    void CalFlashDerivFIM_n();
    /// Calculate Relative Permeability and Capillary and some derivatives for each Bulk
    void CalKrPcDerivFIM();
    /// Update value of last step for FIM.
    void UpdateLastStepFIM();
    /// Allocate Maxmimum memory for internal Matirx for FIM
    void AllocateMatFIM(LinearSystem& myLS) const;
    /// Assemble Matrix for FIM
    void AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const;
    void AssembleMatFIM_n(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Return the Solution to Reservoir Pressure and moles of Components for FIM
    /// Exactly, it's a Newton step.
    void GetSolutionFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
                        const OCP_DBL& dSmax);
    void GetSolutionFIM_n(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
        const OCP_DBL& dSmax);
    /// Calculate the Resiual for FIM, it's also RHS of Linear System
    void CalResFIM(ResFIM& resFIM, const OCP_DBL& dt);
    /// Reset FIM
    void ResetFIM(const OCP_BOOL& flag);
    /// Return NRdPmax
    OCP_DBL GetNRdPmax(){ return bulk.GetNRdPmax(); }
    /// Return NRdSmax
    OCP_DBL GetNRdSmax(OCP_USI& index) { return bulk.CalNRdSmax(index); }
    /// Return NRdNmax
    OCP_DBL GetNRdNmax() { return bulk.GetNRdNmax(); }
    /// Return NRdSmaxP
    OCP_DBL GetNRdSmaxP() { return bulk.GetNRdSmaxP(); }
    void    PrintSolFIM(const string& outfile) const;
    void    ShowRes(const vector<OCP_DBL>& res) const;

    /////////////////////////////////////////////////////////////////////
    // AIMs, AIMt
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for AIMt
    void AllocateAuxAIMt();
    /// Allocate Maxmimum memory for internal Matirx for local FIM
    void AllocateMatAIMt(LinearSystem& myLS) const;
    /// Setup FIMBulk
    void SetupFIMBulk(const OCP_BOOL& NRflag = OCP_FALSE) { conn.SetupFIMBulk(bulk, NRflag); }
    void AddFIMBulk() { conn.AddFIMBulk(bulk); }
    void SetupFIMBulkBoundAIMs() { conn.SetupFIMBulkBoundAIMs(bulk); }
    
    /// Calculate Flash for local FIM, some derivatives are needed
    void CalFlashDerivAIM(const OCP_BOOL& IfAIMs);
    /// Calculate Relative Permeability and Capillary and some derivatives for each Bulk
    void CalKrPcDerivAIM(const OCP_BOOL& IfAIMs);
    /// Calculate the Resiual for local FIM, it's also RHS of Linear System
    void CalResAIMt(ResFIM& resFIM, const OCP_DBL& dt);
    /// Assemble Matrix for AIMt ---- local FIM here
    void AssembleMatAIMt(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Return the Solution to Reservoir Pressure and moles of Components for FIM
    /// Exactly, it's a Newton step.
    void GetSolutionAIMt(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
        const OCP_DBL& dSmax);

    /// Allocate memory for auxiliary variables used for AIMs
    void AllocateAuxAIMs();
    /// Calculate the Resiual for AIMs, it's also RHS of Linear System
    void CalResAIMs(ResFIM& resFIM, const OCP_DBL& dt);
    /// Assemble Matrix for AIMs
    void AssembleMatAIMs(LinearSystem& myLS, vector<OCP_DBL>& res, const OCP_DBL& dt) const;
    void GetSolutionAIMs(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
        const OCP_DBL& dSmax);
    void ResetValAIM();
    OCP_DBL CalCFLAIM(const OCP_DBL& dt);
    /// Update value of last step for IMPEC
    void UpdateLastStepAIM();

    /// Allocate memory for auxiliary variables used for FIM
    void AllocateAuxAIMc();
    /// Assemble Matrix for AIMc
    void AssembleMatAIMc(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Calculate the Resiual for FIM, it's also RHS of Linear System
    void CalResAIMc(ResFIM& resFIM, const OCP_DBL& dt);
    void CalFlashAIMc();
    void CalFlashAIMc01();
    void CalKrPcAIMc();
    /// Calculate Flash for local FIM, some derivatives are needed
    void CalFlashDerivAIMc();
    /// Calculate Relative Permeability and Capillary and some derivatives for each Bulk
    void CalKrPcDerivAIMc();
    void GetSolutionAIMc(const vector<OCP_DBL>& u, const OCP_DBL& dPmax,
        const OCP_DBL& dSmax);
    void InitAIMc();
    void UpdatePj() { bulk.UpdatePj(); }
    
};

#endif /* end if __RESERVOIR_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/