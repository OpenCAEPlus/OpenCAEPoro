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
#include "AllWells.hpp"
#include "Bulk.hpp"
#include "BulkConn.hpp"
#include "Grid.hpp"
#include "ParamRead.hpp"

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
    friend class Out4RPT;
    friend class Out4VTK;

    // temp
    friend class OCP_IMPEC;
    friend class OCP_FIM;
    friend class OCP_FIMn;
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
    void SetupIsoT();
    /// Apply the control of ith critical time point.
    void ApplyControl(const USI& i);
    /// Calculate Well Properties at the beginning of each time step.
    void PrepareWell();
    /// Calculate Flux between Bulk and Wells.
    void CalWellFlux();
    /// Calculate Trans of Wells.
    void CalWellTrans();
    /// Calculate pore of Bulks.
    void CalRock();
    /// Calculate Relative Permeability and Capillary for each Bulk
    void CalKrPc();
    /// Calculate Maximum Change of some reference variables for IMPEC
    void CalMaxChange();
    /// Calculate num of Injection, Production
    void CalIPRT(const OCP_DBL& dt);
    /// Check if abnormal Pressure occurs
    OCP_INT CheckP(const OCP_BOOL& bulkCheck = OCP_TRUE,
                   const OCP_BOOL& wellCheck = OCP_TRUE);
    /// Check if abnormal Pressure occurs
    OCP_BOOL CheckNi();
    /// Check error between Fluids and Pores
    OCP_BOOL CheckVe(const OCP_DBL& Vlim) const;
    /// Return the num of Bulk
    OCP_USI GetBulkNum() const { return bulk.GetBulkNum(); }
    /// Return the num of Well
    USI GetWellNum() const { return allWells.GetWellNum(); }
    /// Return the num of Components
    USI  GetComNum() const { return bulk.GetComNum(); }
    void SetupWellBulk() { allWells.SetupWellBulk(bulk); }


protected:

    Grid     grid;     ///< Init Grid class.
    Bulk     bulk;     ///< Active grid info
    AllWells allWells; ///< AllWells class.
    BulkConn conn;     ///< BulkConn class.

public:
    /// Calculate the CFL number, including bulks and wells for IMPEC
    OCP_DBL CalCFL(const OCP_DBL& dt);
    /// Return NRdPmax
    OCP_DBL GetNRdPmax() { return bulk.GetNRdPmax(); }
    /// Return NRdSmax
    OCP_DBL GetNRdSmax(OCP_USI& index) { return bulk.CalNRdSmax(index); }
    /// Return NRdNmax
    OCP_DBL GetNRdNmax() { return bulk.GetNRdNmax(); }
    /// Return NRdSmaxP
    OCP_DBL GetNRdSmaxP() { return bulk.GetNRdSmaxP(); }
    void    PrintSolFIM(const string& outfile) const;
    void OutInfoFinal() const { bulk.OutMixtureIters(); }

    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for IMPEC
    void AllocateIMPEC_IsoT();
    /// Initialize the properties of Reservoir for IMPEC
    void InitIMPEC();
    /// Allocate maximal memory for internal Matrix for IMPEC
    void AllocateMatIMPEC(LinearSystem& myLS) const;
    /// Assemble matrix for IMPEC
    void AssembleMatIMPEC(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Return the Solution to Reservoir Pressure for IMPEC
    void GetSolutionIMPEC(const vector<OCP_DBL>& u);
    /// Calculate flux between bulks, bulks and wells
    void CalFLuxIMPEC();
    /// Calculate flux between bulks
    void CalConnFluxIMPEC() { conn.CalFluxIMPEC(bulk); }
    /// Calculate Ni according to Flux
    void MassConserveIMPEC(const OCP_DBL& dt);
    /// Calculate Flash For IMPEC
    void CalFlashIMPEC() { bulk.CalFlashIMPEC(); }
    /// Reset IMPEC
    void ResetVal01IMPEC();
    /// Reset IMPEC
    void ResetVal02IMPEC();
    /// Reset IMPEC
    void ResetVal03IMPEC();
    /// Update value of last step for IMPEC
    void UpdateLastStepIMPEC();


    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for FIM.
    void AllocateFIM_IsoT();
    /// Initialize the properties of Reservoir for FIM
    void InitFIM();
    /// Allocate maximal memory for internal Matrix for FIM
    void AllocateMatFIM_IsoT(LinearSystem& myLS) const;
    /// Assemble Matrix for FIM
    void AssembleMatFIM(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Get the Solution -- Reservoir Pressure and moles of Components for FIM
    void GetSolutionFIM(const vector<OCP_DBL>& u, const OCP_DBL& dPmax, const OCP_DBL& dSmax);
    /// Calculate Flash for FIM, some derivatives are needed
    void CalFlashFIM() { bulk.CalFlashFIM(); }
    /// Calculate Relative Permeability and Capillary and some derivatives for each Bulk
    void CalKrPcFIM() { bulk.CalKrPcFIM(); }
    /// Calculate the Residual for FIM, it's also RHS of Linear System
    void CalResFIM(OCPRes& resFIM, const OCP_DBL& dt, const OCP_BOOL& resetRes0);
    /// Reset FIM
    void ResetFIM();
    /// Update value of last step for FIM.
    void UpdateLastStepFIM();

    /////////////////////////////////////////////////////////////////////
    // FIMn
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used for FIMn
    void AllocateFIMn_IsoT();
    /// Initialize the properties of Reservoir for FIMn
    void InitFIMn();
    void AssembleMatFIMn(LinearSystem& myLS, const OCP_DBL& dt) const;

    void GetSolutionFIMn(const vector<OCP_DBL>& u,
        const OCP_DBL& dPmax,
        const OCP_DBL& dSmax);
    void CalFlashFIMn() { bulk.CalFlashFIMn(); }
    void ResetFIMn();
    void UpdateLastStepFIMn();

    /////////////////////////////////////////////////////////////////////
    // AIMc
    /////////////////////////////////////////////////////////////////////

public:

    /// Allocate memory for variables used for AIMc
    void AllocateAIMc_IsoT();
    /// Initialize the properties of Reservoir for AIMc
    void InitAIMc();
    /// Assemble Matrix for AIMc
    void AssembleMatAIMc(LinearSystem& myLS, const OCP_DBL& dt) const;
    /// Perform flash calculation with Ni for Explicit bulk -- Update partial properties
    void CalFlashAIMcEp() { bulk.CalFlashAIMcEp(); }
    /// Perform flash calculation with Ni for Explicit bulk -- Update all properties
    void CalFlashAIMcEa() { bulk.CalFlashAIMcEa(); }
    /// Perform flash calculation with Ni for Implicit bulk
    void CalFlashAIMcI() { bulk.CalFlashAIMcI(); }
    /// Calculate relative permeability and capillary pressure for Explicit bulk
    void CalKrPcAIMcE() { bulk.CalKrPcAIMcE(); }
    /// Calculate relative permeability and capillary pressure for Implicit bulk
    void CalKrPcAIMcI() { bulk.CalKrPcAIMcI(); }
    /// Calculate the Residual for AIMc, it's also RHS of Linear System
    void CalResAIMc(OCPRes& resAIMc, const OCP_DBL& dt, const OCP_BOOL& resetRes0);
    /// Get the solution for AIMc after a Newton iteration.
    void GetSolutionAIMc(const vector<OCP_DBL>& u, const OCP_DBL& dPmax, const OCP_DBL& dSmax);
    /// Reset AIMc
    void ResetAIMc();
    /// Update values of last step for AIMc
    void UpdateLastStepAIMc();
    /// Setup FIMBulk
    void SetupFIMBulk(const OCP_BOOL& NRflag = OCP_FALSE) { conn.SetupFIMBulk(bulk, NRflag); }
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