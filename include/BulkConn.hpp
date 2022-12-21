/*! \file    BulkConn.hpp
 *  \brief   BulkConn class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONN_HEADER__
#define __BULKCONN_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "Bulk.hpp"
#include "DenseMat.hpp"
#include "Grid.hpp"
#include "LinearSystem.hpp"
#include "OCPStructure.hpp"

using namespace std;

/// Connection between two bulks (bId, eId); usually, indices bId > eId.
//  Note: Bulks are the active grid cells.
class BulkPair
{
    friend class BulkConn;

public:
    /// Default constructor.
    BulkPair() = default;

    /// Setup BulkPair with bId and eId.
    BulkPair(const OCP_USI& BId, const OCP_USI& EId, const USI& direct,
        const OCP_DBL& AreaB, const OCP_DBL& AreaE)
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE){};

    OCP_USI GetBId() const { return bId; } ///< Return beginning index.
    OCP_USI GetEId() const { return eId; } ///< Return ending index.

protected:
    OCP_USI bId;        ///< Beginning index of a pair.
    OCP_USI eId;        ///< Ending index of a pair.
    OCP_DBL area;       ///< Effective area
    USI     direction;  ///< 1-x, 2-y, 3-z
    USI     conntype;
    OCP_DBL areaB;      ///< Area of intersecting faces from Begin grid
    OCP_DBL areaE;      ///< Area of intersecting faces from End grid
};

/// Properties and operations on connections between bulks (active grids).
//  Note: BulkConn is a core component of reservoir, it contains all properties and
//  operations on connections between bulks (active grids). You can traverse all the
//  connections through an effective iterator. Flow calculations between active
// bulks and matrix assembling with contributions from bulks only are included.
class BulkConn
{
    // temp
    friend class MyMetisTest;
    friend class Out4VTK;

public:

    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:

    /// Print information of connections on screen.
    void PrintConnectionInfo(const Grid& myGrid) const;
    void PrintConnectionInfoCoor(const Grid& myGrid) const;


    /////////////////////////////////////////////////////////////////////
    // General Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Setup active connections
    void Setup(const GridInitInfo& initGrid);
    /// Setup sparsity pattern of the coefficient matrix.
    void SetupMatSparsity(LinearSystem& myLS) const;
    /// Allocate memory for the coefficient matrix.
    void AllocateMat(LinearSystem& myLS) const;
    /// Calculate the effective area used for flow
    void CalAkd(const Bulk& myBulk);
    /// Return number of bulks.
    OCP_USI GetBulkNum() const { return numBulk; }

protected:

    OCP_USI numBulk; ///< Number of bulks (active grid cells).
    OCP_USI numConn; ///< Number of connections between bulks.

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /// Self-pointer, the indiecs of the i-th bulk in neighbor[i]: numBulk
    vector<USI> selfPtr;

    /// Number of neighbors of the i-th bulk: numBulk, self-included
    vector<USI> neighborNum;

    /// All connections (pair of indices) between bulks: numConn.
    //  Note: In each pair, the index of first bulk is greater than the second. The data
    //  in iteratorConn is generated from neighbor.
    vector<BulkPair> iteratorConn;


    /////////////////////////////////////////////////////////////////////
    // Physical Variables
    /////////////////////////////////////////////////////////////////////

public:
    /// Calculate the CFL number for flow between bulks???
    void CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const;
protected:

    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI> upblock;          ///< Index of upwinding bulk of connections : numConn * numPhase.
    vector<OCP_DBL> upblock_Rho;      ///< Mass density of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Trans;    ///< Transmissibility of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Velocity; ///< Volume flow rate of phase from upblock: numConn * numsPhase.
    vector<OCP_DBL> Adkt;             ///< Thermal conductivity between neighbors

    // Last time step
    vector<OCP_USI> lupblock;          ///< last upblock
    vector<OCP_DBL> lupblock_Rho;      ///< last upblock_Rho
    vector<OCP_DBL> lupblock_Trans;    ///< last upblock_Trans
    vector<OCP_DBL> lupblock_Velocity; ///< last upblock_Velocity
    vector<OCP_DBL> lAdkt;             ///< last Adkt

    // Derivatives 
    vector<OCP_DBL> AdktP;             ///< d Adkt / d P, oreder: connections -> bId.P -> eId.P
    vector<OCP_DBL> AdktT;             ///< d Adkt / d T, oreder: connections -> bId.T -> eId.T
    vector<OCP_DBL> AdktS;             ///< d Adkt / d S, oreder: connections -> bId.phase -> eId.phase

    // Last time step
    vector<OCP_DBL> lAdktP;            ///< last AdktP
    vector<OCP_DBL> lAdktT;            ///< last AdktT
    vector<OCP_DBL> lAdktS;            ///< last AdktS

    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for variables used by the IMPEC method.
    void AllocateIMPEC_IsoT(const USI& np);
    /// Assmeble coefficient matrix for IMPEC, terms related to bulks only.
    void AssembleMatIMPEC(LinearSystem& myLS, const Bulk& myBulk,
                          const OCP_DBL& dt) const;
    /// Calculate flux information about flow between bulks for IMPEC.
    void CalFluxIMPEC(const Bulk& myBulk);
    /// Update mole composition of each bulk according to mass conservation for IMPEC.
    void MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const;
    /// Reset variables needed for IMPEC
    void ResetIMPEC();
    /// Update value of last step for IMPEC.
    void UpdateLastStepIMPEC();

    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used by the FIM method.
    void AllocateFIM_IsoT(const USI& np);

    /// Assmeble coefficient matrix for FIM, terms related to bulks only.
    void AssembleMat_FIM(LinearSystem& myLS, const Bulk& myBulk,
                         const OCP_DBL& dt) const;

    /// Calculate resiual for the Newton iteration in FIM.
    void CalResFIM(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);
    /// rho = (S1*rho1 + S2*rho2)/(S1+S2)
    void CalFluxFIMS(const Grid& myGrid, const Bulk& myBulk);
    void CalResFIMS(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);

    /////////////////////////////////////////////////////////////////////
    // FIM(new) assemble
    /////////////////////////////////////////////////////////////////////

    /// Assmeble coefficient matrix for FIM, terms related to bulks only.
    /// OCP_NEW_FIM
    void AssembleMat_FIM_new(LinearSystem& myLS, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    void AssembleMat_FIM_new1(LinearSystem& myLS, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    /// OCP_NEW_FIM rho = (S1*rho1 + S2*rho2)/(S1+S2)
    void AssembleMat_FIM_newS(LinearSystem& myLS, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    /// OCP_NEW_FIMn
    void AssembleMat_FIM_new_n(LinearSystem& myLS, const Bulk& myBulk,
        const OCP_DBL& dt) const;

    /////////////////////////////////////////////////////////////////////
    // AIMc
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used by the AIMc method.
    void AllocateAIMc_IsoT(const USI& np);
    void AssembleMat_AIMc(LinearSystem& myLS, const Bulk& myBulk, const OCP_DBL& dt) const;
    /// Calculate resiual for the Newton iteration in FIM.
    void CalResAIMc(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);
    /// Determine which bulk are treated Implicit
    void SetupFIMBulk(Bulk& myBulk, const OCP_BOOL& NRflag = OCP_FALSE) const;
    /// Setup k-neighbor for well bulks
    void SetupWellBulk_K(Bulk& myBulk) const;
};


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/17/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/