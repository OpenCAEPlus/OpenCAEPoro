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

/// Connection between two bulks (BId, EId); usually, indices BId > EId.
//  Note: Bulks are the active grid cells.
class BulkPair
{
    friend class BulkConn;

public:
    /// Default constructor.
    BulkPair() = default;

    /// Setup BulkPair with bId and eId.
    BulkPair(const OCP_USI& bId, const OCP_USI& eId)
        : BId(bId)
        , EId(eId){};

    OCP_USI GetBId() const { return BId; } ///< Return beginning index.
    OCP_USI GetEId() const { return EId; } ///< Return ending index.

private:
    OCP_USI BId; ///< Beginning index of a pair.
    OCP_USI EId; ///< Ending index of a pair.
};

/// Properties and operations on connections between bulks (active grids).
//  Note: BulkConn is a core component of reservoir, it contains all properties and
//  operations on connections between bulks (active grids). You can traverse all the
//  connections through an effective iterator. Flow calculations between active
// bulks and matrix assembling with contributions from bulks only are included.
class BulkConn
{

public:
    BulkConn() = default; ///< Default constructor.

    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:
    /// Setup active connections and calculate necessary properties using Grid and Bulk.
    void Setup(const Grid& myGrid, const Bulk& myBulk);

    /// Generate iteratorConn of active connections from neighbor.
    void CalIteratorConn();

    /// Allocate memory for the coefficient matrix.
    void AllocateMat(LinearSystem& myLS) const;

    /// Setup sparsity pattern of the coefficient matrix.
    void SetupMatSparsity(LinearSystem& myLS) const;

    /// Update physcial values of the previous step.
    void UpdateLastStep();

    /// Reset physcial values of the current step with the previous step.
    void Reset();

    /// Check differences between the current and previous steps.
    void CheckDiff() const;

    /// Return number of bulks.
    OCP_USI GetBulkNum() const { return numBulk; }

    /// Print information of connections on screen.
    void PrintConnectionInfo() const;

private:
    OCP_USI numBulk; ///< Number of bulks (active grid cells).
    OCP_USI numConn; ///< Number of connections between bulks.

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /// Self-pointer, the indiecs of the i-th bulk in neighbor[i]: activeGridNum???
    vector<USI> selfPtr;

    /// Number of neighbors of the i-th bulk: activeGridNum.
    vector<USI> neighborNum;

    /// All connections (pair of indices) between bulks: numConn.
    //  Note: In each pair, the index of first bulk is greater than the second. The data
    //  in iteratorConn is generated from neighbor.
    vector<BulkPair> iteratorConn;

    /// Effective area of each connection ordered as iteratorConn: numConn.
    vector<OCP_DBL> area;

    ////// Current Time Step
    /// Index of upwinding bulk of connections: numConn * nums of phase.
    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI> upblock;
    /// Mass density of phase from upblock: numConn * nums of phase.
    vector<OCP_DBL> upblock_Rho;
    /// Transmissibility of phase from upblock: numConn * nums of phase.
    vector<OCP_DBL> upblock_Trans;
    /// Flow volume rate of phase from upblock: numConn * nums of phase.
    vector<OCP_DBL> upblock_Velocity;

    ////// Previous Time Step
    vector<OCP_USI> lastUpblock;          ///< Upwinding cell at the last time step
    vector<OCP_DBL> lastUpblock_Rho;      ///< Density at the last time step
    vector<OCP_DBL> lastUpblock_Trans;    ///< Transmisbility at the last time step
    vector<OCP_DBL> lastUpblock_Velocity; ///< Velocity at the last time step

    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used by the IMPEC method.
    void AllocateAuxIMPEC(const USI& np);

    /// Assmeble coefficient matrix for IMPEC, terms related to bulks only.
    void AssembleMatIMPEC(LinearSystem& myLS, const Bulk& myBulk,
                          const OCP_DBL& dt) const;

    /// Calculate the CFL number for flow between bulks.
    OCP_DBL CalCFLIMPEC(const Bulk& myBulk, const OCP_DBL& dt) const;

    /// Calculate the CFL number for flow between bulks???
    void CalCFL01IMPEC(const Bulk& myBulk, const OCP_DBL& dt) const;

    /// Calculate flux information about flow between bulks for IMPEC.
    void CalFluxIMPEC(const Bulk& myBulk);

    /// Update mole composition of each bulk according to mass conservation for IMPEC.
    void MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const;

    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used by the FIM method.
    void AllocateAuxFIM(const USI& np);

    /// Assmeble coefficient matrix for FIM, terms related to bulks only.
    void AssembleMat_FIM(LinearSystem& myLS, const Bulk& myBulk,
                         const OCP_DBL& dt) const;

    /// Calculate flux for FIM, considering upwinding.
    void CalFluxFIM(const Bulk& myBulk);

    /// Calculate resiual for the Newton iteration in FIM.
    void CalResFIM(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);

    /// Reset upblock for FIM.
    void ResetUpblockFIM() { upblock = lastUpblock; }

    /// Update lastUpblock for FIM.
    void UpdateLastUpblockFIM() { lastUpblock = upblock; }
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