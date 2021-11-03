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
#include "Grid.hpp"
#include "LinearSolver.hpp"
#include "OCPStructure.hpp"
#include "DenseMat.hpp"

using namespace std;

/// Connection between indices of two bulks (BId, EId); usually, BId > EId.
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

    OCP_USI GetBId()const { return BId; }
    OCP_USI GetEId()const { return EId; }

private:
    OCP_USI BId;
    OCP_USI EId;
};

/// BulkConn is a core component of reservoir, it contains all properties and
/// operations about connections between bulks(active grids). due to the activity of
/// bulks, almost all connections are meaningful. you can traverse all the connections
/// through the iterator in it, which is effective. flow calculation between active
/// bulks, matrix assembling contributed only by bulks are included in it.
class BulkConn
{
    friend class LinearSolver;


public:
    BulkConn() = default;


    /////////////////////////////////////////////////////////////////////
    // General
    /////////////////////////////////////////////////////////////////////

public:
    /// Setup active connections and calculate necessary property from Grid and Bulk.
    /// It should be called after Grid and Bulk Setup.
    void Setup(const Grid& myGrid, const Bulk& myBulk);
    /// Initialize the size of variable related to neighbor.
    void InitSize(const Bulk& myBulk);
    /// Setup variable related to neighbor.
    void CalConn(const Grid& myGrid, const USI& np);
    /// Generate iteratorConn of active connections from neighbor.
    void CalIteratorConn();
    /// Calculate all effective area of connections.
    void CalArea(const Grid& myGrid, const Bulk& myBulk);
    /// Calculate effective area of connections.
    OCP_DBL CalAkd(const Grid& myGrid, const Bulk& myBulk, const OCP_USI& bIdb,
        const OCP_USI& eIdb) const;
    /// Allocate memory for Matrix, it should be called only once at the beginning.
    void AllocateMat(LinearSolver& mySolver) const;
    /// Setup sparsity pattern of Matrix begin assembling Matrix.
    void SetupMatSparsity(LinearSolver& mySolver) const;
    /// Update value of last step
    void UpdateLastStep();
    /// Reset current step to last step
    void Reset();
    /// Check the difference from last step -- test
    void CheckDiff() const;
    /// Return num of bulk.
    OCP_USI GetBulkNum() const { return numBulk; }
    /// Print information about connection on screen.
    void GetConnectionInfo() const;

private:

    OCP_USI numBulk; ///< Num of bulks (active grids).
    OCP_USI numConn; ///< Num of connections between Bulks.

    vector<vector<OCP_USI>>
        neighbor; ///< The ith row stores the ith bulk's neighbor, which is sort in
                  ///< increasing order: activeGridNum.
    vector<USI> selfPtr;     ///< The ith row stores the location of the ith bulk in
                             ///< neighbor[i]: activeGridNum.
    vector<USI> neighborNum; ///< The ith row stores num of neighbor of the ith bulk:
                             ///< activeGridNum.
    /// Contains all the connections, in which the index of first bulk is greater than
    /// the ones of second bulk. the iteratorConn is generated from neighbor: numConn.
    vector<BulkPair> iteratorConn;
    vector<OCP_DBL> area; ///< Effective area for each connections, which are ordered
                          ///< the same as iteratorConn: numConn.
    /// Upblock of connections.
    /// Upblock is identified by difference of pressure between phases: numConn * nums
    /// of phase.
    vector<OCP_USI> upblock;
    vector<OCP_DBL>
        upblock_Rho; ///< Mass density of phase from upblock: numConn * nums of phase.
    vector<OCP_DBL> upblock_Trans; ///< Transmissibility of phase from upblock: numConn
                                   ///< * nums of phase.
    vector<OCP_DBL> upblock_Velocity; ///< Flow rate of volume of phase from upblock:
                                      ///< numConn * nums of phase.
    // For last time step
    vector<OCP_USI> lastUpblock;
    vector<OCP_DBL> lastUpblock_Rho;
    vector<OCP_DBL> lastUpblock_Trans;
    vector<OCP_DBL> lastUpblock_Velocity;



    /////////////////////////////////////////////////////////////////////
    // IMPEC
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for IMPEC
    void AllocateAuxIMPEC(const USI& np);
    /// Assmeble Matrix for IMPEC, parts only related to bulks are considered.
    void AssembleMatIMPEC(LinearSolver& mySolver, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    /// calculate the CFL number of flow between bulks.
    OCP_DBL CalCFLIMPEC(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate the CFL number of flow between bulks.
    void CalCFL01IMPEC(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate main information about flow between bulks for IMPEC
    void CalFluxIMPEC(const Bulk& myBulk);
    /// Update moles of component in each bulk according to mass conserve equations at
    /// current timestep for IMPEC
    void MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const;



    /////////////////////////////////////////////////////////////////////
    // FIM
    /////////////////////////////////////////////////////////////////////

public:
    /// Allocate memory for auxiliary variables used for FIM
    void AllocateAuxFIM(const USI& np);
    /// Assmeble Matrix for FIM, parts only related to bulks are considered.
    void AssembleMat_FIM(LinearSolver& mySolver, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    /// calculate Upblock for FIM
    void CalFluxFIM(const Bulk& myBulk);
    /// Calculate Resiual for FIM
    void CalResFIM(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);

};


#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/