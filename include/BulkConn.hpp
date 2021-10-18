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

#ifndef __CONNECTION_BB_HEADER__
#define __CONNECTION_BB_HEADER__

// Standard header files
#include <vector>

// OpenCAEPoro header files
#include "Bulk.hpp"
#include "Grid.hpp"
#include "LinearSolver.hpp"

using namespace std;

/// Connection between indices of two bulks (BId, EId); usually, BId > EId.
class BB_Pair
{
    friend class Connection_BB;

public:
    /// Default constructor.
    BB_Pair() = default;

    /// Setup BB_Pair with bId and eId.
    BB_Pair(const OCP_USI& bId, const OCP_USI& eId)
        : BId(bId)
        , EId(eId){};

private:
    OCP_USI BId;
    OCP_USI EId;
};

/// Connection_BB is a core component of reservoir, it contains all properties and
/// operations about connections between bulks(active grids). due to the activity of
/// bulks, almost all connections are meaningful. you can traverse all the connections
/// through the iterator in it, which is effective. flow calculation between active
/// bulks, matrix assembling contributed only by bulks are included in it.
class Connection_BB
{
    friend class Solver<OCP_DBL>;

public:
    Connection_BB() = default;

    /// Print information about connection on screen.
    void GetConnectionInfo() const;
    /// Return num of bulk.
    OCP_USI GetBulkNum() const { return numBulk; }

    /// Setup active connections and calculate necessary property from Grid and Bulk.
    /// It should be called after Grid and Bulk Setup.
    void Setup(const Grid& myGrid, const Bulk& myBulk);
    /// initialize the size of variable related to neighbor.
    void InitSize(const Bulk& myBulk);
    /// Setup variable related to neighbor.
    void CalConn(const Grid& myGrid, const USI& np);
    /// generate iteratorConn of active connections from neighbor.
    void CalIteratorConn();
    /// calculate all effective area of connections.
    void CalArea(const Grid& myGrid, const Bulk& myBulk);
    /// calculate effective area of connections.
    OCP_DBL CalAkd(const Grid& myGrid, const Bulk& myBulk, const OCP_USI& bIdb,
                   const OCP_USI& eIdb) const;
    /// calculate the CFL number of flow between bulks.
    OCP_DBL CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate the CFL number of flow between bulks.
    void CalCFL01(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate main information about flow between bulks.
    void CalFlux(const Bulk& myBulk);
    /// update moles of component in each bulk according to mass conserve equations at
    /// current timestep.
    void MassConserve(Bulk& myBulk, const OCP_DBL& dt) const;

    // Assemble Mat
    /// Allocate memory for Matrix, it should be called only once at the beginning.
    template <typename T> void AllocateMat(Solver<T>& mySolver) const;
    /// Setup sparsity pattern of Matrix, it should be called before every time the
    /// linear system setups. actually, part from wells is neglect, which is much less
    /// than bulks.
    template <typename T> void InitAssembleMat(Solver<T>& mySolver) const;
    /// assmeble Matrix, parts only related to bulks are considered.
    void AssembleMat_IMPES(Solver<OCP_DBL>& mySolver, const Bulk& myBulk,
                           const OCP_DBL& dt) const;

    void SetLastStep();
    void Reset();

private:
    // Bulk to Bulk
    OCP_USI numBulk; ///< num of bulks (active grids).
    OCP_USI numConn; ///< num of connections between bulks.

    vector<vector<OCP_USI>>
        neighbor; ///< the ith row stores the ith bulk's neighbor, which is sort in
                  ///< increasing order: activeGridNum.
    vector<USI> selfPtr;     ///< the ith row stores the location of the ith bulk in
                             ///< neighbor[i]: activeGridNum.
    vector<USI> neighborNum; ///< the ith row stores num of neighbor of the ith bulk:
                             ///< activeGridNum.
    /// contains all the connections, in which the index of first bulk is greater than
    /// the ones of second bulk. the iteratorConn is generated from neighbor: numConn.
    vector<BB_Pair> iteratorConn;
    vector<OCP_DBL> area; ///< effective area for each connections, which are ordered
                          ///< the same as iteratorConn: numConn.
    /// upblock of connections.
    /// upblock is identified by difference of pressure between phases: numConn * nums
    /// of phase.
    vector<OCP_USI> upblock;
    vector<OCP_DBL>
        upblock_Rho; ///< mass density of phase from upblock: numConn * nums of phase.
    vector<OCP_DBL> upblock_Trans; ///< transmissibility of phase from upblock: numConn
                                   ///< * nums of phase.
    vector<OCP_DBL> upblock_Velocity; ///< flow rate of volume of phase from upblock:
                                      ///< numConn * nums of phase.
    // For last time step
    vector<OCP_USI> lastUpblock;
    vector<OCP_DBL> lastUpblock_Rho;
    vector<OCP_DBL> lastUpblock_Trans;
    vector<OCP_DBL> lastUpblock_Velocity;
                                      
};

template <typename T> void Connection_BB::AllocateMat(Solver<T>& MySolver) const
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        MySolver.rowCapacity[n] += neighborNum[n];
    }
}

template <typename T> void Connection_BB::InitAssembleMat(Solver<T>& mySolver) const
{
    mySolver.dim = numBulk;
    for (OCP_USI n = 0; n < numBulk; n++) {
        mySolver.colId[n].assign(neighbor[n].begin(), neighbor[n].end());
        mySolver.diagPtr[n] = selfPtr[n];
    }
}

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/