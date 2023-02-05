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
    BulkPair(const OCP_USI& BId,
             const OCP_USI& EId,
             const USI&     direct,
             const OCP_DBL& AreaB,
             const OCP_DBL& AreaE)
        : bId(BId)
        , eId(EId)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE){};

    OCP_USI BId() const { return bId; }   ///< Return beginning index.
    OCP_USI EId() const { return eId; }   ///< Return ending index.
    OCP_DBL Area() const { return area; } ///< Return effective area
    OCP_DBL AreaB() const { return areaB; }
    OCP_DBL AreaE() const { return areaE; }
    USI     Direction() const { return direction; }

protected:
    OCP_USI bId;       ///< Beginning index of a pair.
    OCP_USI eId;       ///< Ending index of a pair.
    OCP_DBL area;      ///< Effective area
    USI     direction; ///< 1-x, 2-y, 3-z
    OCP_DBL areaB;     ///< Area of intersecting faces from Begin grid
    OCP_DBL areaE;     ///< Area of intersecting faces from End grid
};

/// Properties and operations on connections between bulks (active grids).
//  Note: BulkConn is a core component of reservoir, it contains all properties and
//  operations on connections between bulks (active grids). You can traverse all the
//  connections through an effective iterator. Flow calculations between active
// bulks and matrix assembling with contributions from bulks only are included.
class BulkConn
{
    friend class Reservoir;
    // temp
    friend class MyMetisTest;
    friend class Out4VTK;
    friend class IsoT_FIM;
    friend class IsoT_IMPEC;
    friend class IsoT_AIMc;
    friend class IsoT_FIMn;
    friend class T_FIM;

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
    void SetupIsoT(const Grid& myGrid, const Bulk& myBulk);
    void Setup(const Grid& myGrid);
    /// Calculate the effective area used for flow
    void CalAkd(const Bulk& myBulk);
    /// Return number of bulks.
    OCP_USI     GetBulkNum() const { return numBulk; }
    const auto& GetNeighborNum() const { return neighborNum; }

protected:
    OCP_USI numBulk; ///< Number of bulks (active grid cells).
    OCP_USI numConn; ///< Number of connections between bulks.

    /// Neighboring information of each bulk: activeGridNum.
    //  Note: The i-th entry stores the i-th bulk's neighbors, which is sorted in an
    //  increasing order.
    vector<vector<OCP_USI>> neighbor;

    /// Self-pointer, the indices of the i-th bulk in neighbor[i]: numBulk
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

protected:
    //  Note: Upblock is identified by difference of pressure between phases.
    vector<OCP_USI>
        upblock; ///< Index of upwinding bulk of connections : numConn * numPhase.
    vector<OCP_DBL>
        upblock_Rho; ///< Mass density of phase from upblock: numConn * numPhase.
    vector<OCP_DBL>
        upblock_Trans; ///< Transmissibility of phase from upblock: numConn * numPhase.
    vector<OCP_DBL> upblock_Velocity; ///< Volume flow rate of phase from upblock:
                                      ///< numConn * numsPhase.
    vector<OCP_DBL> Adkt;             ///< Thermal conductivity between neighbors

    // Last time step
    vector<OCP_USI> lupblock;          ///< last upblock
    vector<OCP_DBL> lupblock_Rho;      ///< last upblock_Rho
    vector<OCP_DBL> lupblock_Trans;    ///< last upblock_Trans
    vector<OCP_DBL> lupblock_Velocity; ///< last upblock_Velocity
    vector<OCP_DBL> lAdkt;             ///< last Adkt

    // Derivatives
    vector<OCP_DBL> AdktP; ///< d Adkt / d P, order: connections -> bId.P -> eId.P
    vector<OCP_DBL> AdktT; ///< d Adkt / d T, order: connections -> bId.T -> eId.T
    vector<OCP_DBL>
        AdktS; ///< d Adkt / d S, order: connections -> bId.phase -> eId.phase

    // Last time step
    vector<OCP_DBL> lAdktP; ///< last AdktP
    vector<OCP_DBL> lAdktT; ///< last AdktT
    vector<OCP_DBL> lAdktS; ///< last AdktS

public:
    /// rho = (S1*rho1 + S2*rho2)/(S1+S2)
    void CalFluxFIMS(const Grid& myGrid, const Bulk& myBulk);
    void CalResFIMS(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt);
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