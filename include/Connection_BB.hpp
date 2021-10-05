/*! \file    Connection_BB.hpp
 *  \brief   Connection_BB class declaration
 *  \author  Shizhe Li
 *  \date    Oct/05/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __CONNECTION_BB_HEADER__
#define __CONNECTION_BB_HEADER__


#include <vector>
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Solver.hxx"


using namespace std;

class BB_Pair
{
	friend class Connection_BB;
public:
	BB_Pair() = default;
	BB_Pair(OCP_USI bId, OCP_USI eId) : BId(bId), EId(eId) {};

private:
	OCP_USI			BId;
	OCP_USI			EId;
};

/// contains information about connectivity between bulks (active grids).
class Connection_BB
{
	friend class Solver<OCP_DBL>;

public:

	Connection_BB() = default;

	/// print information about connection on screen.
	void getConnectionInfo();
	/// return ActiveBulkNum.
	OCP_USI getActiveBulkNum() { return ActiveBulkNum; }

	/// setup active connections and calculate necessary property from Grid and Bulk.
	/// it should be called after Grid and Bulk setup.
	void setup(const Grid& myGrid, const Bulk& myBulk);
	/// initialize the size of variable related to Neighbor.
	void initSize(const Bulk& myBulk);
	/// setup variable related to Neighbor.
	void initActive(const Grid& myGrid, USI np);
	/// generate Iterator of active connections from Neighbor.
	void getIteratorActive();
	/// calculate all effective area of active connections.
	void calAreaActive(const Grid& myGrid, const Bulk& myBulk);
	/// calculate effective area of active connections.
	OCP_DBL calAkd(const Grid&myGrid, const Bulk& myBulk, OCP_USI bIdb, OCP_USI eIdb);
	/// calculate the CFL number of flow between bulks.
	OCP_DBL calCFL(Bulk& myBulk, OCP_DBL dt);
	/// calculate main information about flow between bulks.
	void calFlux(const Bulk& myBulk);
	/// update moles of component in each bulk according to mass conserve equations at current timestep.
	void massConserve(Bulk& myBulk, OCP_DBL dt);

	// Assemble Mat
	/// allocate memory for Matrix, it should be called only once at the beginning.
	template<typename T>
	void allocateMat(Solver<T>& mySolver) const;
	/// setup sparsity pattern of Matrix, it should be called before every time the linear system setups.
	/// actually, part from wells is neglect, which is much less than bulks.
	template<typename T>
	void initAssembleMat(Solver<T>& mySolver) const;
	/// assmeble Matrix, parts only related to bulks are considered.
	void assembleMat_IMPES(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt) const;


private:
	// Bulk to Bulk
	OCP_USI							ActiveBulkNum;	///< num of bulks (active grids).
	OCP_USI							ActiveConnNum;	///< num of connections between bulks.

	std::vector<vector<OCP_USI>>	Neighbor;		///< the ith row stores the ith bulk's neighbor, which is sort in increasing order: ActiveBulkNum.
	std::vector<USI>				SelfPtr;		///< the ith row stores the location of the ith bulk in Neighbor[i]: ActiveBulkNum.
	std::vector<USI>				NeighborNum;	///< the ith row stores num of neighbor of the ith bulk: ActiveBulkNum.
	/// contains all the connections, in which the index of first bulk is greater than the ones of second bulk.
	/// the Iterator is generated from Neighbor: ActiveConnNum.
	std::vector<BB_Pair>			Iterator;		
	std::vector<OCP_DBL>			Area;			///< effective area for each connections, which are ordered the same as Iterator: ActiveConnNum.
	/// upblock of connections.
	/// upblock is identified by difference of pressure between phases: ActiveConnNum * nums of phase.
	std::vector<OCP_USI>			Upblock;
	std::vector<OCP_DBL>			Upblock_Rho;	///< mass density of phase from upblock: ActiveConnNum * nums of phase.
	std::vector<OCP_DBL>			Upblock_Trans;	///< transmissibility of phase from upblock: ActiveConnNum * nums of phase.
	std::vector<OCP_DBL>			Upblock_Velocity;   ///< flow rate of volume of phase from upblock: ActiveConnNum * nums of phase.
};


template<typename T>
void Connection_BB::allocateMat(Solver<T>& MySolver)  const
{
	for (OCP_USI n = 0; n < ActiveBulkNum; n++) {
		MySolver.RowCapacity[n] += NeighborNum[n];
	}
}

template<typename T>
void Connection_BB::initAssembleMat(Solver<T>& mySolver) const
{
	mySolver.Dim = ActiveBulkNum;
	for (OCP_USI n = 0; n < ActiveBulkNum; n++) {
		mySolver.ColId[n].assign(Neighbor[n].begin(), Neighbor[n].end());
		mySolver.DiagPtr[n] = SelfPtr[n];
	}
}

#endif