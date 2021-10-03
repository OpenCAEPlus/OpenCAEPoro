#pragma once
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
	BB_Pair(int bId, int eId) : BId(bId), EId(eId) {};

private:
	int			BId;
	int			EId;
};


class Connection_BB
{
	friend class Solver<OCP_DBL>;

public:

	Connection_BB() = default;

	// info
	void getConnectionInfo();
	int getActiveBulkNum() { return ActiveBulkNum; }

	// Active Conn & Active Bulk
	// init
	void setup(const Grid& myGrid, const Bulk& myBulk);
	void initSize(const Bulk& myBulk);
	void initActive(const Grid& myGrid, int np);
	void getIteratorActive();
	void calAreaActive(const Grid& myGrid, const Bulk& myBulk);
	OCP_DBL calAkd(const Grid&myGrid, const Bulk& myBulk, int bIdb, int eIdb);

	OCP_DBL calCFL(Bulk& myBulk, OCP_DBL dt);
	void calFlux(const Bulk& myBulk);
	void massConserve(Bulk& myBulk, OCP_DBL dt);

	// Assemble Mat
	template<typename T>
	void allocateMat(Solver<T>& mySolver) const;
	template<typename T>
	void initAssembleMat(Solver<T>& mySolver) const;

	void assembleMat(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt) const;


private:
	// Bulk to Bulk
	int								ActiveBulkNum;
	int								ActiveConnNum;


	std::vector<vector<int>>		Neighbor;
	std::vector<int>				SelfPtr;				// ptr for self in every row of Neighbor
	std::vector<int>				NeighborNum;
	std::vector<BB_Pair>			Iterator;
	std::vector<OCP_DBL>				Area;					// effective area for each CONN
	std::vector<int>				Upblock;				// upblock of connection
	std::vector<OCP_DBL>				Upblock_Rho;			// rhoj in flux
	std::vector<OCP_DBL>				Upblock_Trans;
	std::vector<OCP_DBL>				Upblock_Velocity;       // volume rate
};


template<typename T>
void Connection_BB::allocateMat(Solver<T>& MySolver)  const
{
	for (int n = 0; n < ActiveBulkNum; n++) {
		MySolver.RowCapacity[n] += NeighborNum[n];
	}
}

template<typename T>
void Connection_BB::initAssembleMat(Solver<T>& mySolver) const
{
	mySolver.Dim = ActiveBulkNum;
	for (int n = 0; n < ActiveBulkNum; n++) {
		mySolver.ColId[n].assign(Neighbor[n].begin(), Neighbor[n].end());
		mySolver.DiagPtr[n] = SelfPtr[n];
	}
}