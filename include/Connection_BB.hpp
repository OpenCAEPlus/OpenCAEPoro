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
	friend class Solver<double>;

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
	double calAkd(const Grid&myGrid, const Bulk& myBulk, int bIdb, int eIdb);

	void calFlux(const Bulk& myBulk);
	void massConserve(Bulk& myBulk);

	// Assemble Mat
	template<typename T>
	void allocateMat(Solver<T>& mySolver);
	template<typename T>
	void initAssembleMat(Solver<T>& mySolver);

	void assembleMat(Solver<double>& mySolver, const Bulk& myBulk, double dt);


private:
	// Bulk to Bulk
	int								ActiveBulkNum;
	int								ActiveConnNum;


	std::vector<vector<int>>		Neighbor;
	std::vector<int>				SelfPtr;				// ptr for self in every row of Neighbor
	std::vector<int>				NeighborNum;
	std::vector<BB_Pair>			Iterator;
	std::vector<double>				Area;					// effective area for each CONN
	std::vector<int>				Upblock;				// upblock of connection
	std::vector<double>				Upblock_Rho;			// rhoj in flux
	std::vector<double>				Upblock_Trans;
	std::vector<double>				Upblock_Velocity;
};


template<typename T>
void Connection_BB::allocateMat(Solver<T>& MySolver)
{
	for (int n = 0; n < ActiveBulkNum; n++) {
		MySolver.RowCapacity[n] += NeighborNum[n];
	}
}

template<typename T>
void Connection_BB::initAssembleMat(Solver<T>& mySolver)
{
	mySolver.Dim = ActiveBulkNum;
	for (int n = 0; n < ActiveBulkNum; n++) {
		mySolver.ColId[n].assign(Neighbor[n].begin(), Neighbor[n].end());
		mySolver.DiagPtr[n] = SelfPtr[n];
	}
}