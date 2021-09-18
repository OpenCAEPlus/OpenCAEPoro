#pragma once
#include <vector>
#include "Grid.hpp"
#include "Bulk.hpp"
#include "Solver.hpp"


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
	friend class Solver;

public:

	Connection_BB() = default;

	// info
	void getConnectionInfo();
	int getActiveBulkNum() { return ActiveBulkNum; }

	// Active Conn & Active Bulk
	// init
	void setup(const Grid& myGrid, const Bulk& myBulk);
	void initSize(const Grid& myGrid);
	void initActive(const Grid& myGrid);
	void getIteratorActive();
	void calAreaActive(const Grid& myGrid, const Bulk& myBulk);
	double calAkd(const Grid&myGrid, const Bulk& myBulk, int bIdb, int eIdb);

	void calFlux(const Bulk& myBulk);
	void massConserve(Bulk& myBulk);

	// Assemble Mat
	void allocateMat(Solver& mySolver);
	void initAssembleMat(Solver& mySolver);
	void assembleMat(Solver& mySolver, const Bulk& myBulk);


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

	// Well to Bulk
	int								WellNum;
};
