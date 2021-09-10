#pragma once
#include<vector>

// CSRX
class MAT
{
	friend class Solver;
public:
	void clearData();

private:
	int						NRow;
	int						NCol;
	int						Nnz;		
	std::vector<double>		Val;		
	std::vector<int>		ColId;
	std::vector<int>		RowPtr;
	std::vector<int>		DiagPtr;

};

//void MAT::clearData()
//{
//	Val.clear();		
//	ColId.clear();
//	RowPtr.clear();		
//	DiagPtr.clear();
//}
