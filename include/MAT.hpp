#pragma once
#include<vector>

// CSRX
template <typename T>
class MAT
{
public:
	void clearData();

	int						NRow;
	int						NCol;
	int						Nnz;		
	std::vector<T>			Val;		
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
