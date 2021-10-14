#pragma once
#include<vector>

#include "OpenCAEPoro_consts.hpp"

// CSRX
template <typename T>
class MAT_Faspxx
{
public:
	void ClearData();

	OCP_USI						NRow;
	OCP_USI						NCol;
	OCP_USI						Nnz;
	std::vector<T>			val;		
	std::vector<OCP_USI>		colId;
	std::vector<OCP_USI>		RowPtr;
	std::vector<USI>		diagPtr;

};


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/