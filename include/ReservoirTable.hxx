/*! \file    ReservoirTable.hpp
 *  \brief   ReservoirTable class declaration
 *  \author  Shizhe Li
 *  \date    Oct/07/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __RESERVOIRTABLE_HEADER__
#define __RESERVOIRTABLE_HEADER__


#include "OpenCAEPoro_consts.hpp"
#include <iostream>
#include <vector>

using namespace std;

/// A template Table class, which used to deal with everything about table in OpenCAEPoro
/// such as PVT table, saturation table.
template <typename T>
class ReservoirTable
{
public:
	ReservoirTable() = default;
	ReservoirTable(const USI& row, const USI& col);
	/// setup tables from existing data of table.
	void setup(const vector<vector<T>>& src);
	/// judge if table is empty.
	bool isempty() const { return data.empty(); }
	/// return the column num of table.
	USI getCol() const { return NCol; }
	/// push v into the last column of table.
	void pushCol(const vector<T>& v) { data.push_back(v); }
	/// return the jth column in table to modify or use.
	vector<T>& getCol(const USI& j) { return data[j]; }
	/// setup row nums and col nums of tables, initialize the BId.
	void setRowCol() { NRow = data[0].size(); NCol = data.size(); BId = NRow / 2; }
	/// interpolate the specified monotonically increasing column in table to evaluate all columns.
	USI eval_all(const USI& j, const T& val, vector<T>& outdata, vector<T>& slope);
	/// interpolate the specified monotonically increasing column in table to evaluate the target column.
	T eval(const USI& j, const T& val, const USI& destj);
	/// interpolate the specified monotonically decreasing column in table to evaluate the target column.
	T eval_inv(const USI& j, const T& val, const USI& destj);
	/// display the data of table on screen.
	void display() const;

private:

	USI									NRow;		///< number of rows of the table
	USI									NCol;		///< number of columns of the table
	USI									BId;		///< the starting point of rows when interpolating
	vector<vector<T>>					data;		///< stores the data of table, data[i] contains the ith column of table.
};


template <typename T>
ReservoirTable<T>::ReservoirTable(const USI& row, const USI& col)
{
	NRow = row;
	NCol = col;
	BId = 0;
	data.resize(NCol);
	for (USI j = 0; j < NCol; j++) {
		data[j].resize(NRow);
	}
}

template <typename T>
void ReservoirTable<T>::setup(const std::vector<std::vector<T>>& src)
{
	data = src;
	NCol = data.size();
	NRow = data[0].size();
	BId = NRow / 2;
}


template <typename T>
inline USI ReservoirTable<T>::eval_all(const USI& j, const T& val, vector<T>& outdata, vector<T>& slope)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val >= data[j][BId]) {
		for (USI i = BId + 1; i < NRow; i++) {
			if (val < data[j][i]) {
				BId = i - 1;
				for (USI k = 0; k < NCol; k++) {
					slope[k] = (data[k][BId + 1] - data[k][BId]) / (data[j][BId + 1] - data[j][BId]);
					outdata[k] = data[k][BId] + slope[k] * (val - data[j][BId]);
				}
				return BId;
			}	
		}
		for (USI k = 0; k < NCol; k++) {
			slope[k] = 0;
			outdata[k] = data[k].back();
		}
	}
	else {
		for (OCP_INT i = BId - 1; i >= 0; i--) {
			if (val >= data[j][i]) {
				BId = i;
				for (USI k = 0; k < NCol; k++) {
					slope[k] = (data[k][BId + 1] - data[k][BId]) / (data[j][BId + 1] - data[j][BId]);
					outdata[k] = data[k][BId] + slope[k] * (val - data[j][BId]);
				}
				return BId;
			}		
		}
		for (USI k = 0; k < NCol; k++) {
			slope[k] = 0;
			outdata[k] = data[k].front();
		}
	}
	return BId;
}

template <typename T>
inline T ReservoirTable<T>::eval(const USI& j, const T& val, const USI& destj)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val >= data[j][BId]) {
		for (USI i = BId + 1; i < NRow; i++) {
			if (val < data[j][i]) {
				BId = i - 1;
				OCP_DBL k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].back();
	}
	else {
		for (OCP_INT i = BId - 1; i >= 0; i--) {
			if (val >= data[j][i]) {
				BId = i;
				OCP_DBL k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].front();
	}
}

template <typename T>
inline T ReservoirTable<T>::eval_inv(const USI& j, const T& val, const USI& destj)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val > data[j][BId]) {
		for (OCP_INT i = BId - 1; i >= 0; i--) {
			if (val <= data[j][i]) {
				BId = i;
				OCP_DBL k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].front();
	}
	else {
		for (USI i = BId + 1; i < NRow; i++) {
			if (val >= data[j][i]) {
				BId = i;
				OCP_DBL k = (data[destj][BId] - data[destj][BId - 1]) / (data[j][BId] - data[j][BId - 1]);
				return (data[destj][BId - 1] + k * (val - data[j][BId - 1]));
			}
		}
		return data[destj].back();
	}
}


template <typename T>
void ReservoirTable<T>::display() const
{
	for (USI i = 0; i < NRow; i++) {
		for (USI j = 0; j < NCol; j++) {
			cout << data[j][i] << "\t";
		}
		cout << "\n";
	}
}

#endif