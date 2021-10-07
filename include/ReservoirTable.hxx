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

template <typename T>
class ReservoirTable
{
public:
	ReservoirTable() = default;
	ReservoirTable(const USI row, const USI col);
	void setup(const std::vector<std::vector<T>>& src);
	bool isempty()const { return data.empty(); }

	USI getCol() { return NCol; }
	void pushCol(std::vector<T>& v) { data.push_back(v); }
	std::vector<T>& getCol(int j) { return data[j]; }
	

	void setRowCol() { NRow = data[0].size(); NCol = data.size(); BId = NRow / 2; }
	
	USI eval_all(int j, T val, std::vector<T>& outdata, std::vector<T>& slope);
	T eval(int j, T val, int destj);
	T eval_inv(int j, T val, int destj);

	void display();

private:

	USI									NRow;
	USI									NCol;
	USI									BId;		// search from BId
	std::vector<std::vector<T>>			data;
};


template <typename T>
ReservoirTable<T>::ReservoirTable(const USI row, const USI col)
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
inline USI ReservoirTable<T>::eval_all(int j, T val, std::vector<T>& outdata, std::vector<T>& slope)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val >= data[j][BId]) {
		for (int i = BId + 1; i < NRow; i++) {
			if (val < data[j][i]) {
				BId = i - 1;
				for (int k = 0; k < NCol; k++) {
					slope[k] = (data[k][BId + 1] - data[k][BId]) / (data[j][BId + 1] - data[j][BId]);
					outdata[k] = data[k][BId] + slope[k] * (val - data[j][BId]);
				}
				return BId;
			}	
		}
		for (int k = 0; k < NCol; k++) {
			slope[k] = 0;
			outdata[k] = data[k].back();
		}
	}
	else {
		for (int i = BId - 1; i >= 0; i--) {
			if (val >= data[j][i]) {
				BId = i;
				for (int k = 0; k < NCol; k++) {
					slope[k] = (data[k][BId + 1] - data[k][BId]) / (data[j][BId + 1] - data[j][BId]);
					outdata[k] = data[k][BId] + slope[k] * (val - data[j][BId]);
				}
				return BId;
			}		
		}
		for (int k = 0; k < NCol; k++) {
			slope[k] = 0;
			outdata[k] = data[k].front();
		}
	}
	return BId;
}

template <typename T>
inline T ReservoirTable<T>::eval(int j, T val, int destj)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val >= data[j][BId]) {
		for (int i = BId + 1; i < NRow; i++) {
			if (val < data[j][i]) {
				BId = i - 1;
				OCP_DBL k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].back();
	}
	else {
		for (int i = BId - 1; i >= 0; i--) {
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
inline T ReservoirTable<T>::eval_inv(int j, T val, int destj)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val > data[j][BId]) {
		for (int i = BId - 1; i >= 0; i--) {
			if (val <= data[j][i]) {
				BId = i;
				OCP_DBL k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].front();
	}
	else {
		for (int i = BId + 1; i < NRow; i++) {
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
void ReservoirTable<T>::display()
{
	for (int i = 0; i < NRow; i++) {
		for (int j = 0; j < NCol; j++) {
			std::cout << data[j][i] << "\t";
		}
		std::cout << "\n";
	}
}

#endif