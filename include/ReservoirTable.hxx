#pragma once
#include <vector>

template <typename T>
class ReservoirTable
{
public:


	int eval_all(int j, T val, std::vector<T>& outdata, std::vector<T>& slope);
	T eval(int j, T val, int destj);


private:

	short								NRow;
	short								NCol;
	short								BId;		// search from BId
	std::vector<std::vector<T>>			data;
};


template <typename T>
inline int ReservoirTable<T>::eval_all(int j, T val, std::vector<T>& outdata, std::vector<T>& slope)
{
	// becareful when the memory outdata and slope have not be allocated before

	if (val >= data[j][BId]) {
		for (int i = BId + 1; i < NCol; i++) {
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
		for (int i = BId + 1; i < NCol; i++) {
			if (val < data[j][i]) {
				BId = i - 1;
				double k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].back();
	}
	else {
		for (int i = BId - 1; i >= 0; i--) {
			if (val >= data[j][i]) {
				BId = i;
				double k = (data[destj][BId + 1] - data[destj][BId]) / (data[j][BId + 1] - data[j][BId]);
				return (data[destj][BId] + k * (val - data[j][BId]));
			}
		}
		return data[destj].front();
	}
}
