#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define Paramcheck(exp) \
std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in " << __FILE__ ;

using namespace std;


constexpr inline long long Map_str2int(const char* mystr, int len)
{
	long long res = 0;
	long long t = 100;
	for (int i = 0; i < len; i++) {
		res += (int)mystr[len - 1 - i] * t;
		t *= 100;
	}
	return res;
}

bool ReadLine(ifstream& ifs, vector<string>& result);

void DealDefault(vector<string>& result);


template<typename T>
void DealData(vector<string>& vbuf, vector<int>& obj, vector<T>& region)
{
	//  8*1  16*2  8*3  16*4  ->
	// obj <8, 16, 8, 16>   val  <1, 2, 3, 4>

	obj.resize(0);
	region.resize(0);
	for (auto str : vbuf) {
		auto pos = str.find('*');

		if (pos != string::npos) {
			int len = str.size();
			int num = stoi(str.substr(0, pos));
			int val = stoi(str.substr(pos + 1, len - (pos + 1)));
			obj.push_back(num);
			region.push_back(val);
		}
	}
}