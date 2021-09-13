#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define Paramcheck(exp) \
std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in " << __FILE__ ;

using namespace std;

const string DEFAULT = "DEFAULT";

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

