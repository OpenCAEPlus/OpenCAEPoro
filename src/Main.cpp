#include <iostream>
#include <cstdio>
#include <string>
#include "Connection_BB.hpp"
#include "WellGroup.hpp"
#include "Timing.hxx"
#include "OpenCAEPoro_consts.hpp"
#include <algorithm>

using namespace std;



int main()
{
	int a = 10;
	vector<double> v(5, 0);
	v[1] = 2;
	v[2] = 4;
	v[3] = 6;
	v[4] = 8;
	auto temp = find_if(v.begin(), v.end(), [a](auto b) {return b > a; });
	if (temp == v.end())
		cout << "none";
	else
		cout << *temp << endl;


	//int i = distance(v.begin(), find_if(v.begin(), v.end(), [a](auto b) {return b > a; }));
	//cout << v[i];


	return 0;
}
