#include "ReadTool.hpp"

bool ReadLine(ifstream& ifs, vector<string>& result)
{
	result.resize(0);
	string buf;

	while (!ifs.eof())
	{
		getline(ifs, buf);
		if (buf.empty())
			continue;
		while (buf[0] == ' ' || buf[0] == '\t' || buf[0] == '\r')
			buf.erase(0, 1);
		if (buf.empty() || buf[0] == '#')
			continue;
		if (buf.size() > 1 && (buf[0] == '-' && buf[1] == '-'))
			continue;
		
		break;
	}

	// file ends
	if (buf.empty())
		return false;

	// remove the string behind the '/'
	auto pos = buf.find_first_of('/');
	if (pos != string::npos) {
		buf.erase(pos);
		buf.push_back('/');
	}

	// get rid of  '  and  ,
	for (auto& s : buf) {
		if (s == '\'' || s == ',')
			s = ' ';
	}

	istringstream tmp(buf);
	while (tmp >> buf)
		result.push_back(buf);

	return true;
}


void DealDefault(vector<string>& result)
{
	// m*n  -> <n,...,n> size m
	// m*   -> <DEFAULT,..., DEFAULT> size m
	vector<string>  tmp;
	for (auto str : result) {
		auto pos = str.find('*');
		if (pos == string::npos) {
			tmp.push_back(str);
		}
		else {
			int num = atoi(str.substr(0, pos).c_str());  // non number -> 0
			int len = str.size();
			string val = "DEFAULT";
			if (num == 0) {
				tmp.push_back(str);
			}		
			else {
				if (pos != len - 1) {
					val = str.substr(pos + 1, len - (pos + 1));
				}
				for (int i = 0; i < num; i++)
					tmp.push_back(val);
			}			
		}
	}
	swap(result, tmp);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/