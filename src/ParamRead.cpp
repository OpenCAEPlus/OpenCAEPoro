#include "ParamRead.hpp"


void ParamRead::getDirAndName(string file)
{
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64) 
	// for Window file system
	int pos = file.find_last_of('\\') + 1;
	FileDir = file.substr(0, pos);
	FileName = file.substr(pos, file.size() - pos);
#else
	// for Linux and Mac OSX file system
	int pos = file.find_last_of('/') + 1;
	FileDir = file.substr(0, pos);
	FileName = file.substr(pos, file.size() - pos);
#endif
}

void ParamRead::readFile(string file)
{
	ifstream ifs(file, ios::in);
	if (!ifs) {
		cout << "can not open " << file << "\n";
		exit(0);
	}

	while (!ifs.eof())
	{
		vector<string> vbuf;
		if (!ReadLine(ifs, vbuf))
			break;
		string keyword = vbuf[0];
		switch (Map_str2int(&keyword[0], keyword.size()))
		{
		case Map_str2int("DIMENS", 6):
			cout << "DIMENS" << endl;
			Rs_param.inputDIMENS(ifs);
			Rs_param.outputDIMENS();
			break;

		case Map_str2int("EQUALS", 6):
			Rs_param.inputEQUALS(ifs);
			break;

		case Map_str2int("COPY", 4):
			Rs_param.inputCOPY(ifs);
			break;

		case Map_str2int("MULTIPLY", 8):
			Rs_param.inputMULTIPLY(ifs);
			break;

		case Map_str2int("SWOF", 4):
		case Map_str2int("SGOF", 4):
		case Map_str2int("PVCO", 4):
		case Map_str2int("PVDG", 4):
		case Map_str2int("PVTW", 4):
		case Map_str2int("PBVD", 4):
			Rs_param.inputTABLE(ifs, keyword);
			break;

		case Map_str2int("ROCK", 4):
			cout << "ROCK" << endl;
			break;

		case Map_str2int("GRAVITY", 7):
			cout << "GRAVITY" << endl;
			break;

		case Map_str2int("DENSITY", 7):
			cout << "DENSITY" << endl;
			break;

		case Map_str2int("EQUIL", 5):
			cout << "EQUIL" << endl;
			break;

		default:
			break;
		}
	}



	ifs.close();
}
