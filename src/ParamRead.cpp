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
		case Map_str2int("'DX'", 4):
		case Map_str2int("'DY'", 4):
		case Map_str2int("'DZ'", 4):
		case Map_str2int("'PORO'", 6):
		case Map_str2int("'PERMX'", 7):
		case Map_str2int("'PERMY'", 7):
		case Map_str2int("'PERMZ'", 7):
			Rs_param.inputEQUALS(vbuf);
			break;
		case Map_str2int("'TOPS'", 6):
			cout << "'TOPS'" << endl;
			break;
		case Map_str2int("COPY", 4):
			cout << "COPY" << endl;
			break;
		case Map_str2int("SWOF", 4):
			cout << "SWOF" << endl;
			break;
		case Map_str2int("SGOF", 4):
			cout << "SGOF" << endl;
			break;
		case Map_str2int("PVCO", 4):
			cout << "PVCO" << endl;
			break;
		case Map_str2int("PVDG", 4):
			cout << "PVDG" << endl;
			break;
		case Map_str2int("PVTW", 4):
			cout << "PVTW" << endl;
			break;
		case Map_str2int("ROCK", 4):
			cout << "ROCK" << endl;
			break;
		case Map_str2int("GRAVITY", 4):
			cout << "GRAVITY" << endl;
			break;
		case Map_str2int("DENSITY", 4):
			cout << "DENSITY" << endl;
			break;
		case Map_str2int("EQUIL", 4):
			cout << "EQUIL" << endl;
			break;
		case Map_str2int("PBVD", 4):
			cout << "PBVD" << endl;
			break;
		default:
			break;
		}
	}



	ifs.close();
}
