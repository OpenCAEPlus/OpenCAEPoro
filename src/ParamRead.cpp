#include "ParamRead.hpp"

void ParamRead::init()
{
	Rs_param.init();
	Control_param.init();
}

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
			Rs_param.inputDIMENS(ifs);
			Rs_param.outputDIMENS();
			break;

		case Map_str2int("RTEMP", 5):
			Rs_param.inputRTEMP(ifs);
			break;

		case Map_str2int("EQUALS", 6):
			Rs_param.inputEQUALS(ifs);
			break;

		case Map_str2int("DX", 2):
		case Map_str2int("DY", 2):
		case Map_str2int("DZ", 2):
		case Map_str2int("NTG", 3):
		case Map_str2int("PORO", 4):
		case Map_str2int("TOPS", 4):
		case Map_str2int("PERMX", 5):
		case Map_str2int("PERMY", 5):
		case Map_str2int("PERMZ", 5):
			Rs_param.inputGRID(ifs, keyword);
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
			Rs_param.inputROCK(ifs);
			break;

		case Map_str2int("GRAVITY", 7):
			Rs_param.inputGRAVITY(ifs);
			break;

		case Map_str2int("DENSITY", 7):
			Rs_param.inputDENSITY(ifs);
			break;

		case Map_str2int("EQUIL", 5):
			Rs_param.inputEQUIL(ifs);
			break;

		case Map_str2int("INCLUDE", 7):
			inputINCLUDE(ifs);
			break;

		case Map_str2int("METHOD", 6):
			Control_param.inputMETHOD(ifs);
			break;

		case Map_str2int("TUNING", 6):
			Control_param.inputTUNING(ifs);
			break;

		default:
			break;
		}
	}



	ifs.close();
}


void ParamRead::inputINCLUDE(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	DealDefault(vbuf);

	readFile(FileDir + vbuf[0]);

}


// check
void ParamRead::checkParam()
{
	Rs_param.checkParam();
}
