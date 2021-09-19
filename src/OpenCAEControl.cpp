#include "OpenCAEControl.hpp"

ControlTime::ControlTime(vector<double>& src)
{
	TimeInit = src[0];
	TimeMax = src[1];
	TimeMin = src[2];
	TimeMinChop = src[3];
	TimeMaxIncr = src[4];
	TimeMinCut = src[5];
	TimeCut_F = src[6];
	TimeMaxIncre_F = src[7];
}

ControlError::ControlError(vector<double>& src)
{
	ErrorNL_T = src[1];
	ErrorMB_T = src[2];
	ErrorLS_T = src[3];
	ErrorNL_M = src[5];
	ErrorMB_M = src[6];
	ErrorLS_M = src[7];
}

ControlIter::ControlIter(vector<double>& src)
{
	ItMax_NT = src[0];
	ItMin_NT = src[1];
	ItMax_NTL = src[2];
	ItMin_NTL = src[3];
	DPreNT_M = src[6];
	DSatNT_M = src[7];
	DPreNT_T = src[8];
	Dpre_M = src[9];
	if (Dpre_M < 0)
		Dpre_M = DPreNT_M;
}

void CAEControl::inputParam(ParamControl& CtrlParam)
{
	if (CtrlParam.Method == "IMPES") {
		Method = IMPES;
	}
	else if (CtrlParam.Method == "FIM") {
		Method = FIM;
	}
	else {
		ERRORcheck("Wrong Method !");
		exit(0);
	}
	LinearSolve = CtrlParam.LinearSolve;


	int t = CtrlParam.CriticalTime.size();
	CtrlTimeSet.resize(t);
	CtrlErrorSet.resize(t);
	CtrlIterSet.resize(t);

	int n = CtrlParam.Tuning_T.size();
	vector<int>		ctrlCriticalTime(n + 1);
	for (int i = 0; i < n; i++) {
		ctrlCriticalTime[i] = CtrlParam.Tuning_T[i].d;
	}
	ctrlCriticalTime.back() = t;
	for (int i = 0; i < n; i++) {
		for (int d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
			CtrlTimeSet[d] = ControlTime(CtrlParam.Tuning_T[i].Tuning[0]);
			CtrlErrorSet[d] = ControlError(CtrlParam.Tuning_T[i].Tuning[1]);
			CtrlIterSet[d] = ControlIter(CtrlParam.Tuning_T[i].Tuning[2]);
		}
	}

	cout << "CAEControl::input" << endl;
}

void CAEControl::ApplyControl(int i)
{
	CtrlTime = CtrlTimeSet[i];
	CtrlError = CtrlErrorSet[i];
	CtrlIter = CtrlIterSet[i];
}
