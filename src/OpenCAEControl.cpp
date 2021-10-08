#include "OpenCAEControl.hpp"

ControlTime::ControlTime(vector<OCP_DBL>& src)
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

ControlError::ControlError(vector<OCP_DBL>& src)
{
	ErrorNL_T = src[1];
	ErrorMB_T = src[2];
	ErrorLS_T = src[3];
	ErrorNL_M = src[5];
	ErrorMB_M = src[6];
	ErrorLS_M = src[7];
}

ControlIter::ControlIter(vector<OCP_DBL>& src)
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

void OCP_Control::inputParam(ParamControl& CtrlParam)
{
	Dir = CtrlParam.Dir;
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
	SolveFile = CtrlParam.LinearSolve;
	CriticalTime = CtrlParam.CriticalTime;

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

	cout << "OCP_Control::input" << endl;
}

void OCP_Control::ApplyControl(int i)
{
	CtrlTime = CtrlTimeSet[i];
	CtrlError = CtrlErrorSet[i];
	CtrlIter = CtrlIterSet[i];

	End_time = CriticalTime[i + 1];
}


void OCP_Control::initTime(int i)
{
	OCP_DBL dt = CriticalTime[i + 1] - Current_time;
	if (dt < 0) {
		ERRORcheck("Wrong Time Step");
		exit(0);
	}
	Current_dt = min(dt, CtrlTime.TimeInit);
}

void OCP_Control::setNextTstep(Reservoir& reservoir)
{
	Current_time += Current_dt;

	OCP_DBL dPmax = reservoir.bulk.getdPmax();
	OCP_DBL dNmax = reservoir.bulk.getdNmax();
	OCP_DBL dSmax = reservoir.bulk.getdSmax();
	OCP_DBL dVmax = reservoir.bulk.getdVmax();

	OCP_DBL dPlim = 300;
	OCP_DBL dNlim = 0.3;
	OCP_DBL dSlim = 0.3;
	OCP_DBL dVlim = 0.001;

	OCP_DBL c1 = dPmax / dPlim;
	OCP_DBL c2 = dNmax / dNlim;
	OCP_DBL c3 = dSmax / dSlim;
	OCP_DBL c4 = dVmax / dVlim;

	OCP_DBL c = max(max(c1, c2), max(c3, c4));
	c = max(1.0 / 3, c);
	c = min(10.0 / 3, c);

	Current_dt /= c;

	if (Current_dt > CtrlTime.TimeMax)
		Current_dt = CtrlTime.TimeMax;
	if (Current_dt < CtrlTime.TimeMin)
		Current_dt = CtrlTime.TimeMin;

	OCP_DBL dt = End_time - Current_time;
	if (Current_dt > dt)
		Current_dt = dt;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/