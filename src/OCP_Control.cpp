#include "OCP_Control.hpp"

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

ControlIter::ControlIter(const vector<OCP_DBL>& src)
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

void OCP_Control::InputParam(ParamControl& CtrlParam)
{
	Dir = CtrlParam.dir;
	if (CtrlParam.method == "IMPES") {
		Method = IMPES;
	}
	else if (CtrlParam.method == "FIM") {
		Method = FIM;
	}
	else {
		ERRORcheck("Wrong Method !");
		exit(0);
	}
	solveFile = CtrlParam.linearSolve;
	criticalTime = CtrlParam.criticalTime;

	USI t = CtrlParam.criticalTime.size();
	ctrlTimeSet.resize(t);
	ctrlErrorSet.resize(t);
	ctrlIterSet.resize(t);

	USI n = CtrlParam.tuning_T.size();
	vector<USI>		ctrlCriticalTime(n + 1);
	for (USI i = 0; i < n; i++) {
		ctrlCriticalTime[i] = CtrlParam.tuning_T[i].d;
	}
	ctrlCriticalTime.back() = t;
	for (USI i = 0; i < n; i++) {
		for (USI d = ctrlCriticalTime[i]; d < ctrlCriticalTime[i + 1]; d++) {
			ctrlTimeSet[d] = ControlTime(CtrlParam.tuning_T[i].Tuning[0]);
			ctrlErrorSet[d] = ControlError(CtrlParam.tuning_T[i].Tuning[1]);
			ctrlIterSet[d] = ControlIter(CtrlParam.tuning_T[i].Tuning[2]);
		}
	}

	cout << "OCP_Control::input" << endl;
}

void OCP_Control::ApplyControl(const USI& i)
{
	ctrlTime = ctrlTimeSet[i];
	ctrlError = ctrlErrorSet[i];
	ctrlIter = ctrlIterSet[i];

	End_time = criticalTime[i + 1];
}


void OCP_Control::InitTime(const USI& i)
{
	OCP_DBL dt = criticalTime[i + 1] - Current_time;
	if (dt < 0) {
		ERRORcheck("Wrong Time Step");
		exit(0);
	}
	Current_dt = min(dt, ctrlTime.TimeInit);
}

void OCP_Control::setNextTstep(Reservoir& reservoir)
{
	Current_time += Current_dt;

	OCP_DBL dPmax = reservoir.bulk.GetdPmax();
	OCP_DBL dNmax = reservoir.bulk.GetdNmax();
	OCP_DBL dSmax = reservoir.bulk.GetdSmax();
	OCP_DBL dVmax = reservoir.bulk.GetdVmax();

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

	if (Current_dt > ctrlTime.TimeMax)
		Current_dt = ctrlTime.TimeMax;
	if (Current_dt < ctrlTime.TimeMin)
		Current_dt = ctrlTime.TimeMin;

	OCP_DBL dt = End_time - Current_time;
	if (Current_dt > dt)
		Current_dt = dt;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/