#pragma once
#include <vector>
#include "OpenCAEPoro_consts.hpp"
#include "ParamControl.hpp"

using namespace std;

class ControlTime
{
public:
	ControlTime() = default;
	ControlTime(vector<double>& src);
	double		TimeInit;		// Maximum init step length of next timestep
	double		TimeMax;		// Maximum time step during running
	double		TimeMin;		// Minmum time step during running
	double		TimeMinChop;	// Minimum choppable timestep
	double		TimeMaxIncr;	// Maximum timestep increase factor
	double		TimeMinCut;		// Minimum timestep cutback factor
	double		TimeCut_F;		// Minimum timestep cutback factor after convergence failure
	double		TimeMaxIncre_F;	// Maximum increase factor after a convergence failure
};

class ControlError
{
public:
	ControlError() = default;
	ControlError(vector<double>& src);
	double		ErrorNL_T;		// Target non-linear convergence error
	double		ErrorMB_T;		// Target material balance error
	double		ErrorLS_T;		// Target linear convergence error
	double		ErrorNL_M;		// Maximum non-linear convergence error
	double		ErrorMB_M;		// Maximum material balance error
	double		ErrorLS_M;		// Maximum linear convergence error
};

class ControlIter
{
public:
	ControlIter() = default;
	ControlIter(vector<double>& src);
	int			ItMax_NT;		// Maximum number of Newton iterations in a timestep
	int			ItMin_NT;		// Minimum number of Newton iterations in a timestep
	int			ItMax_NTL;		// Maximum number of linear iterations in a Newton iteration
	int			ItMin_NTL;		// Minimum number of linear iterations in a Newton iteration
	double		DPreNT_M;		// Maximum pressure change at last Newton iteration
	double		DSatNT_M;		// Maximum saturation change at last Newton iteration
	double		DPreNT_T;		// Target pressure change at last Newton iteration
	double		Dpre_M;			// Target maximum pressure change in a timestep
};

class CAEControl
{
	friend class OpenCAEPoro;
public:

	void inputParam(ParamControl& CtrlParam);
	void ApplyControl(int i);

	void initTime(int i);

	double getCurTime() { return Current_time; }
	int getLSiter() { return LS_iter; }
	int getNRiter() { return NR_iter; }

private:

	string						Dir;
	int							Method;
	string						SolveFile;
	vector<double>				CriticalTime;
	double						Current_dt{ 1E5 };
	double						Current_time{ 0 };
	int							Tstep{ 0 };
	int							LS_iter{ 0 };
	int							LS_iter_total{ 0 };
	int							NR_iter{ 0 };
	int							NR_iter_total{ 0 };


	ControlTime					CtrlTime;
	vector<ControlTime>			CtrlTimeSet;

	ControlError				CtrlError;
	vector<ControlError>		CtrlErrorSet;

	ControlIter					CtrlIter;
	vector<ControlIter>			CtrlIterSet;

};
