/*! \file    OCP_Output.hpp
 *  \brief   OCP_Output class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCP_OUTPUT_HEADER__
#define __OCP_OUTPUT_HEADER__

// Standard header files
#include <iostream>
#include <iomanip>

// OpenCAEPoro header files
#include "Reservoir.hpp"
#include "ParamOutput.hpp"
#include "OCP_Control.hpp"

using namespace std;

class OCPIJK
{
public:
	OCPIJK() = default;
	OCPIJK(const USI& i, const USI& j, const USI& k) :I(i), J(j), K(k) {};
	OCPIJK(const COOIJK& src) { I = src.I; J = src.J; K = src.K; };
	OCPIJK& operator= (const COOIJK& src) { I = src.I; J = src.J; K = src.K; return*this; }
	USI			I, J, K;
};

class OCPTypeA
{
public:
	OCPTypeA& operator= (const Type_A_o& src) { activity = src.activity; obj = src.obj; return *this; }
	bool				activity{ false };
	vector<string>		obj;	///< 
	vector<USI>			index; ///< Records the index of bulk or well, whose properties will be printed.
};


class OCPTypeB
{
public:
	OCPTypeB& operator= (const Type_B_o& src) { activity = src.activity; obj.assign(src.obj.begin(), src.obj.end()); return *this; }
	bool				activity{ false };
	vector<OCPIJK>		obj;
	vector<USI>			index;	///< Records the index of bulk or well, whose properties will be printed.
};

/// SumPair is an auxiliary structure storing summary data to output.
class SumPair
{
public:
	SumPair(const string& item, const string& obj, const string& unit) :Item(item), Obj(obj), Unit(unit) {};
	string				Item;
	string				Obj;
	string				Unit;
	vector<OCP_DBL>		val;
};


/// Summary manages the output for summary file, it contains the most interested information
/// in each time step. usually these data will be convert to chart for analysing following.
class Summary
{
public:
	void InputParam(const OutputSummary& summary_param);
	void Setup(const Reservoir& reservoir);
	void SetVal(const Reservoir& reservoir, const OCP_Control& ctrl);
	void PrintInfo(const string& dir);


private:

	vector<SumPair>		Sumdata; ///< Contains all information to be printed.

	bool		FPR{ false };	///< Field average Pressure.
	bool		FOPR{ false };	///< Field oil production rate.
	bool		FOPT{ false };	///< Field total oil production.
	bool		FGPR{ false };	///< Field gas production rate.
	bool		FGPt{ false };	///< Field total gas production.
	bool		FWPR{ false };	///< Field water production rate.
	bool		FWPT{ false };	///< Field total water production.
	bool		FGIR{ false };	///< Field gas injection rate.
	bool		FGIT{ false };	///< Field total gas injection.
	bool		FWIR{ false };	///< Field water injection rate.
	bool		FWIT{ false };	///< Field total water injection.

	OCPTypeA		WOPR;	///< Well oil production rate.
	OCPTypeA		WOPT;	///< Well total oil production.
	OCPTypeA		WGPR;	///< Well gas production rate.
	OCPTypeA		WGPT;	///< Well total gas production.
	OCPTypeA		WWPR;	///< Well water production rate.
	OCPTypeA		WWPT;	///< Well total water production.
	OCPTypeA		WGIR;	///< Well gas injection rate.
	OCPTypeA		WGIT;	///< Well total gas injection.
	OCPTypeA		WWIR;	///< Well water injection rate.
	OCPTypeA		WWIT;	///< Well total water injection.
	OCPTypeA		WBHP;	///< Well pressure.

	OCPTypeB		BPR;	///< Bulk pressure.
};


/// OCP_Output manages different kinds of ways to output. the most commonly used is summary file.
/// which usually give the information of bulks and wells in each timestep, such as average bulks pressure,
/// oil production rate of wells. if other information at critical time is interested in, you can chose
/// the PRT file(to do). also, some infomation will be printed on the screen at the critical time to make sure
/// the program is at the right way. 
class OCP_Output
{
	friend class OpenCAEPoro;
public:

	void InputParam(ParamOutput& param_Output);
	void Setup(Reservoir& reservoir, string& dir);
	void SetVal(const Reservoir& reservoir, const OCP_Control& ctrl);
	void PrintInfo();

private:
	string		Dir;
	Summary		summary;

};

#endif /* end if __OCPOUTPUT_HEADER__ */


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/