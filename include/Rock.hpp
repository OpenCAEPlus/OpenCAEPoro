/*! \file    Rock.hpp
 *  \brief   Rock class declaration
 *  \author  Shizhe Li
 *  \date    Nov/15/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ROCK_HEADER__
#define __ROCK_HEADER__

#include <math.h>

 // OpenCAEPoro header files
#include "OCPConst.hpp"
#include "OCPTable.hpp"
#include "ParamReservoir.hpp"

class Rock
{
public:
	Rock() = default;

	// Calculate porosity and d porosity / d P for isothermal model
	virtual void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit, OCP_DBL& poro, OCP_DBL& dPorodP) const = 0;
	// Calculate porosity and d porosity / d P for thermal model
	virtual void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro) const = 0;
};


class RockIso : public Rock
{
public:
	RockIso() = default;
	void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro) const override { OCP_ABORT("Not Used!"); };
};



class Rock_Linear : public RockIso
{
	// poro = poroinit * (1 + phi)
	// poro = poroInit * (1 + cp1 * (P - Pref) + cp2 / 2 * (P - Pref) * (P - Pref))
public:
	Rock_Linear() = default;
	Rock_Linear(const RockParam& param) { Pref = param.Pref; cp1 = param.Cp1; cp2 = param.Cp2; };
	void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit, OCP_DBL& poro, OCP_DBL& dPorodP) const override;


protected:
	OCP_DBL		Pref;
	OCP_DBL		cp1;
	OCP_DBL     cp2;
};



class RockT : public Rock
{
public:
	RockT() = default;
	void Assign(const RockParam& param) {
		Pref = param.Pref;
		Tref = param.Tref;
		cp = param.Cp1;
		ct = param.Ct;
		cpt = param.Cpt;
		ConstRock = param.ConstRock;
	}
	void CalPoro(const OCP_DBL& P, const OCP_DBL& poroInit, OCP_DBL& poro, OCP_DBL& dPorodP) const override { OCP_ABORT("Not Used!"); };

protected:
	OCP_DBL		Pref;
	OCP_DBL     Tref;
	OCP_DBL		cp;
	OCP_DBL     ct;
	OCP_DBL     cpt;
	OCP_BOOL	ConstRock;
};


class RockT_Linear : public RockT
{
	// poro = poroinit * (1 + phi)
	// poro = poroInit * (1 + cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
	RockT_Linear() = default;
	RockT_Linear(const RockParam& param) { Assign(param); };
	void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro) const override;
};


class RockT_Exp : public RockT
{
	// poro = poroinit * (1 + exp(phi))
	// poro = poroInit * exp(cp*(P-Pref) - ct*(T-Tref) + cpt*(P-Pref)*(T-Tref))
public:
	RockT_Exp() = default;
	RockT_Exp(const RockParam& param) { Assign(param); };
	void CalPoroT(const OCP_DBL& P, const OCP_DBL& T, const OCP_DBL& poroInit, OCP_DBL& poro) const override;
};



#endif /* end if __ROCK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
