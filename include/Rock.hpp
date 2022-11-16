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
	virtual void CalPoro(const OCP_DBL& poroInit, const OCP_DBL& P, OCP_DBL& poro, OCP_DBL& dPorodP) const = 0;
	// Calculate porosity and d porosity / d P for thermal model
	virtual void CalPoroT(const OCP_DBL& P, const OCP_DBL& T) const = 0;
};


class RockIso : public Rock
{
public:
	RockIso() = default;
	void CalPoroT(const OCP_DBL& P, const OCP_DBL& T) const override { OCP_ABORT("Not Used!"); };
};



class Rock_Linear01 : public RockIso
{
	// poro = poroinit * (1 + phi)
	// poro = poroInit * (1 + Rc1 * (P - Pref))
public:
	Rock_Linear01() = default;
	Rock_Linear01(const ParamReservoir& rs_param, const USI& i) {}; // temp
	void CalPoro(const OCP_DBL& poroInit, const OCP_DBL& P, OCP_DBL& poro, OCP_DBL& dPorodP) const override;


protected:
	OCP_DBL		Pref;
	OCP_DBL		Rc1;
};



class Rock_Linear02 : public RockIso
{
	// poro = poroinit * (1 + phi)
	// poro = poroInit * (1 + Rc1 * (P - Pref) + Rc2 / 2 * (P - Pref) * (P - Pref))
public:
	Rock_Linear02() = default;
	Rock_Linear02(const ParamReservoir& rs_param, const USI& i) {}; // temp
	void CalPoro(const OCP_DBL& poroInit, const OCP_DBL& P, OCP_DBL& poro, OCP_DBL& dPorodP) const override;


protected:
	OCP_DBL		Pref;
	OCP_DBL		Rc1;
	OCP_DBL     Rc2;
};




#endif /* end if __ROCK_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/15/2022      Create file                          */
/*----------------------------------------------------------------------------*/
