/*! \file    WellGroup.hpp
 *  \brief   WellGroup class declaration
 *  \author  Shizhe Li
 *  \date    Oct/04/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLGROUP_HEADER__
#define __WELLGROUP_HEADER__


#include "Well.hpp"
#include "ParamWell.hpp"

using namespace std;

/// contains all wells, all operations related to the well will be managed here.
class WellGroup
{
public:
	WellGroup() = default;
	/// input param from ParamWell.
	void inputParam(const ParamWell& Well_param);
	/// setup well and mixture in wellgroup.
	void setup(const Grid& myGrid, const Bulk& myBulk);
	/// complete the information of well according to Grid and Bulk.
	void setupWell(const Grid& myGrid, const Bulk& myBulk);
	/// get the mixture from bulk.
	void setupMixture(const Bulk& myBulk);
	/// calculate the CFL number for each perforation and return the maximum one.
	OCP_DBL calCFL(const Bulk& myBulk, const OCP_DBL& dt) const;
	/// calculate volume flow rate and moles flow rate of each perforation.
	void calFlux(const Bulk& myBulk);
	/// update moles of components in bulk which connects to well, according to the well flow rate.
	void massConserve(Bulk& myBulk, OCP_DBL dt);
	/// calculate memory needed to assemble matrix, only parts related to well are considered here.
	template<typename T>
	void allocateMat(Solver<T>& mySolver) const;
	/// guess the initial well pressure, it equals pressure in bulks where topest perforation locates.
	void init(const Bulk& myBulk);
	/// change operation mode of well at the ith critical time point, which decided by user input.
	void applyControl(USI i);
	/// calculate well properties at the beginning of each time step.
	void prepareWell(const Bulk& myBulk);
	/// calculate injection rate, total injection, production rate, total production for each well.
	void calIPRT(const Bulk& myBulk, OCP_DBL dt);
	/// assemble matrix, parts related to well are included. only for IMPES method.
	/// it should be called after parts related to bulks setups.
	void assemblaMat_WB_IMPES(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, const OCP_DBL& dt) const;
	/// return the num of wells.
	USI getWellNum() const { return WellNum; }
	/// return the name of specified well.
	const string& getWellName(const USI& i) const { return WellG[i].Name; }
	/// return the index of specified well.
	USI getIndex(const string& name) const;

	// Field injection / production
	/// return oil production rate in field.
	OCP_DBL getFOPR() const { return FOPR; }
	/// return total oil production in field.
	OCP_DBL getFOPT() const { return FOPT; }
	/// return gas production rate in field.
	OCP_DBL getFGPR() const { return FGPR; }
	/// return total gas production in field.
	OCP_DBL getFGPT() const { return FGPt; }
	/// return water production rate in field.
	OCP_DBL getFWPR() const { return FWPR; }
	/// return total water production in field.
	OCP_DBL getFWPT() const { return FWPT; }
	/// return gas injection rate in field.
	OCP_DBL getFGIR() const { return FGIR; }
	/// return gas water injection in field.
	OCP_DBL getFGIT() const { return FGIT; }
	/// return water injection rate in field.
	OCP_DBL getFWIR() const { return FWIR; }
	/// return total water injection in field.
	OCP_DBL getFWIT() const { return FWIT; }

	// Well injection / production
	/// return oil production rate of the wth well.
	OCP_DBL getWOPR(USI w) const { return WellG[w].WOPR; }
	/// return total oil production of the wth well.
	OCP_DBL getWOPT(USI w) const { return WellG[w].WOPT; }
	/// return gas production rate of the wth well.
	OCP_DBL getWGPR(USI w) const { return WellG[w].WGPR; }
	/// return total gas production of the wth well.
	OCP_DBL getWGPT(USI w) const { return WellG[w].WGPT; }
	/// return water production rate of the wth well.
	OCP_DBL getWWPR(USI w) const { return WellG[w].WWPR; }
	/// return total water production of the wth well.
	OCP_DBL getWWPT(USI w) const { return WellG[w].WWPT; }
	/// return gas injection rate of the wth well.
	OCP_DBL getWGIR(USI w) const { return WellG[w].WGIR; }
	/// return total gas injection of the wth well.
	OCP_DBL getWGIT(USI w) const { return WellG[w].WGIT; }
	/// return water injection rate of the wth well.
	OCP_DBL getWWIR(USI w) const { return WellG[w].WWIR; }
	/// return total water injection of the wth well.
	OCP_DBL getWWIT(USI w) const { return WellG[w].WWIT; }
	// BHP
	/// return the BHP of wth well.
	OCP_DBL getWBHP(USI w) const { return WellG[w].BHP; }
	/// update pressure in well and well perforation with solution of linear system.
	void getSol_IMPES(const vector<OCP_DBL>& u, const OCP_USI& bid);
	/// reset dG to ldG for each well, dG is a array where the pressure difference between well and perforation are stored.
	void setLastStep() { for (auto& w : WellG)	w.ldG = w.dG; }
	/// check if unreasonable well pressure or perforation pressure occurs.
	int checkP(const Bulk& myBulk);

private:

	USI							WellNum;		///< num of wells.
	std::vector<Well>			WellG;			///< well set.
	
	std::vector<Mixture*>		Flashcal;		///< used to flash calculation for well, uesless now.
	
	OCP_DBL						FGIR{ 0 };		///< gas injection rate in field.
	OCP_DBL						FGIT{ 0 };		///< gas total injection in field.
	OCP_DBL						FWIR{ 0 };		///< water injection rate in field.
	OCP_DBL						FWIT{ 0 };		///< water total injection in field.
	OCP_DBL						FOPR{ 0 };		///< oil production rate in field.
	OCP_DBL						FOPT{ 0 };		///< oil total production in field.
	OCP_DBL						FGPR{ 0 };		///< gas production rate in field.
	OCP_DBL						FGPt{ 0 };		///< gas total production in field.
	OCP_DBL						FWPR{ 0 };		///< water production rate in field.
	OCP_DBL						FWPT{ 0 };		///< water total production in field.
};

template<typename T>
void WellGroup::allocateMat(Solver<T>& mySolver) const
{
	for (USI w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}

#endif