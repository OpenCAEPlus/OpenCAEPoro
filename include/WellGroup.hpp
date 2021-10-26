/*! \file    WellGroup.hpp
 *  \brief   WellGroup class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __WELLGROUP_HEADER__
#define __WELLGROUP_HEADER__

// OpenCAEPoro header files
#include "ParamWell.hpp"
#include "Well.hpp"

using namespace std;

/// WellGroups contains all wells now, it's used to manages all wells uniformly in
/// reservoirs. actually, you can regard it as an interface between wells and other
/// modules.
class WellGroup
{
public:
    WellGroup() = default;
    /// input param from ParamWell.
    void InputParam(const ParamWell& paramWell);
    /// Setup well and mixture in wellgroup.
    void Setup(const Grid& myGrid, const Bulk& myBulk);
    /// complete the information of well according to Grid and Bulk.
    void SetupWell(const Grid& myGrid, const Bulk& myBulk);
    /// get the mixture from bulk.
    void SetupMixture(const Bulk& myBulk);
    /// calculate the CFL number for each perforation and return the maximum one.
    OCP_DBL CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate the CFL number for each perforation and return the maximum one.
    void CalCFL01(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// calculate volume flow rate and moles flow rate of each perforation.
    void CalFlux(const Bulk& myBulk);
    /// calculate dG.
    void CaldG(const Bulk& myBulk);
    /// update moles of components in bulk which connects to well, according to the well
    /// flow rate.
    void MassConserve(Bulk& myBulk, OCP_DBL dt);
    /// calculate memory needed to assemble matrix, only parts related to well are
    /// considered here.
    void AllocateMat(LinearSolver& mySolver, const USI& bulknum) const;
    /// guess the initial well pressure, it equals pressure in bulks where topest
    /// perforation locates.
    void Init(const Bulk& myBulk);
    /// change operation mode of well at the ith critical time point, which decided by
    /// user input.
    void ApplyControl(const USI& i);
    /// calculate well properties at the beginning of each time step.
    void PrepareWell(const Bulk& myBulk);
    /// calculate injection rate, total injection, production rate, total production for
    /// each well.
    void CalIPRT(const Bulk& myBulk, OCP_DBL dt);
    /// assemble matrix, parts related to well are included. only for IMPEC method.
    /// it should be called after parts related to bulks setups.
    void AssemblaMat_WB_IMPEC(LinearSolver& mySolver, const Bulk& myBulk,
                              const OCP_DBL& dt) const;
    void AssemblaMat_WB_FIM(LinearSolver& mySolver, const Bulk& myBulk,
        const OCP_DBL& dt) const;
    /// Return the num of wells.
    USI GetWellNum() const { return numWell; }
    /// Return the name of specified well.
    string GetWellName(const USI& i) const { return wellGroup[i].name; }
    /// Return the index of specified well.
    USI GetIndex(const string& name) const;
    /// Return the num of perforations of well i.
    USI GetWellPerfNum(const USI& i)const { return wellGroup[i].numPerf; }
    USI GetMaxWellPerNum() const;
    // Field injection / production
    /// return oil production rate in field.
    OCP_DBL GetFOPR() const { return FOPR; }
    /// return total oil production in field.
    OCP_DBL GetFOPT() const { return FOPT; }
    /// return gas production rate in field.
    OCP_DBL GetFGPR() const { return FGPR; }
    /// return total gas production in field.
    OCP_DBL GetFGPT() const { return FGPt; }
    /// return water production rate in field.
    OCP_DBL GetFWPR() const { return FWPR; }
    /// return total water production in field.
    OCP_DBL GetFWPT() const { return FWPT; }
    /// return gas injection rate in field.
    OCP_DBL GetFGIR() const { return FGIR; }
    /// return gas water injection in field.
    OCP_DBL GetFGIT() const { return FGIT; }
    /// return water injection rate in field.
    OCP_DBL GetFWIR() const { return FWIR; }
    /// return total water injection in field.
    OCP_DBL GetFWIT() const { return FWIT; }

    // Well injection / production
    /// return oil production rate of the wth well.
    OCP_DBL GetWOPR(const USI& w) const { return wellGroup[w].WOPR; }
    /// return total oil production of the wth well.
    OCP_DBL GetWOPT(const USI& w) const { return wellGroup[w].WOPT; }
    /// return gas production rate of the wth well.
    OCP_DBL GetWGPR(const USI& w) const { return wellGroup[w].WGPR; }
    /// return total gas production of the wth well.
    OCP_DBL GetWGPT(const USI& w) const { return wellGroup[w].WGPT; }
    /// return water production rate of the wth well.
    OCP_DBL GetWWPR(const USI& w) const { return wellGroup[w].WWPR; }
    /// return total water production of the wth well.
    OCP_DBL GetWWPT(const USI& w) const { return wellGroup[w].WWPT; }
    /// return gas injection rate of the wth well.
    OCP_DBL GetWGIR(const USI& w) const { return wellGroup[w].WGIR; }
    /// return total gas injection of the wth well.
    OCP_DBL GetWGIT(const USI& w) const { return wellGroup[w].WGIT; }
    /// return water injection rate of the wth well.
    OCP_DBL GetWWIR(const USI& w) const { return wellGroup[w].WWIR; }
    /// return total water injection of the wth well.
    OCP_DBL GetWWIT(const USI& w) const { return wellGroup[w].WWIT; }
    // BHP
    /// return the BHP of wth well.
    OCP_DBL GetWBHP(const USI& w) const { return wellGroup[w].BHP; }
    /// Return the pth dG of wth well.
    OCP_DBL GetWellDg(const USI& w, const USI& p)const { return wellGroup[w].dG[p]; }

    /// update pressure in well and well perforation with solution of linear system.
    void GetSol_IMPEC(const vector<OCP_DBL>& u, const OCP_USI& bId);
    void GetSol_FIM(const vector<OCP_DBL>& u, const OCP_USI& bId, const USI& len);
    void CalResFIM(const Bulk& myBulk, const OCP_DBL& dt) const;
    /// reset dG to ldG for each well, dG is a array where the pressure difference
    /// between well and perforation are stored.
    void UpdateLastStep()
    {
        for (auto& w : wellGroup) w.ldG = w.dG;
    }
    void ResetDg(){ for (auto& w : wellGroup) w.dG = w.ldG; }
    /// check if unreasonable well pressure or perforation pressure occurs.
    OCP_INT CheckP(const Bulk& myBulk);

private:
    USI          numWell;   ///< num of wells.
    vector<Well> wellGroup; ///< well set.

    vector<Mixture*> flashCal; ///< used to flash calculation for well, uesless now.

    OCP_DBL FGIR{0}; ///< gas injection rate in field.
    OCP_DBL FGIT{0}; ///< gas total injection in field.
    OCP_DBL FWIR{0}; ///< water injection rate in field.
    OCP_DBL FWIT{0}; ///< water total injection in field.
    OCP_DBL FOPR{0}; ///< oil production rate in field.
    OCP_DBL FOPT{0}; ///< oil total production in field.
    OCP_DBL FGPR{0}; ///< gas production rate in field.
    OCP_DBL FGPt{0}; ///< gas total production in field.
    OCP_DBL FWPR{0}; ///< water production rate in field.
    OCP_DBL FWPT{0}; ///< water total production in field.
};


#endif /* end if __WELLGROUP_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/