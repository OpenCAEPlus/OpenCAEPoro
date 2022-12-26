/*! \file    AcceleratePVT.hpp
 *  \brief   AcceleratePVT class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ACCELERATEPVT_HEADER__
#define __ACCELERATEPVT_HEADER__

#include "OCPConst.hpp"
#include "DenseMat.hpp"
#include<vector>


using namespace std;

class SkipStaAnaly
{

    friend class MixtureComp;
public:
    void SetUseSkip(const OCP_BOOL& flag) { ifUseSkip = flag; }
    OCP_BOOL IfUseSkip() const { return ifUseSkip; }

    void Allocate(const OCP_USI& numBulk, const USI& np, const USI& nc);
    void SetFlagSkip(const OCP_USI& n, const OCP_BOOL& flagSkip) { flag[n] = flagSkip; }
    void AssignValue(const OCP_USI& n, const OCP_DBL& minEigenSkip,
        const OCP_DBL& PSkip, const OCP_DBL& TSkip, const vector<OCP_DBL>& ziSkip);

    OCP_BOOL IfSkip(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL& Ntin,
        const vector<OCP_DBL>& Niin, const OCP_USI& n) const;

    USI CalFlashTypeIMPEC(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL& Ntin,
        const vector<OCP_DBL>& Niin, const OCP_USI& n);
    USI CalFlashTypeFIM(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL& Ntin,
        const vector<OCP_DBL>& Niin, const OCP_DBL* S, const USI& np, const OCP_USI& n);

    void ResetToLastTimeStep();
    void UpdateLastTimeStep();

protected:
    OCP_BOOL         ifUseSkip{ OCP_TRUE };  ///< If true, then Skip option will be used
    USI              numPhase;               ///< Num of phase used in phase equilibrium calculation
    USI              numCom;                 ///< Num of componnets used in phase equilibrium calculation
    mutable OCP_USI  bulkId;                 ///< Index of work bulk

    vector<OCP_BOOL> flag;
    vector<OCP_DBL> minEigen;
    vector<OCP_DBL> P;
    vector<OCP_DBL> T;
    vector<OCP_DBL> zi;

    vector<OCP_BOOL> lflag;
    vector<OCP_DBL> lminEigen;
    vector<OCP_DBL> lP;
    vector<OCP_DBL> lT;
    vector<OCP_DBL> lzi;
};


#endif /* end if __ACCELERATEPVT_HEADER__ */


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/