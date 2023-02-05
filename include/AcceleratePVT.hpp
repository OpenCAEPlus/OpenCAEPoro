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

#include "DenseMat.hpp"
#include "OCPConst.hpp"
#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////////////
// Skip Stability Analysis
/////////////////////////////////////////////////////////////////////

class SkipStaAnaly
{

public:
    /// Set ifUseSkip to true or false
    void SetUseSkip(const OCP_BOOL& flag) { ifUseSkip = flag; }
    /// Return ifUseSkip
    OCP_BOOL IfUseSkip() const { return ifUseSkip; }
    /// Allocate memory for SkipStaAnaly term
    void Setup(const OCP_USI& numBulk, const USI& np, const USI& nc);
    /// Set flag for skipping
    void SetFlagSkip(const OCP_USI& n, const OCP_BOOL& flagSkip) { flag[n] = flagSkip; }
    /// Update variables used for determine if skipping will happen
    void AssignValue(const OCP_USI&         n,
                     const OCP_DBL&         minEigenSkip,
                     const OCP_DBL&         PSkip,
                     const OCP_DBL&         TSkip,
                     const vector<OCP_DBL>& ziSkip);
    /// Determine if skipping will happen
    OCP_BOOL IfSkip(const OCP_DBL&         Pin,
                    const OCP_DBL&         Tin,
                    const OCP_DBL&         Ntin,
                    const vector<OCP_DBL>& Niin,
                    const OCP_USI&         n) const;
    /// Calculate the ftype for IMPEC
    USI CalFtypeIMPEC(const OCP_DBL&         Pin,
                      const OCP_DBL&         Tin,
                      const OCP_DBL&         Ntin,
                      const vector<OCP_DBL>& Niin,
                      const OCP_USI&         n);
    /// Calculate the ftype for FIM
    USI CalFtypeFIM(const OCP_DBL&         Pin,
                    const OCP_DBL&         Tin,
                    const OCP_DBL&         Ntin,
                    const vector<OCP_DBL>& Niin,
                    const OCP_DBL*         S,
                    const USI&             np,
                    const OCP_USI&         n) const;
    /// Reset SkipStaAnaly term to last time step
    void ResetToLastTimeStep();
    /// Update SkipStaAnaly term at last time step
    void UpdateLastTimeStep();

protected:
    OCP_BOOL ifSetup{OCP_FALSE}; ///< Only one setup is needed.

    OCP_BOOL ifUseSkip{OCP_TRUE}; ///< If true, then Skip option will be used
    USI      numPhase; ///< Num of phase used in phase equilibrium calculation
    USI      numCom;   ///< Num of components used in phase equilibrium calculation

    vector<OCP_BOOL> flag;     ///< If true, skip will be test
    vector<OCP_DBL>  minEigen; ///< minimum eigenvalue used for testing skipping
    vector<OCP_DBL>  P;        ///< Pressure at last step
    vector<OCP_DBL>  T;        ///< Temperature at last step
    vector<OCP_DBL>  zi;       ///< Mole fraction of components(for test) at last step

    vector<OCP_BOOL> lflag;     ///< Last flag
    vector<OCP_DBL>  lminEigen; ///< Last min eigenvalue
    vector<OCP_DBL>  lP;        ///< Last P
    vector<OCP_DBL>  lT;        ///< Last T
    vector<OCP_DBL>  lzi;       ///< Last zi
};

#endif /* end if __ACCELERATEPVT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/