/*! \file    AcceleratePVT.cpp
 *  \brief   AcceleratePVT class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "AcceleratePVT.hpp"

/////////////////////////////////////////////////////////////////////
// Skip Stability Analysis
/////////////////////////////////////////////////////////////////////

void SkipStaAnaly::Setup(const OCP_USI& numBulk, const USI& np, const USI& nc)
{
    if (!ifSetup) {
        ifSetup = OCP_TRUE;

        numPhase = np;
        numCom   = nc;

        flag.resize(numBulk);
        P.resize(numBulk);
        T.resize(numBulk);
        minEigen.resize(numBulk);
        zi.resize(numBulk * numCom);

        lflag.resize(numBulk);
        lP.resize(numBulk);
        lT.resize(numBulk);
        lminEigen.resize(numBulk);
        lzi.resize(numBulk * numCom);
    }
}

void SkipStaAnaly::AssignValue(const OCP_USI&         n,
                               const OCP_DBL&         minEigenSkip,
                               const OCP_DBL&         PSkip,
                               const OCP_DBL&         TSkip,
                               const vector<OCP_DBL>& ziSkip)
{
    minEigen[n] = minEigenSkip;
    P[n]        = PSkip;
    T[n]        = TSkip;
    for (USI i = 0; i < numCom; i++) {
        zi[n * numCom + i] = ziSkip[i];
    }
}

OCP_BOOL SkipStaAnaly::IfSkip(const OCP_DBL&         Pin,
                              const OCP_DBL&         Tin,
                              const OCP_DBL&         Ntin,
                              const vector<OCP_DBL>& Niin,
                              const OCP_USI&         n) const
{
    if (flag[n]) {
        if (fabs(1 - P[n] / Pin) >= minEigen[n] / 10) {
            return OCP_FALSE;
        }
        if (fabs(T[n] - Tin) >= minEigen[n] * 10) {
            return OCP_FALSE;
        }
        // OCP_DBL Nt_w = Ntin - Niin[numCom];
        for (USI i = 0; i < numCom; i++) {
            if (fabs(Niin[i] / Ntin - zi[n * numCom + i]) >= minEigen[n] / 10) {
                return OCP_FALSE;
            }
        }
        return OCP_TRUE;
    } else {
        return OCP_FALSE;
    }
}

USI SkipStaAnaly::CalFtypeIMPEC(const OCP_DBL&         Pin,
                                const OCP_DBL&         Tin,
                                const OCP_DBL&         Ntin,
                                const vector<OCP_DBL>& Niin,
                                const OCP_USI&         n)
{
    if (ifUseSkip) {
        if (IfSkip(Pin, Tin, Ntin, Niin, n)) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}

USI SkipStaAnaly::CalFtypeFIM(const OCP_DBL&         Pin,
                              const OCP_DBL&         Tin,
                              const OCP_DBL&         Ntin,
                              const vector<OCP_DBL>& Niin,
                              const OCP_DBL*         S,
                              const USI&             np,
                              const OCP_USI&         n) const
{
    if (ifUseSkip) {
        if (IfSkip(Pin, Tin, Ntin, Niin, n)) {
            return 1;
        } else if (np >= 2) {
            for (USI j = 0; j < numPhase; j++) {
                if (S[j] < 1E-4) {
                    return 0; // phases change (predicted)
                }
            }
            return 2;
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}

void SkipStaAnaly::ResetToLastTimeStep()
{
    flag     = lflag;
    minEigen = lminEigen;
    P        = lP;
    T        = lT;
    zi       = lzi;
}

void SkipStaAnaly::UpdateLastTimeStep()
{
    lflag     = flag;
    lminEigen = minEigen;
    lP        = P;
    lT        = T;
    lzi       = zi;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/