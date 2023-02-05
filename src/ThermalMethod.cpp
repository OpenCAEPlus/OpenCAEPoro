/*! \file    ThermalMethod.cpp
 *  \brief   ThermalMethod class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ThermalMethod.hpp"

void T_FIM::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    AllocateReservoir(rs);
    AllocateLinearSystem(ls, rs, ctrl);
}

void T_FIM::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalThermalConduct(rs.conn, rs.bulk);

    rs.allWells.InitBHP(rs.bulk);
    UpdateLastTimeStep(rs);
}

void T_FIM::Prepare(Reservoir& rs, const OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalRes(rs, ctrl.GetCurTime() + ctrl.GetCurDt(), ctrl.GetCurDt(), OCP_TRUE);
}

void T_FIM::AssembleMat(LinearSystem&    ls,
                        const Reservoir& rs,
                        const OCP_DBL&   t,
                        const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, t, dt);
    AssembleMatWells(ls, rs, dt);
    ls.AssembleRhsCopy(rs.bulk.res.resAbs);
}

void T_FIM::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl) const
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    // Assemble external linear solver with internal A and b
    ls.AssembleMatLinearSolver();
    // Solve linear system
    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }
    // Record time, iterations
    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

#ifdef DEBUG
    // Output A, b, x
    // ls.OutputLinearSystem("testA_FIMT.out", "testb_FIMT.out");
    // ls.OutputSolution("testx_FIMT.out");
    // Check if inf or nan occurs in solution
    ls.CheckSolution();
#endif // DEBUG

    GetSolution(rs, ls.GetSolution(), ctrl);
    // rs.GetSolution01FIM(ls.GetSolution());
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

OCP_BOOL T_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    if (!ctrl.Check(rs, {"BulkNi", "BulkP", "BulkT"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    CalRock(rs.bulk);

    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalThermalConduct(rs.conn, rs.bulk);
    CalHeatLoss(rs.bulk, ctrl.GetCurTime() + dt, dt);

    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, ctrl.GetCurTime() + dt, dt, OCP_FALSE);

    return OCP_TRUE;
}

OCP_BOOL T_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

    cout << "### NR : " + to_string(ctrl.iterNR) + "    Res:    " << setprecision(2)
         << scientific << rs.bulk.res.maxRelRes0_V << setw(12)
         << rs.bulk.res.maxRelRes_V << setw(12) << rs.bulk.res.maxRelRes_N << setw(12)
         << rs.bulk.res.maxWellRelRes_mol << setw(12) << rs.bulk.res.maxRelRes_E
         << setw(12) << fabs(NRdPmax) << setw(12) << fabs(NRdSmax) << endl;

    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
         rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol &&
         rs.bulk.res.maxRelRes_E <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        if (!ctrl.Check(rs, {"WellP"})) {
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            return OCP_TRUE;
        }

    } else if (ctrl.iterNR >= ctrl.ctrlNR.maxNRiter) {
        ctrl.current_dt *= ctrl.ctrlTime.cutFacNR;
        ResetToLastTimeStep(rs, ctrl);
        cout << "### WARNING: NR not fully converged! Cut time step size and repeat!  "
                "current dt = "
             << fixed << setprecision(3) << ctrl.current_dt << " days\n";
        return OCP_FALSE;
    } else {
        return OCP_FALSE;
    }
}

void T_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
    ctrl.UpdateIters();
}

void T_FIM::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);
    bk.vr.resize(nb);
    bk.Hr.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);
    bk.lvr.resize(nb);
    bk.lHr.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.poroT.resize(nb);
    bk.vrP.resize(nb);
    bk.vrT.resize(nb);
    bk.HrT.resize(nb);

    bk.lporoP.resize(nb);
    bk.lporoT.resize(nb);
    bk.lvrP.resize(nb);
    bk.lvrT.resize(nb);
    bk.lHrT.resize(nb);

    // Fluid
    bk.phaseNum.resize(nb);
    bk.Nt.resize(nb);
    bk.Ni.resize(nb * nc);
    bk.vf.resize(nb);
    bk.T.resize(nb);
    bk.P.resize(nb);
    bk.Pb.resize(nb);
    bk.Pj.resize(nb * np);
    bk.Pc.resize(nb * np);
    bk.phaseExist.resize(nb * np);
    bk.S.resize(nb * np);
    bk.xij.resize(nb * np * nc);
    bk.rho.resize(nb * np);
    bk.xi.resize(nb * np);
    bk.mu.resize(nb * np);
    bk.kr.resize(nb * np);
    bk.Uf.resize(nb);
    bk.H.resize(nb * np);
    bk.kt.resize(nb);

    bk.lphaseNum.resize(nb);
    bk.lNt.resize(nb);
    bk.lNi.resize(nb * nc);
    bk.lvf.resize(nb);
    bk.lT.resize(nb);
    bk.lP.resize(nb);
    bk.lPj.resize(nb * np);
    bk.lPc.resize(nb * np);
    bk.lphaseExist.resize(nb * np);
    bk.lS.resize(nb * np);
    bk.lxij.resize(nb * np * nc);
    bk.lrho.resize(nb * np);
    bk.lxi.resize(nb * np);
    bk.lmu.resize(nb * np);
    bk.lkr.resize(nb * np);
    bk.lUf.resize(nb);
    bk.lH.resize(nb * np);
    bk.lkt.resize(nb);

    // derivatives
    bk.vfP.resize(nb);
    bk.vfT.resize(nb);
    bk.vfi.resize(nb * nc);
    bk.rhoP.resize(nb * np);
    bk.rhoT.resize(nb * np);
    bk.rhox.resize(nb * nc * np);
    bk.xiP.resize(nb * np);
    bk.xiT.resize(nb * np);
    bk.xix.resize(nb * nc * np);
    bk.muP.resize(nb * np);
    bk.muT.resize(nb * np);
    bk.mux.resize(nb * nc * np);
    bk.dPcj_dS.resize(nb * np * np);
    bk.dKr_dS.resize(nb * np * np);
    bk.UfP.resize(nb);
    bk.UfT.resize(nb);
    bk.Ufi.resize(nb * nc);
    bk.HT.resize(nb * np);
    bk.Hx.resize(nb * np * nc);
    bk.ktP.resize(nb);
    bk.ktT.resize(nb);
    bk.ktS.resize(nb * np);

    bk.lvfP.resize(nb);
    bk.lvfT.resize(nb);
    bk.lvfi.resize(nb * nc);
    bk.lrhoP.resize(nb * np);
    bk.lrhoT.resize(nb * np);
    bk.lrhox.resize(nb * nc * np);
    bk.lxiP.resize(nb * np);
    bk.lxiT.resize(nb * np);
    bk.lxix.resize(nb * nc * np);
    bk.lmuP.resize(nb * np);
    bk.lmuT.resize(nb * np);
    bk.lmux.resize(nb * nc * np);
    bk.ldPcj_dS.resize(nb * np * np);
    bk.ldKr_dS.resize(nb * np * np);
    bk.UfP.resize(nb);
    bk.UfT.resize(nb);
    bk.Ufi.resize(nb * nc);
    bk.HT.resize(nb * np);
    bk.Hx.resize(nb * np * nc);
    bk.ktP.resize(nb);
    bk.ktT.resize(nb);
    bk.ktS.resize(nb * np);

    // FIM-Specified
    bk.maxLendSdP = (nc + 2) * (nc + 1) * np;
    bk.dSec_dPri.resize(nb * bk.maxLendSdP);
    bk.bRowSizedSdP.resize(nb);
    bk.pSderExist.resize(nb * np);
    bk.pVnumCom.resize(nb * np);

    bk.ldSec_dPri.resize(nb * bk.maxLendSdP);
    bk.lbRowSizedSdP.resize(nb);
    bk.lpSderExist.resize(nb * np);
    bk.lpVnumCom.resize(nb * np);

    // Allocate Residual
    bk.res.SetupT(nb, rs.allWells.numWell, nc);

    // NR
    bk.NRstep.resize(nb);
    bk.NRphaseNum.resize(nb);
    bk.dSNR.resize(nb * np);
    bk.dSNRP.resize(nb * np);
    bk.dNNR.resize(nb * nc);
    bk.dPNR.resize(nb);
    bk.dTNR.resize(nb);

    // BulkConn
    BulkConn&      conn    = rs.conn;
    const OCP_USI& numConn = conn.numConn;

    conn.upblock.resize(numConn * np);
    conn.upblock_Rho.resize(numConn * np);
    conn.upblock_Velocity.resize(numConn * np);
    conn.Adkt.resize(numConn);
    conn.AdktP.resize(numConn * 2);
    conn.AdktT.resize(numConn * 2);
    conn.AdktS.resize(numConn * 2 * np);
}

void T_FIM::AllocateLinearSystem(LinearSystem&     ls,
                                 const Reservoir&  rs,
                                 const OCPControl& ctrl)
{
    ls.AllocateRowMem(rs.GetBulkNum() + rs.GetWellNum(), rs.GetComNum() + 2);
    ls.AllocateColMem(rs.conn.GetNeighborNum(), rs.allWells.GetWell2Bulk());
    ls.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void T_FIM::InitRock(Bulk& bk) const
{

    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        if (bk.bType[n] == 0) {
            // non fluid bulk
            bk.poroInit[n] = 0;
            bk.poro[n]     = 0;
            bk.rockVp[n]   = 0;
            bk.vr[n]       = bk.v[n];
        } else {
            // fluid bulk
            bk.poroInit[n] *= bk.ntg[n];
        }
    }
}

void T_FIM::CalRock(Bulk& bk) const
{
    const OCP_USI nb = bk.numBulk;
    OCP_DBL       vr, vrP, vrT;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // fluid exists
            bk.rock[bk.ROCKNUM[n]]->CalPoroT(bk.P[n], bk.T[n], bk.poroInit[n],
                                             bk.poro[n], bk.poroP[n], bk.poroT[n], vr,
                                             vrP, vrT);

            bk.rockVp[n] = bk.v[n] * bk.poro[n];
            bk.vr[n]     = bk.v[n] * vr;
            bk.vrP[n]    = bk.v[n] * vrP;
            bk.vrT[n]    = bk.v[n] * vrT;
        }
        bk.rock[bk.ROCKNUM[n]]->CalRockHT(bk.T[n], bk.Hr[n], bk.HrT[n]);
    }
}

void T_FIM::InitFlash(Bulk& bk) const
{
    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;
    const OCP_USI nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            bk.flashCal[bk.PVTNUM[n]]->InitFlashFIM(bk.P[n], bk.Pb[n], bk.T[n],
                                                    &bk.S[n * np], bk.rockVp[n],
                                                    &bk.Ni[n * nc], n);
            for (USI i = 0; i < nc; i++) {
                bk.Ni[n * nc + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
            }
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::CalFlash(Bulk& bk) const
{
    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;
    const OCP_USI nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                &bk.S[n * np], bk.phaseNum[n],
                                                &bk.xij[n * np * nc], n);
            PassFlashValue(bk, n);
        }
    }
}

void T_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();
    bk.Uf[n]       = bk.flashCal[pvtnum]->GetUf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bk.S[bIdp + j]    = bk.flashCal[pvtnum]->GetS(j);
        bk.dSNR[bIdp + j] = bk.S[bIdp + j] - bk.dSNR[bIdp + j];
        if (bk.phaseExist[bIdp + j]) {
            if (fabs(bk.maxNRdSSP) < fabs(bk.dSNR[bIdp + j] - bk.dSNRP[bIdp + j])) {
                bk.maxNRdSSP       = bk.dSNR[bIdp + j] - bk.dSNRP[bIdp + j];
                bk.index_maxNRdSSP = n;
            }
        }
        bk.phaseExist[bIdp + j] = bk.flashCal[pvtnum]->GetPhaseExist(j);
        if (bk.phaseExist[bIdp + j]) {
            bk.phaseNum[n]++;
            bk.rho[bIdp + j] = bk.flashCal[pvtnum]->GetRho(j);
            bk.xi[bIdp + j]  = bk.flashCal[pvtnum]->GetXi(j);
            bk.mu[bIdp + j]  = bk.flashCal[pvtnum]->GetMu(j);
            bk.H[bIdp + j]   = bk.flashCal[pvtnum]->GetH(j);

            // Derivatives
            bk.rhoP[bIdp + j] = bk.flashCal[pvtnum]->GetRhoP(j);
            bk.rhoT[bIdp + j] = bk.flashCal[pvtnum]->GetRhoT(j);
            bk.xiP[bIdp + j]  = bk.flashCal[pvtnum]->GetXiP(j);
            bk.xiT[bIdp + j]  = bk.flashCal[pvtnum]->GetXiT(j);
            bk.muP[bIdp + j]  = bk.flashCal[pvtnum]->GetMuP(j);
            bk.muT[bIdp + j]  = bk.flashCal[pvtnum]->GetMuT(j);
            bk.HT[bIdp + j]   = bk.flashCal[pvtnum]->GetHT(j);

            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXij(j, i);
                bk.rhox[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetRhoX(j, i);
                bk.xix[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXiX(j, i);
                bk.mux[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetMuX(j, i);
                bk.Hx[bIdp * nc + j * nc + i]   = bk.flashCal[pvtnum]->GetHx(j, i);
            }
        }
    }
    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    bk.vfT[n] = bk.flashCal[pvtnum]->GetVfT();
    bk.UfP[n] = bk.flashCal[pvtnum]->GetUfP();
    bk.UfT[n] = bk.flashCal[pvtnum]->GetUfT();

    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
        bk.Ufi[n * nc + i] = bk.flashCal[pvtnum]->GetUfi(i);
    }

    Dcopy(bk.maxLendSdP, &bk.dSec_dPri[n * bk.maxLendSdP],
          &bk.flashCal[pvtnum]->GetDXsDXp()[0]);
}

void T_FIM::CalKrPc(Bulk& bk) const
{
    const USI& np = bk.numPhase;

    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        if (bk.bType[n] > 0) {
            OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPcDeriv(&bk.S[bId], &bk.kr[bId], &bk.Pc[bId],
                                                &bk.dKr_dS[bId * np],
                                                &bk.dPcj_dS[bId * np], n);
            for (USI j = 0; j < np; j++)
                bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
        }
    }
}

void T_FIM::CalThermalConduct(BulkConn& conn, Bulk& bk) const
{
    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // fluid bulk
            OCP_DBL tmp = 0;
            for (USI j = 0; j < np; j++) {
                tmp += bk.S[n * np + j] * bk.thconp[j];
                bk.ktS[n * np + j] = bk.poro[n] * bk.thconp[j];
            }
            bk.kt[n]  = bk.poro[n] * tmp + (1 - bk.poro[n]) * bk.thconr[n];
            bk.ktP[n] = bk.poroP[n] * (tmp - bk.thconr[n]);
            bk.ktT[n] = bk.poroT[n] * (tmp - bk.thconr[n]);
        } else {
            // non fluid bulk
            bk.kt[n] = bk.thconr[n];
        }
    }

    OCP_USI bId, eId;
    OCP_DBL areaB, areaE, T1, T2;
    OCP_DBL tmpB, tmpE;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        if (bk.bType[bId] > 0 && bk.bType[eId] > 0) {
            // fluid bulk connections
            areaB        = conn.iteratorConn[c].AreaB();
            areaE        = conn.iteratorConn[c].AreaE();
            T1           = bk.kt[bId] * areaB;
            T2           = bk.kt[eId] * areaE;
            conn.Adkt[c] = 1 / (1 / T1 + 1 / T2);

            tmpB                  = pow(conn.Adkt[c], 2) / pow(T1, 2) * areaB;
            tmpE                  = pow(conn.Adkt[c], 2) / pow(T2, 2) * areaE;
            conn.AdktP[c * 2 + 0] = tmpB * bk.ktP[bId];
            conn.AdktP[c * 2 + 1] = tmpE * bk.ktP[eId];
            conn.AdktT[c * 2 + 0] = tmpB * bk.ktT[bId];
            conn.AdktT[c * 2 + 1] = tmpE * bk.ktT[eId];
            for (USI j = 0; j < np; j++) {
                conn.AdktS[c * np * 2 + j]      = tmpB * bk.ktS[bId * np + j];
                conn.AdktS[c * np * 2 + np + j] = tmpE * bk.ktS[eId * np + j];
            }
        }
    }
}

void T_FIM::CalHeatLoss(Bulk& bk, const OCP_DBL& t, const OCP_DBL& dt) const
{
    bk.hLoss.CalHeatLoss(bk.bLocation, bk.T, bk.lT, bk.initT, t, dt);
}

void T_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    Bulk& bk = rs.bulk;

    // Rock
    bk.poro   = bk.lporo;
    bk.rockVp = bk.lrockVp;
    bk.vr     = bk.lvr;
    bk.Hr     = bk.lHr;
    // derivatives
    bk.poroP = bk.lporoP;
    bk.poroT = bk.lporoT;
    bk.vrP   = bk.lvrP;
    bk.vrT   = bk.lvrT;
    bk.HrT   = bk.lHrT;

    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
    bk.T          = bk.lT;
    bk.P          = bk.lP;
    bk.Pj         = bk.lPj;
    bk.Pc         = bk.lPc;
    bk.phaseExist = bk.lphaseExist;
    bk.S          = bk.lS;
    bk.xij        = bk.lxij;
    bk.rho        = bk.lrho;
    bk.xi         = bk.lxi;
    bk.mu         = bk.lmu;
    bk.kr         = bk.lkr;
    bk.Uf         = bk.lUf;
    bk.H          = bk.lH;
    bk.kt         = bk.lkt;
    // derivatives
    bk.vfP       = bk.lvfP;
    bk.vfT       = bk.lvfT;
    bk.vfi       = bk.lvfi;
    bk.rhoP      = bk.lrhoP;
    bk.rhoT      = bk.lrhoT;
    bk.rhox      = bk.lrhox;
    bk.xiP       = bk.lxiP;
    bk.xiT       = bk.lxiT;
    bk.xix       = bk.lxix;
    bk.muP       = bk.lmuP;
    bk.muT       = bk.lmuT;
    bk.mux       = bk.lmux;
    bk.dPcj_dS   = bk.ldPcj_dS;
    bk.dKr_dS    = bk.ldKr_dS;
    bk.UfP       = bk.lUfP;
    bk.UfT       = bk.lUfT;
    bk.Ufi       = bk.lUfi;
    bk.HT        = bk.lHT;
    bk.Hx        = bk.lHx;
    bk.ktP       = bk.lktP;
    bk.ktT       = bk.lktT;
    bk.ktS       = bk.lktS;
    bk.dSec_dPri = bk.ldSec_dPri;

    bk.hLoss.ResetToLastTimeStep();

    BulkConn& conn = rs.conn;

    conn.Adkt  = conn.lAdkt;
    conn.AdktP = conn.lAdktP;
    conn.AdktT = conn.lAdktT;
    conn.AdktS = conn.lAdktS;

    // Wells
    rs.allWells.ResetBHP();
    rs.allWells.CalTrans(bk);
    rs.allWells.CaldG(bk);
    rs.allWells.CalFlux(bk);

    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();

    CalRes(rs, ctrl.GetCurTime() + ctrl.GetCurDt(), ctrl.GetCurDt(), OCP_TRUE);
}

void T_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    // Bulk
    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lrockVp = bk.rockVp;
    bk.lvr     = bk.vr;
    bk.lHr     = bk.Hr;
    // derivatives
    bk.lporoP = bk.poroP;
    bk.lporoT = bk.poroT;
    bk.lvrP   = bk.vrP;
    bk.lvrT   = bk.vrT;
    bk.lHrT   = bk.HrT;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
    bk.lT          = bk.T;
    bk.lP          = bk.P;
    bk.lPj         = bk.Pj;
    bk.lPc         = bk.Pc;
    bk.lphaseExist = bk.phaseExist;
    bk.lS          = bk.S;
    bk.lxij        = bk.xij;
    bk.lrho        = bk.rho;
    bk.lxi         = bk.xi;
    bk.lmu         = bk.mu;
    bk.lkr         = bk.kr;
    bk.lUf         = bk.Uf;
    bk.lH          = bk.H;
    bk.lkt         = bk.kt;
    // derivatives
    bk.lvfP       = bk.vfP;
    bk.lvfT       = bk.vfT;
    bk.lvfi       = bk.vfi;
    bk.lrhoP      = bk.rhoP;
    bk.lrhoT      = bk.rhoT;
    bk.lrhox      = bk.rhox;
    bk.lxiP       = bk.xiP;
    bk.lxiT       = bk.xiT;
    bk.lxix       = bk.xix;
    bk.lmuP       = bk.muP;
    bk.lmuT       = bk.muT;
    bk.lmux       = bk.mux;
    bk.ldPcj_dS   = bk.dPcj_dS;
    bk.ldKr_dS    = bk.dKr_dS;
    bk.lUfP       = bk.UfP;
    bk.lUfT       = bk.UfT;
    bk.lUfi       = bk.Ufi;
    bk.lHT        = bk.HT;
    bk.lHx        = bk.Hx;
    bk.lktP       = bk.ktP;
    bk.lktT       = bk.ktT;
    bk.lktS       = bk.ktS;
    bk.ldSec_dPri = bk.dSec_dPri;

    bk.hLoss.UpdateLastTimeStep();

    BulkConn& conn = rs.conn;

    conn.lAdkt  = conn.Adkt;
    conn.lAdktP = conn.AdktP;
    conn.lAdktT = conn.AdktT;
    conn.lAdktS = conn.AdktS;

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

void T_FIM::CalRes(Reservoir&      rs,
                   const OCP_DBL&  t,
                   const OCP_DBL&  dt,
                   const OCP_BOOL& resetRes0)
{
    const Bulk& bk   = rs.bulk;
    const USI   nb   = bk.numBulk;
    const USI   np   = bk.numPhase;
    const USI   nc   = bk.numCom;
    const USI   len  = nc + 2;
    OCPRes&     Res  = bk.res;
    BulkConn&   conn = rs.conn;

    Res.SetZero();

    // Bulk to Bulk

    OCP_USI bId, eId, uId, bIdb;
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // Fluid bulk
            bId  = n * len;
            bIdb = n * nc;
            // Volume Conservation
            Res.resAbs[bId] = bk.rockVp[n] - bk.vf[n];
            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                Res.resAbs[n * len + 1 + i] = bk.Ni[bIdb + i] - bk.lNi[bIdb + i];
            }
            // Energy Conservation
            Res.resAbs[n * len + nc + 1] =
                (bk.vf[n] * bk.Uf[n] + bk.vr[n] * bk.Hr[n]) -
                (bk.lvf[n] * bk.lUf[n] + bk.lvr[n] * bk.lHr[n]);
        } else {
            // Non fluid bulk
            Res.resAbs[n * len + nc + 1] = bk.vr[n] * bk.Hr[n] - bk.lvr[n] * bk.lHr[n];
        }

        // Heat Loss
        if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
            const OCP_DBL lambda = bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
            const OCP_DBL kappa  = bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
            // dT
            Res.resAbs[n * len + nc + 1] +=
                dt * kappa *
                (2 * (bk.T[n] - bk.initT[n]) / sqrt(lambda * t) - bk.hLoss.p[n]) *
                bk.dx[n] * bk.dy[n];
        }
    }

    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL rho, dP, dT, dNi, Akd;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();

        // Thermal conductive term all exist
        dT = bk.T[bId] - bk.T[eId];
        Res.resAbs[bId * len + 1 + nc] += conn.Adkt[c] * dT * dt;
        Res.resAbs[eId * len + 1 + nc] -= conn.Adkt[c] * dT * dt;

        if (bk.bType[bId] > 0 && bk.bType[eId] > 0) {
            // Fluid Bulk Connection
            Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();

            for (USI j = 0; j < np; j++) {
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;

                const OCP_BOOL exbegin = bk.phaseExist[bId_np_j];
                const OCP_BOOL exend   = bk.phaseExist[eId_np_j];

                if ((exbegin) && (exend)) {
                    rho = (bk.rho[bId_np_j] + bk.rho[eId_np_j]) / 2;
                } else if (exbegin && (!exend)) {
                    rho = bk.rho[bId_np_j];
                } else if ((!exbegin) && (exend)) {
                    rho = bk.rho[eId_np_j];
                } else {
                    conn.upblock[c * np + j]     = bId;
                    conn.upblock_Rho[c * np + j] = 0;
                    continue;
                }

                uId = bId;
                dP  = (bk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * bk.depth[bId]) -
                     (bk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * bk.depth[eId]);
                if (dP < 0) {
                    uId = eId;
                }
                conn.upblock_Rho[c * np + j] = rho;
                conn.upblock[c * np + j]     = uId;
                uId_np_j                     = uId * np + j;

                if (bk.phaseExist[uId_np_j]) {
                    conn.upblock_Velocity[c * np + j] =
                        Akd * bk.kr[uId_np_j] / bk.mu[uId_np_j] * dP;
                } else {
                    conn.upblock_Velocity[c * np + j] = 0;
                    continue;
                }

                const OCP_DBL tmpV =
                    dt * conn.upblock_Velocity[c * np + j] * bk.xi[uId_np_j];

                for (USI i = 0; i < nc; i++) {
                    dNi = tmpV * bk.xij[uId_np_j * nc + i];
                    Res.resAbs[bId * len + 1 + i] += dNi;
                    Res.resAbs[eId * len + 1 + i] -= dNi;
                }
                Res.resAbs[bId * len + 1 + nc] += tmpV * bk.H[uId_np_j];
                Res.resAbs[eId * len + 1 + nc] -= tmpV * bk.H[uId_np_j];
            }
        }
    }

    // Well to Bulk
    USI wId = nb * len;
    for (const auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            // Well to Bulk

            for (USI p = 0; p < wl.PerfNum(); p++) {
                const OCP_USI k = wl.PerfLocation(p);
                // Mass Conservation
                for (USI i = 0; i < nc; i++) {
                    Res.resAbs[k * len + 1 + i] += wl.PerfQi_lbmol(p, i) * dt;
                }
            }

            if (wl.WellType() == INJ) {
                // Energy Conservation -- Well to Bulk
                const OCP_DBL Hw =
                    bk.flashCal[0]->CalInjWellEnthalpy(wl.InjTemp(), &wl.InjZi()[0]);
                for (USI p = 0; p < wl.PerfNum(); p++) {
                    const OCP_USI k = wl.PerfLocation(p);
                    Res.resAbs[k * len + 1 + nc] +=
                        wl.PerfInjQt_ft3(p) * wl.PerfXi(p) * Hw * dt;
                }

                // Well Self --  Injection
                switch (wl.OptMode()) {
                    case BHP_MODE:
                        Res.resAbs[wId] = wl.BHP() - wl.MaxBHP();
                        break;
                    case RATE_MODE:
                    case ORATE_MODE:
                    case GRATE_MODE:
                    case WRATE_MODE:
                    case LRATE_MODE:
                        Res.resAbs[wId] = wl.MaxRate();
                        for (USI i = 0; i < nc; i++) {
                            Res.resAbs[wId] += wl.Qi_lbmol(i);
                        }
                        // if (opt.reInj) {
                        //     for (auto& w : opt.connWell) {
                        //         OCP_DBL tmp = 0;
                        //         for (USI i = 0; i < numCom; i++) {
                        //             tmp += allWell[w].qi_lbmol[i];
                        //         }
                        //         tmp *= opt.reInjFactor;
                        //         Res.resAbs[bId] += tmp;
                        //     }
                        // }
                        Res.maxWellRelRes_mol =
                            max(Res.maxWellRelRes_mol,
                                fabs(Res.resAbs[wId] / wl.MaxRate()));
                        break;
                    default:
                        OCP_ABORT("Wrong well opt mode!");
                        break;
                }
            } else {
                // Energy Conservation -- Bulk to Well
                for (USI p = 0; p < wl.PerfNum(); p++) {
                    const OCP_USI k = wl.PerfLocation(p);
                    // Mass Conservation
                    for (USI j = 0; j < np; j++) {
                        Res.resAbs[k * len + 1 + nc] += wl.PerfProdQj_ft3(p, j) *
                                                        bk.xi[k * np + j] *
                                                        bk.H[k * np + j] * dt;
                    }
                }

                // Well Self --  Production
                switch (wl.OptMode()) {
                    case BHP_MODE:
                        Res.resAbs[wId] = wl.BHP() - wl.MinBHP();
                        break;
                    case RATE_MODE:
                    case ORATE_MODE:
                    case GRATE_MODE:
                    case WRATE_MODE:
                    case LRATE_MODE:
                        wl.CalProdWeight(bk);
                        Res.resAbs[wId] = -wl.MaxRate();
                        for (USI i = 0; i < nc; i++) {
                            Res.resAbs[wId] += wl.Qi_lbmol(i) * wl.ProdWeight(i);
                        }
                        Res.maxWellRelRes_mol =
                            max(Res.maxWellRelRes_mol,
                                fabs(Res.resAbs[wId] / wl.MaxRate()));
                        break;
                    default:
                        OCP_ABORT("Wrong well opt mode!");
                        break;
                }
            }
            wId += len;
        }
    }

    // Calculate RelRes
    OCP_DBL tmp;
    for (OCP_USI n = 0; n < nb; n++) {

        // Energy equations always exist
        OCP_DBL eT = bk.vf[n] * bk.Uf[n] + bk.vr[n] * bk.Hr[n];
        tmp        = fabs(Res.resAbs[n * len + nc + 1] / eT);
        if (Res.maxRelRes_E < tmp) {
            Res.maxRelRes_E = tmp;
            Res.maxId_E     = n;
        }

        if (bk.bType[n] > 0) {
            // Fluid Bulk
            for (USI i = 0; i < len - 1; i++) {
                tmp = fabs(Res.resAbs[n * len + i] / bk.rockVp[n]);
                if (Res.maxRelRes_V < tmp) {
                    Res.maxRelRes_V = tmp;
                    Res.maxId_V     = n;
                }
            }

            for (USI i = 1; i < len - 1; i++) {
                tmp = fabs(Res.resAbs[n * len + i] / bk.Nt[n]);
                if (Res.maxRelRes_N < tmp) {
                    Res.maxRelRes_N = tmp;
                    Res.maxId_N     = n;
                }
            }
        }
    }

    Dscalar(Res.resAbs.size(), -1.0, Res.resAbs.data());
    if (resetRes0) Res.SetInitRes();
}

void T_FIM::AssembleMatBulks(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   t,
                             const OCP_DBL&   dt) const
{

    const Bulk&     bk     = rs.bulk;
    const BulkConn& conn   = rs.conn;
    const OCP_USI   nb     = bk.numBulk;
    const USI       np     = bk.numPhase;
    const USI       nc     = bk.numCom;
    const USI       ncol   = nc + 2;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        // Mass consevation
        bmat[i * ncol + i] = 1;
    }
    vector<OCP_DBL> bmatR(bmat);
    bmatR[0] = 1;
    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bType[n] > 0) {
            // Fluid Bulk
            // Volume consevation
            // dP
            bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
            // dNi
            for (USI i = 0; i < nc; i++) {
                bmat[i + 1] = -bk.vfi[n * nc + i];
            }
            // dT
            bmat[nc + 1] = bk.v[n] * bk.poroT[n] - bk.vfT[n];
            // Energy consevation
            // dP
            bmat[ncol * (ncol - 1)] =
                bk.vfP[n] * bk.Uf[n] + bk.vf[n] * bk.UfP[n] + bk.vrP[n] * bk.Hr[n];
            // dNi
            for (USI i = 0; i < nc; i++) {
                bmat[ncol * (ncol - 1) + i + 1] =
                    bk.vfi[n * nc + i] * bk.Uf[n] + bk.vf[n] * bk.Ufi[n * nc + i];
            }
            // dT
            bmat[ncol * ncol - 1] = bk.vfT[n] * bk.Uf[n] + bk.vf[n] * bk.UfT[n] +
                                    bk.vrT[n] * bk.Hr[n] + bk.vr[n] * bk.HrT[n];

            // Heat Loss iterm
            if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
                const OCP_DBL lambda =
                    bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
                const OCP_DBL kappa =
                    bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
                // dT
                bmat[ncol * ncol - 1] += dt * kappa *
                                         (2 / sqrt(lambda * t) - bk.hLoss.pT[n]) *
                                         bk.dx[n] * bk.dy[n];
            }

            ls.NewDiag(n, bmat);
        } else {
            // Non Fluid Bulk
            // Energy consevation
            // dT
            bmatR[ncol * ncol - 1] = bk.vrT[n] * bk.Hr[n] + bk.vr[n] * bk.HrT[n];

            // Heat Loss iterm
            if (bk.hLoss.IfHeatLoss() && bk.bLocation[n] > 0) {
                const OCP_DBL lambda =
                    bk.bLocation[n] == 1 ? bk.hLoss.obD : bk.hLoss.ubD;
                const OCP_DBL kappa =
                    bk.bLocation[n] == 1 ? bk.hLoss.obK : bk.hLoss.ubK;
                // dT
                bmatR[ncol * ncol - 1] += dt * kappa *
                                          (2 / sqrt(lambda * t) - bk.hLoss.pT[n]) *
                                          bk.dx[n] * bk.dy[n];
            }

            ls.NewDiag(n, bmatR);
        }
    }

    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ, transH;
    vector<OCP_DBL> dFdXpB(bsize, 0);  // begin bulk: dF / dXp
    vector<OCP_DBL> dFdXpE(bsize, 0);  // end   bulk: dF / dXp
    vector<OCP_DBL> dFdXsB(bsize2, 0); // begin bulk: dF / dXs
    vector<OCP_DBL> dFdXsE(bsize2, 0); // end   bulk: dF / dXs
    OCP_DBL*        dFdXpU;            // up    bulk: dF / dXp
    OCP_DBL*        dFdXpD;            // down  bulk: dF / dXp
    OCP_DBL*        dFdXsU;            // up    bulk: dF / dXs
    OCP_DBL*        dFdXsD;            // down  bulk: dF / dXs

    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j, uId_np_j, dId_np_j;
    OCP_BOOL phaseExistBj, phaseExistEj, phaseExistDj;
    OCP_DBL  xi, xij, kr, mu, rhox, xiP, xiT, xix, muP, muT, mux, H, HT, Hx;
    OCP_DBL  dP, dT, dGamma;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  tmp;

    for (OCP_USI c = 0; c < conn.numConn; c++) {

        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();

        if (bk.bType[bId] > 0 && bk.bType[eId] > 0) {
            // Fluid Bulk Connection
            Akd    = CONV1 * CONV2 * conn.iteratorConn[c].Area();
            dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

            for (USI j = 0; j < np; j++) {
                uId      = conn.upblock[c * np + j];
                uId_np_j = uId * np + j;
                if (!bk.phaseExist[uId_np_j]) continue;
                bId_np_j     = bId * np + j;
                eId_np_j     = eId * np + j;
                phaseExistBj = bk.phaseExist[bId_np_j];
                phaseExistEj = bk.phaseExist[eId_np_j];

                if (bId == uId) {
                    dFdXpU       = &dFdXpB[0];
                    dFdXpD       = &dFdXpE[0];
                    dFdXsU       = &dFdXsB[0];
                    dFdXsD       = &dFdXsE[0];
                    phaseExistDj = phaseExistEj;
                    dId_np_j     = eId_np_j;
                } else {
                    dFdXpU       = &dFdXpE[0];
                    dFdXpD       = &dFdXpB[0];
                    dFdXsU       = &dFdXsE[0];
                    dFdXsD       = &dFdXsB[0];
                    phaseExistDj = phaseExistBj;
                    dId_np_j     = bId_np_j;
                }
                if (phaseExistDj) {
                    rhoWghtU = 0.5;
                    rhoWghtD = 0.5;
                } else {
                    rhoWghtU = 1;
                    rhoWghtD = 0;
                }

                dP = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                     conn.upblock_Rho[c * np + j] * dGamma;
                dT     = bk.T[bId] - bk.T[eId];
                xi     = bk.xi[uId_np_j];
                kr     = bk.kr[uId_np_j];
                mu     = bk.mu[uId_np_j];
                xiP    = bk.xiP[uId_np_j];
                xiT    = bk.xiT[uId_np_j];
                muP    = bk.muP[uId_np_j];
                muT    = bk.muT[uId_np_j];
                H      = bk.H[uId_np_j];
                HT     = bk.HT[uId_np_j];
                transJ = Akd * kr / mu;

                // Mass Conservation
                for (USI i = 0; i < nc; i++) {
                    xij     = bk.xij[uId_np_j * nc + i];
                    transIJ = transJ * xi * xij;

                    // dP
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    tmp = transJ * xiP * xij * dP;
                    tmp += -transIJ * muP / mu * dP;
                    dFdXpU[(i + 1) * ncol] +=
                        (tmp - transIJ * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
                    dFdXpD[(i + 1) * ncol] +=
                        -transIJ * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

                    // dT
                    tmp = transJ * xiT * xij * dP;
                    tmp += -transIJ * muT / mu * dP;
                    dFdXpU[(i + 2) * ncol - 1] +=
                        (tmp - transIJ * rhoWghtU * bk.rhoT[uId_np_j] * dGamma);
                    dFdXpD[(i + 2) * ncol - 1] +=
                        -transIJ * rhoWghtD * bk.rhoT[dId_np_j] * dGamma;

                    // dS
                    for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * bk.dPcj_dS[bId_np_j * np + k];
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * bk.dPcj_dS[eId_np_j * np + k];
                        dFdXsU[(i + 1) * ncol2 + k] +=
                            Akd * bk.dKr_dS[uId_np_j * np + k] / mu * xi * xij * dP;
                    }
                    // dxij
                    for (USI k = 0; k < nc; k++) {
                        rhox = bk.rhox[uId_np_j * nc + k];
                        xix  = bk.xix[uId_np_j * nc + k];
                        mux  = bk.mux[uId_np_j * nc + k];
                        tmp  = -transIJ * rhoWghtU * rhox * dGamma;
                        tmp += transJ * xix * xij * dP;
                        tmp += -transIJ * mux / mu * dP;
                        dFdXsU[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        dFdXsD[(i + 1) * ncol2 + np + j * nc + k] +=
                            -transIJ * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
                    }
                    dFdXsU[(i + 1) * ncol2 + np + j * nc + i] += transJ * xi * dP;
                }

                // Energy Conservation
                transH = transJ * xi * H;
                // dP
                dFdXpB[(ncol - 1) * ncol] += transH;
                dFdXpE[(ncol - 1) * ncol] -= transH;

                tmp = transJ * xiP * H * dP;
                tmp += -transJ * xi * muP / mu * dP * H;
                dFdXpU[(ncol - 1) * ncol] +=
                    (tmp - transH * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
                dFdXpD[(ncol - 1) * ncol] +=
                    -transH * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

                // dT
                tmp = transJ * xiT * H * dP;
                tmp += transJ * xi * HT * dP;
                tmp += -transH * muT / mu * dP;
                dFdXpU[ncol * ncol - 1] +=
                    (tmp - transH * rhoWghtU * bk.rhoT[uId_np_j] * dGamma);
                dFdXpD[ncol * ncol - 1] +=
                    -transH * rhoWghtD * bk.rhoT[dId_np_j] * dGamma;

                // dS
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(nc + 1) * ncol2 + k] +=
                        transH * bk.dPcj_dS[bId_np_j * np + k];
                    dFdXsE[(nc + 1) * ncol2 + k] -=
                        transH * bk.dPcj_dS[eId_np_j * np + k];
                    dFdXsU[(nc + 1) * ncol2 + k] +=
                        Akd * bk.dKr_dS[uId_np_j * np + k] / mu * xi * H * dP;
                }
                // dxij
                for (USI k = 0; k < nc; k++) {
                    rhox = bk.rhox[uId_np_j * nc + k];
                    xix  = bk.xix[uId_np_j * nc + k];
                    mux  = bk.mux[uId_np_j * nc + k];
                    Hx   = bk.Hx[uId_np_j * nc + k];
                    tmp  = -transH * rhoWghtU * rhox * dGamma;
                    tmp += transJ * xix * H * dP;
                    tmp += transJ * xi * Hx * dP;
                    tmp += -transH * mux / mu * dP;
                    dFdXsU[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
                    dFdXsD[(nc + 1) * ncol2 + np + j * nc + k] +=
                        -transH * rhoWghtD * bk.rhox[dId_np_j * nc + k] * dGamma;
                }
            }
        }

        // Thermal Conduction always exist
        // dP
        dFdXpB[(ncol - 1) * ncol] += conn.AdktP[c * 2 + 0] * dT;
        dFdXpE[(ncol - 1) * ncol] += conn.AdktP[c * 2 + 1] * dT;
        // dT
        dFdXpB[ncol * ncol - 1] += conn.Adkt[c] + conn.AdktT[c * 2 + 0] * dT;
        dFdXpE[ncol * ncol - 1] += -conn.Adkt[c] + conn.AdktT[c * 2 + 1] * dT;
        // dS
        for (OCP_USI j = 0; j < np; j++) {
            dFdXsB[(nc + 1) * ncol2 + j] += conn.AdktS[c * np + j] * dT;
            dFdXsE[(nc + 1) * ncol2 + j] += conn.AdktS[c * np + np + j] * dT;
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &bk.dSec_dPri[bId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        ls.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &bk.dSec_dPri[eId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        ls.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        ls.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}

void T_FIM::AssembleMatWells(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {

            switch (wl.WellType()) {
                case INJ:
                    AssembleMatInjWells(ls, rs.bulk, wl, dt);
                    break;
                case PROD:
                    AssembleMatProdWells(ls, rs.bulk, wl, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}

void T_FIM::AssembleMatInjWells(LinearSystem&  ls,
                                const Bulk&    bk,
                                const Well&    wl,
                                const OCP_DBL& dt) const
{

    const USI nc     = bk.numCom;
    const USI np     = bk.numPhase;
    const USI ncol   = nc + 2;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, muT, dP, Hw;
    OCP_DBL transJ, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    Hw = bk.flashCal[0]->CalInjWellEnthalpy(wl.InjTemp(), &wl.InjZi()[0]);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bk.P[n] - wl.BHP() - wl.DG(p);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bk.phaseExist[n_np_j]) continue;

            mu  = bk.mu[n_np_j];
            muP = bk.muP[n_np_j];
            muT = bk.muT[n_np_j];

            for (USI i = 0; i < nc; i++) {

                // Mass Conservation
                if (!wl.IfUseUnweight()) {
                    transIJ = wl.PerfTransj(p, j) * wl.PerfXi(p) * wl.InjZi(i);
                    // dQ / dP
                    dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                    dQdXpW[(i + 1) * ncol] += -transIJ;

                    // dQ / dT
                    dQdXpB[(i + 2) * ncol - 1] += transIJ * (-dP * muT / mu);
                    dQdXpW[(i + 2) * ncol - 1] += 0;
                } else {
                    // dQ / dP
                    transIJ = wl.PerfTransInj(p) * wl.PerfXi(p) * wl.InjZi(i);
                    dQdXpB[(i + 1) * ncol] += transIJ;
                    dQdXpW[(i + 1) * ncol] += -transIJ;
                }

                if (!wl.IfUseUnweight()) {
                    // dQ / dS
                    for (USI k = 0; k < np; k++) {
                        dQdXsB[(i + 1) * ncol2 + k] +=
                            CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * wl.PerfXi(p) *
                            wl.InjZi(i) * bk.dKr_dS[n_np_j * np + k] * dP / mu;
                    }
                    // dQ / dxij
                    for (USI k = 0; k < nc; k++) {
                        dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                            -transIJ * dP / mu * bk.mux[n_np_j * nc + k];
                    }
                }
            }

            // Energy Conservation
            if (!wl.IfUseUnweight()) {
                transJ = wl.PerfTransj(p, j) * wl.PerfXi(p);
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw * (1 - dP * muP / mu);
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;

                // dQ / dT
                dQdXpB[(nc + 2) * ncol - 1] += transJ * Hw * (-dP * muT / mu);
                dQdXpW[(nc + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(nc + 1) * ncol2 + k] +=
                        CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * wl.PerfXi(p) *
                        bk.dKr_dS[n_np_j * np + k] * dP / mu * Hw;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] +=
                        -transJ * dP / mu * bk.mux[n_np_j * nc + k] * Hw;
                }
            } else {
                transJ = wl.PerfTransInj(p) * wl.PerfXi(p);
                // dQ / dP
                dQdXpB[(nc + 1) * ncol] += transJ * Hw;
                dQdXpW[(nc + 1) * ncol] += -transJ * Hw;
            }

            if (wl.IfUseUnweight()) break;
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);

        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (wl.OptMode()) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                bmat[ncol * ncol - 1] = 1;
                ls.AddDiag(wId, bmat);

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                ls.AddDiag(wId, bmat);
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void T_FIM::AssembleMatProdWells(LinearSystem&  ls,
                                 const Bulk&    bk,
                                 const Well&    wl,
                                 const OCP_DBL& dt) const
{
    const USI np     = bk.numPhase;
    const USI nc     = bk.numCom;
    const USI ncol   = nc + 2;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, xiP, xiT, mu, muP, muT, dP, transIJ, transJ, H, HT, Hx, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Prod Weight
    if (wl.OptMode() != BHP_MODE) wl.CalProdWeight(bk);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!bk.phaseExist[n_np_j]) continue;

            dP  = bk.Pj[n_np_j] - wl.BHP() - wl.DG(p);
            xi  = bk.xi[n_np_j];
            mu  = bk.mu[n_np_j];
            xiP = bk.xiP[n_np_j];
            xiT = bk.xiT[n_np_j];
            muP = bk.muP[n_np_j];
            muT = bk.muT[n_np_j];
            H   = bk.H[n_np_j];
            HT  = bk.HT[n_np_j];

            // Mass Conservation
            for (USI i = 0; i < nc; i++) {
                xij = bk.xij[n_np_j * nc + i];
                Hx  = bk.Hx[n_np_j * nc + i];
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                          dP * wl.PerfTransj(p, j) * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dT
                dQdXpB[(i + 2) * ncol - 1] +=
                    transIJ * (-dP * muT / mu) + dP * wl.PerfTransj(p, j) * xij * xiT;
                dQdXpW[(i + 2) * ncol - 1] += 0;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * dP / mu * xi *
                          xij * bk.dKr_dS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bk.dPcj_dS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * wl.PerfTransj(p, j) * xij *
                          (bk.xix[n_np_j * nc + k] - xi / mu * bk.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    wl.PerfTransj(p, j) * xi * dP;
            }

            // Energy Conservation
            transJ = wl.PerfTransj(p, j) * xi;
            // dQ / dP
            dQdXpB[(nc + 1) * ncol] +=
                transJ * (1 - dP * muP / mu) * H + dP * wl.PerfTransj(p, j) * xiP * H;
            dQdXpW[(nc + 1) * ncol] += -transJ * H;

            // dQ / dT
            dQdXpB[(nc + 2) * ncol - 1] += transJ * (-dP * muT / mu) * H +
                                           dP * wl.PerfTransj(p, j) * xiT * H +
                                           transJ * dP * HT;
            dQdXpW[(nc + 2) * ncol - 1] += 0;

            // dQ / dS
            for (USI k = 0; k < np; k++) {
                tmp = CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * dP / mu * xi *
                      bk.dKr_dS[n_np_j * np + k] * H;
                // capillary pressure
                tmp += transJ * bk.dPcj_dS[n_np_j * np + k] * H;
                dQdXsB[(nc + 1) * ncol2 + k] += tmp;
            }

            // dQ / dxij
            for (USI k = 0; k < nc; k++) {
                tmp =
                    dP * wl.PerfTransj(p, j) *
                        (bk.xix[n_np_j * nc + k] - xi / mu * bk.mux[n_np_j * nc + k]) *
                        H +
                    transJ * dP * Hx;
                dQdXsB[(nc + 1) * ncol2 + np + j * nc + k] += tmp;
            }
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        ls.AddDiag(n, bmat);
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        ls.NewOffDiag(n, wId, bmat);

        // Well
        switch (wl.OptMode()) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * wl.ProdWeight(i);
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                bmat[ncol * ncol - 1] = 1;
                ls.AddDiag(wId, bmat);

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, wl.ProdWeight(i), bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                ls.AddDiag(wId, bmat);
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void T_FIM::GetSolution(Reservoir&             rs,
                        const vector<OCP_DBL>& u,
                        const OCPControl&      ctrl) const
{

    // Bulk
    const OCP_DBL dSmaxlim = ctrl.ctrlNR.NRdSmax;
    // const OCP_DBL dPmaxlim = ctrl.ctrlNR.NRdPmax;

    Bulk&           bk  = rs.bulk;
    const OCP_USI   nb  = bk.numBulk;
    const USI       np  = bk.numPhase;
    const USI       nc  = bk.numCom;
    const USI       row = np * (nc + 1);
    const USI       col = nc + 2;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;
    bk.NRdTmax    = 0;

    for (OCP_USI n = 0; n < nb; n++) {
        // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

        if (bk.bType[n] > 0) {
            // Fluid Bulk

            chopmin = 1;
            // compute the chop
            fill(dtmp.begin(), dtmp.end(), 0.0);
            DaAxpby(np, col, 1, &bk.dSec_dPri[n * bk.maxLendSdP], u.data() + n * col, 1,
                    dtmp.data());

            for (USI j = 0; j < np; j++) {
                choptmp = 1;
                if (fabs(dtmp[j]) > dSmaxlim) {
                    choptmp = dSmaxlim / fabs(dtmp[j]);
                } else if (bk.S[n * np + j] + dtmp[j] < 0.0) {
                    choptmp = 0.9 * bk.S[n * np + j] / fabs(dtmp[j]);
                }
                chopmin = min(chopmin, choptmp);
            }

            // dS
            for (USI j = 0; j < np; j++) {
                bk.dSNRP[n * np + j] = chopmin * dtmp[j];
                bk.S[n * np + j] += bk.dSNRP[n * np + j];
            }

            // dP
            OCP_DBL dP = u[n * col];
            if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;
            bk.P[n] += dP;
            bk.dPNR[n] = dP;

            // dNi
            bk.NRstep[n] = chopmin;
            for (USI i = 0; i < nc; i++) {
                bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
                if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
                    bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];

                bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];
            }
        }

        // dT
        OCP_DBL dT = u[n * col + col - 1];
        if (fabs(bk.NRdTmax) < fabs(dT)) bk.NRdTmax = dT;
        bk.T[n] += dT;
        bk.dTNR[n] = dT;
    }

    // Well
    OCP_USI wId = nb * col;
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            wl.SetBHP(wl.BHP() + u[wId]);
            wId += col;
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/