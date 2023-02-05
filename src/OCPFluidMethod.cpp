/*! \file    OCPFluidMethod.cpp
 *  \brief   Definition of solution methods for fluid part in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Nov/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPFluidMethod.hpp"

////////////////////////////////////////////
// IsothermalMethod
////////////////////////////////////////////

void IsothermalMethod::InitRock(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.poroInit[n] *= bk.ntg[n];
    }
}

void IsothermalMethod::CalRock(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.rock[bk.ROCKNUM[n]]->CalPoro(bk.P[n], bk.poroInit[n], bk.poro[n],
                                        bk.poroP[n]);
        bk.rockVp[n] = bk.v[n] * bk.poro[n];
    }
}

////////////////////////////////////////////
// IsoT_IMPEC
////////////////////////////////////////////

void IsoT_IMPEC::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Memory of auxiliary variables for IMPEC
    AllocateReservoir(rs);
    // Allocate Memory of Matrix for IMPEC
    AllocateLinearSystem(ls, rs, ctrl);
}

/// Initialize reservoir
void IsoT_IMPEC::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);

    CalBulkFlux(rs);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_IMPEC::Prepare(Reservoir& rs, OCPControl& ctrl)
{
    rs.allWells.PrepareWell(rs.bulk);
    rs.CalCFL(ctrl.GetCurDt());
    ctrl.Check(rs, {"CFL"});
}

void IsoT_IMPEC::AssembleMat(LinearSystem&    ls,
                             const Reservoir& rs,
                             const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
}

void IsoT_IMPEC::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    ls.AssembleMatLinearSolver();

#ifdef DEBUG
    // ls.OutputLinearSystem("testA_IMPEC.out", "testb_IMPEC.out");
#endif // DEBUG

    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

#ifdef DEBUG
    // ls.OutputSolution("testx_IMPEC.out");
#endif // DEBUG

    // rs.GetSolutionIMPEC(ls.GetSolution());
    GetSolution(rs, ls.GetSolution());
    ls.ClearData();
}

OCP_BOOL IsoT_IMPEC::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // First check : Pressure check
    if (!ctrl.Check(rs, {"BulkP", "WellP"})) {
        return OCP_FALSE;
    }

    // Calculate Flux between bulks and between bulks and wells
    CalFlux(rs);

    // Second check : CFL check
    rs.CalCFL(dt);
    if (!ctrl.Check(rs, {"CFL"})) {
        ResetToLastTimeStep01(rs, ctrl);
        cout << "CFL is too big" << endl;
        return OCP_FALSE;
    }

    MassConserve(rs, dt);

    // Third check: Ni check
    if (!ctrl.Check(rs, {"BulkNi"})) {
        ResetToLastTimeStep02(rs, ctrl);
        return OCP_FALSE;
    }

    CalRock(rs.bulk);
    CalFlash(rs.bulk);

    // Fouth check: Volume error check
    if (!ctrl.Check(rs, {"BulkVe"})) {
        ResetToLastTimeStep03(rs, ctrl);
        return OCP_FALSE;
    }

    CalKrPc(rs.bulk);
    CalBulkFlux(rs);

    return OCP_TRUE;
}

OCP_BOOL IsoT_IMPEC::FinishNR(const Reservoir& rs) { return OCP_TRUE; }

void IsoT_IMPEC::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    // ctrl.CalNextTstepIMPEC(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dN", "dS", "eV"});
}

void IsoT_IMPEC::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.lporoP.resize(nb);

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
    bk.vj.resize(nb * np);
    bk.xij.resize(nb * np * nc);
    bk.rho.resize(nb * np);
    bk.xi.resize(nb * np);
    bk.mu.resize(nb * np);
    bk.kr.resize(nb * np);

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
    bk.vj.resize(nb * np);
    bk.lxij.resize(nb * np * nc);
    bk.lrho.resize(nb * np);
    bk.lxi.resize(nb * np);
    bk.lmu.resize(nb * np);
    bk.lkr.resize(nb * np);

    // derivatives
    bk.vfP.resize(nb);
    bk.vfi.resize(nb * nc);

    bk.lvfP.resize(nb);
    bk.lvfi.resize(nb * nc);

    // others
    bk.cfl.resize(nb * np);

    BulkConn& conn = rs.conn;

    conn.upblock.resize(conn.numConn * np);
    conn.upblock_Rho.resize(conn.numConn * np);
    conn.upblock_Trans.resize(conn.numConn * np);
    conn.upblock_Velocity.resize(conn.numConn * np);

    conn.lupblock.resize(conn.numConn * np);
    conn.lupblock_Rho.resize(conn.numConn * np);
    conn.lupblock_Trans.resize(conn.numConn * np);
    conn.lupblock_Velocity.resize(conn.numConn * np);
}

void IsoT_IMPEC::AllocateLinearSystem(LinearSystem&     ls,
                                      const Reservoir&  rs,
                                      const OCPControl& ctrl)
{
    ls.AllocateRowMem(rs.GetBulkNum() + rs.GetWellNum(), 1);
    ls.AllocateColMem(rs.conn.GetNeighborNum(), rs.allWells.GetWell2Bulk());
    ls.SetupLinearSolver(SCALARFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void IsoT_IMPEC::InitFlash(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.flashCal[bk.PVTNUM[n]]->InitFlashIMPEC(bk.P[n], bk.Pb[n], bk.T[n],
                                                  &bk.S[n * bk.numPhase], bk.rockVp[n],
                                                  bk.Ni.data() + n * bk.numCom, n);
        for (USI i = 0; i < bk.numCom; i++) {
            bk.Ni[n * bk.numCom + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::CalFlash(Bulk& bk)
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {

        bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * bk.numCom],
                                              bk.phaseNum[n],
                                              &bk.xij[n * bk.numPhase * bk.numCom], n);
        PassFlashValue(bk, n);
    }
}

void IsoT_IMPEC::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();

    for (USI j = 0; j < np; j++) {
        // Important! Saturation must be passed no matter if the phase exists. This is
        // because it will be used to calculate relative permeability and capillary
        // pressure at each time step. Make sure that all saturations are updated at
        // each step!
        bk.phaseExist[bIdp + j] = bk.flashCal[pvtnum]->GetPhaseExist(j);
        bk.S[bIdp + j]          = bk.flashCal[pvtnum]->GetS(j);
        if (bk.phaseExist[bIdp + j]) {
            bk.phaseNum[n]++;
            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetXij(j, i);
            }
            bk.vj[bIdp + j]  = bk.flashCal[pvtnum]->GetVj(j);
            bk.rho[bIdp + j] = bk.flashCal[pvtnum]->GetRho(j);
            bk.xi[bIdp + j]  = bk.flashCal[pvtnum]->GetXi(j);
            bk.mu[bIdp + j]  = bk.flashCal[pvtnum]->GetMu(j);
        }
    }

    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }
}

void IsoT_IMPEC::CalKrPc(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        OCP_USI bId = n * bk.numPhase;
        bk.flow[bk.SATNUM[n]]->CalKrPc(&bk.S[bId], &bk.kr[bId], &bk.Pc[bId], n);
        for (USI j = 0; j < bk.numPhase; j++)
            bk.Pj[n * bk.numPhase + j] = bk.P[n] + bk.Pc[n * bk.numPhase + j];
    }
}

void IsoT_IMPEC::CalFlux(Reservoir& rs) const
{
    CalBulkFlux(rs);
    rs.allWells.CalFlux(rs.bulk);
}

void IsoT_IMPEC::CalBulkFlux(Reservoir& rs) const
{
    const Bulk& bk   = rs.bulk;
    BulkConn&   conn = rs.conn;
    const USI   np   = bk.numPhase;

    // calculate a step flux using iteratorConn
    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j;
    OCP_BOOL exbegin, exend, exup;
    OCP_DBL  rho, dP, Akd;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            exbegin = bk.phaseExist[bId_np_j];
            exend   = bk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho = (bk.rho[bId_np_j] + bk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho = bk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho = bk.rho[eId_np_j];
            } else {
                conn.upblock[c * np + j] = bId;
                continue;
            }

            dP = (bk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * bk.depth[bId]) -
                 (bk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * bk.depth[eId]);
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            } else {
                uId  = bId;
                exup = exbegin;
            }

            conn.upblock_Rho[c * np + j] = rho;
            conn.upblock[c * np + j]     = uId;

            if (exup) {
                conn.upblock_Trans[c * np + j] =
                    Akd * bk.kr[uId * np + j] / bk.mu[uId * np + j];
                conn.upblock_Velocity[c * np + j] = conn.upblock_Trans[c * np + j] * dP;
            } else {
                conn.upblock_Trans[c * np + j]    = 0;
                conn.upblock_Velocity[c * np + j] = 0;
            }
        }
    }
}

void IsoT_IMPEC::MassConserve(Reservoir& rs, const OCP_DBL& dt) const
{

    // Bulk to Bulk
    Bulk&           bk   = rs.bulk;
    const BulkConn& conn = rs.conn;

    const USI np = bk.numPhase;
    const USI nc = bk.numCom;

    OCP_USI bId, eId, uId;
    OCP_USI uId_np_j;
    OCP_DBL phaseVelocity, dNi;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();

        for (USI j = 0; j < np; j++) {
            uId      = conn.upblock[c * np + j];
            uId_np_j = uId * np + j;

            if (!bk.phaseExist[uId_np_j]) continue;

            phaseVelocity = conn.upblock_Velocity[c * np + j];
            for (USI i = 0; i < nc; i++) {
                dNi = dt * phaseVelocity * bk.xi[uId_np_j] * bk.xij[uId_np_j * nc + i];
                bk.Ni[eId * nc + i] += dNi;
                bk.Ni[bId * nc + i] -= dNi;
            }
        }
    }

    // Well to Bulk
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            for (USI p = 0; p < wl.PerfNum(); p++) {
                OCP_USI k = wl.PerfLocation(p);
                for (USI i = 0; i < nc; i++) {
                    bk.Ni[k * nc + i] -= wl.PerfQi_lbmol(p, i) * dt;
                }
            }
        }
    }
}

void IsoT_IMPEC::AssembleMatBulks(LinearSystem&    ls,
                                  const Reservoir& rs,
                                  const OCP_DBL&   dt) const
{

    const Bulk&     bk   = rs.bulk;
    const BulkConn& conn = rs.conn;
    const OCP_USI   nb   = bk.numBulk;
    const USI       np   = bk.numPhase;
    const USI       nc   = bk.numCom;

    ls.AddDim(nb);

    // accumulate term
    OCP_DBL Vpp, Vp, vf, vfP, P;
    for (OCP_USI n = 0; n < nb; n++) {
        vf  = bk.vf[n];
        vfP = bk.vfP[n];
        P   = bk.lP[n];
        Vpp = bk.v[n] * bk.poroP[n];
        Vp  = bk.rockVp[n];

        ls.NewDiag(n, Vpp - vfP);
        ls.AddRhs(n, (Vpp - vfP) * P + dt * (vf - Vp));
    }

    // flux term
    OCP_USI bId, eId, uId_np_j;
    OCP_DBL valupi, valdowni, valup, rhsup, valdown, rhsdown;
    OCP_DBL dD, tmp;

    // Be careful when first bulk has no neighbors!
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        dD  = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

        valup   = 0;
        rhsup   = 0;
        valdown = 0;
        rhsdown = 0;

        for (USI j = 0; j < np; j++) {
            uId_np_j = conn.upblock[c * np + j] * np + j;
            if (!bk.phaseExist[uId_np_j]) continue;

            valupi   = 0;
            valdowni = 0;

            for (USI i = 0; i < nc; i++) {
                valupi += bk.vfi[bId * nc + i] * bk.xij[uId_np_j * nc + i];
                valdowni += bk.vfi[eId * nc + i] * bk.xij[uId_np_j * nc + i];
            }

            tmp = bk.xi[uId_np_j] * conn.upblock_Trans[c * np + j] * dt;
            valup += tmp * valupi;
            valdown += tmp * valdowni;
            tmp *= conn.upblock_Rho[c * np + j] * dD -
                   (bk.Pc[bId * np + j] - bk.Pc[eId * np + j]);
            rhsup += tmp * valupi;
            rhsdown -= tmp * valdowni;
        }
        ls.AddDiag(bId, valup);
        ls.AddDiag(eId, valdown);
        ls.NewOffDiag(bId, eId, -valup);
        ls.NewOffDiag(eId, bId, -valdown);
        ls.AddRhs(bId, rhsup);
        ls.AddRhs(eId, rhsdown);
    }
}

void IsoT_IMPEC::AssembleMatWells(LinearSystem&    ls,
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

    // for Reinjection
    // for (auto& wG : wellGroup) {
    //    if (wG.reInj) {
    //        for (auto& prod : wellGroup[wG.prodGroup].wIdPROD) {
    //            if (wells[prod].IsOpen()) {
    //                wells[prod].AssembleMatReinjection_IMPEC(myBulk, myLS, dt, wells,
    //                    wG.wIdINJ);
    //            }
    //        }
    //    }
    //}
}

void IsoT_IMPEC::AssembleMatInjWells(LinearSystem&  ls,
                                     const Bulk&    bk,
                                     const Well&    wl,
                                     const OCP_DBL& dt) const
{
    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);

    OCP_DBL Vfi_zi, valb, valw, bb, bw;

    const USI nc = bk.numCom;

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);

        Vfi_zi = 0;
        for (USI i = 0; i < nc; i++) {
            Vfi_zi += bk.vfi[n * nc + i] * wl.InjZi(i);
        }

        valw = dt * wl.PerfXi(p) * wl.PerfTransInj(p);
        bw   = valw * wl.DG(p);
        valb = valw * Vfi_zi;
        bb   = valb * wl.DG(p);

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (wl.OptMode()) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                ls.AddDiag(wId, valw);
                ls.NewOffDiag(wId, n, -valw);
                ls.AddRhs(wId, -bw);
                break;
            case BHP_MODE:
                ls.NewOffDiag(wId, n, 0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (wl.OptMode()) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            ls.AddRhs(wId, dt * wl.MaxRate());
            break;
        case BHP_MODE:
            ls.AddDiag(wId, dt);
            ls.AddRhs(wId, dt * wl.MaxBHP());
            ls.AssignGuess(wId, wl.MaxBHP());
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void IsoT_IMPEC::AssembleMatProdWells(LinearSystem&  ls,
                                      const Bulk&    bk,
                                      const Well&    wl,
                                      const OCP_DBL& dt) const
{

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, 0.0);

    // Set Prod Weight
    if (wl.OptMode() != BHP_MODE) wl.CalProdWeight(bk);

    const USI np = bk.numPhase;
    const USI nc = bk.numCom;

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);

        OCP_DBL valb = 0;
        OCP_DBL bb   = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw   = 0;

        for (USI j = 0; j < np; j++) {
            if (!bk.phaseExist[n * np + j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < nc; i++) {
                tempb += bk.vfi[n * nc + i] * bk.xij[n * np * nc + j * nc + i];
                tempw += wl.ProdWeight(i) * bk.xij[n * np * nc + j * nc + i];
            }
            OCP_DBL trans = dt * wl.PerfTransj(p, j) * bk.xi[n * np + j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = wl.DG(p) - bk.Pc[n * np + j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        ls.AddDiag(n, valb);
        ls.NewOffDiag(n, wId, -valb);
        ls.AddRhs(n, bb);

        // Well to Bulk
        switch (wl.OptMode()) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                ls.AddDiag(wId, -valw);
                ls.NewOffDiag(wId, n, valw);
                ls.AddRhs(wId, bw);
                break;
            case BHP_MODE:
                ls.NewOffDiag(wId, n, 0.0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (wl.OptMode()) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            ls.AddRhs(wId, dt * wl.MaxRate());
            break;
        case BHP_MODE:
            ls.AddDiag(wId, dt);
            ls.AddRhs(wId, dt * wl.MinBHP());
            ls.AssignGuess(wId, wl.MinBHP());
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void IsoT_IMPEC::GetSolution(Reservoir& rs, const vector<OCP_DBL>& u)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;

    // Bulk
    for (OCP_USI n = 0; n < nb; n++) {
        bk.P[n] = u[n];
        for (USI j = 0; j < np; j++) {
            bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
        }
    }

    // Well
    USI wId = nb;
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {
            wl.SetBHP(u[wId]);
            wl.CalPerfP();
            wId++;
        }
    }
}

void IsoT_IMPEC::ResetToLastTimeStep01(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    rs.bulk.Pj = rs.bulk.lPj;
    // Bulk Conn
    rs.conn.upblock          = rs.conn.lupblock;
    rs.conn.upblock_Rho      = rs.conn.lupblock_Rho;
    rs.conn.upblock_Trans    = rs.conn.lupblock_Trans;
    rs.conn.upblock_Velocity = rs.conn.lupblock_Velocity;

    // Iters
    ctrl.ResetIterNRLS();
}

void IsoT_IMPEC::ResetToLastTimeStep02(Reservoir& rs, OCPControl& ctrl)
{
    // Bulk
    rs.bulk.Ni = rs.bulk.lNi;
    rs.bulk.Pj = rs.bulk.lPj;
    // Bulk Conn
    rs.conn.upblock          = rs.conn.lupblock;
    rs.conn.upblock_Rho      = rs.conn.lupblock_Rho;
    rs.conn.upblock_Trans    = rs.conn.lupblock_Trans;
    rs.conn.upblock_Velocity = rs.conn.lupblock_Velocity;

    // Iters
    ctrl.ResetIterNRLS();
}

void IsoT_IMPEC::ResetToLastTimeStep03(Reservoir& rs, OCPControl& ctrl)
{
    Bulk& bk = rs.bulk;
    // Rock
    bk.rockVp = bk.lrockVp;
    bk.poro   = bk.lporo;
    bk.poroP  = bk.lporoP;

    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
    bk.Pj         = bk.lPj;
    bk.phaseExist = bk.lphaseExist;
    bk.S          = bk.lS;
    bk.vj         = bk.lvj;
    bk.xij        = bk.lxij;
    bk.rho        = bk.lrho;
    bk.xi         = bk.lxi;
    bk.mu         = bk.lmu;

    // derivatives
    bk.vfP = bk.lvfP;
    bk.vfi = bk.lvfi;

    // Bulk Conn
    rs.conn.upblock          = rs.conn.lupblock;
    rs.conn.upblock_Rho      = rs.conn.lupblock_Rho;
    rs.conn.upblock_Trans    = rs.conn.lupblock_Trans;
    rs.conn.upblock_Velocity = rs.conn.lupblock_Velocity;

    // Optional Features
    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();
}

void IsoT_IMPEC::UpdateLastTimeStep(Reservoir& rs) const
{

    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lporoP  = bk.poroP;
    bk.lrockVp = bk.rockVp;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
    bk.lP          = bk.P;
    bk.lPj         = bk.Pj;
    bk.lPc         = bk.Pc;
    bk.lphaseExist = bk.phaseExist;
    bk.lS          = bk.S;
    bk.lvj         = bk.vj;
    bk.lxij        = bk.xij;
    bk.lrho        = bk.rho;
    bk.lxi         = bk.xi;
    bk.lmu         = bk.mu;
    bk.lkr         = bk.kr;

    // derivatives
    bk.lvfP = bk.vfP;
    bk.lvfi = bk.vfi;

    BulkConn& conn = rs.conn;

    conn.lupblock          = conn.upblock;
    conn.lupblock_Rho      = conn.upblock_Rho;
    conn.lupblock_Trans    = conn.upblock_Trans;
    conn.lupblock_Velocity = conn.upblock_Velocity;

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

////////////////////////////////////////////
// IsoT_FIM
////////////////////////////////////////////

void IsoT_FIM::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate memory for reservoir
    AllocateReservoir(rs);
    // Allocate memory for linear system
    AllocateLinearSystem(ls, rs, ctrl);
}

void IsoT_FIM::InitReservoir(Reservoir& rs) const
{
    // Calculate initial bulk pressure and temperature and water saturation
    rs.bulk.InitPTSw(50);
    // Initialize rock property
    InitRock(rs.bulk);
    CalRock(rs.bulk);
    // Initialize fluid properties
    InitFlash(rs.bulk);
    CalKrPc(rs.bulk);
    // Initialize well pressure
    rs.allWells.InitBHP(rs.bulk);
    // Update variables at last time step
    UpdateLastTimeStep(rs);
}

void IsoT_FIM::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    // Calculate well property at the beginning of next time step
    rs.allWells.PrepareWell(rs.bulk);
    // Calculate initial residual
    CalRes(rs, dt, OCP_TRUE);
}

void IsoT_FIM::AssembleMat(LinearSystem&    ls,
                           const Reservoir& rs,
                           const OCP_DBL&   dt) const
{
    // Assemble matrix
#ifdef OCP_OLD_FIM
    AssembleMatBulks(ls, rs, dt);
    AssembleMatWells(ls, rs, dt);
    // rs.allWells.AssemblaMatFIM(ls, rs.bulk, dt);
#else
    AssembleMatBulksNew(ls, rs, dt);
    AssembleMatWellsNew(ls, rs, dt);
#endif // OCP_OLD_FIM
    // Assemble rhs -- from residual
    ls.AssembleRhsCopy(rs.bulk.res.resAbs);
}

void IsoT_FIM::SolveLinearSystem(LinearSystem& ls,
                                 Reservoir&    rs,
                                 OCPControl&   ctrl) const
{
#ifdef DEBUG
    // Check if inf or nan occurs in A and b
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
    // ls.OutputLinearSystem("testA_FIM.out", "testb_FIM.out");
    // ls.OutputSolution("testx_FIM.out");
    // Check if inf or nan occurs in solution
    ls.CheckSolution();
#endif // DEBUG

    // Get solution from linear system to Reservoir
    GetSolution(rs, ls.GetSolution(), ctrl);
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

OCP_BOOL IsoT_FIM::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    if (!ctrl.Check(rs, {"BulkNi", "BulkP"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update fluid property
    CalFlash(rs.bulk);
    CalKrPc(rs.bulk);
    // Update rock property
    CalRock(rs.bulk);
    // Update well property
    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);
    // Update residual
    CalRes(rs, dt, OCP_FALSE);

    return OCP_TRUE;
}

OCP_BOOL IsoT_FIM::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
         rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
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

void IsoT_FIM::FinishStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
}

void IsoT_FIM::AllocateReservoir(Reservoir& rs)
{
    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    // Rock
    bk.poro.resize(nb);
    bk.rockVp.resize(nb);

    bk.lporo.resize(nb);
    bk.lrockVp.resize(nb);

    // derivatives
    bk.poroP.resize(nb);
    bk.lporoP.resize(nb);

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

    // derivatives
    bk.vfP.resize(nb);
    bk.vfi.resize(nb * nc);
    bk.rhoP.resize(nb * np);
    bk.rhox.resize(nb * nc * np);
    bk.xiP.resize(nb * np);
    bk.xix.resize(nb * nc * np);
    bk.muP.resize(nb * np);
    bk.mux.resize(nb * nc * np);
    bk.dPcj_dS.resize(nb * np * np);
    bk.dKr_dS.resize(nb * np * np);

    bk.lvfP.resize(nb);
    bk.lvfi.resize(nb * nc);
    bk.lrhoP.resize(nb * np);
    bk.lrhox.resize(nb * nc * np);
    bk.lxiP.resize(nb * np);
    bk.lxix.resize(nb * nc * np);
    bk.lmuP.resize(nb * np);
    bk.lmux.resize(nb * nc * np);
    bk.ldPcj_dS.resize(nb * np * np);
    bk.ldKr_dS.resize(nb * np * np);

    // FIM-Specified
    bk.maxLendSdP = (nc + 1) * (nc + 1) * np;
    bk.dSec_dPri.resize(nb * bk.maxLendSdP);
    bk.bRowSizedSdP.resize(nb);
    bk.pSderExist.resize(nb * np);
    bk.pVnumCom.resize(nb * np);

    bk.ldSec_dPri.resize(nb * bk.maxLendSdP);
    bk.lbRowSizedSdP.resize(nb);
    bk.lpSderExist.resize(nb * np);
    bk.lpVnumCom.resize(nb * np);

    // Allocate Residual
    bk.res.Setup_IsoT(nb, rs.allWells.numWell, nc);

    // NR
    bk.NRstep.resize(nb);
    bk.NRphaseNum.resize(nb);
    bk.dSNR.resize(nb * np);
    bk.dSNRP.resize(nb * np);
    bk.dNNR.resize(nb * nc);
    bk.dPNR.resize(nb);

    // BulkConn
    BulkConn& conn = rs.conn;

    conn.upblock.resize(conn.numConn * np);
    conn.upblock_Rho.resize(conn.numConn * np);
    conn.upblock_Velocity.resize(conn.numConn * np);
}

void IsoT_FIM::AllocateLinearSystem(LinearSystem&     ls,
                                    const Reservoir&  rs,
                                    const OCPControl& ctrl)
{
    ls.AllocateRowMem(rs.GetBulkNum() + rs.GetWellNum(), rs.GetComNum() + 1);
    ls.AllocateColMem(rs.conn.GetNeighborNum(), rs.allWells.GetWell2Bulk());
    ls.SetupLinearSolver(VECTORFASP, ctrl.GetWorkDir(), ctrl.GetLsFile());
}

void IsoT_FIM::InitFlash(Bulk& bk) const
{
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        bk.flashCal[bk.PVTNUM[n]]->InitFlashFIM(bk.P[n], bk.Pb[n], bk.T[n],
                                                &bk.S[n * bk.numPhase], bk.rockVp[n],
                                                bk.Ni.data() + n * bk.numCom, n);
        for (USI i = 0; i < bk.numCom; i++) {
            bk.Ni[n * bk.numCom + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
        }
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::CalFlash(Bulk& bk)
{
    bk.maxNRdSSP       = 0;
    bk.index_maxNRdSSP = 0;

    for (OCP_USI n = 0; n < bk.numBulk; n++) {

        bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * bk.numCom],
                                            &bk.S[n * bk.numPhase], bk.phaseNum[n],
                                            &bk.xij[n * bk.numPhase * bk.numCom], n);
        PassFlashValue(bk, n);
    }
}

void IsoT_FIM::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];
    USI           len    = 0;

    bk.phaseNum[n] = 0;
    bk.Nt[n]       = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]       = bk.flashCal[pvtnum]->GetVf();

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

            // Derivatives
            bk.rhoP[bIdp + j] = bk.flashCal[pvtnum]->GetRhoP(j);
            bk.xiP[bIdp + j]  = bk.flashCal[pvtnum]->GetXiP(j);
            bk.muP[bIdp + j]  = bk.flashCal[pvtnum]->GetMuP(j);

            for (USI i = 0; i < nc; i++) {
                bk.xij[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXij(j, i);
                bk.rhox[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetRhoX(j, i);
                bk.xix[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetXiX(j, i);
                bk.mux[bIdp * nc + j * nc + i]  = bk.flashCal[pvtnum]->GetMuX(j, i);
            }
        }

        bk.pSderExist[bIdp + j] = bk.flashCal[pvtnum]->GetPSderExist(j);
        bk.pVnumCom[bIdp + j]   = bk.flashCal[pvtnum]->GetPVnumCom(j);
        if (bk.pSderExist[bIdp + j]) len++;
        len += bk.pVnumCom[bIdp + j];
    }

    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }

#ifdef OCP_OLD_FIM
    Dcopy(bk.maxLendSdP, &bk.dSec_dPri[n * bk.maxLendSdP],
          &bk.flashCal[pvtnum]->GetDXsDXp()[0]);
#else
    bk.bRowSizedSdP[n] = len;
    len *= (nc + 1);
    Dcopy(len, &bk.dSec_dPri[n * bk.maxLendSdP], &bk.flashCal[pvtnum]->GetDXsDXp()[0]);
#endif // OCP_OLD_FIM
}

void IsoT_FIM::CalKrPc(Bulk& bk) const
{
    const USI& np = bk.numPhase;
    for (OCP_USI n = 0; n < bk.numBulk; n++) {
        const OCP_USI bId = n * np;
        bk.flow[bk.SATNUM[n]]->CalKrPcDeriv(&bk.S[bId], &bk.kr[bId], &bk.Pc[bId],
                                            &bk.dKr_dS[bId * np], &bk.dPcj_dS[bId * np],
                                            n);
        for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
    }
}

void IsoT_FIM::CalRes(Reservoir& rs, const OCP_DBL& dt, const OCP_BOOL& resetRes0) const
{
    const Bulk& bk   = rs.bulk;
    const USI   nb   = bk.numBulk;
    const USI   np   = bk.numPhase;
    const USI   nc   = bk.numCom;
    const USI   len  = nc + 1;
    OCPRes&     Res  = bk.res;
    BulkConn&   conn = rs.conn;

    Res.SetZero();

    // Bulk to Bulk

    OCP_USI bId, eId, uId, bIdb;
    // Accumalation Term
    for (OCP_USI n = 0; n < nb; n++) {
        bId             = n * len;
        bIdb            = n * nc;
        Res.resAbs[bId] = bk.rockVp[n] - bk.vf[n];
        for (USI i = 0; i < nc; i++) {
            Res.resAbs[bId + 1 + i] = bk.Ni[bIdb + i] - bk.lNi[bIdb + i];
        }
    }

    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL rho, dP, dNi, Akd;

    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = bk.phaseExist[bId_np_j];
            OCP_BOOL exend   = bk.phaseExist[eId_np_j];

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

            OCP_DBL tmp =
                dt * Akd * bk.xi[uId_np_j] * bk.kr[uId_np_j] / bk.mu[uId_np_j] * dP;

            for (USI i = 0; i < nc; i++) {
                dNi = tmp * bk.xij[uId_np_j * nc + i];
                Res.resAbs[bId * len + 1 + i] += dNi;
                Res.resAbs[eId * len + 1 + i] -= dNi;
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
                for (USI i = 0; i < nc; i++) {
                    Res.resAbs[k * len + 1 + i] += wl.PerfQi_lbmol(p, i) * dt;
                }
            }
            // Well Self
            if (wl.WellType() == INJ) {
                // Injection
                switch (wl.OptMode()) {
                    case BHP_MODE:
                        // bhp = opt.maxBHP;
                        // Res.resAbs[bId] = bhp - opt.maxBHP;
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
                // Production
                switch (wl.OptMode()) {
                    case BHP_MODE:
                        // bhp = opt.minBHP;
                        // Res.resAbs[bId] = bhp - opt.minBHP;
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

        for (USI i = 0; i < len; i++) {
            tmp = fabs(Res.resAbs[n * len + i] / bk.rockVp[n]);
            if (Res.maxRelRes_V < tmp) {
                Res.maxRelRes_V = tmp;
                Res.maxId_V     = n;
            }
            Res.resRelV[n] += tmp * tmp;
        }
        Res.resRelV[n] = sqrt(Res.resRelV[n]);

        for (USI i = 1; i < len; i++) {
            tmp = fabs(Res.resAbs[n * len + i] / bk.Nt[n]);
            if (Res.maxRelRes_N < tmp) {
                Res.maxRelRes_N = tmp;
                Res.maxId_N     = n;
            }
            Res.resRelN[n] += tmp * tmp;
        }
        Res.resRelN[n] = sqrt(Res.resRelN[n]);
    }

    Dscalar(Res.resAbs.size(), -1.0, Res.resAbs.data());
    if (resetRes0) Res.SetInitRes();
}

void IsoT_FIM::AssembleMatBulks(LinearSystem&    ls,
                                const Reservoir& rs,
                                const OCP_DBL&   dt) const
{
    const Bulk&     bk     = rs.bulk;
    const BulkConn& conn   = rs.conn;
    const OCP_USI   nb     = bk.numBulk;
    const USI       np     = bk.numPhase;
    const USI       nc     = bk.numCom;
    const USI       ncol   = nc + 1;
    const USI       ncol2  = np * nc + np;
    const USI       bsize  = ncol * ncol;
    const USI       bsize2 = ncol * ncol2;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
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
    OCP_DBL  kr, mu, xi, xij, xiP, muP, rhox, xix, mux;
    OCP_DBL  dP, dGamma;
    OCP_DBL  rhoWghtU, rhoWghtD;
    OCP_DBL  tmp;

    for (OCP_USI c = 0; c < conn.numConn; c++) {

        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

        bId    = conn.iteratorConn[c].BId();
        eId    = conn.iteratorConn[c].EId();
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
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            muP    = bk.muP[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // dP
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = transJ * xiP * xij * dP;
                tmp += -transIJ * muP / mu * dP;
                dFdXpU[(i + 1) * ncol] +=
                    (tmp - transIJ * rhoWghtU * bk.rhoP[uId_np_j] * dGamma);
                dFdXpD[(i + 1) * ncol] +=
                    -transIJ * rhoWghtD * bk.rhoP[dId_np_j] * dGamma;

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

void IsoT_FIM::AssembleMatBulksNew(LinearSystem&    ls,
                                   const Reservoir& rs,
                                   const OCP_DBL&   dt) const
{
    const Bulk&     bk      = rs.bulk;
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nb      = bk.numBulk;
    const USI       np      = bk.numPhase;
    const USI       nc      = bk.numCom;
    const USI       ncol    = nc + 1;
    const USI       ncol2   = np * nc + np;
    const USI       bsize   = ncol * ncol;
    const USI       bsize2  = ncol * ncol2;
    const USI       lendSdP = bk.maxLendSdP;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

        USI jxB = 0;
        USI jxE = 0;
        ncolB   = 0;
        ncolE   = 0;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[bId * np + j];
            phaseExistE[j] = bk.phaseExist[eId * np + j];
            phasedS_B[j]   = bk.pSderExist[bId * np + j];
            phasedS_E[j]   = bk.pSderExist[eId * np + j];
            if (phasedS_B[j]) jxB++;
            if (phasedS_E[j]) jxE++;
            pVnumComB[j] = bk.pVnumCom[bId * np + j];
            pVnumComE[j] = bk.pVnumCom[eId * np + j];
            ncolB += pVnumComB[j];
            ncolE += pVnumComE[j];
        }
        ncolB += jxB;
        ncolE += jxE;

        for (USI j = 0; j < np; j++) {
            uId = conn.upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                 conn.upblock_Rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            muP    = bk.muP[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            rhoP   = bk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-bk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-bk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = bk.rhox[bId_np_j * nc + k] / 2;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * bk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = bk.rhox[eId_np_j * nc + k] / 2;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * bk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                }
            }
            jxB += pVnumComB[j];
            jxE += pVnumComE[j];
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &bk.dSec_dPri[bId * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        ls.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &bk.dSec_dPri[eId * lendSdP], 1,
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

void IsoT_FIM::AssembleMatBulksNewS(LinearSystem&    ls,
                                    const Reservoir& rs,
                                    const OCP_DBL&   dt) const
{

    const Bulk&     bk      = rs.bulk;
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nb      = bk.numCom;
    const USI       np      = bk.numPhase;
    const USI       nc      = bk.numCom;
    const USI       ncol    = nc + 1;
    const USI       ncol2   = np * nc + np;
    const USI       bsize   = ncol * ncol;
    const USI       bsize2  = ncol * ncol2;
    const USI       lendSdP = bk.maxLendSdP;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<USI>      pEnumComB(np, 0);
    vector<USI>      pEnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;
    OCP_DBL wghtb, wghte;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

        const USI npB = bk.phaseNum[bId];
        ncolB         = npB;
        const USI npE = bk.phaseNum[eId];
        ncolE         = npE;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[bId * np + j];
            phaseExistE[j] = bk.phaseExist[eId * np + j];
            pEnumComB[j]   = bk.pVnumCom[bId * np + j];
            pEnumComE[j]   = bk.pVnumCom[eId * np + j];
            ncolB += pEnumComB[j];
            ncolE += pEnumComE[j];
        }

        USI jxB = npB;
        USI jxE = npE;
        for (USI j = 0; j < np; j++) {
            uId = conn.upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pEnumComB[j];
                jxE += pEnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                 conn.upblock_Rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            muP    = bk.muP[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            rhoP   = bk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] += transIJ * dGamma *
                                              (-bk.rhoP[bId_np_j] * bk.S[bId_np_j]) /
                                              (bk.S[bId_np_j] + bk.S[eId_np_j]);
                    dFdXpE[(i + 1) * ncol] += transIJ * dGamma *
                                              (-bk.rhoP[eId_np_j] * bk.S[eId_np_j]) /
                                              (bk.S[bId_np_j] + bk.S[eId_np_j]);
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0;
                        wghte = 0;
                        if (j1 == j && phaseExistE[j]) {
                            tmp   = -dGamma / ((bk.S[bId_np_j] + bk.S[eId_np_j]) *
                                             (bk.S[bId_np_j] + bk.S[eId_np_j]));
                            wghtb = tmp * bk.rho[bId_np_j] * bk.S[eId_np_j];
                            wghte = tmp * bk.rho[eId_np_j] * bk.S[bId_np_j];
                        }

                        if (phaseExistB[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * (bk.dPcj_dS[bId_np_j * np + j1] + wghtb);
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phaseExistE[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * (bk.dPcj_dS[eId_np_j * np + j1] + wghte);
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pEnumComB[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    } else {
                        wghtb = bk.S[bId_np_j] / (bk.S[bId_np_j] + bk.S[eId_np_j]);
                        wghte = bk.S[eId_np_j] / (bk.S[bId_np_j] + bk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComB[j]; k++) {
                            rhox = bk.rhox[bId_np_j * nc + k] * wghtb;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * dGamma * bk.rhox[eId_np_j * nc + k] * wghte;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0;
                        wghte = 0;
                        if (j1 == j && phaseExistB[j]) {
                            tmp   = -dGamma / ((bk.S[bId_np_j] + bk.S[eId_np_j]) *
                                             (bk.S[bId_np_j] + bk.S[eId_np_j]));
                            wghtb = tmp * bk.rho[bId_np_j] * bk.S[eId_np_j];
                            wghte = tmp * bk.rho[eId_np_j] * bk.S[bId_np_j];
                        }

                        if (phaseExistB[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * (bk.dPcj_dS[bId_np_j * np + j1] + wghtb);
                            j1SB++;
                        }
                        if (phaseExistE[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * (bk.dPcj_dS[eId_np_j * np + j1] + wghte);
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pEnumComE[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    } else {
                        wghtb = bk.S[bId_np_j] / (bk.S[bId_np_j] + bk.S[eId_np_j]);
                        wghte = bk.S[eId_np_j] / (bk.S[bId_np_j] + bk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComE[j]; k++) {
                            rhox = bk.rhox[eId_np_j * nc + k] * wghte;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * dGamma * bk.rhox[bId_np_j * nc + k] * wghtb;
                        }
                        // WARNING !!!
                        if (i < pEnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                }
            }
            jxB += pEnumComB[j];
            jxE += pEnumComE[j];
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &bk.dSec_dPri[bId * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        ls.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &bk.dSec_dPri[eId * lendSdP], 1,
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

void IsoT_FIM::AssembleMatWells(LinearSystem&    ls,
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

    //// for Reinjection
    // for (auto& wG : rs.allWells.wellGroup) {
    //     if (wG.IfReInj()) {
    //         for (auto& prod : rs.allWells.wellGroup[wG.prodGroup].wIdPROD) {
    //             if (rs.allWells.wells[prod].IsOpen()) {
    //                 rs.allWells.wells[prod].AssembleMatReinjection_FIM(rs.bulk, ls,
    //                 dt, rs.allWells.wells,
    //                     wG.wIdINJ);
    //             }
    //         }
    //     }
    // }
}

void IsoT_FIM::AssembleMatInjWells(LinearSystem&  ls,
                                   const Bulk&    bk,
                                   const Well&    wl,
                                   const OCP_DBL& dt) const
{
    const USI nc     = bk.numCom;
    const USI np     = bk.numPhase;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL mu, muP, dP, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

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

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * wl.PerfXi(p) * wl.InjZi(i);
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

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

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Bulk - Bulk -- add
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void IsoT_FIM::AssembleMatProdWells(LinearSystem&  ls,
                                    const Bulk&    bk,
                                    const Well&    wl,
                                    const OCP_DBL& dt) const
{
    const USI nc     = bk.numCom;
    const USI np     = bk.numPhase;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;

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
            muP = bk.muP[n_np_j];
            xiP = bk.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = bk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                          dP * wl.PerfTransj(p, j) * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * dP / mu * xi *
                          xij * bk.dKr_dS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * bk.dPcj_dS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dCij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * wl.PerfTransj(p, j) * xij *
                          (bk.xix[n_np_j * nc + k] - xi / mu * bk.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] +=
                    wl.PerfTransj(p, j) * xi * dP;
            }
        }

        // Bulk - Bulk -- add
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bk.dSec_dPri[n * bsize2], 1,
                bmat.data());

        Dscalar(bsize, dt, bmat.data());
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * wl.ProdWeight(i);
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void IsoT_FIM::AssembleMatWellsNew(LinearSystem&    ls,
                                   const Reservoir& rs,
                                   const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {

            switch (wl.WellType()) {
                case INJ:
                    AssembleMatInjWellsNew(ls, rs.bulk, wl, dt);
                    break;
                case PROD:
                    AssembleMatProdWellsNew(ls, rs.bulk, wl, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}

void IsoT_FIM::AssembleMatInjWellsNew(LinearSystem&  ls,
                                      const Bulk&    bk,
                                      const Well&    wl,
                                      const OCP_DBL& dt) const
{
    const USI nc      = bk.numCom;
    const USI np      = bk.numPhase;
    const USI lendSdP = bk.maxLendSdP;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;

    vector<OCP_DBL>  bmat(bsize, 0);
    vector<OCP_DBL>  bmat2(bsize, 0);
    vector<OCP_DBL>  dQdXpB(bsize, 0);
    vector<OCP_DBL>  dQdXpW(bsize, 0);
    vector<OCP_DBL>  dQdXsB(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    USI              ncolB;
    OCP_USI          n_np_j;

    OCP_DBL mu, muP, dP, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bk.P[n] - wl.BHP() - wl.DG(p);

        USI jxB = 0;
        ncolB   = 0;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[n * np + j];
            phasedS_B[j]   = bk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = bk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            mu     = bk.mu[n_np_j];
            muP    = bk.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * wl.PerfXi(p) * wl.InjZi(i);
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * wl.PerfXi(p) *
                            wl.InjZi(i) * bk.dKr_dS[n_np_j * np + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * bk.mux[n_np_j * nc + k];
                }
            }
            jxB += pVnumComB[j];
        }

        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Bulk - Bulk -- add
        ls.AddDiag(n, bmat);

        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void IsoT_FIM::AssembleMatProdWellsNew(LinearSystem&  ls,
                                       const Bulk&    bk,
                                       const Well&    wl,
                                       const OCP_DBL& dt) const
{
    const USI nc      = bk.numCom;
    const USI np      = bk.numPhase;
    const USI lendSdP = bk.maxLendSdP;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;

    vector<OCP_DBL>  bmat(bsize, 0);
    vector<OCP_DBL>  bmat2(bsize, 0);
    vector<OCP_DBL>  dQdXpB(bsize, 0);
    vector<OCP_DBL>  dQdXpW(bsize, 0);
    vector<OCP_DBL>  dQdXsB(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    USI              ncolB;
    OCP_USI          n_np_j;

    OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Prod Weight
    if (wl.OptMode() != BHP_MODE) wl.CalProdWeight(bk);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        USI jxB = 0;
        ncolB   = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[n * np + j];
            phasedS_B[j]   = bk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = bk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            dP     = bk.Pj[n_np_j] - wl.BHP() - wl.DG(p);
            xi     = bk.xi[n_np_j];
            mu     = bk.mu[n_np_j];
            muP    = bk.muP[n_np_j];
            xiP    = bk.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = bk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                          dP * wl.PerfTransj(p, j) * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * dP / mu *
                              xi * xij * bk.dKr_dS[n_np_j * np + j1];
                        // capillary pressure
                        tmp += transIJ * bk.dPcj_dS[n_np_j * np + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * wl.PerfTransj(p, j) * xij *
                          (bk.xix[n_np_j * nc + k] - xi / mu * bk.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] += wl.PerfTransj(p, j) * xi * dP;
            }
            jxB += pVnumComB[j];
        }
        // Bulk - Bulk -- add
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        ls.AddDiag(n, bmat);
        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * wl.ProdWeight(i);
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, wl.ProdWeight(i), bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void IsoT_FIM::GetSolution(Reservoir&             rs,
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
    const USI       col = nc + 1;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;

    for (OCP_USI n = 0; n < nb; n++) {
        // const vector<OCP_DBL>& scm = satcm[SATNUM[n]];

        chopmin = 1;
        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);

#ifdef OCP_OLD_FIM
        DaAxpby(row, col, 1, &bk.dSec_dPri[n * bk.maxLendSdP], u.data() + n * col, 1,
                dtmp.data());
        const OCP_BOOL newFIM = OCP_FALSE;
#else
        DaAxpby(bk.bRowSizedSdP[n], col, 1, &bk.dSec_dPri[n * bk.maxLendSdP],
                u.data() + n * col, 1, dtmp.data());
        const OCP_BOOL newFIM = OCP_TRUE;
#endif // OCP_OLD_FIM

        USI js = 0;
        for (USI j = 0; j < np; j++) {
            if (!bk.pSderExist[n * np + j] && newFIM) {
                continue;
            }

            choptmp = 1;
            if (fabs(dtmp[js]) > dSmaxlim) {
                choptmp = dSmaxlim / fabs(dtmp[js]);
            } else if (bk.S[n * np + j] + dtmp[js] < 0.0) {
                choptmp = 0.9 * bk.S[n * np + j] / fabs(dtmp[js]);
            }

            // if (fabs(S[n_np_j] - scm[j]) > TINY &&
            //     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
            //     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));

            chopmin = min(chopmin, choptmp);
            js++;
        }

        // dS
        js = 0;
        for (USI j = 0; j < np; j++) {
            if (!bk.pSderExist[n * np + j] && newFIM) {
                bk.dSNRP[n * np + j] = 0;
                continue;
            }
            bk.dSNRP[n * np + j] = chopmin * dtmp[js];
            bk.S[n * np + j] += bk.dSNRP[n * np + j];
            js++;
        }

        // dxij   ---- Compositional model only
        if (bk.IfUseEoS()) {
            if (bk.phaseNum[n] >= 3) {
                // num of Hydroncarbon phase >= 2
                OCP_USI bId = 0;
                for (USI j = 0; j < 2; j++) {
                    bId = n * np * nc + j * nc;
                    for (USI i = 0; i < bk.numComH; i++) {
                        bk.xij[bId + i] += chopmin * dtmp[js];
                        js++;
                    }
#ifdef OCP_OLD_FIM
                    js++;
#endif // OCP_OLD_FIM
                }
            }
        }

        // dP
        OCP_DBL dP = u[n * col];
        // choptmp = dPmaxlim / fabs(dP);
        // chopmin = min(chopmin, choptmp);
        if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;
        bk.P[n] += dP; // seems better
        bk.dPNR[n] = dP;

        // dNi
        bk.NRstep[n] = chopmin;
        for (USI i = 0; i < nc; i++) {
            bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
            if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
                bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];

            bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

            // if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
            //     bk.Ni[n * nc + i] = 1E-20;
            // }
        }
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

void IsoT_FIM::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    Bulk& bk = rs.bulk;

    // Rock
    bk.poro   = bk.lporo;
    bk.poroP  = bk.lporoP;
    bk.rockVp = bk.lrockVp;
    // Fluid
    bk.phaseNum   = bk.lphaseNum;
    bk.Nt         = bk.lNt;
    bk.Ni         = bk.lNi;
    bk.vf         = bk.lvf;
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
    // derivatives
    bk.vfP     = bk.lvfP;
    bk.vfi     = bk.lvfi;
    bk.rhoP    = bk.lrhoP;
    bk.rhox    = bk.lrhox;
    bk.xiP     = bk.lxiP;
    bk.xix     = bk.lxix;
    bk.muP     = bk.lmuP;
    bk.mux     = bk.lmux;
    bk.dPcj_dS = bk.ldPcj_dS;
    bk.dKr_dS  = bk.ldKr_dS;
    // FIM-Specified
    bk.bRowSizedSdP = bk.lbRowSizedSdP;
    bk.dSec_dPri    = bk.ldSec_dPri;
    bk.pSderExist   = bk.lpSderExist;
    bk.pVnumCom     = bk.lpVnumCom;

    // Wells
    rs.allWells.ResetBHP();
    rs.allWells.CalTrans(bk);
    rs.allWells.CaldG(bk);
    rs.allWells.CalFlux(bk);

    // Optional Features
    rs.optFeatures.ResetToLastTimeStep();

    // Iters
    ctrl.ResetIterNRLS();

    // Residual
    CalRes(rs, ctrl.GetCurDt(), OCP_TRUE);
}

void IsoT_FIM::UpdateLastTimeStep(Reservoir& rs) const
{
    Bulk& bk = rs.bulk;

    // Rock
    bk.lporo   = bk.poro;
    bk.lporoP  = bk.poroP;
    bk.lrockVp = bk.rockVp;

    // Fluid
    bk.lphaseNum   = bk.phaseNum;
    bk.lNt         = bk.Nt;
    bk.lNi         = bk.Ni;
    bk.lvf         = bk.vf;
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

    // derivatives
    bk.lvfP     = bk.vfP;
    bk.lvfi     = bk.vfi;
    bk.lrhoP    = bk.rhoP;
    bk.lrhox    = bk.rhox;
    bk.lxiP     = bk.xiP;
    bk.lxix     = bk.xix;
    bk.lmuP     = bk.muP;
    bk.lmux     = bk.mux;
    bk.ldPcj_dS = bk.dPcj_dS;
    bk.ldKr_dS  = bk.dKr_dS;

    // FIM-Specified
    bk.lbRowSizedSdP = bk.bRowSizedSdP;
    bk.ldSec_dPri    = bk.dSec_dPri;
    bk.lpSderExist   = bk.pSderExist;
    bk.lpVnumCom     = bk.pVnumCom;

    rs.allWells.UpdateLastTimeStepBHP();
    rs.optFeatures.UpdateLastTimeStep();
}

////////////////////////////////////////////
// OCP_FIMn
////////////////////////////////////////////

/// Setup FIMn
void IsoT_FIMn::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    AllocateReservoir(rs);
    // Allocate memory for internal matrix structure
    IsoT_FIM::AllocateLinearSystem(ls, rs, ctrl);
}

void IsoT_FIMn::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    InitFlash(rs.bulk);
    IsoT_FIM::CalKrPc(rs.bulk);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

/// Prepare for Assembling matrix.
void IsoT_FIMn::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    rs.allWells.PrepareWell(rs.bulk);
    IsoT_FIM::CalRes(rs, dt, OCP_TRUE);
}

/// Assemble Matrix
void IsoT_FIMn::AssembleMat(LinearSystem&    ls,
                            const Reservoir& rs,
                            const OCP_DBL&   dt) const
{
    AssembleMatBulksNew(ls, rs, dt);
    AssembleMatWellsNew(ls, rs, dt);
    ls.AssembleRhsAccumulate(rs.bulk.res.resAbs);
}

/// Solve the linear system.
void IsoT_FIMn::SolveLinearSystem(LinearSystem& ls,
                                  Reservoir&    rs,
                                  OCPControl&   ctrl) const
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    ls.AssembleMatLinearSolver();

    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }
    // cout << "LS step = " << status << endl;

#ifdef DEBUG
    // ls.OutputLinearSystem("testA_FIMn.out", "testb_FIMn.out");
    // ls.OutputSolution("testx_FIMn.out");
    // ls.CheckSolution();
#endif // DEBUG

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    GetSolution(rs, ls.GetSolution(), ctrl);
    // rs.PrintSolFIM(ctrl.workDir + "testPNi.out");
    ls.ClearData();
}

/// Update properties of fluids.
OCP_BOOL IsoT_FIMn::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // First check: Ni check and bulk Pressure check
    if (!ctrl.Check(rs, {"BulkNi", "BulkP"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    // Update reservoir properties
    CalFlash(rs.bulk);
    IsoT_FIM::CalKrPc(rs.bulk);
    CalRock(rs.bulk);

    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, dt, OCP_FALSE);

    return OCP_TRUE;
}

/// Finish a Newton-Raphson iteration.
OCP_BOOL IsoT_FIMn::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI dSn;

    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
         rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        if (!ctrl.Check(rs, {"WellP"})) {
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            return OCP_TRUE;
        }

    } else if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
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

/// Finish a time step.
void IsoT_FIMn::FinishStep(Reservoir& rs, OCPControl& ctrl) const
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
}

void IsoT_FIMn::AllocateReservoir(Reservoir& rs)
{
    IsoT_FIM::AllocateReservoir(rs);

    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;
    const OCP_USI nc = bk.numCom;

    bk.nj.resize(nb * np);
    bk.res_n.resize(nb * (np + np * nc));
    bk.resPc.resize(nb);

    bk.lnj.resize(nb * np);
    bk.lres_n.resize(nb * (np + np * nc));
    bk.lresPc.resize(nb);
}

void IsoT_FIMn::InitFlash(Bulk& bk) const
{
    OCP_FUNCNAME;

    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;
    const OCP_USI nc = bk.numCom;

    if (bk.ifComps) {
        for (OCP_USI n = 0; n < nb; n++) {
            bk.flashCal[bk.PVTNUM[n]]->InitFlashFIMn(bk.P[n], bk.Pb[n], bk.T[n],
                                                     &bk.S[n * np], bk.rockVp[n],
                                                     bk.Ni.data() + n * nc, n);
            for (USI i = 0; i < nc; i++) {
                bk.Ni[n * nc + i] = bk.flashCal[bk.PVTNUM[n]]->GetNi(i);
            }
            PassFlashValue(bk, n);
        }
    } else {
        OCP_ABORT("Not Completed in BLKOIL MODEL!");
    }
}

void IsoT_FIMn::CalFlash(Bulk& bk)
{

    const OCP_USI nb = bk.numBulk;
    const OCP_USI np = bk.numPhase;
    const OCP_USI nc = bk.numCom;

    if (bk.ifComps) {
        vector<USI> flagB(np, 0);
        bk.maxNRdSSP       = 0;
        bk.index_maxNRdSSP = 0;

        for (OCP_USI n = 0; n < nb; n++) {

            for (USI j = 0; j < np; j++) flagB[j] = bk.phaseExist[n * np + j];

            bk.flashCal[bk.PVTNUM[n]]->FlashFIMn(
                bk.P[n], bk.T[n], &bk.Ni[n * nc], &bk.S[n * np], &bk.xij[n * np * nc],
                &bk.nj[n * np], &flagB[0], bk.phaseNum[n], n);

            PassFlashValue(bk, n);
        }
    } else {
        OCP_ABORT("Not completed!");
    }
}

void IsoT_FIMn::PassFlashValue(Bulk& bk, const OCP_USI& n) const
{
    IsoT_FIM::PassFlashValue(bk, n);

    const USI np     = bk.numPhase;
    const USI nc     = bk.numCom;
    const USI pvtnum = bk.PVTNUM[n];

    for (USI j = 0; j < bk.numPhase; j++) {
        if (bk.phaseExist[n * np + j]) {
            bk.nj[n * np + j] = bk.flashCal[pvtnum]->GetNj(j);
        }
    }

    Dcopy(bk.bRowSizedSdP[n], &bk.res_n[0] + n * (np * (nc + 1)),
          &bk.flashCal[pvtnum]->GetRes()[0]);
    bk.resPc[n] = bk.flashCal[pvtnum]->GetResPc();
}

void IsoT_FIMn::AssembleMatBulksNew(LinearSystem&    ls,
                                    const Reservoir& rs,
                                    const OCP_DBL&   dt) const
{

    const Bulk&     bk      = rs.bulk;
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nb      = bk.numBulk;
    const USI       np      = bk.numPhase;
    const USI       nc      = bk.numCom;
    const USI       ncol    = nc + 1;
    const USI       ncol2   = np * nc + np;
    const USI       bsize   = ncol * ncol;
    const USI       bsize2  = ncol * ncol2;
    const USI       lendSdP = bk.maxLendSdP;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> tmpb(ncol, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }

        ls.NewDiag(n, bmat);
        ls.AddRhs(n * ncol, -bk.resPc[n]);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId = conn.iteratorConn[c].BId();
        eId = conn.iteratorConn[c].EId();
        Akd = CONV1 * CONV2 * conn.iteratorConn[c].Area();
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);

        const USI npB = bk.phaseNum[bId];
        const USI npE = bk.phaseNum[eId];

        USI jxB = 0;
        USI jxE = 0;
        ncolB   = 0;
        ncolE   = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[bId * np + j];
            phaseExistE[j] = bk.phaseExist[eId * np + j];
            phasedS_B[j]   = bk.pSderExist[bId * np + j];
            phasedS_E[j]   = bk.pSderExist[eId * np + j];
            if (phasedS_B[j]) jxB++;
            if (phasedS_E[j]) jxE++;
            pVnumComB[j] = bk.pVnumCom[bId * np + j];
            pVnumComE[j] = bk.pVnumCom[eId * np + j];
            ncolB += pVnumComB[j];
            ncolE += pVnumComE[j];
        }
        ncolB += jxB;
        ncolE += jxE;

        for (USI j = 0; j < np; j++) {
            uId = conn.upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                 conn.upblock_Rho[c * np + j] * dGamma;
            xi     = bk.xi[uId_np_j];
            kr     = bk.kr[uId_np_j];
            mu     = bk.mu[uId_np_j];
            muP    = bk.muP[uId_np_j];
            xiP    = bk.xiP[uId_np_j];
            rhoP   = bk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = bk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-bk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-bk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / bk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] +=
                                xi * transJ * dP / bk.nj[uId_np_j];
                    } else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = bk.rhox[bId_np_j * nc + k] / 2;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / bk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * bk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] +=
                                xi * transJ * dP / bk.nj[uId_np_j];
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu * bk.dKr_dS[uId_np_j * np + j1] *
                                  dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = bk.rhox[uId_np_j * nc + k];
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / bk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] +=
                                xi * transJ * dP / bk.nj[uId_np_j];
                    } else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = bk.rhox[eId_np_j * nc + k] / 2;
                            xix  = bk.xix[uId_np_j * nc + k];
                            mux  = bk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / bk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * bk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] +=
                                xi * transJ * dP / bk.nj[uId_np_j];
                    }
                }
            }
            jxB += pVnumComB[j];
            jxE += pVnumComE[j];
        }

        // Assemble rhs
        // Begin
        if (npB > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolB, -1.0, dFdXsB.data(), &bk.res_n[bId * ncol2], 1.0,
                    &tmpb[0]);
            ls.AddRhs(bId, tmpb);
        }

        // Assemble mat
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &bk.dSec_dPri[bId * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        ls.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // Assemble rhs
        // End
        if (npE > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolE, -1.0, dFdXsE.data(), &bk.res_n[eId * ncol2], 1.0,
                    &tmpb[0]);
            ls.AddRhs(eId, tmpb);
        }

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &bk.dSec_dPri[eId * lendSdP], 1,
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

void IsoT_FIMn::AssembleMatWellsNew(LinearSystem&    ls,
                                    const Reservoir& rs,
                                    const OCP_DBL&   dt) const
{
    for (auto& wl : rs.allWells.wells) {
        if (wl.IsOpen()) {

            switch (wl.WellType()) {
                case INJ:
                    AssembleMatInjWellsNew(ls, rs.bulk, wl, dt);
                    break;
                case PROD:
                    AssembleMatProdWellsNew(ls, rs.bulk, wl, dt);
                    break;
                default:
                    OCP_ABORT("Wrong well type");
            }
        }
    }
}

void IsoT_FIMn::AssembleMatInjWellsNew(LinearSystem&  ls,
                                       const Bulk&    bk,
                                       const Well&    wl,
                                       const OCP_DBL& dt) const
{

    const USI nc      = bk.numCom;
    const USI np      = bk.numPhase;
    const USI lendSdP = bk.maxLendSdP;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> tmpb(ncol, 0);
    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);
    vector<char>    phaseExistB(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<char>    phasedS_B(np, OCP_FALSE);
    USI             ncolB;

    OCP_DBL mu, muP, dP, transIJ;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        const OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = bk.P[n] - wl.BHP() - wl.DG(p);

        const USI npB = bk.phaseNum[n];
        USI       jxB = 0;
        ncolB         = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[n * np + j];
            phasedS_B[j]   = bk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = bk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            mu     = bk.mu[n_np_j];
            muP    = bk.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * wl.PerfXi(p) * wl.InjZi(i);
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * wl.PerfXi(p) *
                            wl.InjZi(i) * bk.dKr_dS[n_np_j * np + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * bk.mux[n_np_j * nc + k];
                }
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(), &bk.res_n[n * ncol2], 1.0,
                    &tmpb[0]);
            ls.AddRhs(n, tmpb);
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Bulk - Bulk -- add
        ls.AddDiag(n, bmat);
        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

void IsoT_FIMn::AssembleMatProdWellsNew(LinearSystem&  ls,
                                        const Bulk&    bk,
                                        const Well&    wl,
                                        const OCP_DBL& dt) const
{
    const USI nc      = bk.numCom;
    const USI np      = bk.numPhase;
    const USI lendSdP = bk.maxLendSdP;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    OCP_USI   n_np_j;

    vector<OCP_DBL> tmpb(ncol, 0);
    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);
    vector<char>    phaseExistB(np, OCP_FALSE);
    vector<char>    phasedS_B(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    USI             ncolB;

    OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;

    const OCP_USI wId = ls.AddDim(1) - 1;
    ls.NewDiag(wId, bmat);

    // Set Prod Weight
    if (wl.OptMode() != BHP_MODE) wl.CalProdWeight(bk);

    for (USI p = 0; p < wl.PerfNum(); p++) {
        OCP_USI n = wl.PerfLocation(p);
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        const USI npB = bk.phaseNum[n];
        USI       jxB = 0;
        ncolB         = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = bk.phaseExist[n * np + j];
            phasedS_B[j]   = bk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = bk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            dP     = bk.Pj[n_np_j] - wl.BHP() - wl.DG(p);
            xi     = bk.xi[n_np_j];
            mu     = bk.mu[n_np_j];
            muP    = bk.muP[n_np_j];
            xiP    = bk.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = bk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = wl.PerfTransj(p, j) * xi * xij;
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                          dP * wl.PerfTransj(p, j) * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * wl.PerfWI(p) * wl.PerfMultiplier(p) * dP / mu *
                              xi * xij * bk.dKr_dS[n_np_j * np + j1];
                        // capillary pressure
                        tmp += transIJ * bk.dPcj_dS[n_np_j * np + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * wl.PerfTransj(p, j) * xij *
                          (bk.xix[n_np_j * nc + k] - xi / mu * bk.mux[n_np_j * nc + k]);
                    tmp -= transIJ * dP / bk.nj[n_np_j];
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] +=
                        wl.PerfTransj(p, j) * xi * dP / bk.nj[n_np_j];
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(), &bk.res_n[n * ncol2], 1.0,
                    &tmpb[0]);
            ls.AddRhs(n, tmpb);
        }

        // Bulk - Bulk -- add
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        ls.AddDiag(n, bmat);
        // Bulk - Well -- insert
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
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * wl.ProdWeight(i);
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &bk.dSec_dPri[n * lendSdP],
                        1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, wl.ProdWeight(i), bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                ls.NewOffDiag(wId, n, bmat2);
                break;

            case BHP_MODE:
                // Well - Well -- add
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < ncol; i++) {
                    bmat[i * ncol + i] = 1;
                }
                ls.AddDiag(wId, bmat);

                // Well - Bulk -- insert
                fill(bmat.begin(), bmat.end(), 0.0);
                ls.NewOffDiag(wId, n, bmat);
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
}

/// Update P, Ni, BHP after linear system is solved
void IsoT_FIMn::GetSolution(Reservoir&             rs,
                            const vector<OCP_DBL>& u,
                            const OCPControl&      ctrl) const
{
    // For saturations changes:
    // 1. maximum changes must be less than dSmaxlim,
    // 2. if phases become mobile/immobile, then set it to crtical point,
    // 3. After saturations are determined, then scale the nij to conserve Volume
    // equations.

    // Bulk
    const OCP_DBL dSmaxlim = ctrl.ctrlNR.NRdSmax;
    // const OCP_DBL dPmaxlim = ctrl.ctrlNR.NRdPmax;

    Bulk&           bk  = rs.bulk;
    const OCP_USI   nb  = bk.numBulk;
    const USI       np  = bk.numPhase;
    const USI       nc  = bk.numCom;
    const USI       row = np * (nc + 1);
    const USI       col = nc + 1;
    vector<OCP_DBL> dtmp(row, 0);
    vector<OCP_DBL> tmpNij(np * nc, 0);
    OCP_USI         n_np_j;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;

    OCP_DBL dSmax;
    OCP_DBL dP;
    OCP_DBL chop = 1;

    for (OCP_USI n = 0; n < nb; n++) {
        const vector<OCP_DBL>& scm = bk.satcm[bk.SATNUM[n]];

        dP         = u[n * col];
        bk.dPNR[n] = dP;
        bk.P[n] += dP; // seems better
        if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;

        // rockVp[n] = rockVntg[n] * (1 + rockC1 * (P[n] - rockPref));
        // cout << scientific << setprecision(6) << dP << "   " << n << endl;

        dSmax = 0;
        chop  = 1;

        const USI cNp = bk.phaseNum[n];
        const USI len = bk.bRowSizedSdP[n];
        Dcopy(len, &dtmp[0], &bk.res_n[n * row]);
        DaAxpby(len, col, 1.0, &bk.dSec_dPri[n * bk.maxLendSdP], u.data() + n * col,
                1.0, dtmp.data());

        // Calculate dSmax
        USI js = 0;
        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (bk.pSderExist[n_np_j]) {
                dSmax = max(fabs(dtmp[js]), dSmax);
                js++;
            }
        }

        // Calculate chop
        if (dSmax > dSmaxlim) {
            chop  = min(chop, dSmaxlim / dSmax);
            dSmax = dSmaxlim;
        }
        js = 0;
        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (bk.pSderExist[n_np_j]) {
                if (fabs(bk.S[n_np_j] - scm[j]) > TINY &&
                    (bk.S[n_np_j] - scm[j]) / (chop * dtmp[js]) < 0)
                    chop *= min(1.0, -((bk.S[n_np_j] - scm[j]) / (chop * dtmp[js])));
                if (bk.S[n_np_j] + chop * dtmp[js] < 0)
                    chop *= min(1.0, -bk.S[n_np_j] / (chop * dtmp[js]));
                js++;
            }
        }

        // Calculate S, phaseExist, xij, nj
        fill(tmpNij.begin(), tmpNij.end(), 0.0);
        // fill(Ni.begin() + n * nc, Ni.begin() + n * nc + nc, 0.0);

        js     = 0;
        USI jx = cNp;
        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (bk.pSderExist[n_np_j]) {
                bk.dSNRP[n_np_j] = chop * dtmp[js];
                js++;

                // S[n_np_j] += chop * dtmp[js];
                // if (S[n_np_j] < TINY) {
                //     pSderExist[n_np_j] = OCP_FALSE;
                // }
                // js++;
                Daxpy(nc, bk.nj[n_np_j], &bk.xij[n_np_j * nc], &tmpNij[j * nc]);
                for (USI i = 0; i < bk.pVnumCom[j]; i++) {
                    tmpNij[j * nc + i] += chop * dtmp[jx + i];
                }
                jx += bk.pVnumCom[j];
                bk.nj[n_np_j] = Dnorm1(nc, &tmpNij[j * nc]);
                for (USI i = 0; i < nc; i++) {
                    bk.xij[n_np_j * nc + i] = tmpNij[j * nc + i] / bk.nj[n_np_j];
                    // Ni[n * nc + i] += tmpNij[j * nc + i];
                }
            }
        }

        // Vf correction
        /*        OCP_DBL dVmin = 100 * rockVp[n];
                USI index = 0;
                for (USI j = 0; j < numPhase; j++) {
                    if (phaseExist[n * numPhase + j]) {
                        tmpxi[j] = flashCal[PVTNUM[n]]->XiPhase(P[n], T, &xij[(n *
           numPhase + j) * nc]); tmpVf[j] = nj[n * numPhase + j] / (S[n * numPhase +
           j] * tmpxi[j]); if (fabs(tmpVf[j] - rockVp[n]) < dVmin) { dVmin =
           fabs(tmpVf[j] - rockVp[n]); index = j;
                        }
                    }
                } */
        // for (USI j = 0; j < numPhase; j++) {
        //     if (phaseExist[n * numPhase + j]) {
        //         nj[n * numPhase + j] *= (tmpVf[index] / tmpVf[j]);
        //         for (USI i = 0; i < nc; i++) {
        //             Ni[n * nc + i] += nj[n * numPhase + j] * xij[(n * numPhase +
        //             j) * nc + i];
        //         }
        //     }
        // }

        bk.NRstep[n] = chop;
        for (USI i = 0; i < nc; i++) {
            bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chop;
            if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
                bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];
            bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];
        }
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

void IsoT_FIMn::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.bulk.nj    = rs.bulk.lnj;
    rs.bulk.res_n = rs.bulk.lres_n;
    rs.bulk.resPc = rs.bulk.lresPc;

    IsoT_FIM::ResetToLastTimeStep(rs, ctrl);
}

void IsoT_FIMn::UpdateLastTimeStep(Reservoir& rs) const
{
    rs.bulk.lnj    = rs.bulk.nj;
    rs.bulk.lres_n = rs.bulk.res_n;
    rs.bulk.lresPc = rs.bulk.resPc;
    IsoT_FIM::UpdateLastTimeStep(rs);
}

////////////////////////////////////////////
// IsoT_AIMc
////////////////////////////////////////////

void IsoT_AIMc::Setup(Reservoir& rs, LinearSystem& ls, const OCPControl& ctrl)
{
    // Allocate Bulk and BulkConn Memory
    AllocateReservoir(rs);
    // Allocate memory for internal matrix structure
    IsoT_FIM::AllocateLinearSystem(ls, rs, ctrl);
}

void IsoT_AIMc::InitReservoir(Reservoir& rs) const
{
    rs.bulk.InitPTSw(50);

    InitRock(rs.bulk);
    CalRock(rs.bulk);

    IsoT_IMPEC::InitFlash(rs.bulk);
    IsoT_IMPEC::CalKrPc(rs.bulk);

    rs.allWells.InitBHP(rs.bulk);

    UpdateLastTimeStep(rs);
}

void IsoT_AIMc::Prepare(Reservoir& rs, const OCP_DBL& dt)
{
    rs.allWells.PrepareWell(rs.bulk);
    CalRes(rs, dt, OCP_TRUE);

    // Set FIM Bulk
    rs.CalCFL(dt);
    rs.allWells.SetupWellBulk(rs.bulk);
    SetFIMBulk(rs);
    //  Calculate FIM Bulk properties
    CalFlashI(rs.bulk);
    CalKrPcI(rs.bulk);

    UpdateLastTimeStep(rs);
    // rs.bulk.ShowFIMBulk(OCP_FALSE);
}

void IsoT_AIMc::AssembleMat(LinearSystem&    ls,
                            const Reservoir& rs,
                            const OCP_DBL&   dt) const
{
    AssembleMatBulks(ls, rs, dt);
    IsoT_FIM::AssembleMatWellsNew(ls, rs, dt);
    ls.AssembleRhsCopy(rs.bulk.res.resAbs);
}

void IsoT_AIMc::SolveLinearSystem(LinearSystem& ls, Reservoir& rs, OCPControl& ctrl)
{
#ifdef DEBUG
    ls.CheckEquation();
#endif // DEBUG

    ls.AssembleMatLinearSolver();

    GetWallTime Timer;
    Timer.Start();
    int status = ls.Solve();
    if (status < 0) {
        status = ls.GetNumIters();
    }

#ifdef DEBUG
    // ls.OutputLinearSystem("testA_AIMc.out", "testb_AIMc.out");
    // ls.OutputSolution("testx_AIMc.out");
    // ls.CheckSolution();
#endif // DEBUG

    ctrl.RecordTimeLS(Timer.Stop() / 1000);
    ctrl.UpdateIterLS(status);
    ctrl.UpdateIterNR();

    GetSolution(rs, ls.GetSolution(), ctrl);
    ls.ClearData();
}

OCP_BOOL IsoT_AIMc::UpdateProperty(Reservoir& rs, OCPControl& ctrl)
{
    OCP_DBL& dt = ctrl.current_dt;

    // First check: Ni check and bulk Pressure check
    if (!ctrl.Check(rs, {"BulkNi", "BulkP"})) {
        ResetToLastTimeStep(rs, ctrl);
        cout << "Cut time step size and repeat! current dt = " << fixed
             << setprecision(3) << dt << " days\n";
        return OCP_FALSE;
    }

    CalFlashI(rs.bulk);
    CalFlashEp(rs.bulk);
    CalKrPcI(rs.bulk);

    CalRock(rs.bulk);

    rs.allWells.CalTrans(rs.bulk);
    rs.allWells.CalFlux(rs.bulk);

    CalRes(rs, dt, OCP_FALSE);
    return OCP_TRUE;
}

OCP_BOOL IsoT_AIMc::FinishNR(Reservoir& rs, OCPControl& ctrl)
{
    OCP_USI       dSn;
    const OCP_DBL NRdSmax = rs.GetNRdSmax(dSn);
    const OCP_DBL NRdPmax = rs.GetNRdPmax();
    // const OCP_DBL NRdNmax = rs.GetNRdNmax();

#ifdef DEBUG
    // cout << "### DEBUG: Residuals = " << setprecision(3) << scientific
    //      << resAIMc.maxRelRes0_V << "  " << resAIMc.maxRelRes_V << "  "
    //      << resAIMc.maxRelRes_N << "  " << NRdPmax << "  " << NRdSmax << endl;
#endif

    if (((rs.bulk.res.maxRelRes_V <= rs.bulk.res.maxRelRes0_V * ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_V <= ctrl.ctrlNR.NRtol ||
          rs.bulk.res.maxRelRes_N <= ctrl.ctrlNR.NRtol) &&
         rs.bulk.res.maxWellRelRes_mol <= ctrl.ctrlNR.NRtol) ||
        (fabs(NRdPmax) <= ctrl.ctrlNR.NRdPmin &&
         fabs(NRdSmax) <= ctrl.ctrlNR.NRdSmin)) {

        if (!ctrl.Check(rs, {"WellP"})) {
            ResetToLastTimeStep(rs, ctrl);
            return OCP_FALSE;
        } else {
            CalFlashEa(rs.bulk);
            CalKrPcE(rs.bulk);
            return OCP_TRUE;
        }

    } else if (ctrl.iterNR > ctrl.ctrlNR.maxNRiter) {
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

/// Finish a time step.
void IsoT_AIMc::FinishStep(Reservoir& rs, OCPControl& ctrl) const
{
    rs.CalIPRT(ctrl.GetCurDt());
    rs.CalMaxChange();
    UpdateLastTimeStep(rs);
    ctrl.CalNextTimeStep(rs, {"dP", "dS", "iter"});
}

/// Allocate memory for reservoir
void IsoT_AIMc::AllocateReservoir(Reservoir& rs)
{
    IsoT_FIM::AllocateReservoir(rs);

    Bulk&         bk = rs.bulk;
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    bk.vj.resize(nb * np);
    bk.lvj.resize(nb * np);

    bk.xijNR.resize(nb * np * nc);
    bk.cfl.resize(nb * np);
    bk.bulkTypeAIM.Setup(nb);
}

void IsoT_AIMc::SetFIMBulk(Reservoir& rs)
{
    Bulk&           bk   = rs.bulk;
    const BulkConn& conn = rs.conn;
    const OCP_USI   nb   = bk.numBulk;
    const USI       np   = bk.numPhase;
    const USI       nc   = bk.numCom;

    bk.bulkTypeAIM.Init();

    OCP_USI  bIdp, bIdc;
    OCP_BOOL flag;

    for (OCP_USI n = 0; n < nb; n++) {
        bIdp = n * np;
        bIdc = n * nc;
        flag = OCP_FALSE;
        // CFL
        for (USI j = 0; j < np; j++) {
            if (bk.cfl[bIdp + j] > 0.8) {
                flag = OCP_TRUE;
                break;
            }
        }
        // Volume error
        if (!flag) {
            if ((fabs(bk.vf[n] - bk.rockVp[n]) / bk.rockVp[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
        }

        // NR Step
        if (!flag && OCP_FALSE) {
            // dP
            if (fabs(bk.dPNR[n] / bk.P[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
            // dNi
            if (!flag) {
                for (USI i = 0; i < bk.numCom; i++) {
                    if (fabs(bk.dNNR[bIdc + i] / bk.Ni[bIdc + i]) > 1E-3) {
                        flag = OCP_TRUE;
                        break;
                    }
                }
            }
        }

        if (flag) {
            // find it
            // bk.map_Bulk2FIM[n] = 1;
            for (auto& v : conn.neighbor[n]) {
                // n is included also
                bk.bulkTypeAIM.SetBulkType(v, 1);
            }
        }
    }

    // add WellBulk's 2-neighbor as Implicit bulk
    for (auto& p : bk.wellBulkId) {
        for (auto& v : conn.neighbor[p]) {
            for (auto& v1 : conn.neighbor[v]) bk.bulkTypeAIM.SetBulkType(v1, 1);
        }
    }
}

void IsoT_AIMc::CalFlashEp(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                  bk.phaseNum[n],
                                                  &bk.xijNR[n * np * nc], n);
            // bk.PassFlashValueAIMcEp(n);
            PassFlashValueEp(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashEa(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashIMPEC(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                  bk.phaseNum[n], &bk.xij[n * np * nc],
                                                  n);
            // bk.PassFlashValueAIMcEa(n);

            IsoT_IMPEC::PassFlashValue(bk, n);
        }
    }
}

void IsoT_AIMc::CalFlashI(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;
    const USI     nc = bk.numCom;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk

            bk.flashCal[bk.PVTNUM[n]]->FlashFIM(bk.P[n], bk.T[n], &bk.Ni[n * nc],
                                                &bk.S[n * np], bk.phaseNum[n],
                                                &bk.xij[n * np * nc], n);
            IsoT_FIM::PassFlashValue(bk, n);
            for (USI j = 0; j < np; j++) {
                bk.vj[n * np + j] = bk.vf[n] * bk.S[n * np + j];
            }
        }
    }
}

void IsoT_AIMc::PassFlashValueEp(Bulk& bk, const OCP_USI& n)
{
    // only var about volume needs, some flash var also
    OCP_FUNCNAME;

    const USI     np     = bk.numPhase;
    const USI     nc     = bk.numCom;
    const OCP_USI bIdp   = n * np;
    const USI     pvtnum = bk.PVTNUM[n];

    bk.Nt[n]  = bk.flashCal[pvtnum]->GetNt();
    bk.vf[n]  = bk.flashCal[pvtnum]->GetVf();
    bk.vfP[n] = bk.flashCal[pvtnum]->GetVfP();
    for (USI i = 0; i < nc; i++) {
        bk.vfi[n * nc + i] = bk.flashCal[pvtnum]->GetVfi(i);
    }

    bk.phaseNum[n] = 0;
    for (USI j = 0; j < np; j++) {
        if (bk.flashCal[pvtnum]->GetPhaseExist(j)) {
            bk.phaseNum[n]++;

            // IMPORTANT -- need for next Flash
            // But xij in nonlinear equations has been modified
            for (USI i = 0; i < nc; i++) {
                bk.xijNR[bIdp * nc + j * nc + i] = bk.flashCal[pvtnum]->GetXij(j, i);
            }
        }
    }
}

void IsoT_AIMc::CalKrPcE(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // Explicit bulk
            const OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPc(&bk.S[bId], &bk.kr[bId], &bk.Pc[bId], n);
            for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::CalKrPcI(Bulk& bk)
{
    const OCP_USI nb = bk.numBulk;
    const USI     np = bk.numPhase;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfFIMbulk(n)) {
            // Implicit bulk
            const OCP_USI bId = n * np;
            bk.flow[bk.SATNUM[n]]->CalKrPcDeriv(&bk.S[bId], &bk.kr[bId], &bk.Pc[bId],
                                                &bk.dKr_dS[bId * np],
                                                &bk.dPcj_dS[bId * np], n);
            for (USI j = 0; j < np; j++) bk.Pj[bId + j] = bk.P[n] + bk.Pc[bId + j];
        }
    }
}

void IsoT_AIMc::AssembleMatBulks(LinearSystem&    ls,
                                 const Reservoir& rs,
                                 const OCP_DBL&   dt) const
{
    const Bulk&     bk      = rs.bulk;
    const BulkConn& conn    = rs.conn;
    const OCP_USI   nb      = bk.numBulk;
    const USI       np      = bk.numPhase;
    const USI       nc      = bk.numCom;
    const USI       ncol    = nc + 1;
    const USI       ncol2   = np * nc + np;
    const USI       bsize   = ncol * ncol;
    const USI       bsize2  = ncol * ncol2;
    const USI       lendSdP = bk.maxLendSdP;

    ls.AddDim(nb);

    vector<OCP_DBL> bmat(bsize, 0);
    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < nb; n++) {
        bmat[0] = bk.v[n] * bk.poroP[n] - bk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -bk.vfi[n * nc + i];
        }
        ls.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    OCP_DBL*         dFdXpI;
    OCP_DBL*         dFdXsI;
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistI(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_I(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    vector<USI>      pVnumComI(np, 0);
    OCP_USI          ncolB, ncolE, ncolI;

    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j, uId_np_j, ibId_np_j;
    OCP_DBL  kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL  dP, dGamma;
    OCP_DBL  tmp;
    OCP_BOOL bIdFIM, eIdFIM, uIdFIM;

    for (OCP_USI c = 0; c < conn.numConn; c++) {
        bId    = conn.iteratorConn[c].BId();
        eId    = conn.iteratorConn[c].EId();
        Akd    = CONV1 * CONV2 * conn.iteratorConn[c].Area();
        dGamma = GRAVITY_FACTOR * (bk.depth[bId] - bk.depth[eId]);
        bIdFIM = eIdFIM = OCP_FALSE;
        if (bk.bulkTypeAIM.IfFIMbulk(bId)) bIdFIM = OCP_TRUE;
        if (bk.bulkTypeAIM.IfFIMbulk(eId)) eIdFIM = OCP_TRUE;

        if (!bIdFIM && !eIdFIM) {
            // both are explicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = bk.phaseExist[bId * np + j];
                phaseExistE[j] = bk.phaseExist[eId * np + j];
                uId            = conn.upblock[c * np + j];
                phaseExistU    = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    continue;
                }
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                xi       = bk.xi[uId_np_j];
                kr       = bk.kr[uId_np_j];
                mu       = bk.mu[uId_np_j];
                transJ   = Akd * xi * kr / mu;
                for (USI i = 0; i < nc; i++) {
                    xij     = bk.xij[uId_np_j * nc + i];
                    transIJ = xij * transJ;
                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;
                }
            }
        } else if ((bIdFIM && !eIdFIM) || (!bIdFIM && eIdFIM)) {

            // one is explicit, one is implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            ncolI   = 0;
            USI jxI = 0;
            if (bIdFIM) {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = bk.phaseExist[bId * np + j];
                    phasedS_I[j]   = bk.pSderExist[bId * np + j];
                    if (phasedS_I[j]) jxI++;
                    pVnumComI[j] = bk.pVnumCom[bId * np + j];
                    ncolI += pVnumComI[j];
                }
                dFdXpI = &dFdXpB[0];
                dFdXsI = &dFdXsB[0];
            } else {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = bk.phaseExist[eId * np + j];
                    phasedS_I[j]   = bk.pSderExist[eId * np + j];
                    if (phasedS_I[j]) jxI++;
                    pVnumComI[j] = bk.pVnumCom[eId * np + j];
                    ncolI += pVnumComI[j];
                }
                dFdXpI = &dFdXpE[0];
                dFdXsI = &dFdXsE[0];
            }
            ncolI += jxI;

            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = bk.phaseExist[bId * np + j];
                phaseExistE[j] = bk.phaseExist[eId * np + j];
                uId            = conn.upblock[c * np + j];
                phaseExistU    = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    continue;
                }

                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                dP       = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                     conn.upblock_Rho[c * np + j] * dGamma;
                xi     = bk.xi[uId_np_j];
                kr     = bk.kr[uId_np_j];
                mu     = bk.mu[uId_np_j];
                transJ = Akd * kr / mu;
                uIdFIM = (uId == bId ? bIdFIM : eIdFIM);

                if (bIdFIM)
                    ibId_np_j = bId_np_j;
                else
                    ibId_np_j = eId_np_j;

                if (uIdFIM) {
                    // upblock is implicit
                    OCP_DBL rhoWgt = 1.0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;

                    muP  = bk.muP[ibId_np_j];
                    xiP  = bk.xiP[ibId_np_j];
                    rhoP = bk.rhoP[ibId_np_j];
                    for (USI i = 0; i < nc; i++) {
                        xij     = bk.xij[ibId_np_j * nc + i];
                        transIJ = xij * xi * transJ;

                        // Pressure -- Primary var
                        dFdXpB[(i + 1) * ncol] += transIJ;
                        dFdXpE[(i + 1) * ncol] -= transIJ;

                        tmp = transIJ * (-rhoP * dGamma) * rhoWgt;
                        tmp += xij * transJ * xiP * dP;
                        tmp += -transIJ * muP / mu * dP;
                        dFdXpI[(i + 1) * ncol] += tmp;

                        USI j1S = 0;
                        // Saturation -- Second var
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_I[j1]) {
                                dFdXsI[(i + 1) * ncolI + j1S] +=
                                    transIJ * bk.dPcj_dS[ibId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      bk.dKr_dS[ibId_np_j * np + j1] * dP;
                                dFdXsI[(i + 1) * ncolI + j1S] += tmp;
                                j1S++;
                            }
                        }
                        // Cij
                        for (USI k = 0; k < pVnumComI[j]; k++) {
                            rhox = bk.rhox[ibId_np_j * nc + k] * rhoWgt;
                            xix  = bk.xix[ibId_np_j * nc + k];
                            mux  = bk.mux[ibId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComI[j])
                            dFdXsI[(i + 1) * ncolI + jxI + i] += xi * transJ * dP;
                    }
                    jxI += pVnumComI[j];
                } else {
                    // downblock is implicit
                    OCP_DBL rhoWgt = 0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;
                    rhoP = bk.rhoP[ibId_np_j];

                    for (USI i = 0; i < nc; i++) {
                        xij     = bk.xij[uId_np_j * nc + i];
                        transIJ = xij * xi * transJ;

                        // Pressure -- Primary var
                        dFdXpB[(i + 1) * ncol] += transIJ;
                        dFdXpE[(i + 1) * ncol] -= transIJ;

                        tmp = transIJ * (-rhoP * dGamma) * rhoWgt;
                        dFdXpI[(i + 1) * ncol] += tmp;

                        USI j1S = 0;
                        // Saturation -- Second var
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_I[j1]) {
                                dFdXsI[(i + 1) * ncolI + j1S] +=
                                    transIJ * bk.dPcj_dS[ibId_np_j * np + j1];
                                j1S++;
                            }
                        }
                        // Cij
                        for (USI k = 0; k < pVnumComI[j]; k++) {
                            rhox = bk.rhox[ibId_np_j * nc + k] * rhoWgt;
                            tmp  = -transIJ * rhox * dGamma;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                    }
                    jxI += pVnumComI[j];
                }
            }
        } else {
            // both are implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            USI jxB = 0;
            USI jxE = 0;
            ncolB   = 0;
            ncolE   = 0;

            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = bk.phaseExist[bId * np + j];
                phaseExistE[j] = bk.phaseExist[eId * np + j];
                phasedS_B[j]   = bk.pSderExist[bId * np + j];
                phasedS_E[j]   = bk.pSderExist[eId * np + j];
                if (phasedS_B[j]) jxB++;
                if (phasedS_E[j]) jxE++;
                pVnumComB[j] = bk.pVnumCom[bId * np + j];
                pVnumComE[j] = bk.pVnumCom[eId * np + j];
                ncolB += pVnumComB[j];
                ncolE += pVnumComE[j];
            }
            ncolB += jxB;
            ncolE += jxE;

            for (USI j = 0; j < np; j++) {
                uId = conn.upblock[c * np + j];

                phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    jxB += pVnumComB[j];
                    jxE += pVnumComE[j];
                    continue;
                }

                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                dP       = bk.Pj[bId_np_j] - bk.Pj[eId_np_j] -
                     conn.upblock_Rho[c * np + j] * dGamma;
                xi     = bk.xi[uId_np_j];
                kr     = bk.kr[uId_np_j];
                mu     = bk.mu[uId_np_j];
                muP    = bk.muP[uId_np_j];
                xiP    = bk.xiP[uId_np_j];
                rhoP   = bk.rhoP[uId_np_j];
                transJ = Akd * kr / mu;

                for (USI i = 0; i < nc; i++) {
                    xij     = bk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    tmp = xij * transJ * xiP * dP;
                    tmp += -transIJ * muP / mu * dP;
                    if (!phaseExistE[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else if (!phaseExistB[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpE[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpB[(i + 1) * ncol] +=
                            transIJ * (-bk.rhoP[bId_np_j] * dGamma) / 2;
                        dFdXpE[(i + 1) * ncol] +=
                            transIJ * (-bk.rhoP[eId_np_j] * dGamma) / 2;
                        if (bId == uId) {
                            dFdXpB[(i + 1) * ncol] += tmp;
                        } else {
                            dFdXpE[(i + 1) * ncol] += tmp;
                        }
                    }

                    // Second var
                    USI j1SB = 0;
                    USI j1SE = 0;
                    if (bId == uId) {
                        // Saturation
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_B[j1]) {
                                dFdXsB[(i + 1) * ncolB + j1SB] +=
                                    transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      bk.dKr_dS[uId_np_j * np + j1] * dP;
                                dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                                j1SB++;
                            }
                            if (phasedS_E[j1]) {
                                dFdXsE[(i + 1) * ncolE + j1SE] -=
                                    transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                                j1SE++;
                            }
                        }
                        // Cij
                        if (!phaseExistE[j]) {
                            for (USI k = 0; k < pVnumComB[j]; k++) {
                                rhox = bk.rhox[uId_np_j * nc + k];
                                xix  = bk.xix[uId_np_j * nc + k];
                                mux  = bk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        } else {
                            for (USI k = 0; k < pVnumComB[j]; k++) {
                                rhox = bk.rhox[bId_np_j * nc + k] / 2;
                                xix  = bk.xix[uId_np_j * nc + k];
                                mux  = bk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                                dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                    -transIJ * bk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        }
                    } else {
                        // Saturation
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_B[j1]) {
                                dFdXsB[(i + 1) * ncolB + j1SB] +=
                                    transIJ * bk.dPcj_dS[bId_np_j * np + j1];
                                j1SB++;
                            }
                            if (phasedS_E[j1]) {
                                dFdXsE[(i + 1) * ncolE + j1SE] -=
                                    transIJ * bk.dPcj_dS[eId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      bk.dKr_dS[uId_np_j * np + j1] * dP;
                                dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                                j1SE++;
                            }
                        }
                        // Cij
                        if (!phaseExistB[j]) {
                            for (USI k = 0; k < pVnumComE[j]; k++) {
                                rhox = bk.rhox[uId_np_j * nc + k];
                                xix  = bk.xix[uId_np_j * nc + k];
                                mux  = bk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComE[j])
                                dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                        } else {
                            for (USI k = 0; k < pVnumComE[j]; k++) {
                                rhox = bk.rhox[eId_np_j * nc + k] / 2;
                                xix  = bk.xix[uId_np_j * nc + k];
                                mux  = bk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                                dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                    -transIJ * bk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                            }
                            // WARNING !!!
                            if (i < pVnumComE[j])
                                dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                        }
                    }
                }
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
            }
        }

        if (bIdFIM && !eIdFIM)
            ncolB = ncolI;
        else if (!bIdFIM && eIdFIM)
            ncolE = ncolI;

        // Assemble
        bmat = dFdXpB;
        if (bIdFIM) {
            DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &bk.dSec_dPri[bId * lendSdP],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        ls.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        ls.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        if (eIdFIM) {
            DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &bk.dSec_dPri[eId * lendSdP],
                    1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        ls.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        ls.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif
    }
}

void IsoT_AIMc::GetSolution(Reservoir&             rs,
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
    const USI       col = nc + 1;
    vector<OCP_DBL> dtmp(row, 0);
    OCP_DBL         chopmin = 1;
    OCP_DBL         choptmp = 0;
    OCP_USI         n_np_j;

    bk.dSNR       = bk.S;
    bk.NRphaseNum = bk.phaseNum;
    bk.NRdPmax    = 0;
    bk.NRdNmax    = 0;

    for (OCP_USI n = 0; n < nb; n++) {
        if (bk.bulkTypeAIM.IfIMPECbulk(n)) {
            // IMPEC Bulk
            // Pressure
            OCP_DBL dP = u[n * col];
            bk.NRdPmax = max(bk.NRdPmax, fabs(dP));
            bk.P[n] += dP; // seems better
            bk.dPNR[n]   = dP;
            bk.NRstep[n] = 1;
            // Ni
            for (USI i = 0; i < nc; i++) {
                bk.dNNR[n * nc + i] = u[n * col + 1 + i];
                bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

                // if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
                //     bk.Ni[n * nc + i] = 1E-20;
                // }
            }
            // Pj
            for (USI j = 0; j < np; j++) {
                bk.Pj[n * np + j] = bk.P[n] + bk.Pc[n * np + j];
            }
            continue;
        }

        chopmin = 1;
        // compute the chop
        fill(dtmp.begin(), dtmp.end(), 0.0);
        DaAxpby(bk.bRowSizedSdP[n], col, 1, &bk.dSec_dPri[n * bk.maxLendSdP],
                u.data() + n * col, 1, dtmp.data());

        USI js = 0;
        for (USI j = 0; j < np; j++) {
            if (!bk.pSderExist[n * np + j]) {
                continue;
            }
            n_np_j = n * np + j;

            choptmp = 1;
            if (fabs(dtmp[js]) > dSmaxlim) {
                choptmp = dSmaxlim / fabs(dtmp[js]);
            } else if (bk.S[n_np_j] + dtmp[js] < 0.0) {
                choptmp = 0.9 * bk.S[n_np_j] / fabs(dtmp[js]);
            }

            // if (fabs(S[n_np_j] - scm[j]) > TINY &&
            //     (S[n_np_j] - scm[j]) / (choptmp * dtmp[js]) < 0)
            //     choptmp *= min(1.0, -((S[n_np_j] - scm[j]) / (choptmp * dtmp[js])));

            chopmin = min(chopmin, choptmp);
            js++;
        }

        // dS
        js = 0;
        for (USI j = 0; j < np; j++) {
            if (!bk.pSderExist[n * np + j]) {
                bk.dSNRP[n * np + j] = 0;
                continue;
            }
            bk.dSNRP[n * np + j] = chopmin * dtmp[js];
            js++;
        }

        // dxij   ---- Compositional model only
        if (bk.IfUseEoS()) {
            if (bk.phaseNum[n] >= 3) {
                OCP_USI bId = 0;
                for (USI j = 0; j < 2; j++) {
                    bId = n * np * nc + j * nc;
                    for (USI i = 0; i < bk.numComH; i++) {
                        bk.xij[bId + i] += chopmin * dtmp[js];
                        js++;
                    }
                }
            }
        }

        // dP
        OCP_DBL dP = u[n * col];
        if (fabs(bk.NRdPmax) < fabs(dP)) bk.NRdPmax = dP;
        bk.P[n] += dP; // seems better
        bk.dPNR[n] = dP;

        // dNi
        bk.NRstep[n] = chopmin;
        for (USI i = 0; i < nc; i++) {
            bk.dNNR[n * nc + i] = u[n * col + 1 + i] * chopmin;
            if (fabs(bk.NRdNmax) < fabs(bk.dNNR[n * nc + i]) / bk.Nt[n])
                bk.NRdNmax = bk.dNNR[n * nc + i] / bk.Nt[n];

            bk.Ni[n * nc + i] += bk.dNNR[n * nc + i];

            // if (bk.Ni[n * nc + i] < 0 && bk.Ni[n * nc + i] > -1E-3) {
            //     bk.Ni[n * nc + i] = 1E-20;
            // }
        }
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

void IsoT_AIMc::ResetToLastTimeStep(Reservoir& rs, OCPControl& ctrl)
{
    rs.bulk.vj    = rs.bulk.lvj;
    rs.bulk.xijNR = rs.bulk.lxij;
    IsoT_FIM::ResetToLastTimeStep(rs, ctrl);
}

void IsoT_AIMc::UpdateLastTimeStep(Reservoir& rs) const
{
    IsoT_FIM::UpdateLastTimeStep(rs);

    rs.bulk.lvj   = rs.bulk.vj;
    rs.bulk.xijNR = rs.bulk.xij;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/01/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update output                        */
/*----------------------------------------------------------------------------*/