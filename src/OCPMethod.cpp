/*! \file    OCPMethod.cpp
 *  \brief   OCPMethod class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPMethod.hpp"


void OCPMethod::CalLsColCapacity(LinearSolver& ls, const Reservoir& rs)
{
	// Bulk to Bulk
	const BulkConn& conn = rs.conn;
	const OCP_USI nb = conn.numBulk;
	for (OCP_USI n = 0; n < nb; n++) {
		ls.rowCapacity[n] = conn.neighborNum[n];
	}
	// Well to Bulk
	OCP_USI k;
	const USI nw = rs.wellgroup.numWell;
	const USI maxNum = rs.wellgroup.GetMaxWellPerNum() + 1;
	const vector<Well>& well = rs.wellgroup.wellGroup;
	for (USI w = 0; w < nw; w++) {
		for (USI p = 0; p < well[w].numPerf; p++) {
			k = well[w].perf[p].location;
			ls.rowCapacity[k]++;
		}
		ls.rowCapacity[nb + w] += maxNum;
	}
}


void OCPMethod::InitMatSparsity(LinearSolver& ls, const BulkConn& conn)
{
	ls.dim = conn.numBulk;
	for (OCP_USI n = 0; n < conn.numBulk; n++) {
		ls.colId[n].assign(conn.neighbor[n].begin(), conn.neighbor[n].end());
		ls.diagPtr[n] = conn.selfPtr[n];
	}
}


void OCPMethod::SetupIMPEC(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl)
{
	AllocateRsIMPEC(rs);
	SetupLsIMPEC(ls, rs, ctrl);
}


void OCPMethod::AllocateRsIMPEC(Reservoir& rs)
{
	Bulk& bulk = rs.bulk;
	const OCP_USI nb = bulk.numBulk;
	const OCP_USI np = bulk.numPhase;
	const USI nc = bulk.numCom;
	// Derivatives
	bulk.vfi.resize(nb * nc);
	bulk.vfp.resize(nb);
	bulk.lvfi.resize(nb * nc);
	bulk.lvfp.resize(nb * nc);

	// Connections
	BulkConn& conn = rs.conn;
	const OCP_USI NC = conn.numConn;
	conn.upblock.resize(NC * np);
	conn.upblock_Rho.resize(NC * np);
	conn.upblock_Trans.resize(NC * np);
	conn.upblock_Velocity.resize(NC * np);
	conn.lastUpblock.resize(NC * np);
	conn.lastUpblock_Rho.resize(NC * np);
	conn.lastUpblock_Trans.resize(NC * np);
	conn.lastUpblock_Velocity.resize(NC * np);
}


void OCPMethod::SetupLsIMPEC(LinearSolver& ls, const Reservoir& rs, const OCPControl& ctrl)
{
	ls.SetupParam(ctrl.workDir, ctrl.lsFile);
	AllocateLsIMPEC(ls, rs);
}


void OCPMethod::AllocateLsIMPEC(LinearSolver& ls, const Reservoir& rs)
{
	const OCP_USI nb = rs.bulk.numBulk;
	const USI nw = rs.wellgroup.numWell;
	ls.AllocateRowMem(nb + nw, 1);
	CalLsColCapacity(ls, rs);
	ls.AllocateColMem();
	ls.AllocateFasp();
}


void OCPMethod::AssembleLsIMPEC(LinearSolver& ls, const Reservoir& rs, const OCP_DBL& dt)
{
	InitMatSparsity(ls, rs.conn);
	AssembleLsConnIMPEC(ls, rs.conn, rs.bulk, dt);
	AssembleLsWellIMPEC(ls, rs.wellgroup, rs.bulk, dt);
}


void OCPMethod::AssembleLsConnIMPEC(LinearSolver& ls, const BulkConn& conn, const Bulk& bulk, const OCP_DBL& dt)
{
	// accumulate term
	OCP_DBL Vp0, Vp, vf, vfp, P;
	const OCP_DBL cr = bulk.rockC1;
	for (OCP_USI n = 0; n < bulk.numBulk; n++) {
		Vp0 = bulk.rockVpInit[n];
		Vp = bulk.rockVp[n];
		vfp = bulk.vfp[n];
		P = bulk.P[n];
		vf = bulk.vf[n];

		OCP_DBL temp = cr * Vp0 - vfp;
		ls.diagVal[n] = temp;
		ls.b[n] = temp * P + dt * (vf - Vp);
	}

	// flux term
	OCP_USI bId, eId, uId;
	const USI     np = bulk.numPhase;
	const USI     nc = bulk.numCom;
	OCP_DBL valupi, valdowni;
	OCP_DBL valup, rhsup, valdown, rhsdown;
	// OCP_USI    lastbId = -1;
	OCP_USI lastbId = conn.iteratorConn[0].GetEId();
	for (OCP_USI c = 0; c < conn.numConn; c++) {
		bId = conn.iteratorConn[c].GetBId();
		eId = conn.iteratorConn[c].GetEId();
		valup = 0;
		rhsup = 0;
		valdown = 0;
		rhsdown = 0;

		for (USI j = 0; j < np; j++) {
			uId = conn.upblock[c * np + j];
			if (!bulk.phaseExist[uId * np + j]) continue;

			valupi = 0;
			valdowni = 0;

			for (USI i = 0; i < nc; i++) {
				valupi +=
					bulk.vfi[bId * nc + i] * bulk.cij[uId * np * nc + j * nc + i];
				valdowni +=
					bulk.vfi[eId * nc + i] * bulk.cij[uId * np * nc + j * nc + i];
			}
			OCP_DBL dD = bulk.depth[bId] - bulk.depth[eId];
			OCP_DBL dPc = bulk.Pc[bId * np + j] - bulk.Pc[eId * np + j];
			OCP_DBL temp = bulk.xi[uId * np + j] * conn.upblock_Trans[c * np + j] * dt;
			valup += temp * valupi;

			//if (!isfinite(valup)) {
			//    cout << "###ERROR   ";
			//    ERRORcheck("NAN or INF in MAT");
			//}

			valdown += temp * valdowni;
			temp *= conn.upblock_Rho[c * np + j] * GRAVITY_FACTOR * dD - dPc;
			rhsup += temp * valupi;
			rhsdown -= temp * valdowni;
		}

		USI diagptr = conn.selfPtr[bId];
		if (bId != lastbId) {
			// new bulk
			assert(ls.val[bId].size() == diagptr);
			ls.val[bId].push_back(ls.diagVal[bId]);
			lastbId = bId;
		}

		ls.val[bId][diagptr] += valup;
		ls.val[bId].push_back(-valup);
		ls.val[eId].push_back(-valdown);
		ls.diagVal[eId] += valdown;
		ls.b[bId] += rhsup;
		ls.b[eId] += rhsdown;
	}

	// Add the rest of diag value
	// important!
	for (OCP_USI n = 0; n < conn.numBulk; n++) {
		if (ls.val[n].size() == conn.selfPtr[n])
			ls.val[n].push_back(ls.diagVal[n]);
	}
}


void OCPMethod::AssembleLsWellIMPEC(LinearSolver& ls, const WellGroup& wellG, const Bulk& bulk, const OCP_DBL& dt)
{
	const USI nw = wellG.numWell;
	const vector<Well>& well = wellG.wellGroup;
	for (USI w = 0; w < nw; w++) {
		if (well[w].WellState()) {
			switch (well[w].WellType())
			{
			case INJ:
				AssembleLsINJWellIMPEC(ls, well[w], bulk, dt);
				break;
			case PROD:
				AssembleLsPRODWellIMPEC(ls, well[w], bulk, dt);
				break;
			default:
				OCP_ABORT("Wrong Well Type!");
				break;
			}
		}
	}
}


void OCPMethod::AssembleLsINJWellIMPEC(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt)
{
	const USI     nc = bulk.numCom;
	const OCP_USI wId = ls.dim;
	// important !
	ls.dim++;

	for (USI p = 0; p < well.numPerf; p++) {
		OCP_USI k = well.perf[p].location;

		OCP_DBL Vfi_zi = 0;
		for (USI i = 0; i < nc; i++) {
			Vfi_zi += bulk.vfi[k * nc + i] * well.opt.zi[i];
		}

		// OCP_DBL valw = dt * perf[p].xi * perf[p].transj[0];
		OCP_DBL valw = dt * well.perf[p].xi * well.perf[p].transINJ;
		OCP_DBL bw = valw * well.dG[p];
		OCP_DBL valb = valw * Vfi_zi;
		OCP_DBL bb = valb * well.dG[p];

		// Bulk to Well

		// diag
		USI ptr = ls.diagPtr[k];
		ls.val[k][ptr] += valb;
		// off diag
		ls.colId[k].push_back(wId);
		ls.val[k].push_back(-valb);
		// b
		ls.b[k] += bb;

		// Well to Bulk
		switch (well.opt.optMode) {
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// diag
			ls.diagVal[wId] += valw;
			// off diag
			ls.colId[wId].push_back(k);
			ls.val[wId].push_back(-valw);
			// b
			ls.b[wId] -= bw;
			break;
		case BHP_MODE:
			ls.colId[wId].push_back(k);
			ls.val[wId].push_back(0);
			break;
		default:
			OCP_ABORT("Wrong Well Opt mode");
			exit(0);
		}
	}

	// Well Self
	OCP_ASSERT(ls.val[wId].size() == numPerf, "Wrong Size!");
	// the order of perforation is not necessarily in order
	switch (well.opt.optMode) {
	case RATE_MODE:
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		ls.colId[wId].push_back(wId);
		ls.diagPtr[wId] = well.numPerf;
		ls.val[wId].push_back(ls.diagVal[wId]);
		// b
		ls.b[wId] += dt * well.opt.maxRate;
		break;
	case BHP_MODE:
		// diag
		ls.colId[wId].push_back(wId);
		ls.diagPtr[wId] = well.numPerf;
		ls.val[wId].push_back(dt);
		// b
		ls.b[wId] += dt * well.opt.maxBHP;
		// u   initial value
		ls.u[wId] = well.opt.maxBHP;
		break;
	default:
		OCP_ABORT("Wrong Well Opt mode in function");
		exit(0);
	}
}
void OCPMethod::AssembleLsPRODWellIMPEC(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt)
{
	USI     np = bulk.numPhase;
	USI     nc = bulk.numCom;
	OCP_USI wId = ls.dim;
	// important !
	ls.dim++;

	for (USI p = 0; p < well.numPerf; p++) {
		OCP_USI k = well.perf[p].location;

		OCP_DBL valb = 0;
		OCP_DBL bb = 0;
		OCP_DBL valw = 0;
		OCP_DBL bw = 0;

		for (USI j = 0; j < np; j++) {
			if (!bulk.phaseExist[k * np + j]) continue;

			OCP_DBL tempb = 0;
			OCP_DBL tempw = 0;

			for (USI i = 0; i < nc; i++) {
				tempb += bulk.vfi[k * nc + i] * bulk.cij[k * np * nc + j * nc + i];
				tempw += well.opt.zi[i] * bulk.cij[k * np * nc + j * nc + i];
			}
			OCP_DBL trans = dt * well.perf[p].transj[j] * bulk.xi[k * np + j];
			valb += tempb * trans;
			valw += tempw * trans;

			OCP_DBL dP = well.dG[p] - bulk.Pc[k * np + j];
			bb += tempb * trans * dP;
			bw += tempw * trans * dP;
		}

		// Bulk to Well
		// diag
		USI ptr = ls.diagPtr[k];
		ls.val[k][ptr] += valb;
		// off diag
		ls.colId[k].push_back(wId);
		ls.val[k].push_back(-valb);
		// b
		ls.b[k] += bb;

		// Well to Bulk
		switch (well.opt.optMode) {
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// diag  !!! attention! sign is -
			ls.diagVal[wId] -= valw;
			// off diag
			ls.colId[wId].push_back(k);
			ls.val[wId].push_back(valw);
			// b
			ls.b[wId] += bw;
			break;
		case BHP_MODE:
			// off diag
			ls.colId[wId].push_back(k);
			ls.val[wId].push_back(0);
			break;
		default:
			OCP_ABORT("Wrong Well Opt mode");
			exit(0);
		}
	}

	// Well Self
	OCP_ASSERT(ls.val[wId].size() == numPerf, "Wrong Size!");
	// the order of perforation is not necessarily in order
	switch (well.opt.optMode) {
	case RATE_MODE:
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		ls.colId[wId].push_back(wId);
		ls.diagPtr[wId] = well.numPerf;
		ls.val[wId].push_back(ls.diagVal[wId]);
		// b
		ls.b[wId] += dt * well.opt.maxRate;
		break;
	case BHP_MODE:
		// diag
		ls.colId[wId].push_back(wId);
		ls.diagPtr[wId] = well.numPerf;
		ls.val[wId].push_back(dt);
		// b
		ls.b[wId] += dt * well.opt.minBHP;
		// u   initial value
		ls.u[wId] = well.opt.minBHP;
		break;
	default:
		OCP_ABORT("Wrong Well Opt mode");
		exit(0);
	}
}


void OCPMethod::SetupFIM(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl)
{
	AllocateRsFIM(rs);
	SetupLsFIM(ls, rs, ctrl);
}

void OCPMethod::AllocateRsFIM(Reservoir& rs)
{
	Bulk& bulk = rs.bulk;
	const OCP_USI nb = bulk.numBulk;
	const USI np = bulk.numPhase;
	const USI nc = bulk.numCom;

	// Derivatives
	bulk.vfi.resize(nb * nc);
	bulk.muP.resize(nb * np * np);
	bulk.xiP.resize(nb * np * np);
	bulk.rhoP.resize(nb * np * np);
	bulk.muC.resize(nb * np * np);
	bulk.xiC.resize(nb * np * np);
	bulk.rhoC.resize(nb * np * np);
	bulk.dSec_dPri.resize(nb * (nc + 1) * np * (nc + 1));
	bulk.dKr_dS.resize(nb * np * np);
	bulk.dPcj_dS.resize(nb * np * np);

	// Connections
	BulkConn& conn = rs.conn;
	const OCP_USI NC = conn.numConn;
	conn.upblock.resize(NC * np);
	conn.upblock_Rho.resize(NC * np);
}


void OCPMethod::SetupLsFIM(LinearSolver& ls, const Reservoir& rs, const OCPControl& ctrl)
{
	ls.SetupParamB(ctrl.workDir, ctrl.lsFile);
	AllocateLsFIM(ls, rs);
}


void OCPMethod::AllocateLsFIM(LinearSolver& ls, const Reservoir& rs)
{
	const OCP_USI nb = rs.bulk.numBulk;
	const USI nw = rs.wellgroup.numWell;
	ls.AllocateRowMem(nb + nw, nb + 1);
	CalLsColCapacity(ls, rs);
	ls.AllocateColMem();
	ls.AllocateBFasp();
}


void OCPMethod::AssembleMatFIM(LinearSolver& ls, const Reservoir& rs, const OCP_DBL& dt)
{
	InitMatSparsity(ls, rs.conn);
	AssembleMatConnFIM(ls, rs.conn, rs.bulk, dt);
	AssembleMatWellFIM(ls, rs.wellgroup, rs.bulk, dt);
}


void OCPMethod::AssembleMatConnFIM(LinearSolver& ls, const BulkConn& conn, const Bulk& bulk, const OCP_DBL& dt)
{
	const USI np = bulk.numPhase;
	const USI nc = bulk.numCom;
	const USI ncol = nc + 1;
	const USI ncol2 = np * nc + np;
	const USI bsize = ncol * ncol;
	const USI bsize2 = ncol * ncol2;

	vector<OCP_DBL> bmat(bsize, 0);

	// Accumulation term
	for (USI i = 1; i < nc + 1; i++) {
		bmat[i * ncol + i] = 1;
	}
	for (OCP_USI n = 0; n < bulk.numBulk; n++) {
		bmat[0] = bulk.rockC1 * bulk.rockVpInit[n] - bulk.vfp[n];
		for (USI i = 0; i < nc; i++) {
			bmat[i + 1] = -bulk.vfi[n * nc + i];
		}
		for (USI i = 0; i < bsize; i++) {
			ls.diagVal[n * bsize + i] = bmat[i];
		}
	}
	// flux term
	OCP_DBL Akd;
	OCP_DBL transJ, transIJ;
	vector<OCP_DBL> dFdXpB(bsize, 0);
	vector<OCP_DBL> dFdXpE(bsize, 0);
	vector<OCP_DBL> dFdXsB(bsize2, 0);
	vector<OCP_DBL> dFdXsE(bsize2, 0);

	OCP_USI bId, eId, uId;
	OCP_USI uId_np_j;
	OCP_DBL kr, mu, xi, cij, rhoP, xiP, muP, rhoC, xiC, muC;
	OCP_DBL dP, dGamma;
	OCP_DBL tmp;

	// Becareful when first bulk has no neighbors!
	OCP_USI lastbId = conn.iteratorConn[0].GetEId();
	for (OCP_USI c = 0; c < conn.numConn; c++) {
		bId = conn.iteratorConn[c].GetBId();
		eId = conn.iteratorConn[c].GetEId();
		dFdXpB.assign(bsize, 0);
		dFdXpE.assign(bsize, 0);
		dFdXsB.assign(bsize2, 0);
		dFdXsE.assign(bsize2, 0);
		Akd = CONV1 * CONV2 * conn.area[c];
		dGamma = GRAVITY_FACTOR * (bulk.depth[bId] - bulk.depth[eId]);

		for (USI j = 0; j < np; j++) {
			uId = conn.upblock[c * np + j];
			uId_np_j = uId * np + j;
			if (!bulk.phaseExist[uId_np_j])    continue;
			dP = bulk.Pj[bId * np + j] - bulk.Pj[eId * np + j]
				- bulk.rho[bId * np + j] * dGamma;
			xi = bulk.xi[uId_np_j];
			kr = bulk.kr[uId_np_j];
			mu = bulk.mu[uId_np_j];
			muP = bulk.muP[uId_np_j];
			xiP = bulk.xiP[uId_np_j];
			rhoP = bulk.rhoP[uId_np_j];
			transJ = Akd * kr / mu;

			for (USI i = 0; i < nc; i++) {
				cij = bulk.cij[uId_np_j * nc + i];

				transIJ = cij * xi * transJ;

				// Pressure -- Primary var
				dFdXpB[(i + 1) * ncol] += transIJ;
				dFdXpE[(i + 1) * ncol] -= transIJ;
				tmp = transIJ * (-rhoP * dGamma);
				tmp += cij * transJ * xiP * dP;
				tmp += -transIJ * muP / mu * dP;
				if (bId == uId) {
					dFdXpB[(i + 1) * ncol] += tmp;
				}
				else {
					dFdXpE[(i + 1) * ncol] += tmp;
				}

				// Saturation -- Second var
				for (USI k = 0; k < np; k++) {
					dFdXsB[(i + 1) * ncol2 + k] += transIJ * bulk.dPcj_dS[bId * np * np + j * np + k];
					dFdXsE[(i + 1) * ncol2 + k] += transIJ * bulk.dPcj_dS[eId * np * np + j * np + k];
					tmp = Akd * cij * xi / mu * bulk.dKr_dS[uId * np * np + j * np + k] * dP;
					if (bId == uId) {
						dFdXsB[(i + 1) * ncol2 + k] += tmp;
					}
					else {
						dFdXsE[(i + 1) * ncol2 + k] += tmp;
					}
				}
				// Cij -- Second var
				for (USI k = 0; k < nc; k++) {
					rhoC = bulk.rhoC[uId_np_j * nc + k];
					xiC = bulk.xiC[uId_np_j * nc + k];
					muC = bulk.muC[uId_np_j * nc + k];
					tmp = -transIJ * rhoC * dGamma;
					tmp += cij * transJ * xiC * dP;
					tmp += -transIJ * muC / mu * dP;
					if (k == i) {
						tmp += xi * transJ * dP;
					}
					if (bId == uId) {
						dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
					}
					else {
						dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
					}
				}
			}
		}

		USI diagptr = conn.selfPtr[bId];

		if (bId != lastbId) {
			// new bulk
			OCP_ASSERT(ls.val[bId].size() == diagptr * bsize, "Wrong Size!");
			OCP_USI id = bId * bsize;
			ls.val[bId].insert(ls.val[bId].end(),
				ls.diagVal.data() + id, ls.diagVal.data() + id + bsize);

			lastbId = bId;
		}

		// Assemble
		bmat = dFdXpB;
		DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &bulk.dSec_dPri[bId * bsize2], 1, bmat.data());
		Dscalar(bsize, dt, bmat.data());
		// Begin
		// Add
		for (USI i = 0; i < bsize; i++) {
			ls.val[bId][diagptr * bsize + i] += bmat[i];
		}
		// End
		// Insert
		Dscalar(bsize, -1, bmat.data());
		ls.val[eId].insert(ls.val[eId].end(), bmat.begin(), bmat.end());

		// End
		bmat = dFdXpE;
		DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &bulk.dSec_dPri[eId * bsize2], 1, bmat.data());
		Dscalar(bsize, dt, bmat.data());
		// Begin
		// Insert
		ls.val[bId].insert(ls.val[bId].end(), bmat.begin(), bmat.end());
		// Add
		Dscalar(bsize, -1, bmat.data());
		for (USI i = 0; i < bsize; i++) {
			ls.diagVal[eId * bsize + i] += bmat[i];
		}
	}
	// Add the rest of diag value
	// important!
	for (OCP_USI n = 0; n < conn.numBulk; n++) {
		if (ls.val[n].size() == conn.selfPtr[n] * bsize)
			ls.val[n].insert(ls.val[n].end(),
				ls.diagVal.data() + n * bsize, ls.diagVal.data() + n * bsize + bsize);
	}
}


void OCPMethod::AssembleMatWellFIM(LinearSolver& ls, const WellGroup& wellG, const Bulk& bulk, const OCP_DBL& dt)
{
	const USI nw = wellG.numWell;
	const vector<Well>& well = wellG.wellGroup;
	for (USI w = 0; w < nw; w++) {
		if (well[w].WellState()) {
			switch (well[w].WellType())
			{
			case INJ:
				AssembleMatINJWellFIM(ls, well[w], bulk, dt);
				break;
			case PROD:
				AssembleMatPRODWellFIM(ls, well[w], bulk, dt);
				break;
			default:
				OCP_ABORT("Wrong Well Type!");
				break;
			}
		}
	}
}


void OCPMethod::AssembleMatINJWellFIM(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt)
{
	const OCP_USI wId = ls.dim;
	// important !
	ls.dim++;

	const USI np = bulk.numPhase;
	const USI nc = bulk.numCom;
	const USI ncol = nc + 1;
	const USI ncol2 = np * nc + np;
	const USI bsize = ncol * ncol;
	const USI bsize2 = ncol * ncol2;

	OCP_DBL kr, mu, muP;
	OCP_DBL dP;
	OCP_DBL transIJ;

	OCP_USI k_np_j;

	vector<OCP_DBL> bmat(bsize, 0);
	vector<OCP_DBL> bmat2(bsize, 0);
	vector<OCP_DBL> dQdXpB(bsize, 0);
	vector<OCP_DBL> dQdXpW(bsize, 0);
	vector<OCP_DBL> dQdXsB(bsize2, 0);

	for (USI p = 0; p < well.numPerf; p++) {
		OCP_USI k = well.perf[p].location;
		dQdXpB.assign(bsize, 0);
		dQdXpW.assign(bsize, 0);
		dQdXsB.assign(bsize2, 0);

		dP = bulk.P[k] - well.perf[p].P;

		for (USI j = 0; j < np; j++) {
			k_np_j = k * np + j;
			if (!bulk.phaseExist[k_np_j]) continue;

			kr = bulk.kr[k_np_j];
			mu = bulk.mu[k_np_j];
			muP = bulk.muP[k_np_j];

			for (USI i = 0; i < nc; i++) {
				// dQ / dP
				transIJ = well.perf[p].transj[j] * well.perf[p].xi * well.opt.zi[i];
				dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
				dQdXpW[(i + 1) * ncol] += -transIJ;

				// dQ / dS
				for (USI k = 0; k < np; k++) {
					dQdXsB[(i + 1) * ncol2 + k] += CONV1 * CONV2 * well.perf[p].WI * well.perf[p].multiplier * well.perf[p].xi
						* well.opt.zi[i] * bulk.dKr_dS[k * np * np + j * np + k] * dP / mu;
				}
				// dQ / dCij
				for (USI k = 0; k < nc; k++) {
					dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += -transIJ * dP / mu * bulk.muC[k * np * nc + j * nc + k];
				}
			}
		}
		// Bulk to Well
		bmat = dQdXpB;
		DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bulk.dSec_dPri[k * bsize2], 1, bmat.data());
		Dscalar(bsize, dt, bmat.data());
		// Add
		USI ptr = ls.diagPtr[k];
		for (USI i = 0; i < bsize; i++) {
			ls.val[k][ptr * bsize + i] += bmat[i];
		}
		// Insert
		bmat = dQdXpW;
		Dscalar(bsize, dt, bmat.data());
		ls.val[k].insert(ls.val[k].end(), bmat.begin(), bmat.end());
		ls.colId[k].push_back(wId);

		// Well
		switch (well.opt.optMode)
		{
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// Diag
			bmat.assign(bsize, 0);
			for (USI i = 0; i < nc; i++) {
				bmat[0] += dQdXpW[(i + 1) * ncol];
				bmat[(i + 1) * ncol + i + 1] = 1;
			}
			for (USI i = 0; i < bsize; i++) {
				ls.diagVal[wId * bsize + i] += bmat[i];
			}

			// OffDiag
			bmat = dQdXpB;
			DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bulk.dSec_dPri[k * bsize2], 1, bmat.data());
			bmat2.assign(bsize, 0);
			for (USI i = 0; i < nc; i++) {
				// Daxpy(ncol, 1, bmat.data() + (i + 1) * ncol, bmat2.data());
				Daxpy(ncol, well.opt.zi[i], bmat.data() + (i + 1) * ncol, bmat2.data());
			}
			ls.val[wId].insert(ls.val[wId].end(), bmat2.begin(), bmat2.end());
			ls.colId[wId].push_back(k);
			break;

		case BHP_MODE:
			// Diag
			bmat.assign(bsize, 0);
			for (USI i = 0; i < nc + 1; i++) {
				bmat[i * ncol + i] = 1;
			}
			// Add
			for (USI i = 0; i < bsize; i++) {
				ls.diagVal[wId * bsize + i] += bmat[i];
			}
			// OffDiag
			bmat.assign(bsize, 0);
			// Insert
			ls.val[wId].insert(ls.val[wId].end(), bmat.begin(), bmat.end());
			ls.colId[wId].push_back(k);
			break;

		default:
			OCP_ABORT("Wrong Well Opt mode!");
			break;
		}
	}
	OCP_ASSERT(ls.val[wId].size() == numPerf * bsize, "Wrong Size!");
	// Well self
	ls.colId[wId].push_back(wId);
	ls.diagPtr[wId] = well.numPerf;
	ls.val[wId].insert(ls.val[wId].end(),
		ls.diagVal.data() + wId * bsize, ls.diagVal.data() + wId * bsize + bsize);
}


void OCPMethod::AssembleMatPRODWellFIM(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt)
{
	const OCP_USI wId = ls.dim;
	// important !
	ls.dim++;

	const USI np = bulk.numPhase;
	const USI nc = bulk.numCom;
	const USI ncol = nc + 1;
	const USI ncol2 = np * nc + np;
	const USI bsize = ncol * ncol;
	const USI bsize2 = ncol * ncol2;

	OCP_DBL cij, xi, kr, mu, muP, xiP;
	OCP_DBL dP;
	OCP_DBL transIJ;
	OCP_DBL tmp;

	OCP_USI k_np_j;

	vector<OCP_DBL> bmat(bsize, 0);
	vector<OCP_DBL> bmat2(bsize, 0);
	vector<OCP_DBL> dQdXpB(bsize, 0);
	vector<OCP_DBL> dQdXpW(bsize, 0);
	vector<OCP_DBL> dQdXsB(bsize2, 0);

	for (USI p = 0; p < well.numPerf; p++) {
		OCP_USI k = well.perf[p].location;
		dQdXpB.assign(bsize, 0);
		dQdXpW.assign(bsize, 0);
		dQdXsB.assign(bsize2, 0);



		for (USI j = 0; j < np; j++) {
			k_np_j = k * np + j;
			if (!bulk.phaseExist[k_np_j]) continue;

			dP = bulk.Pj[k_np_j] - well.perf[p].P;
			xi = bulk.xi[k_np_j];
			kr = bulk.kr[k_np_j];
			mu = bulk.mu[k_np_j];
			muP = bulk.muP[k_np_j];
			xiP = bulk.xiP[k_np_j];

			for (USI i = 0; i < nc; i++) {
				cij = bulk.cij[k_np_j * nc + i];
				// dQ / dP
				transIJ = well.perf[p].transj[j] * xi * cij;
				dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) + dP * well.perf[p].transj[j] * cij * xiP;
				dQdXpW[(i + 1) * ncol] += -transIJ;

				// dQ / dS
				for (USI k = 0; k < np; k++) {
					tmp = CONV1 * CONV2 * well.perf[p].WI * well.perf[p].multiplier * dP / mu * xi * cij * bulk.dKr_dS[k * np * np + j * np + k];
					// capillary pressure
					tmp += transIJ * bulk.dPcj_dS[k * np * np + j * np + k];
					dQdXsB[(i + 1) * ncol2 + k] += tmp;
				}
				// dQ / dCij
				for (USI k = 0; k < nc; k++) {
					tmp = dP * well.perf[p].transj[j] * cij * (bulk.xiC[k_np_j * nc + k] - xi / mu * bulk.muC[k_np_j * nc + k]);
					if (k == i) {
						tmp += well.perf[p].transj[j] * xi * dP;
					}
					dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
				}
			}
		}
		// Bulk to Well
		bmat = dQdXpB;
		DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bulk.dSec_dPri[k * bsize2], 1, bmat.data());
		Dscalar(bsize, dt, bmat.data());
		// Add
		USI ptr = ls.diagPtr[k];
		for (USI i = 0; i < bsize; i++) {
			ls.val[k][ptr * bsize + i] += bmat[i];
		}
		// Insert
		bmat = dQdXpW;
		Dscalar(bsize, dt, bmat.data());
		ls.val[k].insert(ls.val[k].end(), bmat.begin(), bmat.end());
		ls.colId[k].push_back(wId);

		// Well
		switch (well.opt.optMode)
		{
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// Diag
			bmat.assign(bsize, 0);
			for (USI i = 0; i < nc; i++) {
				bmat[0] += dQdXpW[(i + 1) * ncol] * well.opt.zi[i];
				bmat[(i + 1) * ncol + i + 1] = 1;
			}
			for (USI i = 0; i < bsize; i++) {
				ls.diagVal[wId * bsize + i] += bmat[i];
			}

			// OffDiag
			bmat = dQdXpB;
			DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &bulk.dSec_dPri[k * bsize2], 1, bmat.data());
			bmat2.assign(bsize, 0);
			for (USI i = 0; i < nc; i++) {
				Daxpy(ncol, well.opt.zi[i], bmat.data() + (i + 1) * ncol, bmat2.data());
			}
			ls.val[wId].insert(ls.val[wId].end(), bmat2.begin(), bmat2.end());
			ls.colId[wId].push_back(k);
			break;

		case BHP_MODE:
			// Diag
			bmat.assign(bsize, 0);
			for (USI i = 0; i < nc + 1; i++) {
				bmat[i * ncol + i] = 1;
			}
			// Add
			for (USI i = 0; i < bsize; i++) {
				ls.diagVal[wId * bsize + i] += bmat[i];
			}
			// OffDiag
			bmat.assign(bsize, 0);
			// Insert
			ls.val[wId].insert(ls.val[wId].end(), bmat.begin(), bmat.end());
			ls.colId[wId].push_back(k);
			break;

		default:
			OCP_ABORT("Wrong Well Opt mode!");
			break;
		}
	}
	OCP_ASSERT(ls.val[wId].size() == numPerf * bsize, "Wrong Size!");
	// Well self
	ls.colId[wId].push_back(wId);
	ls.diagPtr[wId] = well.numPerf;
	ls.val[wId].insert(ls.val[wId].end(),
		ls.diagVal.data() + wId * bsize, ls.diagVal.data() + wId * bsize + bsize);

}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/