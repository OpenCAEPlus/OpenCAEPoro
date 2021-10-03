#include <cassert>
#include <cmath>
#include <ctime>

#include "Connection_BB.hpp"
#include "OpenCAEPoro_consts.hpp"

// Active Conn & Active Bulk

void Connection_BB::setup(const Grid& myGrid, const Bulk& myBulk)
{
    initSize(myBulk);
    initActive(myGrid, myBulk.Np);
    getIteratorActive();
    calAreaActive(myGrid, myBulk);
}

void Connection_BB::initSize(const Bulk& myBulk)
{
    ActiveConnNum = 0;
    ActiveBulkNum = myBulk.Num;

    Neighbor.resize(ActiveBulkNum);
    SelfPtr.resize(ActiveBulkNum);
    NeighborNum.resize(ActiveBulkNum);
}

void Connection_BB::initActive(const Grid& myGrid, int np)
{
    int nx   = myGrid.Nx;
    int ny   = myGrid.Ny;
    int nz   = myGrid.Nz;
    int nxny = nx * ny;

    // bIdb : begin id in bulk
    // bIdg : begin id in grid
    // eIdb : end id in bulk
    int eIdb;
    for (int bIdb = 0; bIdb < ActiveBulkNum; bIdb++) {
        int bIdg = myGrid.ActiveMap_B2G[bIdb];
        int k    = bIdg / nxny;
        int j    = (bIdg - k * nxny) / nx;
        int i    = bIdg - k * nxny - j * nx;

        int count = 0;
        // up
        if (k > 0) {
            eIdb = myGrid.ActiveMap_G2B[bIdg - nxny];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // back
        if (j > 0) {
            eIdb = myGrid.ActiveMap_G2B[bIdg - nx];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // left
        if (i > 0) {
            eIdb = myGrid.ActiveMap_G2B[bIdg - 1];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        ActiveConnNum += count;
        // self
        Neighbor[bIdb].push_back(bIdb);
        SelfPtr[bIdb] = count;
        count++;

        // right
        if (i < nx - 1) {
            eIdb = myGrid.ActiveMap_G2B[bIdg + 1];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // front
        if (j < ny - 1) {
            eIdb = myGrid.ActiveMap_G2B[bIdg + nx];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // down
        if (k < nz - 1) {
            eIdb = myGrid.ActiveMap_G2B[bIdg + nxny];
            if (eIdb != -1) {
                Neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        NeighborNum[bIdb] = count;
    }

    Upblock.resize(ActiveConnNum * np);
    Upblock_Rho.resize(ActiveConnNum * np);
    Upblock_Trans.resize(ActiveConnNum * np);
    Upblock_Velocity.resize(ActiveConnNum * np);
}

void Connection_BB::getIteratorActive()
{
    Iterator.reserve(ActiveConnNum);
    // generate iterator for BB from Iterator
    for (int bId = 0; bId < ActiveBulkNum; bId++) {
        int beginIt = SelfPtr[bId] + 1;
        int nbc     = NeighborNum[bId];

        for (int c = beginIt; c < nbc; c++) {
            int eId = Neighbor[bId][c];
            Iterator.push_back(BB_Pair(bId, eId));
        }
    }

    assert(Iterator.size() == ActiveConnNum);
}

void Connection_BB::calAreaActive(const Grid& myGrid, const Bulk& myBulk)
{
    // calculate efficient area of interface of bulk to bulk
    // using Iterator
    Area.reserve(ActiveConnNum);

    for (int n = 0; n < ActiveConnNum; n++) {
        int bIdb = Iterator[n].BId;
        int eIdb = Iterator[n].EId;
        // Area.push_back(1);
        Area.push_back(calAkd(myGrid, myBulk, bIdb, eIdb));
    }

    assert(Area.size() == ActiveConnNum);
}

OCP_DBL Connection_BB::calAkd(const Grid& myGrid, const Bulk& myBulk, int bIdb, int eIdb)
{
    int bIdg = myGrid.ActiveMap_B2G[bIdb];
    int eIdg = myGrid.ActiveMap_B2G[eIdb];
    int diff = bIdg - eIdg;
    if (diff < 0) diff *= -1;

    if (diff == 1) {
        // x - direction
        OCP_DBL T1 = myBulk.Rock_Kx[bIdb] * myBulk.Ntg[bIdb] * myBulk.Dy[bIdb] *
                    myBulk.Dz[bIdb] / myBulk.Dx[bIdb];
        OCP_DBL T2 = myBulk.Rock_Kx[eIdb] * myBulk.Ntg[eIdb] * myBulk.Dy[eIdb] *
                    myBulk.Dz[eIdb] / myBulk.Dx[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else if (diff == myGrid.Nx) {
        // y - direction
        OCP_DBL T1 = myBulk.Rock_Ky[bIdb] * myBulk.Ntg[bIdb] * myBulk.Dz[bIdb] *
                    myBulk.Dx[bIdb] / myBulk.Dy[bIdb];
        OCP_DBL T2 = myBulk.Rock_Ky[eIdb] * myBulk.Ntg[eIdb] * myBulk.Dz[eIdb] *
                    myBulk.Dx[eIdb] / myBulk.Dy[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else if (diff == myGrid.Nx * myGrid.Ny) {
        // z - direction  ----  no Ntg
        OCP_DBL T1 =
            myBulk.Rock_Kz[bIdb] * myBulk.Dx[bIdb] * myBulk.Dy[bIdb] / myBulk.Dz[bIdb];
        OCP_DBL T2 =
            myBulk.Rock_Kz[eIdb] * myBulk.Dx[eIdb] * myBulk.Dy[eIdb] / myBulk.Dz[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else {
        // Wrong
        ERRORcheck("Wrong bId and eId in function");
        exit(0);
    }
}

// Connection function, no matter what grid is

OCP_DBL Connection_BB::calCFL(Bulk& myBulk, OCP_DBL dt)
{
    int    np   = myBulk.Np;
    OCP_DBL cfl  = 0;
    OCP_DBL temp = 0;
    for (int c = 0; c < ActiveConnNum; c++) {

        for (int j = 0; j < np; j++) {
            int uId = Upblock[c * np + j];

            if (myBulk.PhaseExist[uId]) {
                temp = fabs(Upblock_Velocity[c * np + j]) * dt;
                temp /= myBulk.Vj[uId * np + j];
                if (cfl < temp) cfl = temp;
            }
        }
    }
    return cfl;
}

void Connection_BB::calFlux(const Bulk& myBulk)
{
    // calculate a step flux using Iterator
    int    bId, eId, uId;
    int    bId_np_j, eId_np_j;
    OCP_DBL Pbegin, Pend, rho;
    int    np = myBulk.Np;

    for (int c = 0; c < ActiveConnNum; c++) {
        bId        = Iterator[c].BId;
        eId        = Iterator[c].EId;
        OCP_DBL Akd = Area[c];

        for (int j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            bool exbegin = myBulk.PhaseExist[bId_np_j];
            bool exend   = myBulk.PhaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = (myBulk.Rho[bId_np_j] + myBulk.Rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.P[eId];
                rho    = myBulk.Rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                Pbegin = myBulk.P[bId];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = myBulk.Rho[eId_np_j];
            } else {
                Upblock[c * np + j] = bId;
                continue;
            }

            uId         = bId;
            bool   exup = exbegin;
            OCP_DBL dP   = (Pbegin - GRAVITY_FACTOR * rho * myBulk.Depth[bId]) -
                        (Pend - GRAVITY_FACTOR * rho * myBulk.Depth[eId]);
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            }
            Upblock_Rho[c * np + j] = rho;
            Upblock[c * np + j]     = uId;
            int    uId_np_j         = uId * np + j;
            OCP_DBL trans =
                CONV1 * CONV2 * Akd * myBulk.Kr[uId_np_j] / myBulk.Mu[uId_np_j];
            Upblock_Trans[c * np + j] = trans;

            if (exup) {
                Upblock_Velocity[c * np + j] = trans * dP;
            } else {
                Upblock_Velocity[c * np + j] = 0;
            }
        }
    }
}

void Connection_BB::massConserve(Bulk& myBulk, OCP_DBL dt)
{
    int np = myBulk.Np;
    int nc = myBulk.Nc;

    for (int c = 0; c < ActiveConnNum; c++) {
        int bId = Iterator[c].BId;
        int eId = Iterator[c].EId;

        for (int j = 0; j < np; j++) {
            int uId = Upblock[c * np + j];
            if (!myBulk.PhaseExist[uId * np + j]) continue;

            int    uId_np_j      = uId * np + j;
            OCP_DBL phaseVelocity = Upblock_Velocity[c * np + j];
            for (int i = 0; i < nc; i++) {
                OCP_DBL dNi = dt * phaseVelocity * myBulk.Xi[uId_np_j] *
                             myBulk.Cij[uId_np_j * nc + i];
                myBulk.Ni[eId * nc + i] += dNi;
                myBulk.Ni[bId * nc + i] -= dNi;
            }
        }
    }
}

void Connection_BB::assembleMat(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt) const
{
    // accumulate term
    OCP_DBL Vp0, Vp, Vf, Vfp, P;
    OCP_DBL cr = myBulk.Rock_C1;
    for (int n = 0; n < ActiveBulkNum; n++) {
        Vp0 = myBulk.Rock_VpInit[n];
        Vp  = myBulk.Rock_Vp[n];
        Vfp = myBulk.Vfp[n];
        P   = myBulk.P[n];
        Vf  = myBulk.Vf[n];


		OCP_DBL temp = cr * Vp0 - Vfp;
		mySolver.DiagVal[n] = temp;
		mySolver.b[n] = temp * P + dt * (Vf - Vp);
	}

	// check 
	//ofstream outb("testb.dat");
	//if (!outb.is_open())
	//	cout << "Can not open " << "testb.dat" << endl;
	//outb << mySolver.Dim << endl;
	//for (int i = 0; i < mySolver.Dim; i++)
	//	outb << mySolver.b[i] << endl;
	//outb.close();

    // flux term
    int    bId, eId, uId;
    int    np = myBulk.Np;
    int    nc = myBulk.Nc;
    OCP_DBL valupi, valdowni;
    OCP_DBL valup, rhsup, valdown, rhsdown;
    int    lastbId = -1;
    for (int c = 0; c < ActiveConnNum; c++) {
        bId = Iterator[c].BId;
        eId = Iterator[c].EId;
        valup = 0;
        rhsup = 0;
        valdown = 0;
        rhsdown = 0;

        for (int j = 0; j < np; j++) {
            uId = Upblock[c * np + j];
            if (!myBulk.PhaseExist[uId * np + j]) continue;

            valupi   = 0;
            valdowni = 0;

            for (int i = 0; i < nc; i++) {
                valupi +=
                    myBulk.Vfi[bId * nc + i] * myBulk.Cij[uId * np * nc + j * nc + i];
                valdowni +=
                    myBulk.Vfi[eId * nc + i] * myBulk.Cij[uId * np * nc + j * nc + i];
            }
            OCP_DBL dD   = myBulk.Depth[bId] - myBulk.Depth[eId];
            OCP_DBL dPc  = myBulk.Pc[bId * np + j] - myBulk.Pc[eId * np + j];
            OCP_DBL temp = myBulk.Xi[uId * np + j] * Upblock_Trans[c * np + j] * dt;
            valup += temp * valupi;
            valdown += temp * valdowni;
            temp *= Upblock_Rho[c * np + j] * GRAVITY_FACTOR * dD - dPc;
            rhsup += temp * valupi;
            rhsdown -= temp * valdowni;
        }

        int diagptr = SelfPtr[bId];
        if (bId != lastbId) {
            // new bulk
            assert(mySolver.Val[bId].size() == diagptr);
            mySolver.Val[bId].push_back(mySolver.DiagVal[bId]);
            lastbId = bId;
        }

        mySolver.Val[bId][diagptr] += valup;
        mySolver.Val[bId].push_back(-valup);
        mySolver.Val[eId].push_back(-valdown);
        mySolver.DiagVal[eId] += valdown;
        mySolver.b[bId] += rhsup;
        mySolver.b[eId] += rhsdown;
    }

    // Add the rest of diag value
    // important!
    for (int n = 0; n < ActiveBulkNum; n++) {
        if (mySolver.Val[n].size() == SelfPtr[n])
            mySolver.Val[n].push_back(mySolver.DiagVal[n]);
    }
}

void Connection_BB::getConnectionInfo()
{
    for (int i = 0; i < ActiveBulkNum; i++) {
        cout << "(" << i << ")"
             << "\t";

        for (auto v : Neighbor[i]) {
            cout << v << "\t";
        }
        cout << "[" << SelfPtr[i] << "]";
        cout << "\t" << NeighborNum[i];
        cout << "\n";
    }

    for (int i = 0; i < ActiveConnNum; i++) {
        cout << Iterator[i].BId << "\t" << Iterator[i].EId << "\n";
    }
}
