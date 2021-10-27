/*! \file    OCPOutput.cpp
 *  \brief   OCPOutput class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPOutput.hpp"

void Summary::InputParam(const OutputSummary& summary_param)
{
    FPR  = summary_param.FPR;
    FOPR = summary_param.FOPR;
    FOPT = summary_param.FOPT;
    FGPR = summary_param.FGPR;
    FGPt = summary_param.FGPt;
    FWPR = summary_param.FWPR;
    FWPT = summary_param.FWPT;
    FGIR = summary_param.FGIR;
    FGIT = summary_param.FGIT;
    FWIR = summary_param.FWIR;
    FWIT = summary_param.FWIT;

    WOPR = summary_param.WOPR;
    WOPT = summary_param.WOPT;
    WGPR = summary_param.WGPR;
    WGPT = summary_param.WGPT;
    WWPR = summary_param.WWPR;
    WWPT = summary_param.WWPT;
    WGIR = summary_param.WGIR;
    WGIT = summary_param.WGIT;
    WWIR = summary_param.WWIR;
    WWIT = summary_param.WWIT;

    WBHP = summary_param.WBHP;
    DG = summary_param.DG;

    BPR = summary_param.BPR;
    SOIL = summary_param.SOIL;
    SGAS = summary_param.SGAS;
    SWAT = summary_param.SWAT;

    cout << "Summary::InputParam" << endl;
}

void Summary::Setup(const Reservoir& reservoir, const OCP_DBL& totalTime)
{
    Sumdata.push_back(SumPair("TIME", "  ", "DAY"));
    Sumdata.push_back(SumPair("NRiter", "  ", "  "));
    Sumdata.push_back(SumPair("LSiter", "  ", "  "));
    if (FPR) Sumdata.push_back(SumPair("FPR", "  ", "PSIA"));
    if (FOPR) Sumdata.push_back(SumPair("FOPR", "  ", "STB/DAY"));
    if (FOPT) Sumdata.push_back(SumPair("FOPT", "  ", "STB"));
    if (FGPR) Sumdata.push_back(SumPair("FGPR", "  ", "MSCF/DAY"));
    if (FGPt) Sumdata.push_back(SumPair("FGPT", "  ", "MSCF"));
    if (FWPR) Sumdata.push_back(SumPair("FWPR", "  ", "STB/DAY"));
    if (FWPT) Sumdata.push_back(SumPair("FWPT", "  ", "STB"));
    if (FGIR) Sumdata.push_back(SumPair("FGIR", "  ", "MSCF/DAY"));
    if (FGIT) Sumdata.push_back(SumPair("FGIT", "  ", "MSCF"));
    if (FWIR) Sumdata.push_back(SumPair("FWIR", "  ", "STB/DAY"));
    if (FWIT) Sumdata.push_back(SumPair("FWIT", "  ", "STB"));

    USI num;
    USI wellnum = reservoir.wellgroup.GetWellNum();
    string wellname;

    if (WOPR.activity) {
        if (WOPR.obj[0] == "All") {
            WOPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WOPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WOPR", wellname, "STB/DAY"));
                WOPR.index.push_back(w);
            }
        } else {
            num = WOPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WOPR", WOPR.obj[w], "STB/DAY"));
                WOPR.index.push_back(reservoir.wellgroup.GetIndex(WOPR.obj[w]));
            }
        }
    }

    if (WOPT.activity) {
        if (WOPT.obj[0] == "All") {
            WOPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WOPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WOPT", wellname, "STB"));
                WOPT.index.push_back(w);
            }
        } else {
            num = WOPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WOPT", WOPT.obj[w], "STB"));
                WOPT.index.push_back(reservoir.wellgroup.GetIndex(WOPT.obj[w]));
            }
        }
    }

    if (WGPR.activity) {
        if (WGPR.obj[0] == "All") {
            WGPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WGPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGPR", wellname, "MSCF/DAY"));
                WGPR.index.push_back(w);
            }
        } else {
            num = WGPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WGPR", WGPR.obj[w], "MSCF/DAY"));
                WGPR.index.push_back(reservoir.wellgroup.GetIndex(WGPR.obj[w]));
            }
        }
    }

    if (WGPT.activity) {
        if (WGPT.obj[0] == "All") {
            WGPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WGPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGPT", wellname, "MSCF"));
                WGPT.index.push_back(w);
            }
        } else {
            num = WGPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WGPT", WGPT.obj[w], "MSCF"));
                WGPT.index.push_back(reservoir.wellgroup.GetIndex(WGPT.obj[w]));
            }
        }
    }

    if (WWPR.activity) {
        if (WWPR.obj[0] == "All") {
            WWPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WWPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWPR", wellname, "STB/DAY"));
                WWPR.index.push_back(w);
            }
        } else {
            num = WWPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WWPR", WWPR.obj[w], "STB/DAY"));
                WWPR.index.push_back(reservoir.wellgroup.GetIndex(WWPR.obj[w]));
            }
        }
    }

    if (WWPT.activity) {
        if (WWPT.obj[0] == "All") {
            WWPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WWPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWPT", wellname, "STB"));
                WWPT.index.push_back(w);
            }
        } else {
            num = WWPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WWPT", WWPT.obj[w], "STB"));
                WWPT.index.push_back(reservoir.wellgroup.GetIndex(WWPT.obj[w]));
            }
        }
    }

    if (WGIR.activity) {
        if (WGIR.obj[0] == "All") {
            WGIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WGIR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGIR", wellname, "MSCF/DAY"));
                WGIR.index.push_back(w);
            }
        } else {
            num = WGIR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WGIR", WGIR.obj[w], "MSCF/DAY"));
                WGIR.index.push_back(reservoir.wellgroup.GetIndex(WGIR.obj[w]));
            }
        }
    }

    if (WGIT.activity) {
        if (WGIT.obj[0] == "All") {
            WGIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WGIT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGIT", wellname, "MSCF"));
                WGIT.index.push_back(w);
            }
        } else {
            num = WGIT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WGIT", WGIT.obj[w], "MSCF"));
                WGIT.index.push_back(reservoir.wellgroup.GetIndex(WGIT.obj[w]));
            }
        }
    }

    if (WWIR.activity) {
        if (WWIR.obj[0] == "All") {
            WWIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WWIR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWIR", wellname, "STB/DAY"));
                WWIR.index.push_back(w);
            }
        } else {
            num = WWIR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WWIR", WWIR.obj[w], "STB/DAY"));
                WWIR.index.push_back(reservoir.wellgroup.GetIndex(WWIR.obj[w]));
            }
        }
    }

    if (WWIT.activity) {
        if (WWIT.obj[0] == "All") {
            WWIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WWIT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWIT", wellname, "STB"));
                WWIT.index.push_back(w);
            }
        } else {
            num = WWIT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WWIT", WWIT.obj[w], "STB"));
                WWIT.index.push_back(reservoir.wellgroup.GetIndex(WWIT.obj[w]));
            }
        }
    }

    if (WBHP.activity) {
        if (WBHP.obj[0] == "All") {
            WBHP.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = reservoir.wellgroup.GetWellName(w);
                WBHP.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WBHP", wellname, "PSIA"));
                WBHP.index.push_back(w);
            }
        } else {
            num = WBHP.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumPair("WBHP", WBHP.obj[w], "PSIA"));
                WBHP.index.push_back(reservoir.wellgroup.GetIndex(WBHP.obj[w]));
            }
        }
    }

    if (DG.activity) {
        if (DG.obj[0] == "All") {
            DG.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                string wellname = reservoir.wellgroup.GetWellName(w);
                DG.obj.push_back(wellname);
                USI perfnum = reservoir.wellgroup.GetWellPerfNum(w);
                for (USI p = 0; p < perfnum; p++) {
                    Sumdata.push_back(SumPair("DG", wellname + " Perf" + to_string(p), "PSIA"));
                    DG.index.push_back(w);
                }
            }
        }
        else {
            num = DG.obj.size();
            for (USI w = 0; w < num; w++) {
                USI wId = reservoir.wellgroup.GetIndex(DG.obj[w]);
                USI perfnum = reservoir.wellgroup.GetWellPerfNum(wId);
                for (USI p = 0; p < perfnum; p++) {
                    Sumdata.push_back(SumPair("DG", DG.obj[w] + " P" + to_string(p), "PSIA"));
                    DG.index.push_back(wId);
                }
            }
        }
    }

    if (BPR.activity) {
        num = BPR.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp = "(" + to_string(BPR.obj[i].I) + "," +
                          to_string(BPR.obj[i].J) + "," + to_string(BPR.obj[i].K) +
                          ")";
            Sumdata.push_back(SumPair("BPR", temp, "PSIA"));
            USI I = BPR.obj[i].I - 1;
            USI J = BPR.obj[i].J - 1;
            USI K = BPR.obj[i].K - 1;
            BPR.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SOIL.activity) {
        num = SOIL.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp = "(" + to_string(SOIL.obj[i].I) + "," +
                to_string(SOIL.obj[i].J) + "," + to_string(SOIL.obj[i].K) +
                ")";
            Sumdata.push_back(SumPair("SOIL", temp, "   "));
            USI I = SOIL.obj[i].I - 1;
            USI J = SOIL.obj[i].J - 1;
            USI K = SOIL.obj[i].K - 1;
            SOIL.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SGAS.activity) {
        num = SGAS.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp = "(" + to_string(SGAS.obj[i].I) + "," +
                to_string(SGAS.obj[i].J) + "," + to_string(SGAS.obj[i].K) +
                ")";
            Sumdata.push_back(SumPair("SGAS", temp, "   "));
            USI I = SGAS.obj[i].I - 1;
            USI J = SGAS.obj[i].J - 1;
            USI K = SGAS.obj[i].K - 1;
            SGAS.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SWAT.activity) {
        num = SWAT.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp = "(" + to_string(SWAT.obj[i].I) + "," +
                to_string(SWAT.obj[i].J) + "," + to_string(SWAT.obj[i].K) +
                ")";
            Sumdata.push_back(SumPair("SWAT", temp, "   "));
            USI I = SWAT.obj[i].I - 1;
            USI J = SWAT.obj[i].J - 1;
            USI K = SWAT.obj[i].K - 1;
            SWAT.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    // Allocate memory
    USI rs = totalTime / 0.1;
    USI cs = Sumdata.size();
    for (USI i = 0; i < cs; i++) {
        Sumdata[i].val.reserve(rs);
    }

    cout << "Summary::Setup" << endl;
}

void Summary::SetVal(const Reservoir& rs, const OCP_Control& ctrl)
{
    USI n = 0;

    // TIME
    Sumdata[n++].val.push_back(ctrl.GetCurTime());
    // NRiter
    Sumdata[n++].val.push_back(ctrl.GetNRiter());
    // LSiter
    Sumdata[n++].val.push_back(ctrl.GetLSiter());

    // FPR
    if (FPR) Sumdata[n++].val.push_back(rs.bulk.CalFPR());
    if (FOPR) Sumdata[n++].val.push_back(rs.wellgroup.GetFOPR());
    if (FOPT) Sumdata[n++].val.push_back(rs.wellgroup.GetFOPT());
    if (FGPR) Sumdata[n++].val.push_back(rs.wellgroup.GetFGPR());
    if (FGPt) Sumdata[n++].val.push_back(rs.wellgroup.GetFGPT());
    if (FWPR) Sumdata[n++].val.push_back(rs.wellgroup.GetFWPR());
    if (FWPT) Sumdata[n++].val.push_back(rs.wellgroup.GetFWPT());
    if (FGIR) Sumdata[n++].val.push_back(rs.wellgroup.GetFGIR());
    if (FGIT) Sumdata[n++].val.push_back(rs.wellgroup.GetFGIT());
    if (FWIR) Sumdata[n++].val.push_back(rs.wellgroup.GetFWIR());
    if (FWIT) Sumdata[n++].val.push_back(rs.wellgroup.GetFWIT());

    USI len = 0;
    // WOPR
    len = WOPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWOPR(WOPR.index[w]));

    // WOPT
    len = WOPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWOPT(WOPT.index[w]));

    // WGPR
    len = WGPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWGPR(WGPR.index[w]));

    // WGPT
    len = WGPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWGPT(WGPT.index[w]));

    // WWPR
    len = WWPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWWPR(WWPR.index[w]));

    // WWPT
    len = WWPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWWPT(WWPT.index[w]));

    // WGIR
    len = WGIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWGIR(WGIR.index[w]));

    // WGIT
    len = WGIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWGIT(WGIT.index[w]));

    // WWIR
    len = WWIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWWIR(WWIR.index[w]));

    // WWIT
    len = WWIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWWIT(WWIT.index[w]));

    // WBHP
    len = WBHP.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.wellgroup.GetWBHP(WBHP.index[w]));
    
    // DG
    len = DG.obj.size();
    for (USI w = 0; w < len; w++) {
        USI numperf = rs.wellgroup.GetWellPerfNum(DG.index[w]);
        for (USI p = 0; p < numperf; p++) {
            Sumdata[n++].val.push_back(rs.wellgroup.GetWellDg(DG.index[w], p));
        }
    }

    // BPR
    len = BPR.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(rs.bulk.GetP(BPR.index[i]));

    // SOIL
    len = SOIL.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(rs.bulk.GetSOIL(SOIL.index[i]));

    // SGAS
    len = SGAS.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(rs.bulk.GetSGAS(SGAS.index[i]));

    // SWAT
    len = SWAT.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(rs.bulk.GetSWAT(SWAT.index[i]));

}

void Summary::PrintInfo(const string& dir) const
{
    string   FileOut = dir + "SUMMARY.out";
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    USI ns  = 10;
    USI col = 10;
    USI row = 0;
    USI num = Sumdata.size();
    USI len = Sumdata[0].val.size();
    USI id  = 0;
    USI ID  = 1;

    while (id != num) {

        outF << "The " << ++row << "th Row\n";

        // Item
        // Time
        outF << "\t" << setw(10) << Sumdata[0].Item;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id].Item;
            id++;
            if (id == num) break;
        }
        outF << "\n";

        // Unit
        // Time
        outF << "\t" << setw(10) << Sumdata[0].Unit;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id].Unit;
            id++;
            if (id == num) break;
        }
        outF << "\n";

        // Obj
        // Time
        outF << "\t" << setw(ns) << Sumdata[0].Obj;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id].Obj;
            id++;
            if (id == num) break;
        }
        outF << "\n";

        // data
        for (USI l = 0; l < len; l++) {

            // Time
            outF << "\t" << setw(ns)  << Sumdata[0].val[l];

            id = ID;
            for (USI i = 1; i < col; i++) {
                outF << "\t" << setw(ns) << Sumdata[id].val[l];
                id++;
                if (id == num) break;
            }
            outF << "\n";
        }

        ID += (col - 1);
    }

    outF.close();
}

void CriticalInfo::Setup(const Reservoir& reservoir, const OCP_DBL& totalTime)
{
    // Allocate memory
    USI rc = totalTime / 0.1;
    time.reserve(rc);
    dt.reserve(rc);
    dPmax.reserve(rc);
    dVmax.reserve(rc);
    dSmax.reserve(rc);
    dNmax.reserve(rc);
    cfl.reserve(rc);
}

void CriticalInfo::SetVal(const Reservoir& reservoir, const OCP_Control& ctrl)
{
    time.push_back(ctrl.GetCurTime());
    dt.push_back(ctrl.GetLastCurDt());
    dPmax.push_back(reservoir.bulk.GetdPmax());
    dVmax.push_back(reservoir.bulk.GetdVmax());
    dSmax.push_back(reservoir.bulk.GetdSmax());
    dNmax.push_back(reservoir.bulk.GetdNmax());
    cfl.push_back(reservoir.cfl);
}

void CriticalInfo::PrintInfo(const string& dir) const
{
    string   FileOut = dir + "FastReview.out";
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    // Item
    outF << "\t" << setw(10) << "Time";
    outF << "\t" << setw(10) << "dt";
    outF << "\t" << setw(10) << "dPmax";
    outF << "\t" << setw(10) << "dVmax";
    outF << "\t" << setw(10) << "dSmax";
    outF << "\t" << setw(10) << "dNmax";
    outF << "\t" << setw(10) << "CFL\n";
    // Unit
    outF << "\t" << setw(10) << "Days";
    outF << "\t" << setw(10) << "Days";
    outF << "\t" << setw(10) << "Psia";
    outF << "\t" << setw(10) << "    ";
    outF << "\t" << setw(10) << "    ";
    outF << "\t" << setw(10) << "    ";
    outF << "\t" << setw(10) << "    \n\n";

    USI n = time.size();
    for (USI i = 0; i < n; i++) {
        outF << "\t" << setw(10) << time[i];
        outF << "\t" << setw(10) << dt[i];
        outF << "\t" << setw(10) << dPmax[i];
        outF << "\t" << setw(10) << dVmax[i];
        outF << "\t" << setw(10) << dSmax[i];
        outF << "\t" << setw(10) << dNmax[i];
        outF << "\t" << setw(10) << cfl[i];
        outF << "\n";
    }

    outF.close();
}

void DetailInfo::InputParam(const OutputDetail& detail_param)
{
    PRE  = detail_param.PRE;
    PGAS = detail_param.PGAS;
    PWAT = detail_param.PWAT;
    SOIL = detail_param.SOIL;
    SGAS = detail_param.SGAS;
    SWAT = detail_param.SWAT;
    DENO = detail_param.DENO;
    DENG = detail_param.DENG;
    DENW = detail_param.DENW;
}

void DetailInfo::Setup(const string& dir)
{
    string   FileOut = dir + "RPT.out";
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }
    outF.close();
}

void DetailInfo::PrintInfo(const string& dir, const Reservoir& rs,
                           const OCP_DBL& days) const
{
    string   FileOut = dir + "RPT.out";
    ofstream outF;
    outF.open(FileOut, ios::app);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    USI     nx  = rs.grid.GetGridNx();
    USI     ny  = rs.grid.GetGridNy();
    OCP_USI num = rs.grid.GetGridNum();
    OCP_USI bId;

    // PRESSURE
    if (PRE) {
        outF << "PRESSURE : psia"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++) {
            if (i % nx == 0) outF << "\n";
            if (i % (nx * ny) == 0) outF << "\n";
            if (rs.grid.MapG2B(i).GetAct()) {
                bId = rs.grid.MapG2B(i).GetId();
                outF << fixed << setprecision(3) << rs.bulk.P[bId] << "   ";
            } else {
                outF << "N   ";
            }
        }
        outF << "\n\n";
    }

    outF.close();
}

void OCP_Output::InputParam(const ParamOutput& paramOutput)
{
    summary.InputParam(paramOutput.summary);
    dtlInfo.InputParam(paramOutput.detailInfo);
}

void OCP_Output::Setup(const Reservoir& reservoir, const OCP_Control& ctrl)
{
    wordDir = ctrl.workDir;
    summary.Setup(reservoir, ctrl.criticalTime.back());
    crtInfo.Setup(reservoir, ctrl.criticalTime.back());
    dtlInfo.Setup(wordDir);
}

void OCP_Output::SetVal(const Reservoir& reservoir, const OCP_Control& ctrl)
{
    summary.SetVal(reservoir, ctrl);
    crtInfo.SetVal(reservoir, ctrl);
}

void OCP_Output::PrintInfo() const
{
    summary.PrintInfo(wordDir);
    crtInfo.PrintInfo(wordDir);
}

void OCP_Output::PrintInfoSched(const Reservoir& rs, const OCP_Control& ctrl, const OCP_DBL& time) const
{   
    OCP_DBL days = ctrl.current_time;
    cout << fixed << setprecision(3) << days << " Days\t";
    cout << ctrl.tstep << "\t";
    cout << time / 1000 << "s";
    cout << "\n";
    dtlInfo.PrintInfo(wordDir, rs, days);
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/