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
    FTR  = summary_param.FTR;
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
    DG   = summary_param.DG;

    BPR  = summary_param.BPR;
    SOIL = summary_param.SOIL;
    SGAS = summary_param.SGAS;
    SWAT = summary_param.SWAT;

    // cout << "Summary::InputParam" << endl;
}

void Summary::Setup(const Reservoir& rs, const OCP_DBL& totalTime)
{
    Sumdata.push_back(SumItem("TIME", "  ", "DAY", "fixed"));
    Sumdata.push_back(SumItem("NRiter", "  ", "  ", "int"));
    Sumdata.push_back(SumItem("LSiter", "  ", "  ", "int"));
    if (FPR) Sumdata.push_back(SumItem("FPR", "  ", "PSIA", "float"));
    if (FTR) Sumdata.push_back(SumItem("FTR", "  ", "F", "float"));
    if (FOPR) Sumdata.push_back(SumItem("FOPR", "  ", "STB/DAY", "float"));
    if (FOPT) Sumdata.push_back(SumItem("FOPT", "  ", "STB", "float"));
    if (FGPR) Sumdata.push_back(SumItem("FGPR", "  ", "MSCF/DAY", "float"));
    if (FGPt) Sumdata.push_back(SumItem("FGPT", "  ", "MSCF", "float"));
    if (FWPR) Sumdata.push_back(SumItem("FWPR", "  ", "STB/DAY", "float"));
    if (FWPT) Sumdata.push_back(SumItem("FWPT", "  ", "STB", "float"));
    if (FGIR) Sumdata.push_back(SumItem("FGIR", "  ", "MSCF/DAY", "float"));
    if (FGIT) Sumdata.push_back(SumItem("FGIT", "  ", "MSCF", "float"));
    if (FWIR) Sumdata.push_back(SumItem("FWIR", "  ", "STB/DAY", "float"));
    if (FWIT) Sumdata.push_back(SumItem("FWIT", "  ", "STB", "float"));

    const Grid&     initGrid = rs.grid;
    const AllWells& wells    = rs.allWells;

    const USI sp = initGrid.GetNumDigitIJK();

    const USI wellnum = wells.GetWellNum();
    string    wellname;
    USI       num;

    if (WOPR.activity) {
        if (WOPR.obj[0] == "All") {
            WOPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WOPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WOPR", wellname, "STB/DAY", "float"));
                WOPR.index.push_back(w);
            }
        } else {
            num = WOPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WOPR", WOPR.obj[w], "STB/DAY", "float"));
                WOPR.index.push_back(wells.GetIndex(WOPR.obj[w]));
            }
        }
    }

    if (WOPT.activity) {
        if (WOPT.obj[0] == "All") {
            WOPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WOPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WOPT", wellname, "STB", "float"));
                WOPT.index.push_back(w);
            }
        } else {
            num = WOPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WOPT", WOPT.obj[w], "STB", "float"));
                WOPT.index.push_back(wells.GetIndex(WOPT.obj[w]));
            }
        }
    }

    if (WGPR.activity) {
        if (WGPR.obj[0] == "All") {
            WGPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGPR", wellname, "MSCF/DAY", "float"));
                WGPR.index.push_back(w);
            }
        } else {
            num = WGPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WGPR", WGPR.obj[w], "MSCF/DAY", "float"));
                WGPR.index.push_back(wells.GetIndex(WGPR.obj[w]));
            }
        }
    }

    if (WGPT.activity) {
        if (WGPT.obj[0] == "All") {
            WGPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGPT", wellname, "MSCF", "float"));
                WGPT.index.push_back(w);
            }
        } else {
            num = WGPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WGPT", WGPT.obj[w], "MSCF", "float"));
                WGPT.index.push_back(wells.GetIndex(WGPT.obj[w]));
            }
        }
    }

    if (WWPR.activity) {
        if (WWPR.obj[0] == "All") {
            WWPR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWPR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWPR", wellname, "STB/DAY", "float"));
                WWPR.index.push_back(w);
            }
        } else {
            num = WWPR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WWPR", WWPR.obj[w], "STB/DAY", "float"));
                WWPR.index.push_back(wells.GetIndex(WWPR.obj[w]));
            }
        }
    }

    if (WWPT.activity) {
        if (WWPT.obj[0] == "All") {
            WWPT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWPT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWPT", wellname, "STB", "float"));
                WWPT.index.push_back(w);
            }
        } else {
            num = WWPT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WWPT", WWPT.obj[w], "STB", "float"));
                WWPT.index.push_back(wells.GetIndex(WWPT.obj[w]));
            }
        }
    }

    if (WGIR.activity) {
        if (WGIR.obj[0] == "All") {
            WGIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGIR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGIR", wellname, "MSCF/DAY", "float"));
                WGIR.index.push_back(w);
            }
        } else {
            num = WGIR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WGIR", WGIR.obj[w], "MSCF/DAY", "float"));
                WGIR.index.push_back(wells.GetIndex(WGIR.obj[w]));
            }
        }
    }

    if (WGIT.activity) {
        if (WGIT.obj[0] == "All") {
            WGIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WGIT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WGIT", wellname, "MSCF", "float"));
                WGIT.index.push_back(w);
            }
        } else {
            num = WGIT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WGIT", WGIT.obj[w], "MSCF", "float"));
                WGIT.index.push_back(wells.GetIndex(WGIT.obj[w]));
            }
        }
    }

    if (WWIR.activity) {
        if (WWIR.obj[0] == "All") {
            WWIR.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWIR.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWIR", wellname, "STB/DAY", "float"));
                WWIR.index.push_back(w);
            }
        } else {
            num = WWIR.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WWIR", WWIR.obj[w], "STB/DAY", "float"));
                WWIR.index.push_back(wells.GetIndex(WWIR.obj[w]));
            }
        }
    }

    if (WWIT.activity) {
        if (WWIT.obj[0] == "All") {
            WWIT.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WWIT.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WWIT", wellname, "STB", "float"));
                WWIT.index.push_back(w);
            }
        } else {
            num = WWIT.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WWIT", WWIT.obj[w], "STB", "float"));
                WWIT.index.push_back(wells.GetIndex(WWIT.obj[w]));
            }
        }
    }

    if (WBHP.activity) {
        if (WBHP.obj[0] == "All") {
            WBHP.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                wellname = wells.GetWellName(w);
                WBHP.obj.push_back(wellname);
                Sumdata.push_back(SumItem("WBHP", wellname, "PSIA", "float"));
                WBHP.index.push_back(w);
            }
        } else {
            num = WBHP.obj.size();
            for (USI w = 0; w < num; w++) {
                Sumdata.push_back(SumItem("WBHP", WBHP.obj[w], "PSIA", "float"));
                WBHP.index.push_back(wells.GetIndex(WBHP.obj[w]));
            }
        }
    }

    if (DG.activity) {
        if (DG.obj[0] == "All") {
            DG.obj.clear();
            for (USI w = 0; w < wellnum; w++) {
                string wellname = wells.GetWellName(w);
                DG.obj.push_back(wellname);
                USI perfnum = wells.GetWellPerfNum(w);
                for (USI p = 0; p < perfnum; p++) {
                    Sumdata.push_back(SumItem("DG", wellname + " Perf" + to_string(p),
                                              "PSIA", "float"));
                    DG.index.push_back(w);
                }
            }
        } else {
            num = DG.obj.size();
            for (USI w = 0; w < num; w++) {
                USI wId     = wells.GetIndex(DG.obj[w]);
                USI perfnum = wells.GetWellPerfNum(wId);
                for (USI p = 0; p < perfnum; p++) {
                    Sumdata.push_back(SumItem("DG", DG.obj[w] + " P" + to_string(p),
                                              "PSIA", "float"));
                    DG.index.push_back(wId);
                }
            }
        }
    }

    if (BPR.activity) {
        num = BPR.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp = GetIJKformat(to_string(BPR.obj[i].I), to_string(BPR.obj[i].J),
                                       to_string(BPR.obj[i].K), sp);
            Sumdata.push_back(SumItem("BPR", temp, "PSIA", "float"));
            USI I = BPR.obj[i].I - 1;
            USI J = BPR.obj[i].J - 1;
            USI K = BPR.obj[i].K - 1;

            const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
            if (tarId < 0)
                OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
                            "in BPR in SUMMARY");
            else
                BPR.index.push_back(tarId);
        }
    }

    if (SOIL.activity) {
        num = SOIL.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp =
                GetIJKformat(to_string(SOIL.obj[i].I), to_string(SOIL.obj[i].J),
                             to_string(SOIL.obj[i].K), sp);
            Sumdata.push_back(SumItem("SOIL", temp, "   ", "float"));
            USI I = SOIL.obj[i].I - 1;
            USI J = SOIL.obj[i].J - 1;
            USI K = SOIL.obj[i].K - 1;

            const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
            if (tarId < 0)
                OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
                            "in SOIL in SUMMARY");
            else
                SOIL.index.push_back(tarId);
        }
    }

    if (SGAS.activity) {
        num = SGAS.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp =
                GetIJKformat(to_string(SGAS.obj[i].I), to_string(SGAS.obj[i].J),
                             to_string(SGAS.obj[i].K), sp);
            Sumdata.push_back(SumItem("SGAS", temp, "   ", "float"));
            USI I = SGAS.obj[i].I - 1;
            USI J = SGAS.obj[i].J - 1;
            USI K = SGAS.obj[i].K - 1;

            const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
            if (tarId < 0)
                OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
                            "in SOIL in SUMMARY");
            else
                SGAS.index.push_back(tarId);
        }
    }

    if (SWAT.activity) {
        num = SWAT.obj.size();
        for (USI i = 0; i < num; i++) {
            string temp =
                GetIJKformat(to_string(SWAT.obj[i].I), to_string(SWAT.obj[i].J),
                             to_string(SWAT.obj[i].K), sp);
            Sumdata.push_back(SumItem("SWAT", temp, "   ", "float"));
            USI I = SWAT.obj[i].I - 1;
            USI J = SWAT.obj[i].J - 1;
            USI K = SWAT.obj[i].K - 1;

            const OCP_INT tarId = initGrid.GetActIndex(I, J, K);
            if (tarId < 0)
                OCP_WARNING("Non Fluid Grid: " + GetIJKformat(I + 1, J + 1, K + 1, sp) +
                            "in SWAT in SUMMARY");
            else
                SWAT.index.push_back(tarId);
        }
    }

    // Allocate memory
    const USI maxRowNum = totalTime / 0.1;
    const USI cs        = Sumdata.size();
    for (USI i = 0; i < cs; i++) {
        Sumdata[i].val.reserve(maxRowNum);
    }

    // cout << "Summary::Setup" << endl;
}

void Summary::SetVal(const Reservoir& rs, const OCPControl& ctrl)
{
    const Bulk&     bulk  = rs.bulk;
    const AllWells& wells = rs.allWells;

    USI n = 0;

    // TIME
    Sumdata[n++].val.push_back(ctrl.GetCurTime());
    // NRiter
    Sumdata[n++].val.push_back(ctrl.GetNRiterT());
    // LSiter
    Sumdata[n++].val.push_back(ctrl.GetLSiterT());

    // FPR
    if (FPR) Sumdata[n++].val.push_back(bulk.CalFPR());
    if (FTR) Sumdata[n++].val.push_back(bulk.CalFTR());
    if (FOPR) Sumdata[n++].val.push_back(wells.GetFOPR());
    if (FOPT) Sumdata[n++].val.push_back(wells.GetFOPT());
    if (FGPR) Sumdata[n++].val.push_back(wells.GetFGPR());
    if (FGPt) Sumdata[n++].val.push_back(wells.GetFGPT());
    if (FWPR) Sumdata[n++].val.push_back(wells.GetFWPR());
    if (FWPT) Sumdata[n++].val.push_back(wells.GetFWPT());
    if (FGIR) Sumdata[n++].val.push_back(wells.GetFGIR());
    if (FGIT) Sumdata[n++].val.push_back(wells.GetFGIT());
    if (FWIR) Sumdata[n++].val.push_back(wells.GetFWIR());
    if (FWIT) Sumdata[n++].val.push_back(wells.GetFWIT());

    USI len = 0;
    // WOPR
    len = WOPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWOPR(WOPR.index[w]));

    // WOPT
    len = WOPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWOPT(WOPT.index[w]));

    // WGPR
    len = WGPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGPR(WGPR.index[w]));

    // WGPT
    len = WGPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGPT(WGPT.index[w]));

    // WWPR
    len = WWPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWPR(WWPR.index[w]));

    // WWPT
    len = WWPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWPT(WWPT.index[w]));

    // WGIR
    len = WGIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGIR(WGIR.index[w]));

    // WGIT
    len = WGIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWGIT(WGIT.index[w]));

    // WWIR
    len = WWIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWIR(WWIR.index[w]));

    // WWIT
    len = WWIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWWIT(WWIT.index[w]));

    // WBHP
    len = WBHP.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(wells.GetWBHP(WBHP.index[w]));

    // DG
    len = DG.obj.size();
    for (USI w = 0; w < len; w++) {
        USI numperf = rs.allWells.GetWellPerfNum(DG.index[w]);
        for (USI p = 0; p < numperf; p++) {
            Sumdata[n++].val.push_back(wells.GetWellDG(DG.index[w], p));
        }
    }

    // BPR
    len = BPR.index.size();
    for (USI i = 0; i < len; i++) Sumdata[n++].val.push_back(bulk.GetP(BPR.index[i]));

    // SOIL
    len = SOIL.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(bulk.GetSOIL(SOIL.index[i]));

    // SGAS
    len = SGAS.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(bulk.GetSGAS(SGAS.index[i]));

    // SWAT
    len = SWAT.index.size();
    for (USI i = 0; i < len; i++)
        Sumdata[n++].val.push_back(bulk.GetSWAT(SWAT.index[i]));
}

/// Write output information in the dir/SUMMARY.out file.
void Summary::PrintInfo(const string& dir) const
{
    string   FileOut = dir + "SUMMARY.out";
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    const USI ns  = 12;
    const USI col = 10;
    const USI num = Sumdata.size();
    const USI len = Sumdata[0].val.size();

    USI row = 0;
    USI id  = 0;
    USI ID  = 1;

    while (id != num) {

        outF << "Row " << ++row << "\n";

        // Item
        outF << "\t" << setw(ns) << Sumdata[0].Item;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Item;
            if (id == num) break;
        }
        outF << "\n";

        // Unit
        outF << "\t" << setw(ns) << Sumdata[0].Unit;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Unit;
            if (id == num) break;
        }
        outF << "\n";

        // Obj Name
        outF << "\t" << setw(ns) << Sumdata[0].Obj;

        id = ID;
        for (USI i = 1; i < col; i++) {
            outF << "\t" << setw(ns) << Sumdata[id++].Obj;
            if (id == num) break;
        }
        outF << "\n";

        // Data
        for (USI l = 0; l < len; l++) {

            // Time
            outF << "\t" << setw(ns) << fixed << setprecision(3) << Sumdata[0].val[l];

            id = ID;
            for (USI i = 1; i < col; i++) {
                if (Sumdata[id].Type == "int") {
                    outF << fixed << setprecision(0);
                } else if (Sumdata[id].Type == "fixed") {
                    outF << fixed << setprecision(3);
                } else if (Sumdata[id].Type == "float") {
                    outF << scientific << setprecision(5);
                }
                outF << "\t" << setw(ns) << Sumdata[id++].val[l];
                if (id == num) break;
            }
            outF << "\n";
        }

        ID += (col - 1);

        outF << "\n";
    }

    outF.close();
}

void CriticalInfo::Setup(const OCP_DBL& totalTime)
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

void CriticalInfo::SetVal(const Reservoir& rs, const OCPControl& ctrl)
{
    const Bulk& bulk = rs.bulk;

    time.push_back(ctrl.GetCurTime());
    dt.push_back(ctrl.GetLastDt());
    dPmax.push_back(bulk.GetdPmax());
    dVmax.push_back(bulk.GeteVmax());
    dSmax.push_back(bulk.GetdSmax());
    dNmax.push_back(bulk.GetdNmax());
    cfl.push_back(bulk.GetMaxCFL());
}

void CriticalInfo::PrintFastReview(const string& dir) const
{
    string    FileOut = dir + "FastReview.out";
    const USI ns      = 12;

    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    // Item
    outF << setw(ns) << "Time";
    outF << setw(ns) << "dt";
    outF << setw(ns) << "dPmax";
    outF << setw(ns) << "dVmax";
    outF << setw(ns) << "dSmax";
    outF << setw(ns) << "dNmax";
    outF << setw(ns) << "CFL" << endl;
    // Unit
    outF << setw(ns) << "Days";
    outF << setw(ns) << "Days";
    outF << setw(ns) << "Psia";
    outF << setw(ns) << "    ";
    outF << setw(ns) << "    ";
    outF << setw(ns) << "    ";
    outF << setw(ns) << "    " << endl;

    USI n = time.size();
    for (USI i = 0; i < n; i++) {
        outF << setw(ns) << fixed << setprecision(3) << time[i];
        outF << setw(ns) << fixed << setprecision(3) << dt[i];
        outF << setw(ns) << scientific << setprecision(3) << dPmax[i];
        outF << setw(ns) << scientific << setprecision(3) << dVmax[i];
        outF << setw(ns) << scientific << setprecision(3) << dSmax[i];
        outF << setw(ns) << scientific << setprecision(3) << dNmax[i];
        outF << setw(ns) << scientific << setprecision(3) << cfl[i];
        outF << "\n";
    }

    outF.close();
}

void BasicGridProperty::SetBasicGridProperty(const BasicGridPropertyParam& param)
{
    PRE  = param.PRE;
    PGAS = param.PGAS;
    PWAT = param.PWAT;
    SOIL = param.SOIL;
    SGAS = param.SGAS;
    SWAT = param.SWAT;
    DENO = param.DENO;
    DENG = param.DENG;
    DENW = param.DENW;
    KRO  = param.KRO;
    KRG  = param.KRG;
    KRW  = param.KRW;
    BOIL = param.BOIL;
    BGAS = param.BGAS;
    BWAT = param.BWAT;
    VOIL = param.VOIL;
    VGAS = param.VGAS;
    VWAT = param.VWAT;
    XMF  = param.XMF;
    YMF  = param.YMF;
    PCW  = param.PCW;
}

void Out4RPT::InputParam(const OutputRPTParam& RPTparam)
{
    useRPT = RPTparam.useRPT;
    if (!useRPT) return;

    bgp.SetBasicGridProperty(RPTparam.bgp);
}

void Out4RPT::Setup(const string& dir, const Reservoir& rs)
{
    if (!useRPT) return;

    string   FileOut = dir + "RPT.out";
    ofstream outF(FileOut);
    if (!outF.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }
    outF.close();

    const Grid& initGrid = rs.grid;

    nx       = initGrid.nx;
    ny       = initGrid.ny;
    numGrid  = initGrid.numGrid;
    IJKspace = initGrid.numDigutIJK;
}

void Out4RPT::PrintRPT(const string&    dir,
                       const Reservoir& rs,
                       const OCP_DBL&   days) const
{

    if (!useRPT) return;

    string   FileOut = dir + "RPT.out";
    ofstream outRPT;
    outRPT.open(FileOut, ios::app);
    if (!outRPT.is_open()) {
        OCP_ABORT("Can not open " + FileOut);
    }

    const Grid& initGrid = rs.grid;
    const Bulk& bulk     = rs.bulk;

    const USI              np     = bulk.numPhase;
    const USI              nc     = bulk.numCom;
    const USI              OIndex = bulk.phase2Index[OIL];
    const USI              GIndex = bulk.phase2Index[GAS];
    const USI              WIndex = bulk.phase2Index[WATER];
    const vector<GB_Pair>& g2bp   = initGrid.map_All2Act;

    outRPT << OCP_SEP02(50) << "\n";

    // Well Info
    USI numWell = rs.allWells.GetWellNum();
    outRPT << "Well Information"
           << "                    ";
    outRPT << fixed << setprecision(3) << days << "  DAYS"
           << "\n";
    // INJ
    for (USI w = 0; w < numWell; w++) {
        if (rs.allWells.wells[w].opt.type == INJ) {
            outRPT << "-------------------------------------"
                   << "\n";
            outRPT << rs.allWells.wells[w].name << "   " << w << "   "
                   << rs.allWells.wells[w].depth << " (feet)     ";
            outRPT << rs.allWells.wells[w].I << "   " << rs.allWells.wells[w].J << "\n";

            if (rs.allWells.wells[w].opt.state == OPEN) {
                outRPT << "OPEN\t" << rs.allWells.wells[w].WGIR << " (MSCF/DAY)\t"
                       << rs.allWells.wells[w].WWIR << " (STB/DAY)"
                       << "\n";
            } else {
                outRPT << "SHUTIN"
                       << "\n";
            }
            // perf
            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
                outRPT << "perf" << p << "   " << rs.allWells.wells[w].perf[p].I
                       << "   " << rs.allWells.wells[w].perf[p].J << "   "
                       << rs.allWells.wells[w].perf[p].K << "   "
                       << rs.allWells.wells[w].perf[p].depth << "   ";
                if (rs.allWells.wells[w].perf[p].state == OPEN) {
                    outRPT << "OPEN";
                } else {
                    outRPT << "SHUTIN";
                }
                outRPT << "   " << rs.allWells.wells[w].perf[p].location << "\n";
            }
        }
    }
    // PROD
    for (USI w = 0; w < numWell; w++) {
        if (rs.allWells.wells[w].opt.type == PROD) {
            outRPT << "-------------------------------------"
                   << "\n";
            outRPT << rs.allWells.wells[w].name << "   " << w << "   "
                   << rs.allWells.wells[w].depth << " (feet)     ";
            outRPT << rs.allWells.wells[w].I << "   " << rs.allWells.wells[w].J << "\n";

            if (rs.allWells.wells[w].opt.state == OPEN) {
                outRPT << "OPEN\t" << rs.allWells.wells[w].WOPR << " (STB/DAY)\t"
                       << rs.allWells.wells[w].WGPR << " (MSCF/DAY)\t"
                       << rs.allWells.wells[w].WWPR << " (STB/DAY)"
                       << "\n";
            } else {
                outRPT << "SHUTIN"
                       << "\n";
            }
            // perf
            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
                outRPT << "perf" << p << "   " << rs.allWells.wells[w].perf[p].I
                       << "   " << rs.allWells.wells[w].perf[p].J << "   "
                       << rs.allWells.wells[w].perf[p].K << "   "
                       << rs.allWells.wells[w].perf[p].depth << "   ";
                if (rs.allWells.wells[w].perf[p].state == OPEN) {
                    outRPT << "OPEN";
                } else {
                    outRPT << "SHUTIN";
                }
                outRPT << "   " << rs.allWells.wells[w].perf[p].location << "\n";
            }
        }
    }

    outRPT << "\n\n";

    static OCP_BOOL flag = OCP_FALSE;
    // Print once
    if (flag) {
        PrintRPT_Scalar(outRPT, "DX : feet", days, &initGrid.dx[0], 1, g2bp, OCP_FALSE);
        PrintRPT_Scalar(outRPT, "DY : feet", days, &initGrid.dy[0], 1, g2bp, OCP_FALSE);
        PrintRPT_Scalar(outRPT, "DZ : feet", days, &initGrid.dz[0], 1, g2bp, OCP_FALSE);
        PrintRPT_Scalar(outRPT, "Depth : feet", days, &initGrid.depth[0], 1, g2bp,
                        OCP_FALSE);
        PrintRPT_Scalar(outRPT, "PERMX : MDarcy", days, &initGrid.kx[0], 1, g2bp,
                        OCP_FALSE);
        PrintRPT_Scalar(outRPT, "PERMY : MDarcy", days, &initGrid.ky[0], 1, g2bp,
                        OCP_FALSE);
        PrintRPT_Scalar(outRPT, "PERMZ : MDarcy", days, &initGrid.kz[0], 1, g2bp,
                        OCP_FALSE);
        flag = OCP_FALSE;
    }

    // PRESSURE
    if (bgp.PRE) {
        PrintRPT_Scalar(outRPT, "PRESSURE : psia", days, &bulk.P[0], 1, g2bp, OCP_TRUE);
    }

    // DENSITY of OIL
    if (bgp.DENO && bulk.oil) {
        PrintRPT_Scalar(outRPT, "DENO : lb/ft3", days, &bulk.rho[OIndex], np, g2bp,
                        OCP_TRUE);
    }
    outRPT << endl;

    // DENSITY of GAS
    if (bgp.DENG && bulk.gas) {
        PrintRPT_Scalar(outRPT, "DENG : lb/ft3", days, &bulk.rho[GIndex], np, g2bp,
                        OCP_TRUE);
    }

    // DENSITY of WATER
    if (bgp.DENW && bulk.water) {
        PrintRPT_Scalar(outRPT, "DENW : lb/ft3", days, &bulk.rho[WIndex], np, g2bp,
                        OCP_TRUE);
    }

    // SATURATION of OIL
    if (bgp.SOIL && bulk.oil) {
        PrintRPT_Scalar(outRPT, "SOIL         ", days, &bulk.S[OIndex], np, g2bp,
                        OCP_TRUE);
    }

    // SATURATION of GAS
    if (bgp.SGAS && bulk.gas) {
        PrintRPT_Scalar(outRPT, "SGAS         ", days, &bulk.S[GIndex], np, g2bp,
                        OCP_TRUE);
    }

    // SATURATION of WATER
    if (bgp.SWAT && bulk.water) {
        PrintRPT_Scalar(outRPT, "SWAT         ", days, &bulk.S[WIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Relative Permeability of OIL
    if (bgp.KRO && bulk.oil) {
        PrintRPT_Scalar(outRPT, "KRO          ", days, &bulk.kr[OIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Relative Permeability of GAS
    if (bgp.KRG && bulk.gas) {
        PrintRPT_Scalar(outRPT, "KRG          ", days, &bulk.kr[GIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Relative Permeability of WATER
    if (bgp.KRW && bulk.water) {
        PrintRPT_Scalar(outRPT, "KRW          ", days, &bulk.kr[WIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Molar Density of OIL
    if (bgp.BOIL && bulk.oil && bulk.IfUseEoS()) {
        PrintRPT_Scalar(outRPT, "BOIL : lb-M/rb", days, &bulk.xi[OIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Molar Density of GAS
    if (bgp.BGAS && bulk.gas && bulk.IfUseEoS()) {
        PrintRPT_Scalar(outRPT, "BGAS : lb-M/rb", days, &bulk.xi[GIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Molar Density of WATER
    if (bgp.BWAT && bulk.water) {
        PrintRPT_Scalar(outRPT, "BWAT : lb-M/rb", days, &bulk.xi[WIndex], np, g2bp,
                        OCP_TRUE, (CONV1 * 19.437216));
    }

    // Viscosity of OIL
    if (bgp.VOIL && bulk.oil) {
        PrintRPT_Scalar(outRPT, "VOIL : lb-M/rb", days, &bulk.mu[OIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Viscosity of GAS
    if (bgp.VGAS && bulk.gas) {
        PrintRPT_Scalar(outRPT, "VGAS : lb-M/rb", days, &bulk.mu[GIndex], np, g2bp,
                        OCP_TRUE);
    }

    // Viscosity of WATER
    if (bgp.VWAT && bulk.water) {
        PrintRPT_Scalar(outRPT, "VWAT : lb-M/rb", days, &bulk.mu[WIndex], np, g2bp,
                        OCP_TRUE);
    }

    // liquid component mole fractions.
    if (bgp.XMF && bulk.IfUseEoS()) {
        for (USI i = 0; i < nc - 1; i++) {
            PrintRPT_Scalar(outRPT, "XMF : Oil  " + to_string(i + 1) + "th Component",
                            days, &bulk.xij[OIndex * nc + i], np * nc, g2bp, OCP_TRUE);
        }
    }

    // gas component mole fractions.
    if (bgp.YMF && bulk.IfUseEoS()) {
        for (USI i = 0; i < nc - 1; i++) {
            PrintRPT_Scalar(outRPT, "YMF : Gas  " + to_string(i + 1) + "th Component",
                            days, &bulk.xij[GIndex * nc + i], np * nc, g2bp, OCP_TRUE);
        }
    }

    // Po - Pw
    if (bgp.PCW) {
        PrintRPT_Scalar(outRPT, "PCW : psia  ", days, &bulk.Pc[WIndex], np, g2bp,
                        OCP_TRUE);

        // PrintRPT_Scalar(outRPT, "PPCW : psia ", days,
        // &rs.optFeatures.scalePcow.scaleVal[0],
        //                 1, g2bp, OCP_TRUE, 103.053 - 0.116);
    }

    outRPT.close();
}

void Out4RPT::GetIJKGrid(USI& i, USI& j, USI& k, const OCP_USI& n) const
{
    // i,j,k begin from 1
    // n must be the index of grids instead bulks
    k = n / (nx * ny) + 1;
    j = (n - (k - 1) * nx * ny) / nx + 1;
    i = n - (k - 1) * nx * ny - (j - 1) * nx + 1;
}

void Out4VTK::InputParam(const OutputVTKParam& VTKParam)
{
    useVTK = VTKParam.useVTK;
    if (!useVTK) return;

    bgp.SetBasicGridProperty(VTKParam.bgp);
}

void Out4VTK::Setup(const string& dir, const Reservoir& rs, const USI& ndates)
{
    if (!useVTK) return;

    string file = dir + "grid" + to_string(index) + ".vtk";
    string newfile;
    string title = "test";

    const Grid& initGrid = rs.grid;

    out4vtk.Init(file, title, VTK_ASCII, VTK_UNSTRUCTURED_GRID,
                 initGrid.polyhedronGrid.size(), rs.allWells.polyhedronWell.size());
    out4vtk.OutputPOINTS(file, initGrid.polyhedronGrid, rs.allWells.polyhedronWell,
                         VTK_FLOAT);
    out4vtk.OutputCELLS(file, initGrid.polyhedronGrid, rs.allWells.polyhedronWell);
    out4vtk.OutputCELL_TYPES(file, initGrid.polyhedronGrid, rs.allWells.polyhedronWell);
    out4vtk.BeginCellData();
    // output dead grid, live grid, well
    vector<USI> tmpW(rs.allWells.numWell, 10);
    out4vtk.OutputCELL_DATA_SCALARS(file, "CellType", VTK_UNSIGNED_INT,
                                    &initGrid.gridTag[0], 1, initGrid.map_All2Act,
                                    OCP_FALSE, &tmpW[0]);

    for (USI i = 1; i < ndates; i++) {
        index++;
        newfile = dir + "grid" + to_string(index) + ".vtk";
        ;
        ifstream source(file, ios::binary);
        ofstream dest(newfile, ios::binary);
        dest << source.rdbuf();
        source.close();
        dest.close();
    }
    index = 0;

#ifdef USE_METIS
    metisTest.Setup(rs);
#endif // USE_METIS
}

void Out4VTK::PrintVTK(const string&    dir,
                       const Reservoir& rs,
                       const OCP_DBL&   days) const
{
    if (!useVTK) return;

    string file = dir + "grid" + to_string(index) + ".vtk";
    // Calulcate Well val for output
    rs.allWells.SetWellVal();

    const Grid&            initGrid = rs.grid;
    const Bulk&            bulk     = rs.bulk;
    const vector<GB_Pair>& g2bp     = initGrid.map_All2Act;
    const vector<OCP_DBL>& well     = rs.allWells.wellVal;
    const USI              np       = bulk.numPhase;
    const USI              OIndex   = bulk.phase2Index[OIL];
    const USI              GIndex   = bulk.phase2Index[GAS];
    const USI              WIndex   = bulk.phase2Index[WATER];

    // output
    if (bgp.PRE)
        out4vtk.OutputCELL_DATA_SCALARS(file, "PRESSURE", VTK_FLOAT, &bulk.P[0], 1,
                                        g2bp, OCP_TRUE, &well[0]);
    if (bgp.SOIL)
        out4vtk.OutputCELL_DATA_SCALARS(file, "SOIL", VTK_FLOAT, &bulk.S[OIndex], np,
                                        g2bp, OCP_TRUE, &well[0]);
    if (bgp.SGAS)
        out4vtk.OutputCELL_DATA_SCALARS(file, "SGAS", VTK_FLOAT, &bulk.S[GIndex], np,
                                        g2bp, OCP_TRUE, &well[0]);
    if (bgp.SWAT)
        out4vtk.OutputCELL_DATA_SCALARS(file, "SWAT", VTK_FLOAT, &bulk.S[WIndex], np,
                                        g2bp, OCP_TRUE, &well[0]);

#ifdef USE_METIS
    if (metisTest.useMetis) {
        // partition and print
        // vertex weights set to 1 now
        metisTest.vwgt.resize(metisTest.nvtxs, 1);
        // metisTest.MyPartitionFunc(METIS_PartGraphRecursive);
        metisTest.MyPartitionFunc(METIS_PartGraphKway);
        metisTest.SetPartitions(initGrid.map_Act2All);
        out4vtk.OutputCELL_DATA_SCALARS(file, "PARTIONS", VTK_UNSIGNED_INT,
                                        &metisTest.partitions[0], 1, g2bp, OCP_FALSE,
                                        &metisTest.partitions[metisTest.ng]);
    }
#endif // USE_METIS

    index++;
}

void OCPOutput::InputParam(const ParamOutput& paramOutput)
{
    summary.InputParam(paramOutput.summary);
    out4RPT.InputParam(paramOutput.outRPTParam);
    out4VTK.InputParam(paramOutput.outVTKParam);
}

void OCPOutput::Setup(const Reservoir& reservoir, const OCPControl& ctrl)
{
    workDir = ctrl.workDir;
    summary.Setup(reservoir, ctrl.criticalTime.back());
    crtInfo.Setup(ctrl.criticalTime.back());
    out4RPT.Setup(workDir, reservoir);
    out4VTK.Setup(workDir, reservoir, ctrl.criticalTime.size());
}

void OCPOutput::SetVal(const Reservoir& reservoir, const OCPControl& ctrl)
{
    summary.SetVal(reservoir, ctrl);
    crtInfo.SetVal(reservoir, ctrl);
}

void OCPOutput::PrintInfo() const
{
    summary.PrintInfo(workDir);
    crtInfo.PrintFastReview(workDir);
}

void OCPOutput::PrintInfoSched(const Reservoir&  rs,
                               const OCPControl& ctrl,
                               const OCP_DBL&    time) const
{
    OCP_DBL days = ctrl.current_time;

    // print timing info on the screen
    if (ctrl.printLevel >= PRINT_MIN) {
        cout << "Timestep " << setw(6) << left << ctrl.numTstep << ": " << fixed
             << setw(10) << setprecision(3) << right << days << " Days"
             << "    Wall time: " << time / 1000 << " Sec" << endl;
    }

    // print to output file
    // TODO: Add a control flag to enable or disable --zcs
    GetWallTime timer;
    timer.Start();
    out4RPT.PrintRPT(workDir, rs, days);
    out4VTK.PrintVTK(workDir, rs, days);
    outputTime += timer.Stop() / 1000;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Jan/09/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/