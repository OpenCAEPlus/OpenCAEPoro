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

void Summary::InputParam(const OutputSummary &summary_param)
{
    FPR = summary_param.FPR;
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

    // cout << "Summary::InputParam" << endl;
}

void Summary::Setup(const Reservoir &reservoir, const OCP_DBL &totalTime)
{
    Sumdata.push_back(SumPair("TIME", "  ", "DAY"));
    Sumdata.push_back(SumPair("NRiter", "  ", "  "));
    Sumdata.push_back(SumPair("LSiter", "  ", "  "));
    if (FPR)
        Sumdata.push_back(SumPair("FPR", "  ", "PSIA"));
    if (FOPR)
        Sumdata.push_back(SumPair("FOPR", "  ", "STB/DAY"));
    if (FOPT)
        Sumdata.push_back(SumPair("FOPT", "  ", "STB"));
    if (FGPR)
        Sumdata.push_back(SumPair("FGPR", "  ", "MSCF/DAY"));
    if (FGPt)
        Sumdata.push_back(SumPair("FGPT", "  ", "MSCF"));
    if (FWPR)
        Sumdata.push_back(SumPair("FWPR", "  ", "STB/DAY"));
    if (FWPT)
        Sumdata.push_back(SumPair("FWPT", "  ", "STB"));
    if (FGIR)
        Sumdata.push_back(SumPair("FGIR", "  ", "MSCF/DAY"));
    if (FGIT)
        Sumdata.push_back(SumPair("FGIT", "  ", "MSCF"));
    if (FWIR)
        Sumdata.push_back(SumPair("FWIR", "  ", "STB/DAY"));
    if (FWIT)
        Sumdata.push_back(SumPair("FWIT", "  ", "STB"));

    const USI sp = reservoir.grid.GetNumDigitIJK();

    USI num;
    USI wellnum = reservoir.allWells.GetWellNum();
    string wellname;

    if (WOPR.activity)
    {
        if (WOPR.obj[0] == "All")
        {
            WOPR.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WOPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WOPR", wellname, "STB/DAY"));
                WOPR.index.push_back(w);
            }
        }
        else
        {
            num = WOPR.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WOPR", WOPR.obj[w], "STB/DAY"));
                WOPR.index.push_back(reservoir.allWells.GetIndex(WOPR.obj[w]));
            }
        }
    }

    if (WOPT.activity)
    {
        if (WOPT.obj[0] == "All")
        {
            WOPT.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WOPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WOPT", wellname, "STB"));
                WOPT.index.push_back(w);
            }
        }
        else
        {
            num = WOPT.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WOPT", WOPT.obj[w], "STB"));
                WOPT.index.push_back(reservoir.allWells.GetIndex(WOPT.obj[w]));
            }
        }
    }

    if (WGPR.activity)
    {
        if (WGPR.obj[0] == "All")
        {
            WGPR.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WGPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGPR", wellname, "MSCF/DAY"));
                WGPR.index.push_back(w);
            }
        }
        else
        {
            num = WGPR.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WGPR", WGPR.obj[w], "MSCF/DAY"));
                WGPR.index.push_back(reservoir.allWells.GetIndex(WGPR.obj[w]));
            }
        }
    }

    if (WGPT.activity)
    {
        if (WGPT.obj[0] == "All")
        {
            WGPT.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WGPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGPT", wellname, "MSCF"));
                WGPT.index.push_back(w);
            }
        }
        else
        {
            num = WGPT.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WGPT", WGPT.obj[w], "MSCF"));
                WGPT.index.push_back(reservoir.allWells.GetIndex(WGPT.obj[w]));
            }
        }
    }

    if (WWPR.activity)
    {
        if (WWPR.obj[0] == "All")
        {
            WWPR.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WWPR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWPR", wellname, "STB/DAY"));
                WWPR.index.push_back(w);
            }
        }
        else
        {
            num = WWPR.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WWPR", WWPR.obj[w], "STB/DAY"));
                WWPR.index.push_back(reservoir.allWells.GetIndex(WWPR.obj[w]));
            }
        }
    }

    if (WWPT.activity)
    {
        if (WWPT.obj[0] == "All")
        {
            WWPT.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WWPT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWPT", wellname, "STB"));
                WWPT.index.push_back(w);
            }
        }
        else
        {
            num = WWPT.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WWPT", WWPT.obj[w], "STB"));
                WWPT.index.push_back(reservoir.allWells.GetIndex(WWPT.obj[w]));
            }
        }
    }

    if (WGIR.activity)
    {
        if (WGIR.obj[0] == "All")
        {
            WGIR.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WGIR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGIR", wellname, "MSCF/DAY"));
                WGIR.index.push_back(w);
            }
        }
        else
        {
            num = WGIR.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WGIR", WGIR.obj[w], "MSCF/DAY"));
                WGIR.index.push_back(reservoir.allWells.GetIndex(WGIR.obj[w]));
            }
        }
    }

    if (WGIT.activity)
    {
        if (WGIT.obj[0] == "All")
        {
            WGIT.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WGIT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WGIT", wellname, "MSCF"));
                WGIT.index.push_back(w);
            }
        }
        else
        {
            num = WGIT.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WGIT", WGIT.obj[w], "MSCF"));
                WGIT.index.push_back(reservoir.allWells.GetIndex(WGIT.obj[w]));
            }
        }
    }

    if (WWIR.activity)
    {
        if (WWIR.obj[0] == "All")
        {
            WWIR.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WWIR.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWIR", wellname, "STB/DAY"));
                WWIR.index.push_back(w);
            }
        }
        else
        {
            num = WWIR.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WWIR", WWIR.obj[w], "STB/DAY"));
                WWIR.index.push_back(reservoir.allWells.GetIndex(WWIR.obj[w]));
            }
        }
    }

    if (WWIT.activity)
    {
        if (WWIT.obj[0] == "All")
        {
            WWIT.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WWIT.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WWIT", wellname, "STB"));
                WWIT.index.push_back(w);
            }
        }
        else
        {
            num = WWIT.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WWIT", WWIT.obj[w], "STB"));
                WWIT.index.push_back(reservoir.allWells.GetIndex(WWIT.obj[w]));
            }
        }
    }

    if (WBHP.activity)
    {
        if (WBHP.obj[0] == "All")
        {
            WBHP.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                wellname = reservoir.allWells.GetWellName(w);
                WBHP.obj.push_back(wellname);
                Sumdata.push_back(SumPair("WBHP", wellname, "PSIA"));
                WBHP.index.push_back(w);
            }
        }
        else
        {
            num = WBHP.obj.size();
            for (USI w = 0; w < num; w++)
            {
                Sumdata.push_back(SumPair("WBHP", WBHP.obj[w], "PSIA"));
                WBHP.index.push_back(reservoir.allWells.GetIndex(WBHP.obj[w]));
            }
        }
    }

    if (DG.activity)
    {
        if (DG.obj[0] == "All")
        {
            DG.obj.clear();
            for (USI w = 0; w < wellnum; w++)
            {
                string wellname = reservoir.allWells.GetWellName(w);
                DG.obj.push_back(wellname);
                USI perfnum = reservoir.allWells.GetWellPerfNum(w);
                for (USI p = 0; p < perfnum; p++)
                {
                    Sumdata.push_back(
                        SumPair("DG", wellname + " Perf" + to_string(p), "PSIA"));
                    DG.index.push_back(w);
                }
            }
        }
        else
        {
            num = DG.obj.size();
            for (USI w = 0; w < num; w++)
            {
                USI wId = reservoir.allWells.GetIndex(DG.obj[w]);
                USI perfnum = reservoir.allWells.GetWellPerfNum(wId);
                for (USI p = 0; p < perfnum; p++)
                {
                    Sumdata.push_back(
                        SumPair("DG", DG.obj[w] + " P" + to_string(p), "PSIA"));
                    DG.index.push_back(wId);
                }
            }
        }
    }

    if (BPR.activity)
    {
        num = BPR.obj.size();
        for (USI i = 0; i < num; i++)
        {
            string temp = GetIJKformat(to_string(BPR.obj[i].I), to_string(BPR.obj[i].J), 
                to_string(BPR.obj[i].K), sp);
            Sumdata.push_back(SumPair("BPR", temp, "PSIA"));
            USI I = BPR.obj[i].I - 1;
            USI J = BPR.obj[i].J - 1;
            USI K = BPR.obj[i].K - 1;
            BPR.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SOIL.activity)
    {
        num = SOIL.obj.size();
        for (USI i = 0; i < num; i++)
        {
            string temp = GetIJKformat(to_string(SOIL.obj[i].I), to_string(SOIL.obj[i].J),
                to_string(SOIL.obj[i].K), sp);
            Sumdata.push_back(SumPair("SOIL", temp, "   "));
            USI I = SOIL.obj[i].I - 1;
            USI J = SOIL.obj[i].J - 1;
            USI K = SOIL.obj[i].K - 1;
            SOIL.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SGAS.activity)
    {
        num = SGAS.obj.size();
        for (USI i = 0; i < num; i++)
        {
            string temp = GetIJKformat(to_string(SGAS.obj[i].I), to_string(SGAS.obj[i].J),
                to_string(SGAS.obj[i].K), sp);
            Sumdata.push_back(SumPair("SGAS", temp, "   "));
            USI I = SGAS.obj[i].I - 1;
            USI J = SGAS.obj[i].J - 1;
            USI K = SGAS.obj[i].K - 1;
            SGAS.index.push_back(reservoir.grid.GetActIndex(I, J, K));
        }
    }

    if (SWAT.activity)
    {
        num = SWAT.obj.size();
        for (USI i = 0; i < num; i++)
        {
            string temp = GetIJKformat(to_string(SWAT.obj[i].I), to_string(SWAT.obj[i].J),
                to_string(SWAT.obj[i].K), sp);
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
    for (USI i = 0; i < cs; i++)
    {
        Sumdata[i].val.reserve(rs);
    }

    // cout << "Summary::Setup" << endl;
}

void Summary::SetVal(const Reservoir &rs, const OCPControl &ctrl)
{
    USI n = 0;

    // TIME
    Sumdata[n++].val.push_back(ctrl.GetCurTime());
    // NRiter
    Sumdata[n++].val.push_back(ctrl.GetNRiterT());
    // LSiter
    Sumdata[n++].val.push_back(ctrl.GetLSiterT());

    // FPR
    if (FPR)
        Sumdata[n++].val.push_back(rs.bulk.CalFPR());
    if (FOPR)
        Sumdata[n++].val.push_back(rs.allWells.GetFOPR());
    if (FOPT)
        Sumdata[n++].val.push_back(rs.allWells.GetFOPT());
    if (FGPR)
        Sumdata[n++].val.push_back(rs.allWells.GetFGPR());
    if (FGPt)
        Sumdata[n++].val.push_back(rs.allWells.GetFGPT());
    if (FWPR)
        Sumdata[n++].val.push_back(rs.allWells.GetFWPR());
    if (FWPT)
        Sumdata[n++].val.push_back(rs.allWells.GetFWPT());
    if (FGIR)
        Sumdata[n++].val.push_back(rs.allWells.GetFGIR());
    if (FGIT)
        Sumdata[n++].val.push_back(rs.allWells.GetFGIT());
    if (FWIR)
        Sumdata[n++].val.push_back(rs.allWells.GetFWIR());
    if (FWIT)
        Sumdata[n++].val.push_back(rs.allWells.GetFWIT());

    USI len = 0;
    // WOPR
    len = WOPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWOPR(WOPR.index[w]));

    // WOPT
    len = WOPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWOPT(WOPT.index[w]));

    // WGPR
    len = WGPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWGPR(WGPR.index[w]));

    // WGPT
    len = WGPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWGPT(WGPT.index[w]));

    // WWPR
    len = WWPR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWWPR(WWPR.index[w]));

    // WWPT
    len = WWPT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWWPT(WWPT.index[w]));

    // WGIR
    len = WGIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWGIR(WGIR.index[w]));

    // WGIT
    len = WGIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWGIT(WGIT.index[w]));

    // WWIR
    len = WWIR.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWWIR(WWIR.index[w]));

    // WWIT
    len = WWIT.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWWIT(WWIT.index[w]));

    // WBHP
    len = WBHP.index.size();
    for (USI w = 0; w < len; w++)
        Sumdata[n++].val.push_back(rs.allWells.GetWBHP(WBHP.index[w]));

    // DG
    len = DG.obj.size();
    for (USI w = 0; w < len; w++)
    {
        USI numperf = rs.allWells.GetWellPerfNum(DG.index[w]);
        for (USI p = 0; p < numperf; p++)
        {
            Sumdata[n++].val.push_back(rs.allWells.GetWellDg(DG.index[w], p));
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

/// Write output information in the dir/SUMMARY.out file.
void Summary::PrintInfo(const string &dir) const
{
    string FileOut = dir + "SUMMARY.out";
    ofstream outF(FileOut);
    if (!outF.is_open())
    {
        OCP_ABORT("Can not open " + FileOut);
    }

    USI ns = 12;
    USI col = 10;
    USI row = 0;
    USI num = Sumdata.size();
    USI len = Sumdata[0].val.size();
    USI id = 0;
    USI ID = 1;

    while (id != num)
    {

        outF << "Row " << ++row << endl;

        // Item
        // Time
        outF << "\t" << setw(ns) << Sumdata[0].Item;

        id = ID;
        for (USI i = 1; i < col; i++)
        {
            outF << "\t" << setw(ns) << Sumdata[id++].Item;
            if (id == num)
                break;
        }
        outF << endl;

        // Unit
        // Time
        outF << "\t" << setw(ns) << Sumdata[0].Unit;

        id = ID;
        for (USI i = 1; i < col; i++)
        {
            outF << "\t" << setw(ns) << Sumdata[id++].Unit;
            if (id == num)
                break;
        }
        outF << endl;

        // Obj
        // Time
        outF << "\t" << setw(ns) << Sumdata[0].Obj;

        id = ID;
        for (USI i = 1; i < col; i++)
        {
            outF << "\t" << setw(ns) << Sumdata[id++].Obj;
            if (id == num)
                break;
        }
        outF << endl;

        // Data
        for (USI l = 0; l < len; l++)
        {

            // Time
            outF << "\t" << setw(ns) << Sumdata[0].val[l];

            id = ID;
            for (USI i = 1; i < col; i++)
            {
                outF << "\t" << setw(ns) << Sumdata[id++].val[l];
                if (id == num)
                    break;
            }
            outF << endl;
        }

        ID += (col - 1);

        outF << endl;
    }

    outF.close();
}

void CriticalInfo::Setup(const OCP_DBL &totalTime)
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

void CriticalInfo::SetVal(const Reservoir &reservoir, const OCPControl &ctrl)
{
    time.push_back(ctrl.GetCurTime());
    dt.push_back(ctrl.GetLastDt());
    dPmax.push_back(reservoir.bulk.GetdPmax());
    dVmax.push_back(reservoir.bulk.GetdVmax());
    dSmax.push_back(reservoir.bulk.GetdSmax());
    dNmax.push_back(reservoir.bulk.GetdNmax());
    cfl.push_back(reservoir.cfl);
}

void CriticalInfo::PrintInfo(const string &dir) const
{
    string FileOut = dir + "FastReview.out";
    ofstream outF(FileOut);
    if (!outF.is_open())
    {
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
    for (USI i = 0; i < n; i++)
    {
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

void DetailInfo::InputParam(const OutputDetail &detail_param)
{
    PRE = detail_param.PRE;
    PGAS = detail_param.PGAS;
    PWAT = detail_param.PWAT;
    SOIL = detail_param.SOIL;
    SGAS = detail_param.SGAS;
    SWAT = detail_param.SWAT;
    DENO = detail_param.DENO;
    DENG = detail_param.DENG;
    DENW = detail_param.DENW;
    KRO = detail_param.KRO;
    KRG = detail_param.KRG;
    KRW = detail_param.KRW;
    BOIL = detail_param.BOIL;
    BGAS = detail_param.BGAS;
    BWAT = detail_param.BWAT;
    VOIL = detail_param.VOIL;
    VGAS = detail_param.VGAS;
    VWAT = detail_param.VWAT;
    XMF = detail_param.XMF;
    YMF = detail_param.YMF;
    PCW = detail_param.PCW;
}

void DetailInfo::Setup(const string &dir)
{
    string FileOut = dir + "RPT.out";
    ofstream outF(FileOut);
    if (!outF.is_open())
    {
        OCP_ABORT("Can not open " + FileOut);
    }
    outF.close();
}

void DetailInfo::PrintInfo(const string &dir, const Reservoir &rs,
                           const OCP_DBL &days) const
{
    string FileOut = dir + "RPT.out";
    ofstream outF;
    outF.open(FileOut, ios::app);
    if (!outF.is_open())
    {
        OCP_ABORT("Can not open " + FileOut);
    }

    const USI     np = rs.bulk.numPhase;
    const USI     nc = rs.bulk.numCom;
    const USI OIndex = rs.bulk.phase2Index[OIL];
    const USI GIndex = rs.bulk.phase2Index[GAS];
    const USI WIndex = rs.bulk.phase2Index[WATER];
    const USI     nx  = rs.grid.GetGridNx();
    const USI     ny  = rs.grid.GetGridNy();
    const OCP_USI num = rs.grid.GetGridNum();
    const USI tmpsp = rs.grid.GetNumDigitIJK();
    OCP_USI bId;
    OCP_USI tmpId;
    USI I, J, K;

    const string sep01(50, '=');
    const string sep02(50, '-');

    outF << sep01 << "\n";

    static OCP_BOOL flag = OCP_FALSE;
    // Print once
    if (flag)
    {
        outF << "DX : feet";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.dx[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "DY : feet";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.dy[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "DZ : feet";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.dz[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "Depth : feet";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(1) << rs.bulk.depth[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "PERMX : MDarcy";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.rockKxInit[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "PERMY : MDarcy";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.rockKyInit[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    if (flag)
    {
        outF << "PERMZ : MDarcy";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(5) << rs.bulk.rockKzInit[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    flag = OCP_FALSE;


    // Well infor
    USI numWell = rs.allWells.GetWellNum();
    outF << "Well Information" << "                    ";
    outF << fixed << setprecision(3) << days << "  DAYS" << endl;
    // INJ
    for (USI w = 0; w < numWell; w++) {
        if (rs.allWells.wells[w].opt.type == INJ) {
            outF << "-------------------------------------" << endl;
            outF << rs.allWells.wells[w].name << "   " << w << "   "
                << rs.allWells.wells[w].depth << " (feet)     ";
            outF << rs.allWells.wells[w].I << "   "
                << rs.allWells.wells[w].J << endl;

            if (rs.allWells.wells[w].opt.state == OPEN) {
                outF << "OPEN\t" << rs.allWells.wells[w].WGIR << " (MSCF/DAY)\t" << rs.allWells.wells[w].WWIR << " (STB/DAY)" << endl;
            }
            else {
                outF << "SHUTIN" << endl;
            }
            // perf
            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
                outF << "perf" << p << "   "
                    << rs.allWells.wells[w].perf[p].I << "   "
                    << rs.allWells.wells[w].perf[p].J << "   "
                    << rs.allWells.wells[w].perf[p].K << "   "
                    << rs.allWells.wells[w].perf[p].depth << "   ";
                if (rs.allWells.wells[w].perf[p].state == OPEN) {
                    outF << "OPEN";
                }
                else {
                    outF << "SHUTIN";
                }
                outF << "   " << rs.allWells.wells[w].perf[p].location << endl;
            }
        }
    }
    // PROD
    for (USI w = 0; w < numWell; w++) {
        if (rs.allWells.wells[w].opt.type == PROD) {
            outF << "-------------------------------------" << endl;
            outF << rs.allWells.wells[w].name << "   " << w << "   "
                << rs.allWells.wells[w].depth << " (feet)     ";
            outF << rs.allWells.wells[w].I << "   "
                << rs.allWells.wells[w].J << endl;

            if (rs.allWells.wells[w].opt.state == OPEN) {
                outF << "OPEN\t" << rs.allWells.wells[w].WOPR << " (STB/DAY)\t" << rs.allWells.wells[w].WGPR << " (MSCF/DAY)\t" << rs.allWells.wells[w].WWPR << " (STB/DAY)" << endl;
            }
            else {
                outF << "SHUTIN" << endl;
            }
            // perf
            for (USI p = 0; p < rs.allWells.wells[w].numPerf; p++) {
                outF << "perf" << p << "   "
                    << rs.allWells.wells[w].perf[p].I << "   "
                    << rs.allWells.wells[w].perf[p].J << "   "
                    << rs.allWells.wells[w].perf[p].K << "   "
                    << rs.allWells.wells[w].perf[p].depth << "   ";
                if (rs.allWells.wells[w].perf[p].state == OPEN) {
                    outF << "OPEN";
                }
                else {
                    outF << "SHUTIN";
                }
                outF << "   " << rs.allWells.wells[w].perf[p].location << endl;
            }
        }
    }

    outF << endl << endl;



    // PRESSURE
    if (PRE)
    {
        outF << "PRESSURE : psia"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(3) << rs.bulk.P[bId];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }

    // DENSITY of OIL
    if (DENO && rs.bulk.oil)
    {
        outF << sep02 << "\n";
        outF << "DENO : lb/ft3"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + OIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(3) << rs.bulk.rho[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(2) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // DENSITY of GAS
    if (DENG && rs.bulk.gas)
    {
        outF << sep02 << "\n";
        outF << "DENG : lb/ft3"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + GIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(3) << rs.bulk.rho[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(2) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // DENSITY of WATER
    if (DENW && rs.bulk.water)
    {
        outF << sep02 << "\n";
        outF << "DENW : lb/ft3"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + WIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(3) << rs.bulk.rho[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(2) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // SATURATION of OIL
    if (SOIL && rs.bulk.oil)
    {
        outF << sep02 << "\n";
        outF << "SOIL"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + OIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.S[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // SATURATION of GAS
    if (SGAS && rs.bulk.gas)
    {
        outF << sep02 << "\n";
        outF << "SGAS"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + GIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.S[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // SATURATION of WATER
    if (SWAT && rs.bulk.water)
    {
        outF << sep02 << "\n";
        outF << "SWAT"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + WIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.S[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Relative Permeability of OIL
    if (KRO && rs.bulk.oil)
    {
        outF << sep02 << "\n";
        outF << "KRO"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + OIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.kr[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Relative Permeability of GAS
    if (KRG && rs.bulk.gas)
    {
        outF << sep02 << "\n";
        outF << "KRG"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + GIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.kr[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Relative Permeability of WATER
    if (KRW && rs.bulk.water)
    {
        outF << sep02 << "\n";
        outF << "KRW"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + WIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.kr[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Molar Density of OIL
    if (BOIL && rs.bulk.oil && rs.bulk.comps)
    {
        outF << sep02 << "\n";
        outF << "BOIL : lb-M/rb"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + OIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.xi[tmpId] * CONV1;
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Molar Density of GAS
    if (BGAS && rs.bulk.gas && rs.bulk.comps)
    {
        outF << sep02 << "\n";
        outF << "BGAS : lb-M/rb"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + GIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.xi[tmpId] * CONV1;
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Molar Density of WATER
    if (BWAT && rs.bulk.water)
    {
        outF << sep02 << "\n";
        outF << "BWAT : lb-M/rb"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + WIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.xi[tmpId] * (CONV1 * 19.437216);
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Viscosity of OIL
    if (VOIL && rs.bulk.oil)
    {
        outF << sep02 << "\n";
        outF << "VOIL : cp"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + OIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.mu[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Viscosity of GAS
    if (VGAS && rs.bulk.gas)
    {
        outF << sep02 << "\n";
        outF << "VGAS : cp"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + GIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.mu[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Viscosity of WATER
    if (VWAT && rs.bulk.water)
    {
        outF << sep02 << "\n";
        outF << "VWAT : cp"
             << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                tmpId = bId * np + WIndex;
                if (rs.bulk.phaseExist[tmpId])
                {
                    outF << setw(10) << fixed << setprecision(5) << rs.bulk.mu[tmpId];
                }
                else
                {
                    outF << setw(9) << fixed << setprecision(4) << 0.0 << "N";
                }
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // liquid component mole fractions.
    if (XMF && rs.bulk.comps) {     
        for (USI i = 0; i < nc - 1; i++) {
            // the ith component
            outF << sep02 << "\n";
            outF << "XMF : Oil  " << i+1 << "th Component"
                << "                   ";
            outF << fixed << setprecision(3) << days << "  DAYS";

            for (OCP_USI n = 0; n < num; n++) {
                if (n % nx == 0) outF << "\n";
                if (n % (nx * ny) == 0) outF << "\n";

                if (n % nx == 0) {
                    rs.grid.GetIJKGrid(I, J, K, n);
                    outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                    // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
                }

                if (rs.grid.MapG2B(n).IsAct()) {
                    bId = rs.grid.MapG2B(n).GetId();
                    tmpId = bId * np + OIndex;
                    if (rs.bulk.phaseExist[tmpId]) {
                        tmpId = tmpId * nc + i;
                        outF << setw(10) << fixed << setprecision(6) << rs.bulk.xij[tmpId];
                    }
                    else {
                        outF << setw(9) << fixed << setprecision(5) << 0.0 << "N";
                    }
                }
                else {
                    outF << setw(10) << " --- ";
                }
            }
            outF << "\n\n";
        }
    }

    // gas component mole fractions.
    if (YMF && rs.bulk.comps) {
        for (USI i = 0; i < nc - 1; i++) {
            // the ith component
            outF << sep02 << "\n";
            outF << "YMF : Gas  " << i << "th Component"
                << "                   ";
            outF << fixed << setprecision(3) << days << "  DAYS";

            for (OCP_USI n = 0; n < num; n++) {
                if (n % nx == 0) outF << "\n";
                if (n % (nx * ny) == 0) outF << "\n";

                if (n % nx == 0) {
                    rs.grid.GetIJKGrid(I, J, K, n);
                    outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                    // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
                }

                if (rs.grid.MapG2B(n).IsAct()) {
                    bId = rs.grid.MapG2B(n).GetId();
                    tmpId = bId * np + GIndex;
                    if (rs.bulk.phaseExist[tmpId]) {
                        tmpId = tmpId * nc + i;
                        outF << setw(10) << fixed << setprecision(6) << rs.bulk.xij[tmpId];
                    }
                    else {
                        outF << setw(9) << fixed << setprecision(5) << 0.0 << "N";
                    }
                }
                else {
                    outF << setw(10) << " --- ";
                }
            }
            outF << "\n\n";
        }
    }

    // surface tension
    if (rs.bulk.miscible && OCP_FALSE)
    {
        outF << sep02 << "\n";
        outF << "STEN"
            << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(10) << fixed << setprecision(5) << rs.bulk.surTen[bId];
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Fk
    if (rs.bulk.miscible && OCP_FALSE)
    {
        outF << sep02 << "\n";
        outF << "FMISC"
            << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(10) << fixed << setprecision(5) << rs.bulk.Fk[bId];
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Fp
    if (rs.bulk.miscible && OCP_FALSE)
    {
        outF << sep02 << "\n";
        outF << "FPC"
            << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(10) << fixed << setprecision(5) << rs.bulk.Fp[bId];
            }
            else
            {
                outF << setw(10) << " --- ";
            }
        }
        outF << "\n\n";
    }

    // Po - Pw
    if (PCW)
    {
        outF << "PCW : psia"
            << "                   ";
        outF << fixed << setprecision(3) << days << "  DAYS";
        for (OCP_USI i = 0; i < num; i++)
        {
            if (i % nx == 0)
                outF << "\n";
            if (i % (nx * ny) == 0)
                outF << "\n\n";

            if (i % nx == 0)
            {
                rs.grid.GetIJKGrid(I, J, K, i);
                outF << GetIJKformat("*", to_string(J), to_string(K), tmpsp);
                // outF << "(*," << setw(3) << J << "," << setw(3) << K << ")";
            }

            if (rs.grid.MapG2B(i).IsAct())
            {                
                bId = rs.grid.MapG2B(i).GetId();
                outF << setw(12) << fixed << setprecision(3) << -rs.bulk.Pc[bId * np + WIndex];
            }
            else
            {
                outF << setw(12) << "-----  ";
            }
        }
        outF << "\n\n\n";
    }


    outF.close();
}

void OCPOutput::InputParam(const ParamOutput &paramOutput)
{
    summary.InputParam(paramOutput.summary);
    dtlInfo.InputParam(paramOutput.detailInfo);
}

void OCPOutput::Setup(const Reservoir &reservoir, const OCPControl &ctrl)
{
    wordDir = ctrl.workDir;
    summary.Setup(reservoir, ctrl.criticalTime.back());
    crtInfo.Setup(ctrl.criticalTime.back());
    dtlInfo.Setup(wordDir);
}

void OCPOutput::SetVal(const Reservoir &reservoir, const OCPControl &ctrl)
{
    summary.SetVal(reservoir, ctrl);
    crtInfo.SetVal(reservoir, ctrl);
}

void OCPOutput::PrintInfo() const
{
    summary.PrintInfo(wordDir);
    crtInfo.PrintInfo(wordDir);
}

void OCPOutput::PrintInfoSched(const Reservoir &rs, const OCPControl &ctrl,
                               const OCP_DBL &time) const
{
    OCP_DBL days = ctrl.current_time;
    cout << "Timestep " << setw(6) << left << ctrl.numTstep
         << ": " << fixed << setw(10) << setprecision(3) << right << days << " Days"
         << "    Wall time: " << time / 1000 << " Sec" << endl;
    dtlInfo.PrintInfo(wordDir, rs, days);
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