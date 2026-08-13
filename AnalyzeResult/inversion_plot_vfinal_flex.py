#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cleaned inversion plotting script for Python 3.7.6
Compatible with:
    numpy  1.18.1
    pandas 1.0.1
Usage
-----
python inversion_plot_vfinal_flex.py STA MCMC_DIR

Example
-------
python inversion_plot_vfinal_flex.py VR03 /path/to/MonteCarlo
Note:
All the input file included in "build_paths" function
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter


# =============================================================================
# USER STYLE / PLOT PARAMETERS
# =============================================================================
FIGSIZE_MAIN = (18, 20)
FIGSIZE_MISFIT = (5, 5)
DPI = 300

VS_XLIM = (0.0, 7.0)
VS_YLIM_FULL = (0.0, 10.0)
VS_YLIM_SHALLOW = (0.0, 2.0)

RF_YLIM = (-0.25, 0.6)
HV_YLIM = (0.0, 2.0)
PH_YLIM = (0.0, 5.0)

POSTERIOR_VS_COLOR = "gray"
POSTERIOR_VS_ALPHA = 0.10
POSTERIOR_RF_ALPHA = 0.20
POSTERIOR_DISP_ALPHA = 0.50

AX_FACE = "lightyellow"

FONT_LABEL = {
    "family": "serif",
    "color": "darkblue",
    "size": 20,
    "weight": "bold",
}
FONT_TITLE_BIG = {
    "family": "serif",
    "color": "darkgreen",
    "size": 30,
    "weight": "bold",
}
FONT_TITLE_PANEL = {
    "family": "serif",
    "color": "black",
    "size": 25,
    "weight": "bold",
}
FONT_LABEL_SMALL = {
    "family": "serif",
    "color": "darkblue",
    "size": 11,
    "weight": "bold",
}
FONT_TITLE_SMALL = {
    "family": "serif",
    "color": "darkgreen",
    "size": 11,
    "weight": "bold",
}

velmodfilename="NTW_1D_Liu21_modified.dat" # name of 1D model located in Project/Vel_mod
# =============================================================================
# HELPERS
# =============================================================================
def die(msg):
    print(msg)
    sys.exit(1)


def file_len(fname):
    with open(fname, "r") as fobj:
        for i, _ in enumerate(fobj):
            pass
    return i + 1 if "i" in locals() else 0


def loadtxt_1d(path, usecols=None):
    arr = np.loadtxt(path, usecols=usecols)
    return np.atleast_1d(arr)


def loadtxt_columns(path, usecols):
    data = np.loadtxt(path, usecols=usecols, unpack=True)
    if isinstance(data, np.ndarray):
        data = np.atleast_1d(data)
    else:
        data = tuple(np.atleast_1d(x) for x in data)
    return data


def safe_read_csv(path, **kwargs):
    if not os.path.isfile(path):
        die("File does not exist: {}".format(path))
    return pd.read_csv(path, **kwargs)


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def velmod(df):
    """
    Convert layer-bottom style model into step-style profile.
    Input columns expected:
        thick, vs, id
    """
    b = df["id"].ne(df["id"].shift())
    if len(b) > 0:
        b.iloc[0] = False

    ins = df.loc[b, ["thick", "vs", "id"]].copy()
    ins["thick"] = df["thick"].shift().loc[b].values

    # use mergesort for stable behavior in older pandas
    df2 = pd.concat([ins, df]).sort_index(kind="mergesort").reset_index(drop=True)
    df2 = pd.concat([df2, df2.tail(1)], ignore_index=True)

    df2["dep"] = df2["thick"].shift(fill_value=0.0)
    df2["thickness"] = df2["thick"] - df2["dep"]
    return df2


def read_minmisfit_keep34(minmisfitfile, depth_round=6):
    """
    Read MC.<sta>.minmisfit.mod (depth, vs, flag) and remove duplicate depths
    by keeping the later row, EXCEPT keep both rows when duplicate depth
    represents a flag transition 3 -> 4.
    """
    df = pd.read_csv(
        minmisfitfile,
        sep=r"\s+",
        header=None,
        names=["depth", "vs", "flag"]
    )

    if df.empty:
        return df

    df["flag"] = df["flag"].astype(int)
    df["depth_key"] = df["depth"].round(depth_round)

    keep_idx = []
    for idx, row in df.iterrows():
        if not keep_idx:
            keep_idx.append(idx)
            continue

        prev_idx = keep_idx[-1]
        same_depth = (df.at[prev_idx, "depth_key"] == row["depth_key"])

        if not same_depth:
            keep_idx.append(idx)
            continue

        prev_flag = df.at[prev_idx, "flag"]
        cur_flag = row["flag"]

        if prev_flag == 3 and cur_flag == 4:
            keep_idx.append(idx)
        else:
            keep_idx[-1] = idx

    out = df.loc[keep_idx, ["depth", "vs", "flag"]].reset_index(drop=True)
    return out


def parse_inversion_flags(indatafile):
    if not os.path.exists(indatafile):
        die("File does not exist: {}".format(indatafile))

    data = np.genfromtxt(indatafile, dtype=int)
    data = np.atleast_1d(data)

    if len(data) < 4:
        die("Bad format in {}. Need at least 4 integer values.".format(indatafile))

    iph = int(data[0])
    igv = int(data[1])
    ihv = int(data[2])
    irf = int(data[3])
    return iph, igv, ihv, irf


def parse_plot_type(connectorfile, line_number=8):
    if not os.path.isfile(connectorfile):
        die("Connector file does not exist: {}".format(connectorfile))

    with open(connectorfile, "r") as fobj:
        for i, line in enumerate(fobj, start=1):
            if i == line_number:
                parts = line.strip().split()
                if not parts:
                    die("Empty plot_type line in {}".format(connectorfile))
                try:
                    return int(parts[0])
                except ValueError:
                    die("Invalid plot_type in line {} of {}".format(line_number, connectorfile))

    die("Could not read line {} from {}".format(line_number, connectorfile))


def parse_report_files(reportfiletmp, reportfile):
    """
    Replace shell grep/awk logic with pure Python parsing.

    From output_*_tmp.txt:
        lines containing "test MC.C"
        use token 3 as nsamplesvisit
        use token 7 as misfit

    From output_*.txt:
        line containing "kai_cri"
        use token 1 as kaimin
        lines containing "posterior:"
        use token 1 as posterior indices
    """
    nsamplesvisit = []
    misfitfn = []
    kaimin = None
    posterior_idx = []

    if os.path.isfile(reportfiletmp):
        with open(reportfiletmp, "r") as fobj:
            for line in fobj:
                if "test MC.C" in line:
                    parts = line.strip().split()
                    if len(parts) > 7:
                        try:
                            nsamplesvisit.append(float(parts[3]))
                            misfitfn.append(float(parts[7]))
                        except ValueError:
                            pass

    if os.path.isfile(reportfile):
        with open(reportfile, "r") as fobj:
            for line in fobj:
                parts = line.strip().split()
                if not parts:
                    continue

                if "kai_cri" in line and len(parts) > 1:
                    try:
                        kaimin = float(parts[1])
                    except ValueError:
                        pass

                if "posterior:" in line and len(parts) > 1:
                    try:
                        posterior_idx.append(int(float(parts[1])))
                    except ValueError:
                        pass

    if len(nsamplesvisit) == 0 or len(misfitfn) == 0:
        die("Could not parse misfit history from {}".format(reportfiletmp))

    if kaimin is None:
        die("Could not parse kai_cri from {}".format(reportfile))

    return np.array(nsamplesvisit), np.array(misfitfn), float(kaimin), np.array(posterior_idx, dtype=int)


def read_posterior_periods(firstline_file, nobs_file):
    """
    Infer number of periods from posterior file first line and
    corresponding Initial.* file.
    """
    with open(firstline_file, "r") as ff:
        firstline = ff.readline().rstrip()

    vals = firstline.split()
    vals = [float(v) for v in vals[2:]]

    indata = pd.read_csv(nobs_file, delim_whitespace=True, names=["per", "vel"])

    tmp = 0.0
    countper = 0
    periods = []

    for ii in range(len(indata)):
        if vals[ii] >= tmp or vals[ii] > 10:
            tmp = vals[ii]
            countper += 1
            periods.append(vals[ii])
        else:
            break

    return countper, periods


def read_posterior_rf_blocks(posteriorfilerf):
    posteriorrftime = []
    posteriorfval = []

    startread = 0
    starttmp = []
    endtmp = []

    with open(posteriorfilerf, "r") as fobj:
        for num, line in enumerate(fobj, 1):
            if ">" in line and startread == 0:
                startread = 1
                starttmp.append(num)
                continue
            elif ">" in line and startread == 1:
                endtmp.append(num - 1)
                starttmp.append(num)
                continue
        endtmp.append(num)

    for ii in range(len(endtmp)):
        timetmp, rftmp = np.genfromtxt(
            posteriorfilerf,
            usecols=(0, 1),
            skip_header=starttmp[ii],
            max_rows=(endtmp[ii] - starttmp[ii]),
            unpack=True
        )
        posteriorrftime.append(np.atleast_1d(timetmp))
        posteriorfval.append(np.atleast_1d(rftmp))

    return posteriorrftime, posteriorfval


def read_posterior_vs_and_depth(posteriorfile, plot_type):
    with open(posteriorfile, "r") as ff:
        posterior_vs_tmp = [x.strip().split() for x in ff if x.strip()]

    with open(posteriorfile + "_depth", "r") as ff:
        posterior_depth_tmp = [x.strip().split() for x in ff if x.strip()]

    if len(posterior_vs_tmp) != len(posterior_depth_tmp):
        die(
            "Inconsistent posterior models in {}: Vs rows={} depth rows={}".format(
                posteriorfile, len(posterior_vs_tmp), len(posterior_depth_tmp)
            )
        )

    posterior_vs = []
    posterior_depth = []

    for i, (vs_row, dep_row) in enumerate(zip(posterior_vs_tmp, posterior_depth_tmp)):
        if len(vs_row) != len(dep_row):
            die(
                "Inconsistent point count in model #{} of {}: Vs={} depth={}".format(
                    i, posteriorfile, len(vs_row), len(dep_row)
                )
            )

        try:
            vs = np.array([float(x) for x in vs_row], dtype=float)
            dep = np.array([float(x) for x in dep_row], dtype=float)
        except Exception as exc:
            die("Non-numeric entry in model #{}: {}".format(i, exc))

        mask = np.isfinite(vs) & np.isfinite(dep)
        vs = vs[mask]
        dep = dep[mask]

        order = np.argsort(dep, kind="mergesort")
        dep = dep[order]
        vs = vs[order]

        if plot_type == 1:
            keep = np.r_[True, dep[1:] != dep[:-1]]
            dep = dep[keep]
            vs = vs[keep]
            if dep.size == 0:
                die("Model #{} empty after de-duplicating depths".format(i))

        posterior_vs.append(vs)
        posterior_depth.append(dep)

    return posterior_vs, posterior_depth


def read_model_space(modelspacemaxfile, modelspaceminfile, plot_type):
    if plot_type == 0:
        mspacemindepth, mspacemin = np.genfromtxt(
            modelspacemaxfile, skip_footer=6, usecols=(0, 1), unpack=True
        )
        mspacemaxdepth, mspacemax = np.genfromtxt(
            modelspaceminfile, skip_footer=6, usecols=(0, 1), unpack=True
        )
        return mspacemindepth, mspacemin, mspacemaxdepth, mspacemax

    mspacemindepth, mspacemin = np.genfromtxt(modelspacemaxfile, usecols=(0, 1), unpack=True)
    mspacemaxdepth, mspacemax = np.genfromtxt(modelspaceminfile, usecols=(0, 1), unpack=True)

    bad_idx = np.where(mspacemin == -1)[0]
    if bad_idx.size:
        cut = bad_idx[0]
        mspacemindepth = mspacemindepth[:cut]
        mspacemin = mspacemin[:cut]
        mspacemaxdepth = mspacemaxdepth[:cut]
        mspacemax = mspacemax[:cut]

    if mspacemindepth.size:
        keep = np.r_[True, mspacemindepth[1:] != mspacemindepth[:-1]]
        mspacemindepth = mspacemindepth[keep]
        mspacemin = mspacemin[keep]

    if mspacemaxdepth.size:
        keep = np.r_[True, mspacemaxdepth[1:] != mspacemaxdepth[:-1]]
        mspacemaxdepth = mspacemaxdepth[keep]
        mspacemax = mspacemax[keep]

    return mspacemindepth, mspacemin, mspacemaxdepth, mspacemax


def build_paths(sta, mdir, pwd):
    project_root = os.path.dirname(pwd)
    data_root = os.path.join(project_root, "data")
    sta_dir = os.path.join(mdir, sta)

    paths = {
        "sta_dir": sta_dir,
        "indatafile": os.path.join(data_root, "{}_data".format(sta), "in.data_{}".format(sta)),
        "phdatafile": os.path.join(data_root, "{}_data".format(sta), "{}.ph".format(sta)),
        "gvdatafile": os.path.join(data_root, "{}_data".format(sta), "{}.gv".format(sta)),
        "hvdatafile": os.path.join(data_root, "{}_data".format(sta), "{}.HV".format(sta)),
        "rfdatafile": os.path.join(data_root, "{}_data".format(sta), "{}.RF".format(sta)),
        "modfile": os.path.join(data_root, "{}_data".format(sta), "mod.{}".format(sta)),
        "connectorfile": os.path.join(data_root, "{}_data".format(sta), "in.connector"),
        "phdataprintfile": os.path.join(pwd, sta, "Initial.ph"),
        "hvdataprintfile": os.path.join(pwd, sta, "Initial.hv"),
        "rfdataprintfile": os.path.join(pwd, sta, "Initial.rf"),
        "velmodfile": os.path.join(project_root, "Vel_mod", velmodfilename),
        "posteriorfile": os.path.join(sta_dir, "MC.{}.all.value".format(sta)),
        "posteriorfileph": os.path.join(sta_dir, "MC.{}.all.ph".format(sta)),
        "posteriorfilerf": os.path.join(sta_dir, "MC.{}.all.rf".format(sta)),
        "posteriorfilee": os.path.join(sta_dir, "MC.{}.all.e".format(sta)),
        "inputmodelfile": os.path.join(sta_dir, "Initial.mod"),
        "inputmodelfileleft": os.path.join(sta_dir, "Initial_left.mod"),
        "inputmodelfileright": os.path.join(sta_dir, "Initial_right.mod"),
        "all_parameters": os.path.join(sta_dir, "MC.{}.all.parameters".format(sta)),
        "modelspacemaxfile": os.path.join(sta_dir, "MC.{}.maxmodelvalue".format(sta)),
        "modelspaceminfile": os.path.join(sta_dir, "MC.{}.minmodelvalue".format(sta)),
        "reportfiletmp": os.path.join(sta_dir, "output_{}_tmp.txt".format(sta)),
        "reportfile": os.path.join(sta_dir, "output_{}.txt".format(sta)),
        "meanpara": os.path.join(sta_dir, "MC.{}.ave.paras".format(sta)),
        "avgmodel_file": os.path.join(sta_dir, "MC.{}.acc.average.mod".format(sta)),
        "mod_q_file": os.path.join(sta_dir, "MC.{}.mod.q".format(sta)),
        "startfileph": os.path.join(sta_dir, "Initial.ph"),
        "startfilehv": os.path.join(sta_dir, "Initial.hv"),
        "startfilerf": os.path.join(sta_dir, "Initial.rf"),
        "avemod_file": os.path.join(sta_dir, "MC.{}.ave.mod".format(sta)),
        "accavg_info": os.path.join(sta_dir, "MC.{}.acc.average.info".format(sta)),
        "minmisfit_file": os.path.join(sta_dir, "MC.{}.minmisfit.mod".format(sta)),
        "minmisfit_info": os.path.join(sta_dir, "MC.{}.minmisfit.info".format(sta)),
        "minmisfit_group": os.path.join(sta_dir, "MC.{}.minmisfit.mod.group".format(sta)),
        "final_ph": os.path.join(sta_dir, "MC.{}.final.ph".format(sta)),
        "mean_ph": os.path.join(sta_dir, "MC.{}.mean.ph".format(sta)),
        "final_e": os.path.join(sta_dir, "MC.{}.final.e".format(sta)),
        "mean_e": os.path.join(sta_dir, "MC.{}.mean.e".format(sta)),
        "final_rf": os.path.join(sta_dir, "MC.{}.final.rf".format(sta)),
        "mean_rf": os.path.join(sta_dir, "MC.{}.mean.rf".format(sta)),
        "minmisfit_p_disp": os.path.join(sta_dir, "MC.{}.minmisfit.p.disp".format(sta)),
        "minmisfit_e_disp": os.path.join(sta_dir, "MC.{}.minmisfit.e.disp".format(sta)),
        "minmisfit_rf": os.path.join(sta_dir, "MC.{}.minmisfit.rf".format(sta)),
        "out_main_png": os.path.join(sta_dir, "{}_MCMC.png".format(sta)),
        "out_misfit_png": os.path.join(sta_dir, "{}_IterVsMisfit.png".format(sta)),
        "out_avg_txt": "{}_Vsv_average.txt".format(sta),
        "out_posterior_txt": "{}_Vsv.txt".format(sta),
    }
    return paths


# =============================================================================
# MAIN
# =============================================================================
def main():
    if len(sys.argv) != 3:
        print("proper usage:")
        print("python inversion_plot_clean.py Sta LocationOfStaDir")
        sys.exit(1)

    sta = str(sys.argv[1])
    mdir = str(sys.argv[2])
    pwd = os.getcwd()

    print("current directory: {}".format(pwd))
    paths = build_paths(sta, mdir, pwd)
    print("Input velocity model: {}".format(paths["velmodfile"]))

    iph, igv, ihv, irf = parse_inversion_flags(paths["indatafile"])
    print("making figure for station: {}".format(sta))

    plot_type = parse_plot_type(paths["connectorfile"])
    if plot_type == 1:
        print("Plot the model in gradient style")
    elif plot_type == 0:
        print("Plot the model in layer-cake style")
    else:
        die("Error! unknown type of model plot")

    nmodelsprint = file_len(paths["all_parameters"])
    print("number of model: {}".format(nmodelsprint))
    print("posteriorfile: {}".format(paths["posteriorfile"]))

    # -------------------------------------------------------------------------
    # Read report / misfit info
    # -------------------------------------------------------------------------
    if not (os.path.exists(paths["reportfile"]) or os.path.exists(paths["reportfiletmp"])):
        die("Missing report files: {} or {}".format(paths["reportfile"], paths["reportfiletmp"]))

    nsamplesvisit, misfitfn, kaimin, indices_to_plot = parse_report_files(
        paths["reportfiletmp"], paths["reportfile"]
    )

    print("kaicri...: {}".format(kaimin))

    # -------------------------------------------------------------------------
    # Read initial / input models
    # -------------------------------------------------------------------------
    inputmodel0 = pd.read_csv(
        paths["velmodfile"], sep=r"\s+", usecols=[0, 1, 7],
        names=["dep", "vs", "vp"], header=None
    )
    inputmodel = pd.read_csv(
        paths["inputmodelfile"], sep=r"\s+", usecols=[0, 1, 7],
        names=["thick", "vs", "id"], header=None
    )
    inputmodelleft = pd.read_csv(
        paths["inputmodelfileleft"], sep=r"\s+", usecols=[0, 1, 7],
        names=["thick", "vs", "id"], header=None
    )
    inputmodelright = pd.read_csv(
        paths["inputmodelfileright"], sep=r"\s+", usecols=[0, 1, 7],
        names=["thick", "vs", "id"], header=None
    )

    inputmodel = velmod(inputmodel)
    inputmodelleft = velmod(inputmodelleft)
    inputmodelright = velmod(inputmodelright)
    mask = inputmodelright["id"] != 4 # none mantle layers
    inputmodelright.loc[mask & (inputmodelright["vs"] > 4.9), "vs"] = 4.9 # for none mantle layer boundary <= 4.9

    inputmodel["dvs"] = inputmodel["vs"].diff()
    moho_idx = inputmodel["dvs"].idxmax() - 1
    moho_depth = inputmodel.loc[moho_idx, "dep"]

    est_moho_depth = np.atleast_1d(np.loadtxt(paths["meanpara"], usecols=(-1,)))

    # -------------------------------------------------------------------------
    # Read posterior dispersion periods
    # -------------------------------------------------------------------------
    countper_ph, posteriorphper = read_posterior_periods(
        paths["posteriorfileph"], paths["phdataprintfile"]
    )
    posteriorph = np.loadtxt(
        paths["posteriorfileph"],
        usecols=range(countper_ph + 2, len(open(paths["posteriorfileph"]).readline().split()))
    )

    inhvdata = pd.read_csv(paths["hvdataprintfile"], delim_whitespace=True, names=["per", "vel"])
    countper_e, posterioreper = read_posterior_periods(
        paths["posteriorfilee"], paths["hvdataprintfile"]
    )
    posteriore = np.loadtxt(
        paths["posteriorfilee"],
        usecols=range(countper_e + 2, len(open(paths["posteriorfilee"]).readline().split()))
    )

    posteriorph = np.atleast_2d(posteriorph)
    posteriore = np.atleast_2d(posteriore)

    # -------------------------------------------------------------------------
    # Read posterior RF / Vs
    # -------------------------------------------------------------------------
    posteriorrftime, posteriorfval = read_posterior_rf_blocks(paths["posteriorfilerf"])
    posteriorVs, posteriordepth = read_posterior_vs_and_depth(paths["posteriorfile"], plot_type)

    with open(paths["out_posterior_txt"], "w") as fobj:
        for i in range(len(posteriordepth)):
            for j in range(len(posteriordepth[i])):
                fobj.write("{} {}\n".format(posteriordepth[i][j], posteriorVs[i][j]))

    # -------------------------------------------------------------------------
    # Read model space
    # -------------------------------------------------------------------------
    mspacemindepth, mspacemin, mspacemaxdepth, mspacemax = read_model_space(
        paths["modelspacemaxfile"], paths["modelspaceminfile"], plot_type
    )

    # -------------------------------------------------------------------------
    # Final average and start model
    # -------------------------------------------------------------------------
    avgdepth, avgvs = np.loadtxt(paths["avgmodel_file"], usecols=(0, 1), unpack=True)
    rmdepthtmp, rmvstmp, rmvptmp, rmrhotmp = np.loadtxt(
        paths["mod_q_file"], usecols=(0, 1, 2, 3), unpack=True
    )

    rmdepth = []
    rmvs = []
    rmvpt = []
    rmrhot = []

    for ii in range(len(rmdepthtmp)):
        if plot_type == 1:
            if (ii < len(rmdepthtmp) - 1) and (rmdepthtmp[ii + 1] == rmdepthtmp[ii]):
                continue
        rmdepth.append(rmdepthtmp[ii])
        rmvs.append(rmvstmp[ii])
        rmvpt.append(rmvptmp[ii])
        rmrhot.append(rmrhotmp[ii])

    rmthick = np.diff(rmdepth)
    rmvsf = []
    rmvp = []
    rmrho = []
    for ii in range(len(rmdepth) - 1):
        rmvsf.append((rmvs[ii + 1] + rmvs[ii]) * 0.5)
        rmvp.append((rmvpt[ii + 1] + rmvpt[ii]) * 0.5)
        rmrho.append((rmrhot[ii + 1] + rmrhot[ii]) * 0.5)

    rmvsucvm = rmvs
    rmdepthucvm = rmdepth

    avgdepthf = avgdepth
    avgvsf = avgvs

    # -------------------------------------------------------------------------
    # Starting model response
    # -------------------------------------------------------------------------
    fcheckper, fcheckRC = np.loadtxt(paths["startfileph"], usecols=(0, 1), unpack=True)
    fcheckperhv, fcheckHV = np.loadtxt(paths["startfilehv"], usecols=(0, 1), unpack=True)
    timerf0, rf0 = np.genfromtxt(paths["startfilerf"], usecols=(0, 1), unpack=True)

    # -------------------------------------------------------------------------
    # Final average model / min misfit model
    # -------------------------------------------------------------------------
    avgdepthnew, avgvsnew, avgnewstd, avgnewstdn, avgnewstdp = np.loadtxt(
        paths["avemod_file"], usecols=(0, 1, 2, 5, 6), unpack=True
    )

    with open(paths["out_avg_txt"], "w") as fobj:
        print("write in!")
        for i in range(len(avgdepthnew)):
            fobj.write("{} {}\n".format(avgdepthnew[i], avgvsnew[i]))
        print("write in end!")

    mmfdf = read_minmisfit_keep34(paths["minmisfit_file"])
    mindepthnew = mmfdf["depth"].to_numpy()
    minvsnew = mmfdf["vs"].to_numpy()

    with open(paths["accavg_info"], "r") as misff:
        firstline = misff.readline().rstrip()
    x = firstline.split()
    misfittpm = np.float64(x[2])
    misfitprint = str(round(misfittpm, 2))

    with open(paths["minmisfit_info"], "r") as misff:
        firstline = misff.readline().rstrip()
    x = firstline.split()
    minmisfit = np.float64(x[2])

    # -------------------------------------------------------------------------
    # Final inversion predictions
    # -------------------------------------------------------------------------
    if iph == 1:
        mcper, mcpf, mcpo, mcpounc = np.loadtxt(paths["final_ph"], usecols=(0, 1, 2, 3), unpack=True)
        mcper0, mcpo0 = np.loadtxt(paths["mean_ph"], usecols=(0, 1), unpack=True)
        mcper_data, mcpo_data, mchvounc_data = np.loadtxt(
            paths["phdatafile"], usecols=(0, 1, 2), unpack=True
        )

    if ihv == 1:
        mceper, mchvf, mchvo, mchvounc = np.loadtxt(paths["final_e"], usecols=(0, 1, 2, 3), unpack=True)
        mceper0, mchvo0 = np.loadtxt(paths["mean_e"], usecols=(0, 1), unpack=True)
        mchvounc = np.loadtxt(paths["hvdatafile"], usecols=(2), unpack=True)

    mcrftime, mcrfval, mcrotime, mcrfo, mcrfounc = np.loadtxt(
        paths["final_rf"], usecols=(0, 1, 2, 3, 4), unpack=True
    )
    mcrftime0, mcrfval0, mcrotime0, mcrfo0, mcrfounc0 = np.loadtxt(
        paths["mean_rf"], usecols=(0, 1, 2, 3, 4), unpack=True
    )

    mincper, mincpf, mincpo, mincpounc = np.loadtxt(
        paths["minmisfit_p_disp"], usecols=(0, 1, 2, 3), unpack=True
    )
    minceper, minchvf, minchvo, minchvounc = np.loadtxt(
        paths["minmisfit_e_disp"], usecols=(0, 1, 2, 3), unpack=True
    )
    mincrftime, mincrfval, mincrotime, mcinrfo, mincrfounc = np.loadtxt(
        paths["minmisfit_rf"], usecols=(0, 1, 2, 3, 4), unpack=True
    )

    # -------------------------------------------------------------------------
    # Extra info
    # -------------------------------------------------------------------------
    misfitavgmodel = None
    if os.path.isfile(paths["accavg_info"]):
        with open(paths["accavg_info"], "r") as ff:
            firstline = ff.readline().rstrip()
        x = firstline.split()
        misfitavgmodel = float(x[2])

    with open(paths["meanpara"], "r") as ff:
        firstline = ff.readline().rstrip()
    x = firstline.split()
    sthickf = float(x[len(x) - 2])
    cthickf = sthickf + float(x[len(x) - 1])

    # -------------------------------------------------------------------------
    # PLOTTING MAIN FIGURE
    # -------------------------------------------------------------------------
    print("making figure...")

    fig = plt.figure(int(nmodelsprint), figsize=FIGSIZE_MAIN)

    ax1 = plt.subplot2grid((3, 5), (1, 0), rowspan=2, colspan=2)
    ax2 = plt.subplot2grid((3, 5), (0, 2), colspan=3)
    ax3 = plt.subplot2grid((3, 5), (1, 2), colspan=3)
    ax4 = plt.subplot2grid((3, 5), (0, 0), colspan=2)
    ax5 = plt.subplot2grid((3, 5), (2, 2), colspan=3)

    for ax in (ax1, ax2, ax3, ax4, ax5):
        ax.set_facecolor(AX_FACE)

    # ---- Vs full depth
    for jj in range(len(posteriorVs) - 1):
        ax1.plot(
            posteriorVs[jj], posteriordepth[jj],
            c=POSTERIOR_VS_COLOR, alpha=POSTERIOR_VS_ALPHA, lw=2
        )

    ax1.plot(minvsnew, mindepthnew, "bs--", lw=2, ms=10, label="MinMisfit Vs", zorder=7, mec="k")
    ax1.plot(avgvsnew, avgdepthnew, "mo--", mfc="w", lw=2, ms=10, label="Final Vs", zorder=7, mec="k")
    ax1.errorbar(avgvsnew, avgdepthnew, xerr=avgnewstd, fmt="-o", ecolor="m", elinewidth=2, capsize=4, alpha=1.0)
    # Model spacing plot 
    # split crust / mantle
    right_non_mantle = inputmodelright[inputmodelright["id"] != 4].copy()
    right_mantle     = inputmodelright[inputmodelright["id"] == 4].copy()

    left_non_mantle  = inputmodelleft[inputmodelleft["id"] != 4].copy()
    left_mantle      = inputmodelleft[inputmodelleft["id"] == 4].copy()

    # non-mantle: normal dashed line
    ax1.plot(
        right_non_mantle["vs"], right_non_mantle["dep"],
        dashes=[8, 8, 16, 8], lw=2, color="g", label="Model Space"
    )
    ax1.plot(
        left_non_mantle["vs"], left_non_mantle["dep"],
        dashes=[8, 8, 16, 8], lw=2, color="g"
    )

    # connect non-mantle to mantle
    if (not right_non_mantle.empty) and (not right_mantle.empty):
        ax1.plot(
            [right_non_mantle["vs"].iloc[-1], right_mantle["vs"].iloc[0]],
            [right_non_mantle["dep"].iloc[-1], right_mantle["dep"].iloc[0]],
            dashes=[8, 8, 16, 8], lw=2, color="g"
        )

    if (not left_non_mantle.empty) and (not left_mantle.empty):
        ax1.plot(
            [left_non_mantle["vs"].iloc[-1], left_mantle["vs"].iloc[0]],
            [left_non_mantle["dep"].iloc[-1], left_mantle["dep"].iloc[0]],
            dashes=[8, 8, 16, 8], lw=2, color="g"
        )

    # mantle: dotted line
    ax1.plot(right_mantle["vs"], right_mantle["dep"], ":", lw=2, color="g")
    ax1.plot(left_mantle["vs"],  left_mantle["dep"],  ":", lw=2, color="g")

    # ax1.plot(inputmodelleft["vs"], inputmodelleft["dep"], "g", label="ModelSpace", dashes=[8, 8, 16, 8])
    # ax1.plot(inputmodelright["vs"], inputmodelright["dep"], "g", dashes=[8, 8, 16, 8])

    ax1.plot(inputmodel["vs"], inputmodel["dep"], "r^-", lw=2.5, ms=10, label="Start Vs")

    ax1.legend(loc="lower left", fontsize=15)
    ticks = np.r_[np.arange(0, 5.0, 1.0), np.arange(5, 51, 5)]
    ax1.set_yticks(ticks)
    ax1.set_xlabel("Vs (km/s)", fontdict=FONT_LABEL)
    ax1.set_ylabel("Depth (km)", fontdict=FONT_LABEL)
    ax1.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax1.tick_params(labeltop=False, labelsize=20)
    ax1.set_ylim(VS_YLIM_FULL)
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.set_xlim(VS_XLIM)
    ax1.grid(color="k", axis="both", linestyle="-", linewidth=2, alpha=0.1, zorder=2)

    title = "Vs of station: [{}] - MinimumMisfit: {} - #posteriorModels:{}".format(
        sta, round(minmisfit, 2), len(posteriorfval)
    )
    fig.suptitle(
        title,
        x=0.5, y=0.985,
        ha="center",
        va="top",
        fontfamily="serif",
        color="darkgreen",
        fontsize=30,
        fontweight="bold"
    )

    # ---- Vs shallow
    ax4.plot(inputmodelleft["vs"], inputmodelleft["dep"], "g", label="ModelSpace", dashes=[8, 8, 16, 8])
    ax4.plot(inputmodelright["vs"], inputmodelright["dep"], "g", dashes=[8, 8, 16, 8])

    for jj in range(len(posteriorVs) - 1):
        ax4.plot(
            posteriorVs[jj], posteriordepth[jj],
            c=POSTERIOR_VS_COLOR, alpha=POSTERIOR_VS_ALPHA, lw=2
        )

    ax4.plot(minvsnew, mindepthnew, "bs--", lw=2, ms=10, label="MinMisfit Vs", zorder=7, mec="k")

    ax4.plot(avgvsnew, avgdepthnew, "mo--", mfc="w", lw=2, ms=10, label="Final Vs", zorder=7, mec="k")
    ax4.errorbar(avgvsnew, avgdepthnew, xerr=avgnewstd, fmt="-o", ecolor="m", elinewidth=2, capsize=4, alpha=1.0)
    ax4.plot(inputmodel["vs"], inputmodel["dep"], "r^-", lw=2.5, ms=10, label="Start Vs")

    # ax4.plot(inputmodel0["vs"], inputmodel0["dep"], "ko-", lw=2.5, ms=10, alpha=0.5, label="Input model")

    ax4.set_yticks(np.arange(0, int(round(np.max(rmdepthucvm), 0)), 5.0))
    ax4.set_xlabel("Vs (km/s)", fontdict=FONT_LABEL)
    ax4.set_ylabel("Depth (km)", fontdict=FONT_LABEL)
    ax4.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax4.tick_params(labeltop=True, labelsize=20)
    ax4.set_ylim(VS_YLIM_SHALLOW)
    ax4.set_yticks(np.arange(0, 2.0, 0.25))
    ax4.set_ylim(ax4.get_ylim()[::-1])
    ax4.set_xlim(VS_XLIM)
    ax4.grid(color="k", axis="both", linestyle="-", linewidth=2, alpha=0.1, zorder=2)

    # ---- Phase velocity
    ax3.set_title("Phase Velocity", loc="left", fontdict=FONT_TITLE_PANEL)
    if iph == 1:
        for ii in range(len(posteriorph[:, 0]) - 1):
            ax3.plot(posteriorphper, posteriorph[ii], c="gray", alpha=POSTERIOR_DISP_ALPHA, lw=2)

        ax3.errorbar(
            mcper_data, mcpo_data, yerr=mchvounc_data,
            fmt="o", color="k", ecolor="k", ms=7,
            elinewidth=1.5, capthick=1.5,
            label="Data Vph", zorder=6
        )
        ax3.plot(mincper, mincpf, "bs", lw=1, ms=6, label="Minmisfit Vph", zorder=7, mec="k")
        ax3.plot(mcper0, mcpo0, "wo", lw=1, ms=6, label="Final Vph", zorder=7, mec="k")

    ax3.plot(fcheckper, fcheckRC, "r^", label="Starting Vph", ms=10, zorder=5)
    ax3.legend(loc="best", fontsize=15)
    ax3.set_xlabel("Period (s)", fontdict=FONT_LABEL)
    ax3.set_ylabel("Vph (km/s)", fontdict=FONT_LABEL)
    ax3.set_ylim(PH_YLIM)
    ax3.set_xlim([posteriorphper[0] - 0.5, posteriorphper[-1] + 0.5])
    ax3.grid(color="k", axis="both", linestyle="-", linewidth=2, alpha=0.1, zorder=2)
    ax3.tick_params(labelsize=20)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")

    # ---- H/V
    ax5.set_title("H/V", loc="left", fontdict=FONT_TITLE_PANEL)
    if ihv == 1:
        for ii in range(len(posteriore[:, 0]) - 1):
            ax5.plot(posterioreper, posteriore[ii], c="gray", alpha=POSTERIOR_DISP_ALPHA, lw=2)

        ax5.errorbar(
            mceper, mchvo, yerr=mchvounc,
            fmt="o", color="k", ms=7, label="Data H/V",
            ecolor="k", elinewidth=1.5, capthick=1.5, zorder=6
        )
        ax5.plot(minceper, minchvf, "bs", lw=1, ms=7, label="minmisfit H/V", zorder=7, mec="k")
        ax5.plot(mceper0, mchvo0, "wo", lw=1, ms=7, label="Final H/V", zorder=7, mec="k")

    ax5.plot(fcheckperhv[:len(fcheckHV)], fcheckHV, "r^", ms=10, label="Starting H/V", zorder=5)
    ax5.legend(loc="upper right", fontsize=15)
    ax5.set_xlabel("Period (s)", fontdict=FONT_LABEL)
    ax5.set_ylabel("H/V", fontdict=FONT_LABEL)
    ax5.set_ylim(HV_YLIM)
    ax5.tick_params(labelsize=20)
    ax5.grid(color="k", axis="both", linestyle="-", linewidth=2, alpha=0.1, zorder=2)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")

    # ---- RF
    ax2.set_title("Receiver function", loc="left", fontdict=FONT_TITLE_PANEL)
    if irf == 1:
        ax2.errorbar(
            mcrotime, mcrfo, yerr=mcrfounc,
            fmt=".:", color="k", ecolor="k",
            elinewidth=1.5, capthick=0.0,
            label="Data RFs", zorder=5, alpha=0.5
        )
        ax2.plot(mincrftime, mincrfval, "b-", lw=1, ms=3, mec="k", label="Minmisfit RF", zorder=7)

        for ii in range(len(posteriorfval) - 1):
            ax2.plot(posteriorrftime[ii], posteriorfval[ii], c="gray", alpha=POSTERIOR_RF_ALPHA, lw=2)

        ax2.plot(mcrftime0, mcrfval0, "wo", lw=1, ms=3, label="Final RF", zorder=7, mec="k")

    ax2.plot(timerf0, rf0, "r-", ms=10, label="Starting RFs", zorder=5)
    ax2.legend(loc="upper right", fontsize=15)
    ax2.set_ylim(RF_YLIM)
    ax2.set_xlabel("Time (s)", fontdict=FONT_LABEL)
    ax2.set_ylabel("RF Amp", fontdict=FONT_LABEL)
    ax2.grid(color="k", axis="both", linestyle="-", linewidth=2, alpha=0.1, zorder=2)
    ax2.tick_params(labelsize=20)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    plt.tight_layout(rect=[0, 0, 1, 0.965])

    plt.savefig(paths["out_main_png"], bbox_inches="tight", transparent=False, pad_inches=0.1, dpi=DPI)
    plt.close()

    # -------------------------------------------------------------------------
    # PLOT MISFIT FIGURE
    # -------------------------------------------------------------------------
    print("Plot misfit toward iterations...")

    fig = plt.figure(2, figsize=FIGSIZE_MISFIT)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    title = "Vs of station: [{}] - MinimumMisfit: {} - #posteriorModels:{}".format(
            sta, round(minmisfit, 2), len(posteriorfval)
    )
    fig.suptitle(
        title,
        x=0.5, y=0.985,
        ha="center",
        va="top",
        fontfamily="serif",
        color="darkgreen",
        fontsize=7.5,
        fontweight="bold"
    )
    ax1.axhline(kaimin, color="k", linestyle="--", lw=2)
    ax1.plot(nsamplesvisit, misfitfn, "o", alpha=0.2, ms=4, mec="k", mfc="none")

    for ii in range(len(misfitfn)):
        if ii in indices_to_plot:
            ax1.plot(nsamplesvisit[ii], misfitfn[ii], "ro", alpha=0.2, ms=5, mec="r", mfc="r")

    ax1.set_ylabel("Zoom Misfit", fontdict=FONT_LABEL_SMALL)
    ax1.set_ylim([kaimin * 2, 0])

    ax2.axhline(kaimin, color="k", linestyle="--", lw=2, label="critical misfit")
    ax2.plot(nsamplesvisit, misfitfn, "o", alpha=0.2, ms=4, mec="k", mfc="none")

    for ii in range(len(misfitfn)):
        if ii in indices_to_plot:
            ax2.plot(nsamplesvisit[ii], misfitfn[ii], "ro", alpha=0.2, ms=5, mec="r", mfc="r")

    ax2.set_xlabel("Iteration number", fontdict=FONT_LABEL_SMALL)
    ax2.set_ylabel("Misfit", fontdict=FONT_LABEL_SMALL)
    ax2.set_ylim(ax2.get_ylim()[::-1])

    plt.savefig(paths["out_misfit_png"], bbox_inches="tight", transparent=False, pad_inches=0.1, dpi=DPI)
    plt.close()
    plt.show()

    print("Done.")


if __name__ == "__main__":
    main()