#!/usr/bin/env python
"""
plot_hpmode_check.py  --  quick check of the higher-mode (.hph) fit.

Overlays observed vs fitted Rayleigh phase velocity for BOTH the fundamental
(*.p.disp) and the 1st higher mode (*.hp.disp), using the post-process outputs:
    MC.{sta}.acc.average.p.disp / .hp.disp   (posterior mean model)
    MC.{sta}.minmisfit.p.disp   / .hp.disp   (min-misfit model)
Each *.disp file has columns:  period  predicted  observed  uncertainty

Usage:  python plot_hpmode_check.py [STA] [MonteCarlo_dir]
Default: STA=00990, dir=../CHT/MonteCarlo (relative to this script)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sta = sys.argv[1] if len(sys.argv) > 1 else "00990"
here = os.path.dirname(os.path.abspath(__file__))
mcdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "CHT", "MonteCarlo")
sdir = os.path.join(mcdir, sta)


def load(tag):
    """Return (per, pred, obs, unc) or None if the file is missing/empty."""
    f = os.path.join(sdir, "MC.{}.{}".format(sta, tag))
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return None
    a = np.atleast_2d(np.loadtxt(f))
    if a.size == 0:
        return None
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


acc_p = load("acc.average.p.disp")
acc_h = load("acc.average.hp.disp")
mmf_p = load("minmisfit.p.disp")
mmf_h = load("minmisfit.hp.disp")

if acc_h is None and mmf_h is None:
    sys.exit("No higher-mode (.hp.disp) output found in {} -- run with .hph enabled first.".format(sdir))

fig, ax = plt.subplots(figsize=(8, 5.5))

# ---- observed (with error bars): fundamental black, higher mode red ----
if acc_p is not None:
    per, _, obs, unc = acc_p
    ax.errorbar(per, obs, yerr=unc, fmt="o", ms=6, color="k", ecolor="k",
                elinewidth=1.4, capsize=3, zorder=6, label="Fundamental — observed")
if acc_h is not None:
    per, _, obs, unc = acc_h
    ax.errorbar(per, obs, yerr=unc, fmt="s", ms=6, color="tab:red", ecolor="tab:red",
                elinewidth=1.4, capsize=3, zorder=6, label="Higher mode — observed")

# ---- fitted curves: acc.average (solid), minmisfit (dashed) ----
if acc_p is not None:
    per, pred, _, _ = acc_p
    ax.plot(per, pred, "-", color="0.35", lw=2, zorder=5, label="Fundamental — acc.average")
if mmf_p is not None:
    per, pred, _, _ = mmf_p
    ax.plot(per, pred, "--", color="0.35", lw=1.6, zorder=5, label="Fundamental — minmisfit")
if acc_h is not None:
    per, pred, _, _ = acc_h
    ax.plot(per, pred, "-", color="tab:red", lw=2, zorder=5, label="Higher mode — acc.average")
if mmf_h is not None:
    per, pred, _, _ = mmf_h
    ax.plot(per, pred, "--", color="darkred", lw=1.6, zorder=5, label="Higher mode — minmisfit")

ax.set_xlabel("Period (s)")
ax.set_ylabel("Rayleigh phase velocity (km/s)")
ax.set_title("Station {} — fundamental + 1st higher-mode phase velocity".format(sta))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, ncol=2)

out = os.path.join(sdir, "{}_hpmode_check.png".format(sta))
fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
