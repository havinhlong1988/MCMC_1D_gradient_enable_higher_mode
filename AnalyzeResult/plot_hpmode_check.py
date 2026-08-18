#!/usr/bin/env python
"""
plot_hpmode_check.py  --  ZOOM-IN on the phase-velocity fit (observed vs predicted).

Shows, on one tightly-zoomed axis, how well the inverted models reproduce the
OBSERVED Vph for BOTH the fundamental (*.p.disp) and the 1st higher mode
(*.hp.disp).  Uses the post-process outputs:
    MC.{sta}.acc.average.p.disp / .hp.disp   (posterior mean model)
    MC.{sta}.minmisfit.p.disp   / .hp.disp   (min-misfit model)
Each *.disp file has columns:  period  predicted  observed  uncertainty

Colours match the main figure (data=black, minmisfit=blue, average=white);
the higher mode uses the SAME colours, only a different marker.
The y-axis auto-zooms to the data+prediction range (the far-off starting model
is deliberately excluded so the fit is easy to read).

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


def load_ensemble(tag):
    """Parse an all.ph/all.hp file (rows 'M idx <N per> <N val>').
    Return (periods, curves[2D]) or None."""
    f = os.path.join(sdir, "MC.{}.{}".format(sta, tag))
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return None
    periods, rows = None, []
    with open(f) as fh:
        for line in fh:
            t = line.split()
            if len(t) < 4:
                continue
            v = [float(x) for x in t[2:]]   # skip 'M' and model index
            n = len(v) // 2
            if n == 0:
                continue
            periods = v[:n]
            rows.append(v[n:])
    if not rows:
        return None
    return np.array(periods), np.array(rows)


acc_p = load("acc.average.p.disp")
acc_h = load("acc.average.hp.disp")
mmf_p = load("minmisfit.p.disp")
mmf_h = load("minmisfit.hp.disp")
post_p = load_ensemble("all.ph")   # fundamental posterior ensemble
post_h = load_ensemble("all.hp")   # higher-mode posterior ensemble

if acc_p is None and mmf_p is None and acc_h is None and mmf_h is None:
    sys.exit("No .p.disp/.hp.disp output found in {}".format(sdir))

fig, ax = plt.subplots(figsize=(8.5, 5.5))

# track y-range over data + predictions only (exclude the starting model)
yvals = []


def track(y):
    if y is not None:
        yvals.extend(np.asarray(y, dtype=float).ravel().tolist())


# ---- POSTERIOR ensemble in the BACKGROUND ----
#   fundamental: gray lines (kept as current);
#   higher mode: FILLED band in magenta (CNan style).
_fund_post = False
if post_p is not None:
    pper, curves = post_p
    for ii in range(curves.shape[0]):
        ax.plot(pper, curves[ii], c="gray", alpha=0.04, lw=1.0, zorder=1)
    _fund_post = True
    track(np.percentile(curves, 2, axis=0))
    track(np.percentile(curves, 98, axis=0))
if post_h is not None:
    pper, curves = post_h
    order = np.argsort(pper)
    lo = np.percentile(curves, 2, axis=0)[order]
    hi = np.percentile(curves, 98, axis=0)[order]
    ax.fill_between(pper[order], lo, hi, color="magenta", alpha=0.30, lw=0,
                    zorder=1, label="Posterior (higher)")
    track(lo)
    track(hi)
if _fund_post:
    ax.plot([], [], c="gray", alpha=0.6, lw=4, label="Posterior (fund)")  # legend proxy

# ---- OBSERVED (black, with error bars): fundamental=circle, higher=diamond ----
if acc_p is not None:
    per, _, obs, unc = acc_p
    ax.errorbar(per, obs, yerr=unc, fmt="o", ms=7, color="k", ecolor="k",
                elinewidth=1.4, capsize=3, zorder=6, label="Data Vph (fund)")
    track(obs)
if acc_h is not None:
    per, _, obs, unc = acc_h
    ax.errorbar(per, obs, yerr=unc, fmt="D", ms=7, color="k", ecolor="k",
                elinewidth=1.4, capsize=3, zorder=6, label="Data Vph (higher)")
    track(obs)

# ---- MIN-MISFIT prediction (blue): fundamental=square, higher=down-triangle ----
if mmf_p is not None:
    per, pred, _, _ = mmf_p
    ax.plot(per, pred, "-s", color="b", ms=7, lw=1.6, zorder=5, mec="k",
            label="Minmisfit Vph (fund)")
    track(pred)
if mmf_h is not None:
    per, pred, _, _ = mmf_h
    ax.plot(per, pred, "-v", color="b", ms=8, lw=1.6, zorder=5, mec="k",
            label="Minmisfit Vph (higher)")
    track(pred)

# ---- AVERAGE / Final prediction (white face): fundamental=circle, higher=diamond ----
if acc_p is not None:
    per, pred, _, _ = acc_p
    ax.plot(per, pred, "--", color="0.4", lw=1.5, zorder=4)
    ax.plot(per, pred, "o", mfc="w", mec="k", ms=7, zorder=4, label="Final Vph (fund)")
    track(pred)
if acc_h is not None:
    per, pred, _, _ = acc_h
    ax.plot(per, pred, "--", color="0.4", lw=1.5, zorder=4)
    ax.plot(per, pred, "D", mfc="w", mec="k", ms=7, zorder=4, label="Final Vph (higher)")
    track(pred)

# ---- tight zoom on data + prediction ----
if yvals:
    ymin, ymax = min(yvals), max(yvals)
    pad = 0.08 * (ymax - ymin) if ymax > ymin else 0.05
    ax.set_ylim(ymin - pad, ymax + pad)

ax.set_xlabel("Period (s)")
ax.set_ylabel("Rayleigh phase velocity Vph (km/s)")
ax.set_title("Station {} - Vph fit (zoom): observed vs predicted, "
             "fundamental + higher mode".format(sta))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8, ncol=2, loc="best")

out = os.path.join(sdir, "{}_hpmode_check.png".format(sta))
fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
