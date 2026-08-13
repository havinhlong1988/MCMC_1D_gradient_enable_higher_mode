# MCMC-CHT Analysis Checkpoint

> Portable session checkpoint for continuing analysis elsewhere (e.g. Claude on Mac).
> Generated: 2026-08-13 · Machine: Linux (vinh-longha) · Project root: `~/Research/MCMC_CHT/20260801_MCMC_CHT_monotonic_vpvs_H14`

This document captures two things:
1. **Checkpoint 1** — state of *our* version (H14 monotonic_vpvs, fundamental-mode only).
2. **Current** — findings on the *colleague's* higher-mode version (CNan), for comparison.

The code itself is **not** in this file — only the analysis/context. To do deep code work on another machine you also need the source files (see "Portability" at the end).

---

## What the code does (one paragraph)

Per-station **Bayesian MCMC (Metropolis) inversion for 1-D shear-velocity (Vs) structure** (crust + mantle), jointly fitting surface-wave dispersion (Rayleigh phase/group), Rayleigh ellipticity (H/V), and receiver functions (RF). C++/Fortran engine in `Codes/MCMC_flex/`, driven by bash scripts (`01`–`06`) that call Python setup (`SetupData/`) and plotting (`AnalyzeResult/`). Three stages: (1) `01_prepare_data.sh` → `SetupData/main.py` builds per-station `mod.{sta}`, `in.para_{sta}`, `in.data_{sta}`, `in.connector`; (2) `02_do_MCMC.sh` runs `do_MC_Para` (inversion) then `do_MC_Para_Post_process_v3` (select+average) then plots; (3) `03`–`06` re-run/re-plot.

---

## CHECKPOINT 1 — H14 monotonic_vpvs version (our version)

### Environment (critical to run)
- Run inside conda env **`MCMC`** = Python 3.7.6, numpy 1.18.1, pandas 1.0.1, matplotlib 3.4.3, gfortran added.
- `conda activate MCMC` **before** running scripts — they call `python` (not `python3`), which only resolves in that env.
- **Install files in project root:** `environment.yml` (portable conda spec), `environment.linux-lock.yml` (exact build-pinned), `requirements.txt` (pip), `INSTALL.md` (system requirements + build steps + macOS caveats). Recreate with `conda env create -f environment.yml`.
- Original documented toolchain was gcc/g++/gfortran 4.8.5. This machine: system g++ 13.3 + conda gfortran 15 → required the fixes below.

### Fixes applied to build + run on this machine
| File | Fix | Why |
|---|---|---|
| `Codes/MCMC_flex/head_c++/gen_random.C` | added `#include <random>` | gcc15 needs it for `std::mt19937` |
| `CALmodel.C`, `CALgroup.C`, `MC.C` | added 7 missing `return` statements | **Root cause of a SIGILL crash** — gcc13 traps (ud2) when a non-void function falls off the end; gcc 4.8.5 didn't. Every freshly-built binary crashed instantly before this fix. |
| `Codes/MCMC_flex/Makefile` | `LDLIBS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib` | so `make` (recompile branch) finds conda libgfortran |
| `02_do_MCMC.sh`, `CHT/MonteCarlo/00_make_file_control.sh` | converted bashisms to POSIX; `sh …`→`bash …`; `rm`→`rm -f`; `stop`(not a cmd)→`break`/`continue`; removed stray double-`break` | script failed under `sh` (dash); false "control not generated" msg |
| `SetupData/parameters.py` | crust gwStep used sediment steps → fixed to `gwStepCrustVel`/`gwStepCrustThick`; set `MC_number_of_cores=1` | latent step-size bug (crust thickness step was 10× too small) |
| `SetupData/src/setup_layercake.py` | `os.exit`→`sys.exit` + real error msg | crash-on-missing-file was itself broken |
| (data) | seeded `CHT/data/` by copying from `CHT/real_data/` | `main.py` writes to `data/` but never creates it; `real_data/` holds raw `.ph` inputs |

### Verified working
- Full chain runs: `do_MC_Para` → `do_MC_Para_Post_process_v3` → plot. Both checkpoints reached: `do_MC_Para done!!!` and `post process finish!!! ready to plot!`.
- Binaries exit code **1 = success** (this codebase's convention; scripts key off stdout strings, not exit code).
- Test station `01390` is **phase-only** (`in.data` = `1 0 0 0`).

### Known OPEN bugs (not yet fixed)
- **Plot crash for phase-only stations**: `AnalyzeResult/inversion_plot_vfinal_flex.py:646` does `np.loadtxt(Initial.hv)` on an empty HV file → `ValueError`. Needs guarding absent HV/gv/RF.
- **Makefile has no header deps** → won't rebuild `do_MC_Para_Post_process_v3` when included `.C`/`.h` change (must force `rm *.o`).

### Forward modeling (surface-wave dispersion)
- Uses **Rayleigh, fundamental mode only**. Gated by:
  - `Codes/MCMC_flex/DISP2/init.f:58` → `mode=1` (compute 1 mode).
  - `Codes/MCMC_flex/DISP2/fast_surf.f:238` → extracts only `cR(i,1)`.
  - `CALforward.C:114` → `cflag=2` forces Rayleigh (Love `kind=1` coded but unused).
- Engine (`calcul.f`, `surfa.f`) IS higher-mode-capable (mode loop, arrays `cR(200,2)`, `nmod=20`) but disabled.

### Post-process: how the final average model is built
- Selects the best **~10%** (`MC_percentage_post_process_select`) of **accepted** posterior models, using percentile misfit cutoffs on joint + per-data-type misfit (`kai_cri`, `minmisfit_disp_o`, `minmisfit_rf_o`; `signs1>0` = accepted).
- Two averages produced:
  - `.acc.average` = **parameter-space mean** (mean of each inversion parameter → forward to one layer-cake). Keeps a **sharp Moho** (single mean crustal thickness).
  - `.ave.mod` = **depth-domain blocky mean** via `get_vs_at_depth_common_grid` (piecewise-constant Vs per layer at 0.05 km grid), with `FORCE_MONOTONIC_MEAN` clamp + ±1σ/min/max. **Smears the Moho jump** because Moho depth varies across the ensemble.
- Blocky-vs-gradient sampler choice is NOT the cause of Moho smearing — the ensemble spread of Moho depth is.

---

## CURRENT — Higher-mode version (colleague "CNan"), for comparison

**Location:** `~/Research/MCMC_CHT/MCMC_joint_inversion_CNan/`
- Code: `Codes_ori/MCMC_PhOnly_Higher_Mode_ori/` (chosen variant; also `_no_mono`, `_0.5`, `_0.9`).
- Full project w/ data + outputs: `MCMC_CHT_hph_misfit_0.7_0.3_ori/`.
- "PhOnly" = phase velocity only (no group); computes **fundamental + 1st higher mode**.

### How it enables the 1st higher mode (5 coordinated changes vs. ours)
1. `DISP2/init.f:58` → **`mode=2`** (ours = 1).
2. `DISP2/fast_surf.f` → **extended subroutine**: 6 extra output args `uR1,uL1,cR1,cL1,rR1,rL1`; extracts `cR0=cR(i,1)` (fundamental) **and `cR1=cR(i,2)`** (higher). NOTE: the *root-level* `fast_surf.f` there is a STALE unused copy; the Makefile links `DISP2/fast_surf.o` (the extended one).
3. `CALforward.C compute_disp` → extended `fast_surf_` prototype + call; period list built from `pper` **+ `hpper`**; stores fundamental→`pvel` (cR0), higher→**`hpvel`** (cR1). One forward call returns both modes.
4. `INITstructure.h` disp struct → adds a parallel higher-mode stream: **`nhpper, hpper, hpvelo, hpvel, unhpvelo`** (mirrors fundamental `npper/pper/pvelo/pvel/unpvelo`).
5. `CALmodel.C` misfit → higher-mode residual `tempv1hmode = Σ((hpvelo−hpvel)/unhpvelo)²`; `phmisfit=√(tempv1hmode/nhpper)`; joint disp misfit combines fundamental+ellipticity+higher: `√((tempv1+tempv3+tempv1hmode)/(npper+neper+nhpper))`. This `_0.7_0.3` variant reweights fundamental/higher 0.7/0.3.

### How input data is organized (fundamental vs higher mode)
- **Two separate files per station, same `{sta}_data/` folder, same 3-column format** `period[s]  phase_velocity[km/s]  uncertainty`:
  - `{sta}.ph`  → **fundamental** Rayleigh phase → routed to `pper…`
  - `{sta}.hph` → **1st higher-mode** Rayleigh phase → routed to `hpper…`
  - (higher-mode phase velocity is faster at a given period, as expected)
- **Control-file `disp` line** pairs them with flags: `disp <surtype> <surn> <flag1> <file1> <flag2> <file2> …`
  - `readdisp1` flag routing: **1 = fundamental phase, 2 = higher-mode phase, 3 = ellipticity, 4 = group**.
  - e.g. `disp 1 2 1 {sta}.ph 2 {sta}.hph` (Rayleigh, 2 inputs). `read_in.C` validates `vargv.size()==3+2*surn`.
- **`.hph` generation:** `make_hph.sh` copies forward-modeled `*_mode1.ph` ("mode1" = 1st higher) → `{sta}.hph`; `Vel_mod_data/make_file.sh` builds them from `merged_higher_weighted.csv`.
- **Caveat:** in this checkout only **97** `.hph` files exist; `make_hph.sh` source path points at Utah CHPC (`/uufs/chpc.utah.edu/…`), so higher-mode picks were generated there — not regenerable here without that source data.

### Build caveat for the CNan version
It's the colleague's **original** code → expect the **same gcc-13 issues we fixed in ours** (missing `return`→SIGILL, boost/`std::mt19937`, `-lgfortran` linking). To run it: use the `MCMC` env + apply the same return fixes + Makefile LDLIBS.

---

## Side-by-side summary

| Aspect | Our version (H14) | CNan higher-mode |
|---|---|---|
| Modes forwarded | Rayleigh fundamental only | Rayleigh fundamental + 1st higher |
| `DISP2/init.f:58` | `mode=1` | `mode=2` |
| `fast_surf` outputs | `cR0,uR0,rR0` | + `cR1,uR1,rR1` |
| Disp struct | `pper/pvelo/pvel…` | + `hpper/hpvelo/hpvel…` |
| Input files/station | `{sta}.ph` (+`.HV`/`.RF`) | `{sta}.ph` **and** `{sta}.hph` |
| `disp` flag for higher | (n/a) | **2** |
| Misfit | fundamental + ellip (+RF) | + higher-mode phase term (0.7/0.3 wt) |

## Open questions / next steps
- Decide whether to port higher-mode support into our version, or run CNan's version directly for comparison (needs `.hph` data + the same build fixes).
- Fix our plotting bug (empty `Initial.hv`) for phase-only stations.
- Consider Moho-aligned averaging if a sharp Moho is wanted in `.ave.mod`.

## Portability note (importing to Claude on Mac)
This `.md` carries **context/findings only**, not the source. On the Mac, Claude can read/paste this file to resume the *analysis thread*, but for code-level work it also needs the relevant source files (or keep code-specific tasks on the machine that has the repo). The two Claude apps do **not** share memory automatically — this file is the manual bridge.
