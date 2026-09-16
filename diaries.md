# Work diary — 2026-09-15

Record of every file modified on this machine during the session, in the order
the work happened.

| | |
|---|---|
| Machine | `Longs-MacBook-Pro-148`, macOS (Darwin 25.6.0), Apple Silicon |
| Repository | `/Users/vinhlongha/Research/MCMC/MCMC_1D_gradient_enable_higher_mode` |
| Branch | `fix/macos-portability` |
| Starting commit | `9b3fce4` "Snapshot current working tree: rerun results for 00740/00990" |
| Ending commit | `6c6db0e` |
| Commits added | 6 |
| Pushed? | **No.** Both `main` and `fix/macos-portability` are 6 commits ahead of `origin`. |
| Reference project used | `~/Research/MCMC/20260510_FM_vpvs_pertubation_huang14_vpvs_model` — **read only, never modified** (verified: no file in it has a modification time from this session) |

At the end, `main` was fast-forwarded to `6c6db0e` with `git branch -f main
fix/macos-portability` rather than `git checkout main && git merge`. Reason:
`main` had no `FM/` directory, so checking it out would have deleted the whole
`FM/` tree from disk and restored it on merge — and `in.connector` was open in
vim at the time (a `.in.connector.swp` was present). Moving the ref is
equivalent because `main` was a strict ancestor, and it leaves the working tree
untouched.

---

## Commit history

```
6c6db0e  in.connector: allow "#" comments, and write it self-documenting
323e13c  Make FM runnable: per-project config, parameterised drivers, first station
75070c2  Add a higher-mode master switch in in.connector (force off even with .hph)
20df7b0  CHT: regenerate station inputs with Vp/Vs perturbation (PARTIAL, 17 of 135)
f1d5cbf  Add FM project skeleton, cloned from CHT
85e5502  Port per-sub-layer Vp/Vs perturbation from the huang14 project
9b3fce4  (starting point, not part of this session)
```

---

## 1. `85e5502` — Port per-sub-layer Vp/Vs perturbation

**7 files, +422 / −32.** Ported from the huang14 reference project. Higher-mode
support left untouched.

### Why

The `in.para` types `-1`/`-2` already existed in this tree but were unusable.
`gen_newpara` computed the thickness-ratio block as *everything after velocity +
total thickness* (`npara_tr` was counted and never used), so any `-1`/`-2` row
was swept into the ratio simplex. Measured by appending one `-2` row to station
00740:

- Vp/Vs came out **0.0196 – 0.1810** instead of ~1.73
- the five real thickness ratios summed to **0.20 – 0.75** instead of 1
- **no error was reported.** `check_para()` would have caught it but is dead
  code, never called from anywhere in the codebase.

### Files

| file | change |
|---|---|
| `Codes/MCMC_flex/INITstructure.h` | `groupdef.value1vpvs`, `paradef.npara_vpvs` |
| `Codes/MCMC_flex/CALpara.C` | `readpara` counts type `-1`; `gen_newpara` rewritten to four contiguous blocks `[velocities \| total thickness \| ratios \| Vp/Vs]` with an `idx_vpvs0+n_vpvs==npara` assert; `mod2para`/`para2mod` read and write `value1vpvs[sub]`; `ibad` initialised in the gaussian loop (was used uninitialised, so `ibad%5000` could trip on garbage and `exit(1)`) |
| `Codes/MCMC_flex/CALgroup.C` | `get_vp` flag 5; `updategroup1` takes the per-sub-layer override |
| `Codes/MCMC_flex/CALmodel.C` | `readmod` makes the expected column count `p_flag`-aware and loads the extra Vp/Vs columns; `write_model` emits them |
| `SetupData/parameters.py` | `sublay_vpvschange`, `sub{sed,crust,mantle}vpvschange`, `uniform_vpvs`, `vpvs_{sed,crust,mantle}`, and the `PertType/Style/Range/gwStep *VpVs` columns |
| `SetupData/src/setup_layercake.py` | sets `p_flag=5` automatically, appends one Vp/Vs per sub-layer to `mod.{sta}`, writes the `-1` rows **last** in `in.para_{sta}` |
| `AnalyzeResult/inversion_plot_vfinal_flex.py` | reads `MC.{sta}.all.value_vp`, draws the Vp posterior cloud and an ensemble Vp/Vs curve with propagated error bars, writes `{sta}_Vp_average.txt` |

### Three deliberate deviations from the reference

1. `updategroup1` does **not** `push_back` into `value1vpvs`. The reference
   appends `nnlay` entries on every call while only `[0,np)` is ever read —
   results identical, without the growth.
2. The `mod.{sta}` Vp/Vs columns are written **only when that group's `p_flag`
   is 5**. The reference writes them unconditionally, which would make `readmod`
   reject the line with the feature off; it never hit this because it is always
   on there.
3. `write_model`'s ratio loop was fixed. It read `ratio[np]` out of bounds after
   the value loop, so `.mod.group` files carried one bogus ratio (`0.000000`)
   instead of all of them.

**Not ported:** `setup_bspline.py`'s Vp/Vs section. Only `updategroup1`
(layer-cake) applies per-sub-layer Vp/Vs in the forward — true in the reference
too — so bspline would write columns the forward ignores. This project is
`mod_type_flag=1`.

### Verified

On a station with Vph + HV + RF + higher mode: `readmod` loads `value1vpvs` for
both groups, `readpara` reports `npara vpvs: 6`, thickness ratios still sum to
`1.000000`, all six Vp/Vs stay inside their priors, and the sampled values
appear exactly as the Vp/Vs of the layer blocks in the forward model
(`2.340 / 1.797 / 2.084 / 1.485 / 1.048 / 1.095`). Higher-mode outputs
(`all.hp`, `mean.hp`, `*.hp.disp`) still produced. Regression: an old
`p_flag=3` station still runs, reporting `npara vpvs: 0` with 12 parameters.

---

## 2. `f1d5cbf` — FM project skeleton

**90 files, +2425.** A second project directory beside `CHT`, cloned from it.

- `FM/station_cor.lst`, `FM/all_station_cor.lst`, `FM/sta_now` — one station
- `FM/data/` — one empty station folder, ready for data
- `FM/MonteCarlo/00_make_file_control.sh` + `FM/MonteCarlo/indata/`
  (`Qmodel.dat`, `in.ph`/`in.gv`/`in.hv`/`in.rf` fallbacks)
- `FM/Vel_mod/` — copied whole from CHT
- `FM/query_data/`, `FM/real_data/` — empty, so git does not track them

Verified end to end before being emptied again: SetupData → control file →
`do_MC_Para` → post-process → figure, with Vph + HV + RF.

---

## 3. `20df7b0` — CHT station inputs regenerated **(PARTIAL)**

**51 files, +136 / −36.**

A SetupData run on CHT with `sublay_vpvschange=1` **stopped partway**. Only
stations **00620–00780 (17 of 135)** were regenerated; the other 118 still have
the old `p_flag=3` models and no Vp/Vs rows. Committed as-is, with the partial
state stated in the commit message rather than quietly baked in.

```
mod:   ... <5 Vs> <5 ratios> 2.050000 1.728000 1.727000 1.730000 1.740000  1 3 5 1.70 1.73 ...
para:  -1 -1 40 0.05 0 0   ...   -1 -1 40 0.05 1 0
```

Each station is internally consistent (its `mod` and `in.para` were written
together), so both kinds run correctly — but the array is currently a mix: 17
stations would invert Vp/Vs and 118 would not.

Also in this commit: `in.data_00740` lost its higher-mode lines 10–11, because
SetupData writes only the 9 data/weight lines. Harmless —
`00_make_file_control.sh` re-appends them whenever `{sta}.hph` exists.
`query_data/*/StartingModel_*.png` are the regenerated starting-model figures
(16, not 17: 00780's had not been written when the run stopped).

**Open item.** Either finish the conversion (`conda activate MCMC && bash
01_prepare_data.sh CHT`) and commit the remaining 118, or `git revert 20df7b0`
to put CHT back. The code port in `85e5502` is independent either way.

---

## 4. `75070c2` — Higher-mode master switch

**5 files, +74 / −14.**

Higher mode used to be driven purely by the existence of a non-empty
`{sta}.hph`, with its weight hard-coded to `1.0` inside
`00_make_file_control.sh`. There was no way to turn it off except deleting the
data file — editing `in.data` lines 10–11 did not work, because the script
overwrote them on every run.

`in.connector` gains two values:

| position | name | meaning |
|---|---|---|
| 11 | `use_higher_mode` | `1` = use `.hph` when present, `0` = **force off everywhere** |
| 12 | `hpw` | weight of the higher mode in the misfit |

| file | change |
|---|---|
| `SetupData/parameters.py` | `use_higher_mode`, `hpw` |
| `SetupData/src/setup_layercake.py` | writes the two values |
| `SetupData/src/setup_bspline.py` | writes them, **and restores its missing line 10** (`MC_QC_thres`) — a bspline-generated `in.connector` was one line short, so every later position was off by one. `setup_layercake.py` always wrote it. |
| `CHT/MonteCarlo/00_make_file_control.sh` | reads them, gates the higher-mode branch |
| `FM/MonteCarlo/00_make_file_control.sh` | same (separate copy of the same script) |

When off, the script writes `in.data` lines 10–11 explicitly as `0` and `0.0`
rather than dropping them: the off state is visible, and the
`is_equal_weight=1` check in `CALmodel.C` (`phw+gvw+hvw+rfw+hpw == 1`) stays
satisfiable — the hard-coded `1.0` used to push that sum to 2.0 and abort with
`Total weighting less than 1, stop!`.

The `sedmonocheck`/`crustmonocheck` reads moved from positions 11/12 to 13/14.
They are dead: never used anywhere in the script, never written by SetupData.

### Verified on a station that **has** a `.hph` file

| switch | disp line | `in.data` 10/11 | result |
|---|---|---|---|
| `1` | `... 3 HV 5 hph` | `1  1.0` | 8 higher-mode outputs, `.mean.hp` from 45 models |
| `0` | `... 3 HV` | `0  0.0` | `No higher-mode data in ensemble; .mean.hp not written` |
| legacy 10-line file | `... 3 HV 5 hph` | — | defaults to ON, exactly as before |

Both states run `do_MC_Para` → post-process → figure with no error. The plotting
script is safe in the forced-off case because it gates on the *content* of
`MC.{sta}.all.hp`, not merely on the file existing.

---

## 5. `323e13c` — FM made runnable

**21 files, +412 / −7.**

### Drivers parameterised, default switched to FM

`01_prepare_data.sh`, `02_do_MCMC.sh`, `04_re_do_MCMC_post_process.sh`,
`05_do_replot_MCMC.sh` took the project from a hard-coded `CHT` in four places.
It now comes from `$1`, defaulting to **FM**. So `bash 02_do_MCMC.sh` runs FM
and `bash 02_do_MCMC.sh CHT` goes back.

### Per-project configuration — the real blocker

`SetupData/parameters.py` is shared by every project, but CHT and FM need
different settings, and editing it for FM would have silently changed CHT.
`SetupData/main.py` now loads `SetupData/parameters_{project}.py` when it exists
and overrides only the names defined there, **in memory**. CHT has no such file,
so it is completely unaffected.

New `SetupData/parameters_FM.py`:

| setting | CHT | FM | why |
|---|---|---|---|
| `modfile` | `NTW_1D_Liu21_modified.dat` | `NTW1d_H14_intp_50km` | FM data are crustal-scale (Vph 3–40 s, RF to 7.5 s) and need the 50 km starting model |
| `ddc` | 0.05 | 0.5 | 0.05 km over a ~30 km crust = 600 sub-layers, past the forward model's limit |
| `phw/gvw/hvw/rfw` | 1.0 / 0 / 0 / 0 | 0.335 / 0 / 0.33 / 0.335 | the shared file would read and plot HV/RF but give them **zero weight** in the misfit |
| `RF_gaussian_width` | 2.5 | 3.0 | must match the gaussian the observed RFs were computed with |
| `MC_number_of_cores` | 5 | 1 | the Fortran forward is called without the OpenMP critical section, so >1 thread can corrupt or crash on macOS |

### Station `121.475-24.9`

`.ph` (41 periods), `.HV` (49 periods) and `.RF` (151 samples, dt 0.05) copied
from the reference project. **Only the data** — `mod.`/`in.para`/`in.data`/
`in.connector` are SetupData's output and were regenerated here. They come out
identical to the reference's, including the per-sub-layer Vp/Vs columns and
`p_flag 5`. `FM/Vel_mod/NTW1d_H14_intp_50km{,.dat}` copied in. Station lists
updated; the `00740` placeholder folder removed.

Smoke-tested in a scratch copy so `FM/MonteCarlo` stayed clean:

```
disp check n-phase, n-group, p-ellip: 41 0 49
RF!!!! 20 151
[readmod] group=0 p_flag=5 read value1vpvs size=5
readpara > npara: 18 | velocity: 6 | total thicness: 1 | ratio: 5 | vpvs: 6
do_MC_Para done!!!
```

---

## 6. `6c6db0e` — `in.connector` comments

**9 files, +226 / −105.**

`in.connector` had grown to 12 bare numbers with nothing to say which was which,
so editing it by hand meant counting lines. Every reader now skips lines whose
first non-blank character is `#`, plus blank lines, and addresses values by
their **position among the data lines** instead of by physical line number.

| file | change |
|---|---|
| `CHT/MonteCarlo/00_make_file_control.sh`, `FM/MonteCarlo/00_make_file_control.sh` | the twelve `awk 'NR==n'` calls replaced by one `conn_get` helper, keeping each value's `printf` format |
| `AnalyzeResult/inversion_plot_vfinal_flex.py` | new `read_connector_values()`; `parse_plot_type()` counts data lines |
| `AnalyzeResult/inversion_plot_vfinal_flex_only_final.py`, `AnalyzeResult/inversion_plot_forward.py` | same, replacing their inline `if i == 9` loops |
| `SetupData/src/setup_layercake.py`, `SetupData/src/setup_bspline.py` | build the file from one list of `(label, value)` pairs and emit `# N. label` above each number, so the two cannot drift apart again |

Resulting file:

```
# in.connector - generated by SetupData, edit the numbers only.
# Lines starting with '#' and blank lines are ignored by every reader.
# The values are read IN ORDER, so do not reorder or delete any of them.
# 1. inversion type: 1 = layer-cake, 2 = B-spline
1
# 2. gaussian width the OBSERVED receiver functions were computed with
3.0
...
```

### Verified

- The control file generated from the new commented `in.connector` is
  **byte-identical** to the one generated from the old plain file.
- `read_connector_values()` returns the 12 values in order; `parse_plot_type()`
  returns 1.
- Three shapes parse correctly: a legacy plain 10-line file (higher mode
  defaults ON), a commented file with the switch set to 0 (forced off), and a
  commented file with extra hand-written notes and blank lines inserted
  mid-file.
- Full chain on the last of those: `do_MC_Para` → post-process (2/2) → figure.

Backwards compatible: a file with no comments is filtered to itself.

---

## 7. `<next commit>` — Plot: Vp model, and a depth range that follows the model

**1 file** (`AnalyzeResult/inversion_plot_vfinal_flex.py`).

### Vp model drawn on the figure

The Vp/Vs curve was already there from `85e5502`; the **Vp model itself** is now
drawn too, as the reference project does — cyan `Final Vp` with error bars, on
both model panels, from the ensemble mean of `MC.{sta}.all.value_vp`. The
velocity axis widens automatically when Vp is present, because mantle Vp reaches
~7.7 km/s and the default `VS_XLIM` stopped at 7.0, which was silently clipping
the deep half of the curve off the panel.

### Depth range taken from the model

`VS_YLIM_FULL = (0, 10)` and `VS_YLIM_SHALLOW = (0, 2)` were hard-coded, so the
50 km FM model was being drawn on a 10 km axis — 80 % of it off the panel.

`read_model_total_thickness()` sums the group thicknesses in `mod.{sta}`, and the
panels follow:

| model | full panel | zoom panel (`ZOOM_FRACTION = 0.10`) |
|---|---|---|
| FM `121.475-24.9` — 30 + 20 | 0–50 km, 5 km ticks | 0–5 km, 0.5 km ticks |
| CHT `00740` — 5 + 10 | 0–15 km, 2 km ticks | 0–1.5 km, 0.2 km ticks |

Ticks come from `depth_ticks()`, which snaps to a 1 / 2 / 2.5 / 5 × 10^k step, so
any model depth gets readable labels. The old constants remain as fallbacks if
`mod.{sta}` cannot be read.

### Bug found and fixed while doing this

The ensemble Vp/Vs profile from `85e5502` averaged the posterior **node by
node**, taking `nnode = min(len(model))` across models. But posterior models do
not share a depth grid — the layer thicknesses are inverted for, so this
station's 56 models carry **116 to 194 nodes each** even though every one of them
reaches 50 km. The profile was therefore being truncated to the shortest model:
**28.85 km instead of 50 km**, with no warning. Every model is now interpolated
onto one common 400-point grid before averaging.

Markers and error bars are thinned to ~25 per panel (`decimate_profile()`); the
curve itself is still drawn at full resolution. Thinning is done per panel, so
the zoom keeps its own markers. Without this, 400 markers per curve rendered as
a solid blob.

### Verified

- FM (50 km, `p_flag=5`): `Model depth 50 km -> full panel 0-50 km, zoom panel
  0-5 km`, `Vp/Vs profile written ... (400 points, 0-50 km)`.
- CHT-style (15 km, `p_flag=3`, no Vp/Vs perturbation): `0-15 km / 0-1.5 km`.
  The Vp curve is still drawn — `value_vp` is always written by the
  post-process, and there Vp is just Vs × 1.70, so the Vp/Vs curve reads flat at
  ~1.7 as it should.
- Station with no `MC.{sta}.all.value_vp` at all: plots fine, no Vp curve.
- Unreadable `mod.{sta}`: falls back to the old fixed limits with a warning.

## 8. `<next commit>` — Posterior cloud for Vp/Vs

**1 file** (`AnalyzeResult/inversion_plot_vfinal_flex.py`).

Vs and Vp already had a posterior cloud behind their mean curve; Vp/Vs only had
the mean ± sigma error bars. The per-model Vp/Vs profiles are now drawn as a
cloud too, so the actual spread of the ratio is visible rather than just a
symmetric error bar.

The curves are the ratio of the two already-interpolated profiles, model by
model, on the same common grid (`vp_arr / vs_arr`, guarded against Vs = 0), so
they cost nothing extra to compute.

Colour: `mediumorchid` at alpha 0.14. It needs its own colour because Vp/Vs
(~1.5–2.3) overlaps the Vs range on the same axis, so a grey cloud there would
read as Vs. The three cloud colours are now constants —
`POSTERIOR_VS_COLOR`/`_ALPHA`, `POSTERIOR_VP_COLOR`/`_ALPHA`,
`POSTERIOR_VPVS_COLOR`/`_ALPHA` — and `SHOW_VPVS_POSTERIOR = False` turns the
new one off.

Verified: FM 50 km and a CHT-style 15 km station both render it on the full and
the zoom panel; a station with no `all.value_vp` skips it cleanly; the
`SHOW_VPVS_POSTERIOR = False` path runs.

## 9. `9933217` — First FM inversion result

**86 files, +267314.** The run itself, committed as a result rather than as code.

```
tt 12, jt 10000, 1 thread, 50 km model
56 posterior models, 18 parameters each
minimum misfit 2.182663
post process finish checkpoint matched
```

Vph + H/V + RF jointly, per-sub-layer Vp/Vs inverted for, higher mode off (this
station has no `.hph`). It confirms the Vp/Vs port end to end on real data: the
five thickness ratios still sum to `1.000000` (0.999999–1.000001) and the crust
Vp/Vs varies over **1.656–2.248** without being dragged into the ratio simplex.

Figures were regenerated after sections 7 and 8, so `121.475-24.9_MCMC.png`
carries the Vp model, the Vp/Vs mean and both posterior clouds, on a 0–50 km
full panel and a 0–5 km zoom.

**Size:** 32 MB, of which 26 MB is `MC.121.475-24.9.out`, a log rewritten on
every run. Committed because the repo already tracks six such logs for CHT, so
excluding this one would be an inconsistent convention rather than a decision.
If repository size becomes a problem the fix is `MC.*.out` in `.gitignore` plus
`git rm --cached` on all seven — that applies to the CHT logs too.

---

## Files changed on disk but **not** in git

| path | note |
|---|---|
| `Codes/MCMC_flex/do_MC_Para`, `do_MC_Para_Post_process_v3`, `*.o` | rebuilt (last at 23:38). Git-ignored on purpose — they are platform-specific. |
| `Codes/MCMC_flex/macbuild/{fast_surf,calcul}.f` | generated by the Makefile's macOS COMMON-block normalisation. Git-ignored. |
| `SetupData/**/__pycache__/*.pyc` | tracked from before `.gitignore` existed, rewritten by every SetupData run. Deliberately kept out of all commits. `git rm --cached -r SetupData/__pycache__ SetupData/src/__pycache__` would stop the churn. |
| ~~`FM/MonteCarlo/...` run outputs~~ | **now committed** (`9933217`), see section 9. |
| `FM/data/121.475-24.9_data/.in.connector.swp` | vim swap file, untracked. Consider adding `*.swp` to `.gitignore`. |
| `~/.claude/projects/.../memory/` | assistant memory notes: `mcmc-mac-build-and-env.md` (updated), `mcmc-vpvs-perturbation.md` (new), `MEMORY.md` |
| scratchpad under `/private/tmp/claude-501/...` | all verification runs, so the project tree stayed clean. Temporary. |

---

## Findings recorded but **not** fixed

| finding | where | status |
|---|---|---|
| `check_para()` is dead code — never called anywhere, though it would have caught exactly the out-of-bounds parameter bug | `Codes/MCMC_flex/CALpara.C:265` | left as-is |
| Vp/Vs prior is percentage-style ±40 % with **no physical floor** (`mod2para` skips the clamp for type `-1`), so bounds reach ~1.04 and a test run sampled 1.045. Anything below √2 ≈ 1.414 is a **negative Poisson's ratio**. `goodmodel` checks Vp monotonicity but never the ratio. | `Codes/MCMC_flex/CALpara.C` | warning block + the absolute-style fix documented beside the setting in `parameters.py` |
| `readQmodel` reports a missing Qmodel file as `#########rf file ... does not exist!` — copy-pasted string | `Codes/MCMC_flex/CALmodel.C:104` | left as-is |
| `get_ave_para` return value ignored at the call site | `Codes/MCMC_flex/do_MC_Para_Post_process_v3.C:1258` | left as-is |
| SetupData crashes with `KeyError: 'inparacrustVs0'` on the default python 3.13 — it writes into `locals()`, which PEP 667 made a snapshot. **Must run under the `MCMC` env (python 3.7.6).** `INSTALL.md`'s claim that the python code runs on the modern stack is wrong for SetupData. | `SetupData/src/setup_*.py` | documented, not rewritten |
| Post-analysis utilities still hard-coded to CHT: `96/97/98/99_*.py`, `00_interpolate_model_1D.py`, `000_copy_data_phasevel.py`. `03_re_prepare_data.sh` has a pre-existing `project_dir="CHTH3"`. | project root | left on purpose — they are not part of the run chain and changing them would move where figures are copied |

### About the reference project

`~/Research/MCMC/20260510_FM_vpvs_pertubation_huang14_vpvs_model` was read only.
Two problems were found in it while verifying, and **fixed only in a scratch
copy**, never in place:

- `FM/MonteCarlo/indata/` is empty — `Qmodel.dat` is missing, so `do_MC_Para`
  dies on startup.
- Its Makefile has no `macbuild/` COMMON-block normalisation and no
  `%.o : %.f` rule, so on macOS its post-process dies with `Bus error: 10` in
  `init_` ← `fast_surf_`. Confirmed by decoding the crash report, and confirmed
  fixed by `sed 's/nper=200,/nper=2000,/' DISP2/fast_surf.f` +
  `sed 's/ndep=20,/ndep=100,/' DISP2/calcul.f`, recompiling the eight Fortran
  objects by hand, then `make`.
- Its stored results are dated Jun 13 while its sources and binaries are Aug 4,
  so those outputs were not produced by its current code.

---

## Open items

1. **CHT is half-converted** — 17 of 135 stations in the new `p_flag=5` format.
   Finish with `bash 01_prepare_data.sh CHT` and commit, or `git revert 20df7b0`.
2. **Nothing is pushed.** `git push origin main fix/macos-portability`.
3. `main` is behind the working branch again — it was fast-forwarded once, at
   `6c6db0e`, and three commits have landed since.
