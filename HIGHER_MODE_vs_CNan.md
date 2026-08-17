# Higher-mode surface wave: where H14 was modified vs CNan reference

**Goal:** forward-model and fit the **fundamental + 1st higher mode** Rayleigh phase
velocity, mirroring CNan's proven higher-mode code, adapted to the H14 code.

- **CNan (reference):** `~/Research/MCMC_CHT/MCMC_joint_inversion_CNan/Codes_ori/MCMC_PhOnly_Higher_Mode_ori/`
- **H14 (ours):**       `Codes/MCMC_flex/`

The higher mode is enabled by 5 coordinated edits, in the **same places** CNan edited,
plus a few adaptations forced by differences between the two codebases.

---

## 1. Map of the 5 change sites (same locations as CNan)

| # | Layer | File | CNan does | H14 does (this repo) | Same place? |
|---|-------|------|-----------|----------------------|-------------|
| A1 | engine: #modes | `DISP2/init.f` | `mode=1` → `mode=2` | keep `mode=1` default; overridden per call (see A3) | ✅ same file |
| A2 | engine: extract mode 2 | `DISP2/fast_surf.f` | `FAST_SURF` +6 out args `uR1/uL1/cR1/cL1/rR1/rL1`; `cR1(i)=cR(i,2)` etc; clear mode-2 arrays | **identical** | ✅ same |
| B  | forward wrapper | `CALforward.C` | extern proto +6 args; alloc `cR1/uR1/rR1…`; period list = pper+gper+eper+hpper; store `cR1`→`disp.hpvel` when `nhpper>0` | same, but keeps H14 `invtype` arg and H14 `disp.period1` (set_union) then **appends `hpper`** | ✅ same |
| C  | data struct | `INITstructure.h` | `dispdef` +`nhpper/hpper/hpvelo/hpvel/unhpvelo` (+`fhphase/phmisfit/hpL`) | same `dispdef` fields; **plus** `indatadef` +`hpflag/hpw` | ✅ same |
| D1 | read `.hph` | `CALmodel.C readdisp1()` | `tflag==2` → higher phase | `tflag==`**`5`** → higher phase | ✅ same, **different flag #** |
| D2 | misfit | `CALmodel.C compute_misfit()` | add `tempv1hmode`, `phmisfit`, fold into equal-weight joint disp sum | add `tempv1h`, `phmisfit`, `hpL`, `tmfhp`; fold into H14 **weighted** misfit (`hpw`,`ihp`) | ✅ same |
| A3 | engine: flexible #modes | `fast_surf.f` + `CALforward.C` | — (always mode=2) | **H14 addition:** `FAST_SURF` gains `nmode` arg; `mode=nmode`; caller passes `(nhpper>0)?2:1` | H14 extra |

---

## 2. What each edit looks like

### A1/A3 — `DISP2/init.f`
```fortran
        mode=1   ! default; FAST_SURF overrides via nmode (mode=2 if .hph)
```
CNan hard-codes `mode=2`. H14 keeps 1 and lets the caller choose per station.

### A2/A3 — `DISP2/fast_surf.f`
```fortran
        subroutine FAST_SURF(... uR0,uL0,cR0,cL0,rR0,rL0,
     &        uR1,uL1,cR1,cL1,rR1,rL1,nmode)     ! +6 mode-2 outputs, +nmode
        ...
        mode=nmode                 ! H14: pick #modes per call (clamp 1..2)
        if (mode.gt.2) mode=2
        if (mode.lt.1) mode=1
        ...
        cR1(i)=cR(i,2)  uR1(i)=uR(i,2)  rR1(i)=rR(i,2)   ! extract 2nd mode
```

### B — `CALforward.C compute_disp()`
```cpp
// extern proto extended with the 6 *1 arrays + int *nmode
int nmode = (model.data.disp.nhpper>0) ? 2 : 1;             // H14: per-station
// period1 = disp.period1 (union of pper/gper/eper) then append hpper
fast_surf_(..., uR1,uL1,cR1,cL1,rR1,rL1, &nmode);
if (model.data.disp.nhpper>0) {                             // store higher mode
  for (...) model.data.disp.hpvel.push_back(cR1[...]);     // cR1 -> disp.hpvel
}
```

### C — `INITstructure.h`
```cpp
struct dispdef  { ... int nhpper,fhphase; double phmisfit,hpL;
                  vector<double> hpper,hpvelo,hpvel,unhpvelo; };
struct indatadef{ ... int hpflag; double hpw; };   // H14 weighted-misfit needs these
```

### D1 — `CALmodel.C readdisp1()`
```cpp
else if (tflag ==5 ) {         // H14 flag 5 (CNan uses 2). KEEPS 2=group,3=ellip.
    model.data.disp.hpper=cv1; hpvelo=cv2; unhpvelo=cv3; nhpper=i; fhphase=1; }
```

### D2 — `CALmodel.C compute_misfit()` + `readindata()`
```cpp
tempv1h = Σ((hpvelo-hpvel)/unhpvelo)^2 ;  phmisfit=√(tempv1h/nhpper);  hpL=exp(-…)
tmfhp   = √((tempv1h*ihp)/nhpper);  Sph=tempv1h;
// folded (with weight hpw, flag ihp) into disp.misfit, model.data.misfit, tS, tL1,
// and EVERY normalisation denominator + the weight-sum check.
// readindata parses OPTIONAL in.data line 10 (hpflag) & 11 (hpw)  -> back-compatible
```

---

## 3. Three deliberate differences from CNan (not a blind copy)

1. **Flag number 5, not 2.** In H14 the `disp` control line already uses `2 = group,
   3 = ellipticity`. CNan moved group to 4 to free up 2 for higher phase. We instead
   gave higher phase a **new flag 5**, leaving group/ellipticity untouched. Must stay
   consistent across `readdisp1`, `00_make_file_control.sh`, and SetupData.

2. **Weighted misfit, not equal-weight sum.** H14 uses per-data weights
   (`phw/gvw/hvw/rfw`); we added a matching `hpw` weight and `hpflag` availability flag
   (via `in.data` lines 10–11) instead of CNan's plain `(…+tempv1hmode)/(…+nhpper)` sum.

3. **Higher mode is computed only where needed.** `nmode=(nhpper>0)?2:1` — stations
   without `.hph` run fundamental-only (faster / no needless higher-mode root search).
   CNan always computes 2 modes.

---

## 4. Verification

- Both binaries (`do_MC_Para`, `do_MC_Para_Post_process_v3`) build clean.
- **Backward compatible:** with no `.hph` (`hpflag=0, hpw=0, nhpper=0`) every higher-mode
  term is zero → existing fundamental-only runs are bit-identical.
- **HV/ellipticity unaffected:** `nmode=1` vs `nmode=2` on the same model gives identical
  fundamental phase vel, **ellipticity (`rR0`)**, and group vel to ~3×10⁻⁷ (float32
  rounding), ~5 orders below data uncertainty. HV fits the fundamental `rR(i,1)`, which
  `calcul.f` stores in column 1 independently of the mode-2 search.

## 5. Still to wire (to actually activate — currently `nhpper` stays 0)
- **F** `00_make_file_control.sh` add `5 {sta}.hph` to the disp line; SetupData write
  `in.data` lines 10–11 (`hpflag`,`hpw`) for stations that have a `.hph`.
- **E** post-process output of the higher-mode curve (CNan: `MC.{sta}.phmode.disp`).
- **G** plotting.
