# Install / System Requirements — MCMC-CHT

How to set up the environment and build the code on a fresh machine. See
[CHECKPOINT.md](CHECKPOINT.md) for the analysis state and applied fixes.

---

## 1. System requirements

| Component | Documented (original) | Verified working here | Notes |
|---|---|---|---|
| OS | Linux (RHEL) | Ubuntu 24.04 (Linux 6.17) | macOS untested for the *build*; fine for Python-only analysis |
| Python | 3.7.6 | 3.7.6 | scripts call `python` (not `python3`) |
| numpy | 1.18.1 | 1.18.1 | `np.float` used in setup → needs numpy < 1.24 |
| pandas | 1.0.1 | 1.0.1 | `read_csv(delim_whitespace=...)` → needs pandas < 3.0 |
| matplotlib | (any) | 3.4.3 | plotting only |
| C++ compiler | g++ 4.8.5 | system g++ 13.3 | needs `-std=c++0x`; see fixes in CHECKPOINT.md |
| Fortran compiler | gfortran 4.8.5 | conda gfortran 15.2 | builds DISP2/ + RF/ objects |
| Boost | 1.55 | 1.55 (bundled) | headers vendored in `Codes/boost-1.55.0-gcc630/` — no separate install |

The pinned Python stack matters: the setup code uses APIs removed in modern
numpy/pandas (`np.float`, `delim_whitespace`). Do **not** substitute newer
versions without patching the code.

---

## 2. Create the Python environment

### Option A — conda (recommended)
```bash
conda env create -f environment.yml     # portable spec
conda activate MCMC
```
For a byte-exact Linux reproduction instead, use the lock file:
```bash
conda env create -f environment.linux-lock.yml
```

### Option B — pip (inside an existing Python 3.7 env)
```bash
pip install -r requirements.txt          # numpy/pandas/matplotlib only
# gfortran must still be installed separately (conda/apt/brew)
```

> **IMPORTANT:** the driver scripts call `python`, which only resolves inside
> the activated env. Always `conda activate MCMC` before running `01_*.sh`,
> `02_*.sh`, etc. (or prepend the env bin to PATH).

---

## 3. Build the C++/Fortran engine

```bash
conda activate MCMC
cd Codes/MCMC_flex
make clean
make do_MC_Para
make do_MC_Para_Post_process_v3
```

The `Makefile` links the conda `libgfortran` via:
```
LDLIBS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
```
so the binaries find `libgfortran.so.5` at runtime without `LD_LIBRARY_PATH`.
This only works while the `MCMC` env is active (uses `$CONDA_PREFIX`).

**If building on a machine with a matching-era gcc (≈4.8–7)**, none of the
gcc-13 fixes are needed. On a modern gcc you need the fixes already applied in
this repo (missing `return` statements, `#include <random>` in
`head_c++/gen_random.C`). See CHECKPOINT.md → "Fixes applied".

---

## 4. macOS caveats (for Claude-on-Mac / analysis)

- **Python-only analysis** (reading data, plotting, running `SetupData`) works on
  macOS with the same conda env.
- **Apple Silicon (arm64):** these *old* pinned packages (Python 3.7.6, numpy
  1.18.1, pandas 1.0.1) may have **no arm64 conda builds**. Force the Intel
  channel under Rosetta:
  ```bash
  CONDA_SUBDIR=osx-64 conda env create -f environment.yml
  conda activate MCMC
  conda config --env --set subdir osx-64
  ```
- **Building the C++/Fortran on macOS is untested.** The Makefile assumes GNU
  tooling and a conda `libgfortran` layout; on macOS you'd use clang++ + a
  gfortran from conda/brew and likely adjust `LDLIBS`/rpath. Recommended:
  keep the build + inversion runs on the Linux box; use macOS for analysis.

---

## 5. Quick verification

```bash
conda activate MCMC
python -c "import sys,numpy,pandas,matplotlib as m; \
  print(sys.version.split()[0], numpy.__version__, pandas.__version__, m.__version__)"
# expect: 3.7.6 1.18.1 1.0.1 3.4.3

cd Codes/MCMC_flex && ./do_MC_Para
# expect usage message: "input [control.file]"  (means the binary runs)
```

## Files in this directory
- `environment.yml` — portable conda spec (recreate the `MCMC` env)
- `environment.linux-lock.yml` — exact build-pinned Linux export
- `requirements.txt` — pip fallback for the Python libraries
- `CHECKPOINT.md` — analysis state, fixes, higher-mode comparison
