# Install / System Requirements — MCMC-CHT

How to set up the environment and build the code on a fresh machine. See
[CHECKPOINT.md](CHECKPOINT.md) for the analysis state and applied fixes.

---

## 1. System requirements

| Component | Documented (original) | Verified working here | Notes |
|---|---|---|---|
| OS | Linux (RHEL) | Ubuntu 24.04 (Linux 6.17) + macOS 15 (arm64) | macOS build **verified working** — see §4 |
| Python | 3.7.6 | 3.7.6 and 3.13.12 | scripts call `python` (not `python3`) — must be on PATH |
| numpy | 1.18.1 | 1.18.1 and 2.4.6 | any version — `np.float` removed from the code |
| pandas | 1.0.1 | 1.0.1 and 3.0.3 | any version — `delim_whitespace` removed from the code |
| matplotlib | (any) | 3.4.3 and 3.10.8 | plotting only |
| C++ compiler | g++ 4.8.5 | g++ 13.3 (Linux), g++-15 (macOS/brew) | needs `-std=c++0x`; see fixes in CHECKPOINT.md |
| Fortran compiler | gfortran 4.8.5 | conda gfortran 15.2 (Linux), gfortran-15 (macOS/brew) | builds DISP2/ + RF/ objects |
| Boost | 1.55 | 1.55 (bundled) | headers vendored in `Codes/boost-1.55.0-gcc630/` — no separate install |

**The Python code now runs on both the old and the modern stack.** It
previously used APIs removed in newer releases; these have been replaced with
spellings that are valid in *every* version, so the pinned env is no longer
mandatory (it is still the reference for byte-exact reproduction):

| Was | Now | Why |
|---|---|---|
| `np.float(x)` | `float(x)` | `np.float` was only an alias for the builtin; removed in numpy 1.24. NumPy's own advice is to use `float`. |
| `delim_whitespace=True` | `sep=r"\s+"` | Removed in pandas 3.0; `sep=r"\s+"` is the documented replacement and works back to pandas 0.x. |

Verified running on macOS with numpy 2.4 / pandas 3.0 / matplotlib 3.10 as well
as on the pinned 1.18/1.0.1 stack.

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
- **Building the C++/Fortran on macOS now works** (verified on Apple Silicon:
  builds clean and runs a full inversion). Requirements:
  ```bash
  brew install gcc          # provides g++-15 / gfortran-15
  ```
  Apple's `g++` is clang, which has **no OpenMP and no libgfortran**, so the
  `Makefile` auto-switches to the Homebrew GNU toolchain on Darwin and sets an
  rpath to `$(brew --prefix gcc)/lib/gcc/current` — the binaries then run with
  no `DYLD_LIBRARY_PATH` and no conda env. Override the compiler version with
  `make CC=g++-14 FC=gfortran-14` if yours differs.
- **macOS COMMON-block fix (automatic).** `DISP2/fast_surf.f` declares
  `nper=200` and `calcul.f` declares `ndep=20`, while `init.f`/`surfa.f` declare
  `nper=2000`/`ndep=100` for the *same* COMMON blocks. GNU ld sizes a common
  block to the **largest** definition (harmless on Linux); the macOS Mach-O
  linker takes the **smallest**, so writes land outside the block and the run
  dies with `Bus error: 10`. The Makefile generates normalised copies under
  `macbuild/` and compiles those. **The repo sources are left untouched**, so
  the CentOS build uses exactly the same files as before.

---

## 4b. Moving the code between macOS and CentOS (avoiding conflicts)

**The rule: ship source, never compiled objects.**

The `.o` objects and the `do_MC_Para*` binaries are architecture-specific —
CentOS builds ELF x86-64, macOS builds Mach-O arm64. They used to be
**committed to git**, which caused two failures when hopping between machines:

1. Every switch rewrote them, and git cannot merge binaries → conflicts on
   every pull.
2. A CentOS checkout could receive macOS objects. Because git sets checkout
   timestamps, `make` often considers them newer than the `.f` sources and
   **skips rebuilding**, so the link fails with `unknown file type` /
   incompatible-format errors.

They are now **untracked and git-ignored**, and are rebuilt on each machine from
the tracked `.f`/`.C` sources. Nothing platform-specific travels with the repo,
so a plain `git pull` on CentOS is safe.

If you move the folder by **copy/rsync instead of git**, clear the artifacts
first so the other machine cannot reuse them:

```bash
cd Codes/MCMC_flex && make clean && make clean-fortran
```

`make clean-fortran` also removes the generated `macbuild/` directory. On the
destination machine just rebuild (§3); `02_do_MCMC.sh` also offers to recompile
when you answer `y` at its prompt.

> Note: the Linux branch of the Makefile links against `$(CONDA_PREFIX)/lib` and
> falls back to `$HOME/anaconda3/envs/MCMC`. On CentOS, `conda activate MCMC`
> before building, or pass the right path:
> `make CONDA_PREFIX=/path/to/envs/MCMC do_MC_Para`.

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
