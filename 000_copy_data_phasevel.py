#!/usr/bin/env python3
"""
000_copy_data_phasevel.py

Copy new phase-velocity dispersion data into the per-station inversion input dirs.

Source (per station), in ~/Research/MCMC_CHT/Vel_mod_data/ :
    {station}.ph     fundamental-mode Rayleigh phase velocity   (period  vel  unc)
    {station}.hph    1st higher-mode  Rayleigh phase velocity   (period  vel  unc)  [optional]

Destination (per station):
    CHT/data/{station}_data/{station}.ph     <- REPLACED with the new file
    CHT/data/{station}_data/{station}.hph    <- copied when the source has one

Notes
-----
- The run is driven by the SOURCE files. Every {station}.ph (and {station}.hph)
  found in SRC_DIR is a candidate.
- Only stations that ALREADY have a CHT/data/{station}_data/ directory are
  updated. Creating that directory + the other control files (in.data, mod.*,
  *.control, in.connector, ...) is SetupData's job, not this script's. Flip
  CREATE_MISSING_DIRS=True if you really want this script to make new dirs.
- Not every station has a .hph; those stay fundamental-only. That is expected.
"""

from pathlib import Path
import shutil
import os
import re


# =======================
# Parameters
# =======================
PWD = os.getcwd()
SRC_DIR = Path("~/Research/MCMC_CHT/Vel_mod_data").expanduser()
DEST_BASE = Path(PWD) / "CHT" / "data"

COPY_PH = True      # copy/replace {station}.ph   (fundamental mode)
COPY_HPH = True     # copy        {station}.hph   (1st higher mode)

OVERWRITE = True    # replace existing destination files
DRY_RUN = False     # if True, only print what WOULD happen (no writes)
VERBOSE = True

# Create CHT/data/{station}_data/ when it does not exist yet.
# The other control files (in.data, mod.*, *.control, in.connector) are still
# written afterwards by SetupData; this only stages the raw .ph/.hph + dir.
CREATE_MISSING_DIRS = True

# Only treat purely-numeric stems as stations (matches the 00000 convention).
# Filters out stray files like "6_8_sec_data.ph".
STATION_ONLY_NUMERIC = True


# =======================
# helpers
# =======================
def safe_copy(src: Path, dst: Path, overwrite: bool = True, dry_run: bool = False) -> str:
    """Copy src -> dst. Returns one of: 'copied', 'skip', 'no_dir'."""
    if not dst.parent.exists():
        if CREATE_MISSING_DIRS:
            if not dry_run:
                dst.parent.mkdir(parents=True, exist_ok=True)
        else:
            return "no_dir"

    if dst.exists() and not overwrite:
        if VERBOSE:
            print(f"[SKIP] exists: {dst}")
        return "skip"

    if dry_run:
        print(f"[DRY]  {src}  ->  {dst}")
        return "copied"

    shutil.copy2(src, dst)
    if VERBOSE:
        print(f"[OK]   {src.name} -> {dst.parent.name}/{dst.name}")
    return "copied"


def source_station_names(src_dir: Path):
    """All station stems that have a .ph and/or .hph in src_dir (sorted)."""
    names = set()
    for p in list(src_dir.glob("*.ph")) + list(src_dir.glob("*.hph")):
        stem = p.stem
        if STATION_ONLY_NUMERIC and not re.fullmatch(r"\d+", stem):
            if VERBOSE:
                print(f"[STRAY] non-station file ignored: {p.name}")
            continue
        names.add(stem)
    return sorted(names)


def main():
    if not SRC_DIR.exists():
        raise FileNotFoundError(f"Source directory not found: {SRC_DIR}")
    if not DEST_BASE.exists():
        raise FileNotFoundError(f"Destination base not found: {DEST_BASE}")

    names = source_station_names(SRC_DIR)

    n_sta = 0
    n_ph_copied = 0
    n_hph_copied = 0
    n_ph_missing_src = 0      # dest dir exists but no source .ph
    ph_only = []              # stations updated with .ph but no .hph in source
    no_dest_dir = []          # source station has no CHT/data/{sta}_data dir

    for name in names:
        dest_dir = DEST_BASE / f"{name}_data"
        src_ph = SRC_DIR / f"{name}.ph"
        src_hph = SRC_DIR / f"{name}.hph"
        dst_ph = dest_dir / f"{name}.ph"
        dst_hph = dest_dir / f"{name}.hph"

        # Skip stations that are not set up for inversion yet.
        if not dest_dir.exists() and not CREATE_MISSING_DIRS:
            no_dest_dir.append(name)
            if VERBOSE:
                print(f"[NODIR] {dest_dir} does not exist -> skipped {name}")
            continue

        n_sta += 1

        # ---- fundamental .ph (replace) ----
        if COPY_PH:
            if src_ph.exists():
                res = safe_copy(src_ph, dst_ph, overwrite=OVERWRITE, dry_run=DRY_RUN)
                if res == "copied":
                    n_ph_copied += 1
            else:
                n_ph_missing_src += 1
                if VERBOSE:
                    print(f"[MISS] no source .ph for {name}")

        # ---- 1st higher mode .hph (optional) ----
        if COPY_HPH:
            if src_hph.exists():
                res = safe_copy(src_hph, dst_hph, overwrite=OVERWRITE, dry_run=DRY_RUN)
                if res == "copied":
                    n_hph_copied += 1
            else:
                ph_only.append(name)

    # dest station dirs that got no source data at all (kept unchanged)
    dest_names = {d.name[:-len("_data")] for d in DEST_BASE.glob("*_data") if d.is_dir()}
    src_names = set(names)
    dest_no_source = sorted(dest_names - src_names)

    print("\n==== Summary ====")
    print(f"Source dir              : {SRC_DIR}")
    print(f"Destination base        : {DEST_BASE}")
    print(f"Mode                    : {'DRY-RUN (no writes)' if DRY_RUN else 'WRITE'}")
    print(f"Source stations (.ph/.hph): {len(names)}")
    print(f"Stations updated        : {n_sta}")
    print(f".ph  copied/replaced    : {n_ph_copied}")
    print(f".hph copied             : {n_hph_copied}")
    print(f"Updated but no src .ph  : {n_ph_missing_src}")
    print(f"Updated, fundamental-only (no .hph): {len(ph_only)}")
    print(f"Source stations w/o dest dir (skipped): {len(no_dest_dir)}")
    if no_dest_dir:
        print(f"    {', '.join(no_dest_dir)}")
    print(f"Dest dirs w/o any source (left as-is) : {len(dest_no_source)}")
    if dest_no_source:
        print(f"    {', '.join(dest_no_source)}")


if __name__ == "__main__":
    main()
