#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combine station profile-map PNG and MCMC PNG side-by-side.

Current filename settings:
- Map figure:
    {station}_{FN}_station_profile_map.png
- MCMC figure:
    MCMC_{station}.png

Output:
- QC_{station}.png

This version matches pairs by STATION NAME, not by lon-lat tag.
It keeps names like 00000 as strings.
"""

import os
from pathlib import Path
from PIL import Image
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

# =======================
# Parameters
# =======================
PWD = os.getcwd()
FN = "CHT"   # change to "CHT" if needed

MAP_DIR = Path(PWD) / "output_figures" / "station_profile_maps"
MCMC_DIR = Path(PWD) / "output_figures" / "1D_model"
OUT_DIR = Path(PWD) / "output_figures" / "QC_map_mcmc"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAP_SUFFIX = f"_{FN}_station_profile_map.png"
MCMC_PREFIX = "MCMC_"

BACKGROUND = (255, 255, 255)   # white
GAP_PX = 20
ALIGN_HEIGHT = True
DPI = 300

REQUESTED_WORKERS = 8
CPU_COUNT = os.cpu_count() or 1
N_WORKERS = min(REQUESTED_WORKERS, CPU_COUNT)


# =======================
# functions
# =======================
def extract_station_from_map(fname: str) -> Optional[str]:
    """
    Map filename format:
        {station}_{FN}_station_profile_map.png

    Example:
        00000_FM_station_profile_map.png
        ABC01_FM_station_profile_map.png
    """
    if not fname.endswith(MAP_SUFFIX):
        return None
    return fname[: -len(MAP_SUFFIX)]


def extract_station_from_mcmc(fname: str) -> Optional[str]:
    """
    MCMC filename format:
        MCMC_{station}.png

    Example:
        MCMC_00000.png
        MCMC_ABC01.png
    """
    prefix = MCMC_PREFIX
    suffix = ".png"
    if not fname.startswith(prefix) or not fname.endswith(suffix):
        return None
    return fname[len(prefix):-len(suffix)]


def combine_pair(map_png: Path, mcmc_png: Path, out_png: Path):
    im_map = Image.open(map_png).convert("RGBA")
    im_mcmc = Image.open(mcmc_png).convert("RGBA")

    if ALIGN_HEIGHT:
        target_h = im_map.height
        if im_mcmc.height != target_h:
            new_w = int(im_mcmc.width * (target_h / im_mcmc.height))
            im_mcmc = im_mcmc.resize((new_w, target_h), resample=Image.LANCZOS)

    out_w = im_map.width + GAP_PX + im_mcmc.width
    out_h = max(im_map.height, im_mcmc.height)

    canvas = Image.new("RGBA", (out_w, out_h), BACKGROUND + (255,))
    canvas.paste(im_map, (0, 0), im_map)
    canvas.paste(im_mcmc, (im_map.width + GAP_PX, 0), im_mcmc)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(out_png, dpi=(DPI, DPI))


def _worker(payload):
    map_png, mcmc_png, out_png = payload
    combine_pair(map_png, mcmc_png, out_png)
    return str(out_png)


def main():
    if not MAP_DIR.exists():
        raise FileNotFoundError(f"MAP_DIR not found: {MAP_DIR}")
    if not MCMC_DIR.exists():
        raise FileNotFoundError(f"MCMC_DIR not found: {MCMC_DIR}")

    maps = {}
    for p in MAP_DIR.glob(f"*{MAP_SUFFIX}"):
        sta = extract_station_from_map(p.name)
        if sta:
            maps[sta] = p

    mcmcs = {}
    for p in MCMC_DIR.glob(f"{MCMC_PREFIX}*.png"):
        sta = extract_station_from_mcmc(p.name)
        if sta:
            mcmcs[sta] = p

    common_stations = sorted(set(maps).intersection(mcmcs))
    map_only = sorted(set(maps) - set(mcmcs))
    mcmc_only = sorted(set(mcmcs) - set(maps))

    print("Found maps     :", len(maps))
    print("Found MCMC     :", len(mcmcs))
    print("Matched pairs  :", len(common_stations))
    print("Map only       :", len(map_only))
    print("MCMC only      :", len(mcmc_only))
    print("OUT_DIR        :", OUT_DIR)
    print("Workers        :", N_WORKERS)

    if map_only:
        print("\nStations with map only:")
        for sta in map_only[:20]:
            print("  ", sta)
        if len(map_only) > 20:
            print(f"  ... ({len(map_only) - 20} more)")

    if mcmc_only:
        print("\nStations with MCMC only:")
        for sta in mcmc_only[:20]:
            print("  ", sta)
        if len(mcmc_only) > 20:
            print(f"  ... ({len(mcmc_only) - 20} more)")

    jobs = []
    for sta in common_stations:
        out_png = OUT_DIR / f"QC_{sta}.png"
        jobs.append((maps[sta], mcmcs[sta], out_png))

    if not jobs:
        print("\nNo matched pairs found.")
        return

    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for fut in as_completed(futs):
            done += 1
            try:
                out = fut.result()
                print(f"[{done}/{len(jobs)}] {out}")
            except Exception as e:
                print(f"[{done}/{len(jobs)}] ERROR: {e}")

    print("Finished.")


if __name__ == "__main__":
    main()