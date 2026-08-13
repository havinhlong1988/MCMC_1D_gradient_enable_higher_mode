#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combine 3 PNG figures side-by-side:

    Map | Previous | Current

Inputs
------
1) Map figure:
   output_figures/station_profile_maps/{station}_{FN}_station_profile_map.png

2) Previous MCMC figure:
   output_figures/1D_model_previous/MCMC_{station}.png

3) Current MCMC figure:
   output_figures/1D_model/MCMC_{station}.png

Output
------
output_figures/QC_map_previous_current/QC_{station}.png

Notes
-----
- Pairing is by station name string, preserving names like 00000.
- Adds text labels above the Previous and Current panels.
- Also adds a Map label for consistency.
"""

import os
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

# =======================
# Parameters
# =======================
PWD = os.getcwd()
FN = "CHT"   # example: "CHT"

MAP_DIR = Path(PWD) / "output_figures" / "station_profile_maps"
PREV_DIR = Path(PWD) / "output_figures" / "1D_model_previous"
CURR_DIR = Path(PWD) / "output_figures" / "1D_model"
OUT_DIR = Path(PWD) / "output_figures" / "QC_map_previous_current"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAP_SUFFIX = f"_{FN}_station_profile_map.png"
MCMC_PREFIX = "MCMC_"
PNG_SUFFIX = ".png"

BACKGROUND = (255, 255, 255)   # white
TEXT_COLOR = (0, 0, 0)         # black
GAP_PX = 20
TOP_PAD_PX = 70                # top band for labels
SIDE_PAD_PX = 0
BOTTOM_PAD_PX = 0
ALIGN_HEIGHT = True
DPI = 300

MAP_LABEL = "Map"
PREVIOUS_LABEL = "Previous"
CURRENT_LABEL = "Current"

REQUESTED_WORKERS = 8
CPU_COUNT = os.cpu_count() or 1
N_WORKERS = min(REQUESTED_WORKERS, CPU_COUNT)


# =======================
# helpers
# =======================
def extract_station_from_map(fname: str) -> Optional[str]:
    """
    Map filename format:
        {station}_{FN}_station_profile_map.png
    """
    if not fname.endswith(MAP_SUFFIX):
        return None
    return fname[:-len(MAP_SUFFIX)]


def extract_station_from_mcmc(fname: str) -> Optional[str]:
    """
    MCMC filename format:
        MCMC_{station}.png
    """
    if not fname.startswith(MCMC_PREFIX) or not fname.endswith(PNG_SUFFIX):
        return None
    return fname[len(MCMC_PREFIX):-len(PNG_SUFFIX)]


def load_font(size: int = 28):
    """
    Try some common TrueType fonts. Fall back to PIL default if not found.
    """
    candidates = [
        "DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
        "/Library/Fonts/Arial Bold.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "arialbd.ttf",
    ]
    for f in candidates:
        try:
            return ImageFont.truetype(f, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def resize_to_target_height(img: Image.Image, target_h: int) -> Image.Image:
    if img.height == target_h:
        return img
    new_w = int(round(img.width * (target_h / img.height)))
    return img.resize((new_w, target_h), resample=Image.LANCZOS)


def text_bbox(draw: ImageDraw.ImageDraw, text: str, font) -> tuple:
    try:
        return draw.textbbox((0, 0), text, font=font)
    except Exception:
        w, h = draw.textsize(text, font=font)
        return (0, 0, w, h)


def draw_centered_label(draw: ImageDraw.ImageDraw, x0: int, x1: int, y0: int, text: str, font):
    bbox = text_bbox(draw, text, font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    x = int(round((x0 + x1 - tw) / 2))
    y = int(round(y0 + max(0, (TOP_PAD_PX - th) / 2)))
    draw.text((x, y), text, fill=TEXT_COLOR, font=font)


def combine_triplet(map_png: Path, prev_png: Path, curr_png: Path, out_png: Path):
    im_map = Image.open(map_png).convert("RGBA")
    im_prev = Image.open(prev_png).convert("RGBA")
    im_curr = Image.open(curr_png).convert("RGBA")

    if ALIGN_HEIGHT:
        target_h = max(im_map.height, im_prev.height, im_curr.height)
        im_map = resize_to_target_height(im_map, target_h)
        im_prev = resize_to_target_height(im_prev, target_h)
        im_curr = resize_to_target_height(im_curr, target_h)

    panel_h = max(im_map.height, im_prev.height, im_curr.height)
    out_h = TOP_PAD_PX + panel_h + BOTTOM_PAD_PX
    out_w = (
        SIDE_PAD_PX
        + im_map.width
        + GAP_PX
        + im_prev.width
        + GAP_PX
        + im_curr.width
        + SIDE_PAD_PX
    )

    canvas = Image.new("RGBA", (out_w, out_h), BACKGROUND + (255,))

    x_map = SIDE_PAD_PX
    x_prev = x_map + im_map.width + GAP_PX
    x_curr = x_prev + im_prev.width + GAP_PX
    y_img = TOP_PAD_PX

    canvas.paste(im_map, (x_map, y_img), im_map)
    canvas.paste(im_prev, (x_prev, y_img), im_prev)
    canvas.paste(im_curr, (x_curr, y_img), im_curr)

    draw = ImageDraw.Draw(canvas)
    font = load_font(size=28)

    draw_centered_label(draw, x_map, x_map + im_map.width, 0, MAP_LABEL, font)
    draw_centered_label(draw, x_prev, x_prev + im_prev.width, 0, PREVIOUS_LABEL, font)
    draw_centered_label(draw, x_curr, x_curr + im_curr.width, 0, CURRENT_LABEL, font)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(out_png, dpi=(DPI, DPI))


def _worker(payload):
    map_png, prev_png, curr_png, out_png = payload
    combine_triplet(map_png, prev_png, curr_png, out_png)
    return str(out_png)


def main():
    if not MAP_DIR.exists():
        raise FileNotFoundError(f"MAP_DIR not found: {MAP_DIR}")
    if not PREV_DIR.exists():
        raise FileNotFoundError(f"PREV_DIR not found: {PREV_DIR}")
    if not CURR_DIR.exists():
        raise FileNotFoundError(f"CURR_DIR not found: {CURR_DIR}")

    maps = {}
    for p in MAP_DIR.glob(f"*{MAP_SUFFIX}"):
        sta = extract_station_from_map(p.name)
        if sta:
            maps[sta] = p

    prevs = {}
    for p in PREV_DIR.glob(f"{MCMC_PREFIX}*{PNG_SUFFIX}"):
        sta = extract_station_from_mcmc(p.name)
        if sta:
            prevs[sta] = p

    currs = {}
    for p in CURR_DIR.glob(f"{MCMC_PREFIX}*{PNG_SUFFIX}"):
        sta = extract_station_from_mcmc(p.name)
        if sta:
            currs[sta] = p

    common_stations = set(maps) & set(prevs) & set(currs)
    common_stations_sorted = sorted(common_stations)

    map_only = sorted(set(maps) - common_stations)
    prev_only = sorted(set(prevs) - common_stations)
    curr_only = sorted(set(currs) - common_stations)

    print("Found maps        :", len(maps))
    print("Found previous    :", len(prevs))
    print("Found current     :", len(currs))
    print("Matched triplets  :", len(common_stations_sorted))
    print("Map only          :", len(map_only))
    print("Previous only     :", len(prev_only))
    print("Current only      :", len(curr_only))
    print("OUT_DIR           :", OUT_DIR)
    print("Workers           :", N_WORKERS)

    if map_only:
        print("\nStations with map only:")
        for sta in map_only[:20]:
            print("  ", sta)
        if len(map_only) > 20:
            print(f"  ... ({len(map_only) - 20} more)")

    if prev_only:
        print("\nStations with previous only:")
        for sta in prev_only[:20]:
            print("  ", sta)
        if len(prev_only) > 20:
            print(f"  ... ({len(prev_only) - 20} more)")

    if curr_only:
        print("\nStations with current only:")
        for sta in curr_only[:20]:
            print("  ", sta)
        if len(curr_only) > 20:
            print(f"  ... ({len(curr_only) - 20} more)")

    jobs = []
    for sta in common_stations_sorted:
        out_png = OUT_DIR / f"QC_{sta}.png"
        jobs.append((maps[sta], prevs[sta], currs[sta], out_png))

    if not jobs:
        print("\nNo matched triplets found.")
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
