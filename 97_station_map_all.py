#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
os.environ["GMT_COMPATIBILITY"] = "6"

from pathlib import Path
import math
import re
import numpy as np
import pandas as pd
import pygmt
from pyproj import Geod

# =========================
# PARAMETERS
# =========================
PWD = os.getcwd()
FN = "CHT"

# station file format:
# name   long   lat   dist
STATION_FILE = os.path.join(PWD, "plot_data", "station_use.lst")
SEP = r"\s+"
HEADER = None

NAME_COL = "name"
LON_COL = "long"
LAT_COL = "lat"
DIST_COL = "dist"

OUT_DIR = Path(PWD) / "output_figures" / "station_profile_maps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------
# Region / map settings
# -------------------------
REGION_PAD_DEG = 0.05
MIN_REGION_SPAN_DEG = 0.10
RELIEF_RESOLUTION = "01s"
RELIEF_INC_DEG = 1.0 / 3600.0
MAP_CPT = "geo"
CPT_SERIES = [0,300,10]
DPI = 300

# -------------------------
# Profile settings
# -------------------------
PROFILE_SAMPLE_SPACING_KM = 0.10
PROFILE_Y_PAD = 0.8
FLIP_PROFILE_X = True   # False: left->right, True: right->left

# -------------------------
# Panel sizes (cm)
# -------------------------
TOP_WIDTH = 15
TOP_HEIGHT = 3.0
PROFILE_WIDTH = 15
PROFILE_HEIGHT = 2.3
MAP_WIDTH = 15

# vertical gaps between panels (cm)
GAP_TOP_PROFILE = 0.6
GAP_PROFILE_MAP = 0.8

# -------------------------
# Plot styles
# -------------------------
LINE_PEN = "1.2p,red"

OTHER_MAP_SYMBOL = "t0.22c"
CURRENT_MAP_SYMBOL = "a0.34c"

OTHER_PROFILE_SYMBOL = "t0.24c"
CURRENT_PROFILE_SYMBOL = "a0.34c"

OTHER_TOPO_SYMBOL = "t0.15c"
CURRENT_TOPO_SYMBOL = "a0.25c"

FONT_LABEL = "9p,Helvetica-Bold,black"


def enforce_modern_gmt():
    os.environ["GMT_COMPATIBILITY"] = "6"
    pygmt.config(GMT_COMPATIBILITY="6")


def normalize_station_name(value) -> str:
    """
    Keep station names as strings and preserve / enforce zero padding like 00000.
    If the field is purely numeric (or numeric with .0), convert to zero-padded 5-char string.
    Otherwise keep the original string.
    """
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return s

    if re.fullmatch(r"[+-]?\d+(\.0+)?", s):
        s = str(int(float(s)))
        return s.zfill(5)

    if s.isdigit():
        return s.zfill(5)

    return s


def read_stations() -> pd.DataFrame:
    station_path = Path(STATION_FILE).expanduser().resolve()
    if not station_path.exists():
        raise FileNotFoundError(f"Station file not found: {station_path}")

    if HEADER is None:
        df = pd.read_csv(
            station_path,
            sep=SEP,
            header=None,
            names=[NAME_COL, LON_COL, LAT_COL, DIST_COL],
            dtype={NAME_COL: "string"},
            engine="python",
        )
    else:
        df = pd.read_csv(
            station_path,
            sep=SEP,
            header=HEADER,
            dtype={NAME_COL: "string"},
            engine="python",
        )

    df[NAME_COL] = df[NAME_COL].map(normalize_station_name)

    for col in [LON_COL, LAT_COL, DIST_COL]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[NAME_COL, LON_COL, LAT_COL, DIST_COL]).reset_index(drop=True)
    df = df.sort_values(DIST_COL).reset_index(drop=True)

    if len(df) < 2:
        raise ValueError("Need at least 2 stations to build the profile and station line.")

    return df


def outward_round_2dec(vmin: float, vmax: float):
    vmin2 = math.floor(vmin * 100.0) / 100.0
    vmax2 = math.ceil(vmax * 100.0) / 100.0
    return vmin2, vmax2


def enforce_min_span(vmin: float, vmax: float, min_span: float):
    span = vmax - vmin
    if span >= min_span:
        return vmin, vmax
    center = 0.5 * (vmin + vmax)
    half = 0.5 * min_span
    return center - half, center + half


def compute_regions(df: pd.DataFrame):
    minlon = float(df[LON_COL].min()) - REGION_PAD_DEG
    maxlon = float(df[LON_COL].max()) + REGION_PAD_DEG
    minlat = float(df[LAT_COL].min()) - REGION_PAD_DEG
    maxlat = float(df[LAT_COL].max()) + REGION_PAD_DEG

    minlon, maxlon = enforce_min_span(minlon, maxlon, MIN_REGION_SPAN_DEG)
    minlat, maxlat = enforce_min_span(minlat, maxlat, MIN_REGION_SPAN_DEG)

    minlon, maxlon = outward_round_2dec(minlon, maxlon)
    minlat, maxlat = outward_round_2dec(minlat, maxlat)
    plot_region = [minlon, maxlon, minlat, maxlat]

    grid_region = [
        math.floor(minlon / RELIEF_INC_DEG) * RELIEF_INC_DEG,
        math.ceil(maxlon / RELIEF_INC_DEG) * RELIEF_INC_DEG,
        math.floor(minlat / RELIEF_INC_DEG) * RELIEF_INC_DEG,
        math.ceil(maxlat / RELIEF_INC_DEG) * RELIEF_INC_DEG,
    ]

    return plot_region, grid_region


def build_dense_track(df: pd.DataFrame, spacing_km: float = 0.1) -> pd.DataFrame:
    geod = Geod(ellps="WGS84")
    rows = []

    for i in range(len(df) - 1):
        lon1 = float(df.loc[i, LON_COL])
        lat1 = float(df.loc[i, LAT_COL])
        lon2 = float(df.loc[i + 1, LON_COL])
        lat2 = float(df.loc[i + 1, LAT_COL])
        d1 = float(df.loc[i, DIST_COL])
        d2 = float(df.loc[i + 1, DIST_COL])

        _, _, seg_m = geod.inv(lon1, lat1, lon2, lat2)
        seg_km = max(seg_m / 1000.0, 1e-6)
        npts = max(2, int(np.ceil(seg_km / spacing_km)) + 1)

        lons = np.linspace(lon1, lon2, npts)
        lats = np.linspace(lat1, lat2, npts)
        dists = np.linspace(d1, d2, npts)

        if i > 0:
            lons = lons[1:]
            lats = lats[1:]
            dists = dists[1:]

        rows.extend(zip(lons, lats, dists))

    return pd.DataFrame(rows, columns=[LON_COL, LAT_COL, DIST_COL])


def sample_topography(track_df: pd.DataFrame, grid_region):
    grid = pygmt.datasets.load_earth_relief(
        resolution=RELIEF_RESOLUTION,
        region=grid_region,
    )

    topo_track = pygmt.grdtrack(
        points=track_df[[LON_COL, LAT_COL]],
        grid=grid,
        newcolname="elev_m",
    )
    topo_track[DIST_COL] = track_df[DIST_COL].to_numpy()
    return topo_track, grid


def station_topography(df: pd.DataFrame, grid) -> pd.DataFrame:
    sta_topo = pygmt.grdtrack(
        points=df[[LON_COL, LAT_COL]],
        grid=grid,
        newcolname="elev_m",
    )
    sta_topo[NAME_COL] = df[NAME_COL].to_numpy()
    sta_topo[DIST_COL] = df[DIST_COL].to_numpy()
    return sta_topo


def estimate_map_height_cm(region, map_width_cm):
    lonmin, lonmax, latmin, latmax = region
    geod = Geod(ellps="WGS84")

    clon = 0.5 * (lonmin + lonmax)
    clat = 0.5 * (latmin + latmax)

    _, _, width_m = geod.inv(lonmin, clat, lonmax, clat)
    _, _, height_m = geod.inv(clon, latmin, clon, latmax)

    width_km = max(width_m / 1000.0, 1e-6)
    height_km = max(height_m / 1000.0, 1e-6)

    return map_width_cm * (height_km / width_km)


def profile_projection(width_cm, height_cm, flip=False):
    sign = "-" if flip else ""
    return f"X{sign}{width_cm}c/{height_cm}c"


def get_profile_xlabel():
    if FLIP_PROFILE_X:
        return "Distance from 1st station (km, reversed)"
    return "Distance from 1st station (km)"


def make_one_figure(
    df: pd.DataFrame,
    topo_track: pd.DataFrame,
    sta_topo: pd.DataFrame,
    grid,
    plot_region,
    current_idx: int,
):
    cur = df.iloc[current_idx]
    cur_topo = sta_topo.iloc[current_idx]

    cur_name = str(cur[NAME_COL])
    cur_lon = float(cur[LON_COL])
    cur_lat = float(cur[LAT_COL])
    cur_dist = float(cur[DIST_COL])
    cur_elev = float(cur_topo["elev_m"])

    dist_min = float(df[DIST_COL].min())
    dist_max = float(df[DIST_COL].max())

    topo_min = float(np.nanmin(topo_track["elev_m"]))
    topo_max = float(np.nanmax(topo_track["elev_m"]))
    topo_range = max(topo_max - topo_min, 300.0)
    topo_pad = 0.08 * topo_range
    baseline = topo_min - topo_pad

    topo_region = [dist_min, dist_max, baseline, topo_max + topo_pad]
    profile_region = [dist_min, dist_max, -PROFILE_Y_PAD, PROFILE_Y_PAD]

    map_height_cm = estimate_map_height_cm(plot_region, MAP_WIDTH)
    shift_to_profile = f"-{PROFILE_HEIGHT + GAP_TOP_PROFILE:.2f}c"
    shift_to_map = f"-{map_height_cm + GAP_PROFILE_MAP:.2f}c"

    topo_poly_x = np.r_[topo_track[DIST_COL].to_numpy(), dist_max, dist_min]
    topo_poly_y = np.r_[topo_track["elev_m"].to_numpy(), baseline, baseline]

    shade = pygmt.grdgradient(grid=grid, azimuth=315, normalize="t0.8")
    xlab = get_profile_xlabel()

    fig = pygmt.Figure()
    pygmt.config(
        MAP_FRAME_TYPE="plain",
        FONT_LABEL="15p,Times-Bold,black",
        FONT_TITLE="15p,Times-Bold,black",
        FONT_ANNOT_PRIMARY="15p,Times-Bold,black",
        FONT_ANNOT_SECONDARY="15p,Times-Roman,black"
    )

    # ---------------------------------
    # 1) topography profile
    # ---------------------------------
    fig.basemap(
        region=topo_region,
        projection=profile_projection(TOP_WIDTH, TOP_HEIGHT, flip=FLIP_PROFILE_X),
        frame=[
            f"xaf+l{xlab}",
            "yaf+lElevation (m)",
            f"+t{FN}: station {cur_name}",
        ],
    )

    fig.plot(
        x=topo_poly_x,
        y=topo_poly_y,
        fill="gray85",
        pen="0.8p,black",
    )

    fig.plot(
        x=topo_track[DIST_COL],
        y=topo_track["elev_m"],
        pen="1.0p,black",
    )

    fig.plot(
        x=sta_topo[DIST_COL],
        y=sta_topo["elev_m"],
        style=OTHER_TOPO_SYMBOL,
        pen="0.25p,black",
        fill="black",
    )

    fig.plot(
        x=[cur_dist],
        y=[cur_elev],
        style=CURRENT_TOPO_SYMBOL,
        pen="0.5p,black",
        fill="red",
    )

    # ---------------------------------
    # 2) stations along profile
    # ---------------------------------
    fig.shift_origin(yshift=shift_to_profile)

    # fig.basemap(
    #     region=profile_region,
    #     projection=profile_projection(PROFILE_WIDTH, PROFILE_HEIGHT, flip=FLIP_PROFILE_X),
    #     frame=[
    #         f"xaf+l{xlab}",
    #         "yaf+lRelative position",
    #         "+tStations along the profile",
    #     ],
    # )

    # fig.plot(
    #     x=[dist_min, dist_max],
    #     y=[0, 0],
    #     pen="1.0p,black",
    # )

    # fig.plot(
    #     x=df[DIST_COL],
    #     y=np.zeros(len(df)),
    #     style=OTHER_PROFILE_SYMBOL,
    #     pen="0.35p,black",
    #     fill="black",
    # )

    # fig.plot(
    #     x=[cur_dist],
    #     y=[0.0],
    #     style=CURRENT_PROFILE_SYMBOL,
    #     pen="0.55p,black",
    #     fill="red",
    # )

    # fig.text(
    #     x=[cur_dist],
    #     y=[0.34],
    #     text=[cur_name],
    #     font=FONT_LABEL,
    #     justify="BC",
    #     fill="white@60",
    #     pen="0.15p,black",
    # )

    # ---------------------------------
    # 3) lower map
    # ---------------------------------
    fig.shift_origin(yshift=shift_to_map)

    if CPT_SERIES is None:
        pygmt.makecpt(cmap=MAP_CPT, series=[float(grid.min()), float(grid.max())],continuous=True)
    else:
        pygmt.makecpt(cmap=MAP_CPT, series=CPT_SERIES,continuous=True)

    fig.grdimage(
        grid=grid,
        region=plot_region,
        projection=f"M{MAP_WIDTH}c",
        shading=shade,
        cmap=True,
        frame=[
            "xaf+lLongitude",
            "yaf+lLatitude",
            "+tHorizontal map",
        ],
    )

    fig.coast(
        region=plot_region,
        projection=f"M{MAP_WIDTH}c",
        shorelines="0.4p,black",
        borders="1/0.4p,black",
        rivers="a/0.2p,blue",
    )

    fig.plot(
        x=df[LON_COL],
        y=df[LAT_COL],
        pen=LINE_PEN,
    )

    fig.plot(
        x=df[LON_COL],
        y=df[LAT_COL],
        style=OTHER_MAP_SYMBOL,
        pen="0.35p,black",
        fill="black",
    )

    fig.plot(
        x=[cur_lon],
        y=[cur_lat],
        style=CURRENT_MAP_SYMBOL,
        pen="0.6p,black",
        fill="red",
    )

    fig.text(
        x=[cur_lon],
        y=[cur_lat],
        text=[cur_name],
        font=FONT_LABEL,
        justify="LM",
        offset="0.12c/0.12c",
        fill="white@60",
        pen="0.15p,black",
    )

    fig.plot(
        data="plot_data/faults/CGS_2021_SFT.txt",
        pen="1.0p,red",
        label="CGS_2021"
    )

    fig.plot(
        data="plot_data/faults/TEM_2016_SFT.txt",
        pen="1.0p,black",
        label="TEM_2016"
    )

    fig.colorbar(
        position="JMR+o0.9c/0c+w8c/0.35c",
        frame='af+l"Elevation (m)"',
    )

    fig.legend()

    out_png = OUT_DIR / f"{cur_name}_{FN}_station_profile_map.png"
    fig.savefig(str(out_png), dpi=DPI)
    return out_png


def main():
    enforce_modern_gmt()

    df = read_stations()
    plot_region, grid_region = compute_regions(df)

    print("Plot region (2 decimals):", [f"{v:.2f}" for v in plot_region])
    print("Grid region (1-sec snapped):", [f"{v:.6f}" for v in grid_region])
    print(f"Stations: {len(df)}")
    print(f"Output  : {OUT_DIR}")

    dense_track = build_dense_track(df, spacing_km=PROFILE_SAMPLE_SPACING_KM)
    topo_track, grid = sample_topography(dense_track, grid_region)
    sta_topo = station_topography(df, grid)

    for i in range(len(df)):
        out = make_one_figure(df, topo_track, sta_topo, grid, plot_region, i)
        print(f"[{i+1}/{len(df)}] Saved: {out}")


if __name__ == "__main__":
    main()