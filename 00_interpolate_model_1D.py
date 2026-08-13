import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# ================= USER SETTINGS =================
pwd = os.getcwd()
maindir = os.path.join(pwd, "CHT", "Vel_mod")

inmodfile = os.path.join(maindir, "NTW_Liu2021_5km_modified")

# Only sample at these depths (km)
base_depths = np.array([0., 0.5, 1., 2., 3., 5., 15.])

extrapolation = 0  # 1=yes, 0=no  (controls mantle_ref if outside mantle range)
out_file = os.path.join(maindir, "NTW_1D_Liu21_modified.dat")

MOHO = 5.0       # km
MANTLE_REF = 15.0  # km (use mantle velocity at this depth for the 2nd 5-km layer)
# =================================================


def make_interp(x, y, extrapolate: bool):
    """Safe 1D linear interpolator. If extrapolate==False, clamp to endpoints."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    if len(x) == 0:
        return None
    if len(x) == 1:
        const = float(y[0])
        return lambda z: np.full_like(np.asarray(z, float), const, dtype=float)

    fillv = "extrapolate" if extrapolate else (y[0], y[-1])
    f = interp1d(x, y, kind="linear", bounds_error=False, fill_value=fillv)
    return lambda z: np.asarray(f(z), float)


# ---- read + sort ----
indata = pd.read_csv(inmodfile, sep=r"\s+", names=["dep", "vs", "vp"], comment="#")
indata = indata.sort_values("dep").reset_index(drop=True)

# ---- split crust / mantle ----
crust_df = indata[indata["dep"] <= MOHO].copy()   # include Moho for crust-side value
mantle_df = indata[indata["dep"] > MOHO].copy()   # mantle strictly below Moho

if mantle_df.empty:
    raise ValueError("No mantle points found (dep > MOHO). Need at least one mantle row.")

# Build interpolators
do_extrap = (extrapolation == 1)
f_vs_cr = make_interp(crust_df["dep"], crust_df["vs"], do_extrap)
f_vp_cr = make_interp(crust_df["dep"], crust_df["vp"], do_extrap)

f_vs_ma = make_interp(mantle_df["dep"], mantle_df["vs"], do_extrap)
f_vp_ma = make_interp(mantle_df["dep"], mantle_df["vp"], do_extrap)

# ---- compute mantle reference values (at 15 km) ----
vs_ref = float(f_vs_ma(MANTLE_REF))
vp_ref = float(f_vp_ma(MANTLE_REF))

# ---- build output with TWO 5.0-km layers ----
depths_in = np.unique(base_depths.astype(float))
# ensure MOHO present in requested depths
if not np.any(np.isclose(depths_in, MOHO)):
    depths_in = np.sort(np.append(depths_in, MOHO))

depths_out, vs_out, vp_out = [], [], []

for d in np.sort(depths_in):
    if np.isclose(d, MOHO):
        # 1) crust-side 5 km
        depths_out.append(MOHO)
        vs_out.append(float(f_vs_cr(MOHO)))
        vp_out.append(float(f_vp_cr(MOHO)))

        # 2) mantle-side 5 km = mantle at 15 km
        depths_out.append(MOHO)
        vs_out.append(vs_ref)
        vp_out.append(vp_ref)

    elif d < MOHO:
        depths_out.append(d)
        vs_out.append(float(f_vs_cr(d)))
        vp_out.append(float(f_vp_cr(d)))
    else:
        # regular mantle sampling for depths > MOHO
        depths_out.append(d)
        vs_out.append(float(f_vs_ma(d)))
        vp_out.append(float(f_vp_ma(d)))

new_model = pd.DataFrame({"dep": depths_out, "vs": vs_out, "vp": vp_out})

# ---- save ----
new_model = new_model.round(5)
new_model.to_csv(out_file, sep=" ", header=False, index=False, float_format="%.5f")

print("Saved:", out_file)
print(new_model)
print("\nCheck Moho jump rows (dep==5.0):")
print(new_model[np.isclose(new_model["dep"], MOHO)])
