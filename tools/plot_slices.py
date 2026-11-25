#!/usr/bin/env python3
"""Plots of the slices written by tools/slices.py (.mat or .npz), as tools/plot_slices.m (python3 with numpy and
matplotlib; scipy for .mat files).

    python3 tools/plot_slices.py slices.mat                       # all fields: slices_<field>.png
    python3 tools/plot_slices.py slices.npz --fields u s1 --clim -0.3 0.3 --plane xz

One figure per field: the xz-planes (one per stored y) in the top row, the xy-plane (at z = z_xy) and the yz-plane
(at x = x_yz) below. Coordinates in units of the half-height h. Fields with both signs use a symmetric blue-white-red
colour scale, the others viridis.
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("file")
ap.add_argument("--fields", nargs="+", default=None)
ap.add_argument("--clim", nargs=2, type=float, default=None)
ap.add_argument("--plane", choices=["all", "xz", "xy", "yz"], default="all")
ap.add_argument("--dpi", type=int, default=150)
a = ap.parse_args()

if a.file.endswith(".npz"):
    S = dict(np.load(a.file))
else:
    from scipy.io import loadmat
    S = loadmat(a.file, squeeze_me=True)
sc = lambda k: float(np.ravel(S[k])[0])
x, y, z = np.ravel(S["x"]), np.ravel(S["y"]), np.ravel(S["z"])
lx, lz = sc("lx"), sc("lz")
fields = a.fields or [str(f).strip() for f in np.ravel(S["fields"])]
yxz = np.ravel(S["y_xz"])
fl = int(sc("fluct")) if "fluct" in S else 0
xp, zp = np.append(x, lx), np.append(z, lz)  # periodic: close the box
base = os.path.splitext(a.file)[0]
for f in fields:
    if "xy_" + f not in S:
        print("field %s not in %s" % (f, a.file)); continue
    XZ = np.asarray(S["xz_" + f], dtype=float)
    if XZ.ndim == 2:
        XZ = XZ[:, :, None]
    XY, YZ = np.asarray(S["xy_" + f], dtype=float), np.asarray(S["yz_" + f], dtype=float)
    XZ = np.concatenate([XZ, XZ[:, :1]], axis=1); XZ = np.concatenate([XZ, XZ[:1]], axis=0)
    XY = np.concatenate([XY, XY[:, :1]], axis=1); YZ = np.concatenate([YZ, YZ[:, :1]], axis=1)
    if a.clim:
        cl = a.clim
    else:
        lo, hi = min(XZ.min(), XY.min(), YZ.min()), max(XZ.max(), XY.max(), YZ.max())
        cl = [-max(-lo, hi), max(-lo, hi)] if lo < 0 < hi else [lo, hi]
    cmap = "RdBu_r" if cl[0] < 0 < cl[1] else "viridis"
    lab = f + ("'" if fl else "")
    doxz, doxy, doyz = a.plane in ("all", "xz"), a.plane in ("all", "xy"), a.plane in ("all", "yz")
    nxz = XZ.shape[2]
    nrow = int(doxz) + int(doxy or doyz)
    ncol = max(nxz if doxz else 0, int(doxy) + int(doyz), 1)
    fig, ax = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow), squeeze=False)
    for axx in ax.flat:
        axx.set_visible(False)
    r = 0
    def show(axx, X, Y, F, xl, yl, title):
        axx.set_visible(True)
        m = axx.pcolormesh(X, Y, F, cmap=cmap, vmin=cl[0], vmax=cl[1], shading="nearest")
        axx.set_aspect("equal"); axx.set_xlabel(xl); axx.set_ylabel(yl); axx.set_title(title, fontsize=9)
        fig.colorbar(m, ax=axx, shrink=0.8)
    if doxz:
        for k in range(nxz):
            show(ax[r, k], xp, zp, XZ[:, :, k], "x/h", "z/h",
                 "%s, xz-plane at y = %.4f (%.4f from wall)" % (lab, yxz[k], 1 - abs(yxz[k])))
        r += 1
    if doxy or doyz:
        c = 0
        if doxy:
            show(ax[r, c], xp, y, XY, "x/h", "y/h", "%s, xy-plane at z = %.4f" % (lab, sc("z_xy"))); c += 1
        if doyz:
            show(ax[r, c], zp, y, YZ, "z/h", "y/h", "%s, yz-plane at x = %.4f" % (lab, sc("x_yz")))
    fig.suptitle("%s   t = %.6g" % (lab, sc("time")))
    fig.tight_layout()
    png = "%s_%s.png" % (base, f)
    fig.savefig(png, dpi=a.dpi)
    plt.close(fig)
    print("wrote", png)
