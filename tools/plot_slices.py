#!/usr/bin/env python3
"""Pictures of the slices written by tools/slices.py (.mat or .npz): every plane as a clean image, without colour
bar, axes, title or margins (the aspect ratio of the plane), coloured with the Scientific Colour Maps of Fabio
Crameri (https://www.fabiocrameri.ch/colourmaps, doi:10.5281/zenodo.1243862). Python 3 with numpy and matplotlib;
scipy for .mat files.

    python3 tools/plot_slices.py s.mat                          # all fields and planes
    python3 tools/plot_slices.py s.npz --fields u s1 --clim -0.3 0.3 --plane xz --size 3000
    python3 tools/plot_slices.py s.mat --annotated              # one figure per field with axes and colour bar

Files: <base>_<field>_xz1.png, _xz2.png, ... (the stored wall distances, in order), <base>_<field>_xy.png,
<base>_<field>_yz.png; the colour range of every field is written to <base>_clim.txt (use --clim for the same range
in a series of pictures, e.g. frames of a film).

Colour maps: --cmap-div (default vik) for fields with both signs (fluctuations, v, w), --cmap (default batlow)
otherwise; any map of the collection (e.g. roma, berlin, broc, cork, lajolla, oslo, ...; reversed: name + '_r').
They are read from a directory with the maps as text files (256 RGB lines): --cmap-dir, or the environment variable
SCIENTIFIC_COLOURMAPS, or $HOME/scientific_colourmaps or /scratch/geert/scientific_colourmaps; the layouts of the
official collection (Zenodo: <name>/<name>.txt) and of the Python package cmcrameri (cmcrameri/cmaps/<name>.txt)
both work (tools/install_scientific_colourmaps.sh installs them); else from an installed cmcrameri package.
"""
import argparse, glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("file")
ap.add_argument("--fields", nargs="+", default=None)
ap.add_argument("--clim", nargs=2, type=float, default=None, help="colour range (default: per field, symmetric if signed)")
ap.add_argument("--plane", choices=["all", "xz", "xy", "yz"], default="all")
ap.add_argument("--cmap", default="batlow", help="map for fields of one sign (default batlow)")
ap.add_argument("--cmap-div", default="vik", help="map for fields with both signs (default vik)")
ap.add_argument("--cmap-dir", default=None, help="directory with the Scientific Colour Maps")
ap.add_argument("--size", type=int, default=2000, help="pixels of the longer side of each image (default 2000)")
ap.add_argument("--shading", choices=["gouraud", "nearest"], default="gouraud",
                help="gouraud: smooth between the grid points (default); nearest: one colour per grid point")
ap.add_argument("--annotated", action="store_true", help="one figure per field with titles, axes and colour bars")
a = ap.parse_args()


# ------------------------------------------------------------------------------------------- colour maps
def find_map_file(name):
    dirs = [a.cmap_dir, os.environ.get("SCIENTIFIC_COLOURMAPS"), os.path.expanduser("~/scientific_colourmaps"),
            "/scratch/geert/scientific_colourmaps"]
    for d in dirs:
        if not d or not os.path.isdir(d):
            continue
        for pat in (f"{name}.txt", f"{name}/{name}.txt", f"cmcrameri/cmaps/{name}.txt", f"**/{name}/{name}.txt",
                    f"**/cmaps/{name}.txt", f"**/{name}.txt"):
            hits = sorted(glob.glob(os.path.join(d, pat), recursive=True))
            if hits:
                return hits[0]
    return None


def colourmap(name):
    rev = name.endswith("_r")
    base = name[:-2] if rev else name
    f = find_map_file(base)
    if f:
        rgb = np.loadtxt(f)[:, :3]
        if rgb.max() > 1.0:
            rgb = rgb / 255.0
        cm = ListedColormap(rgb, name=base)
    else:
        try:
            import cmcrameri.cm as cmc
            cm = getattr(cmc, base)
        except (ImportError, AttributeError):
            if base in plt.colormaps():
                return plt.get_cmap(name)
            sys.exit(f"Scientific Colour Map '{base}' not found: install the maps "
                     f"(tools/install_scientific_colourmaps.sh) and give --cmap-dir, or set SCIENTIFIC_COLOURMAPS")
    return cm.reversed() if rev else cm


cm_div, cm_seq = colourmap(a.cmap_div), colourmap(a.cmap)

# ------------------------------------------------------------------------------------------- data
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
doxz, doxy, doyz = a.plane in ("all", "xz"), a.plane in ("all", "xy"), a.plane in ("all", "yz")


def image(png, X, Y, F, w, h, cmap, cl):
    """F on the grid X (columns), Y (rows), physical size w x h, nothing but the field"""
    s = a.size / max(w, h)
    fig = plt.figure(figsize=(w * s / 100, h * s / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.pcolormesh(X, Y, F, cmap=cmap, vmin=cl[0], vmax=cl[1], shading=a.shading, rasterized=True)
    ax.set_xlim(X[0], X[-1]); ax.set_ylim(Y[0], Y[-1])
    fig.savefig(png, dpi=100)
    plt.close(fig)
    print("wrote", png)


clims = []
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
        parts = ([XZ] if doxz else []) + ([XY] if doxy else []) + ([YZ] if doyz else [])
        lo, hi = min(p.min() for p in parts), max(p.max() for p in parts)
        cl = [-max(-lo, hi), max(-lo, hi)] if lo < 0 < hi else [lo, hi]
    cmap = cm_div if cl[0] < 0 < cl[1] else cm_seq
    clims.append((f, cl[0], cl[1], cmap.name))
    lab = f + ("'" if fl else "")
    if not a.annotated:
        if doxz:
            for k in range(XZ.shape[2]):
                image("%s_%s_xz%d.png" % (base, f, k + 1), xp, zp, XZ[:, :, k], lx, lz, cmap, cl)
        if doxy:
            image("%s_%s_xy.png" % (base, f), xp, y, XY, lx, 2.0, cmap, cl)
        if doyz:
            image("%s_%s_yz.png" % (base, f), zp, y, YZ, lz, 2.0, cmap, cl)
        continue
    # annotated: one figure per field, as tools/plot_slices.m
    nxz = XZ.shape[2]
    nrow = int(doxz) + int(doxy or doyz)
    ncol = max(nxz if doxz else 0, int(doxy) + int(doyz), 1)
    fig, ax = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow), squeeze=False)
    for axx in ax.flat:
        axx.set_visible(False)

    def show(axx, X, Y, F, xl, yl, title):
        axx.set_visible(True)
        m = axx.pcolormesh(X, Y, F, cmap=cmap, vmin=cl[0], vmax=cl[1], shading=a.shading)
        axx.set_aspect("equal"); axx.set_xlabel(xl); axx.set_ylabel(yl); axx.set_title(title, fontsize=9)
        fig.colorbar(m, ax=axx, shrink=0.8)
    r = 0
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
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print("wrote", png)

with open(base + "_clim.txt", "w") as o:
    o.write("# colour range of the pictures of %s (t = %.8g): field, min, max, colour map\n" % (os.path.basename(a.file), sc("time")))
    for f, lo, hi, n in clims:
        o.write("%-4s %15.7e %15.7e  %s\n" % (f, lo, hi, n))
for f, lo, hi, n in clims:
    print("  %-4s colour range %12.5e .. %12.5e (%s)" % (f, lo, hi, n))
