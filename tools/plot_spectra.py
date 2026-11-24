#!/usr/bin/env python3
"""Plots of the spectra written by rsspec (python3 with numpy and matplotlib).

    python3 tools/plot_spectra.py <rsspec output dir> [--retau R --utau U] [--fields u,v,w,s1]

Writes into the directory:
  spec_1d.png   premultiplied streamwise kx Ex(kx) and spanwise kz Ez(kz) against the wavelength, a panel per field,
                a line per y-position
  spec_2d.png   premultiplied 2D spectra kx kz E2(kx, kz) against lambda_x, lambda_z, a panel per (position, field):
                contours at 0.125, 0.25, 0.5, 0.75 of the maximum, as usual for these spectra
--retau R: wavelengths in wall units (lambda+ = lambda Re_tau); --utau U: velocity spectra divided by U^2 (in these
codes u_tau = Re_tau / Re for velocities in units of the bulk or centreline velocity).
"""
import argparse, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--retau", type=float, default=0)
ap.add_argument("--utau", type=float, default=0)
ap.add_argument("--fields", default="")
a = ap.parse_args()


def info(fn):
    with open(fn) as f:
        lines = [l for l in f if l.startswith("#")]
    names = lines[-1].split()[3:] if "2d" in fn else lines[-1].split()[3:]
    m = re.search(r"y\+ = ([0-9.eE+-]+)", "".join(lines))
    y = re.search(r"y = ([0-9.eE+-]+)", "".join(lines))
    lab = f"y+ = {float(m.group(1)):.0f}" if m else (f"y = {float(y.group(1)):.3f}" if y else os.path.basename(fn))
    return [n.split("_", 1)[1] for n in names], lab


js = sorted(int(re.search(r"j(\d+)", f).group(1)) for f in glob.glob(os.path.join(a.dir, "spec_x_j*.dat")))
if not js:
    raise SystemExit("no spec_x_j*.dat in " + a.dir)
names, _ = info(os.path.join(a.dir, f"spec_x_j{js[0]:04d}.dat"))
sel = [n for n in names if not a.fields or n in a.fields.split(",")]
L = a.retau if a.retau > 0 else 1.0
unit = "+" if a.retau > 0 else ""


def scale(n):
    return 1.0 / a.utau ** 2 if (a.utau > 0 and n in ("u", "v", "w")) else 1.0


# 1D
fig, ax = plt.subplots(2, len(sel), figsize=(3.2 * len(sel), 5.6), squeeze=False)
for j in js:
    X = np.loadtxt(os.path.join(a.dir, f"spec_x_j{j:04d}.dat"))
    Z = np.loadtxt(os.path.join(a.dir, f"spec_z_j{j:04d}.dat"))
    _, lab = info(os.path.join(a.dir, f"spec_x_j{j:04d}.dat"))
    for c, n in enumerate(sel):
        q = names.index(n)
        ax[0, c].semilogx(X[1:, 1] * L, X[1:, 0] * X[1:, 2 + q] * scale(n), label=lab)
        ax[1, c].semilogx(Z[1:, 1] * L, Z[1:, 0] * Z[1:, 2 + q] * scale(n), label=lab)
for c, n in enumerate(sel):
    s = "+" if (a.utau > 0 and n in ("u", "v", "w")) else ""
    ax[0, c].set_title(n)
    ax[0, c].set_xlabel(f"$\\lambda_x^{{{unit}}}$")
    ax[1, c].set_xlabel(f"$\\lambda_z^{{{unit}}}$")
    ax[0, c].set_ylabel(f"$k_x E_x^{{{s}}}$")
    ax[1, c].set_ylabel(f"$k_z E_z^{{{s}}}$")
for x in ax.flat:
    x.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(a.dir, "spec_1d.png"), dpi=150)

# 2D
fig, ax = plt.subplots(len(js), len(sel), figsize=(3.0 * len(sel), 2.8 * len(js)), squeeze=False)
for r, j in enumerate(js):
    T = np.loadtxt(os.path.join(a.dir, f"spec_2d_j{j:04d}.dat"))
    _, lab = info(os.path.join(a.dir, f"spec_2d_j{j:04d}.dat"))
    nkx = len(np.unique(T[:, 0]))
    T = T.reshape(nkx, -1, T.shape[1])[1:, 1:]  # without kx = 0 and kz = 0 (infinite wavelength)
    for c, n in enumerate(sel):
        q = names.index(n)
        P = T[:, :, 0] * T[:, :, 1] * T[:, :, 4 + q]
        if P.max() > 0:
            ax[r, c].contourf(T[:, :, 2] * L, T[:, :, 3] * L, P / P.max(), levels=[0.125, 0.25, 0.5, 0.75, 1.0],
                              cmap="Blues")
            ax[r, c].contour(T[:, :, 2] * L, T[:, :, 3] * L, P / P.max(), levels=[0.125, 0.25, 0.5, 0.75],
                             colors="k", linewidths=0.5)
        ax[r, c].set_xscale("log")
        ax[r, c].set_yscale("log")
        ax[r, c].set_title(f"{n}, {lab}", fontsize=9)
        ax[r, c].set_xlabel(f"$\\lambda_x^{{{unit}}}$")
        ax[r, c].set_ylabel(f"$\\lambda_z^{{{unit}}}$")
fig.tight_layout()
for x in ax.flat:
    x.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    x.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.tight_layout()
fig.savefig(os.path.join(a.dir, "spec_2d.png"), dpi=150)
print("written", os.path.join(a.dir, "spec_1d.png"), "and", os.path.join(a.dir, "spec_2d.png"))
