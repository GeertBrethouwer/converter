#!/usr/bin/env python3
"""
2D slices of u, v, w and the scalars from a HiCHAN or an FDCHAN restart file (single or double precision) for
visualisation, saved as a MATLAB .mat file (as HiCHAN's tools/restart_slices.py: read with load() or plot with
tools/plot_slices.m) or, with -o <name>.npz or without scipy, as a numpy .npz file (tools/plot_slices.py).

Both codes hold Fourier coefficients in x; HiCHAN also in z, FDCHAN the values at the points z_k = k lz/nz:
    f(x, z, y_j) = sum_{kx, kz} f^(kx, kz, y_j) exp(i (kx x + kz z)),  kx >= 0 (Hermitian symmetry).
The fields are evaluated exactly (inverse FFTs) at the collocation points y_j of the wall-normal grid, optionally on
a finer x-z grid by zero padding of the Fourier modes (--upsample; for FDCHAN the z-points are first transformed
to their Fourier modes; at --upsample 1 the FDCHAN z-points are reproduced exactly).

Slices:
    xz-planes at given wall-normal positions (nearest grid points)    -> arrays (Nz, Nx, n_y)
    xy-plane  at spanwise position z                                   -> array  (ny, Nx)
    yz-plane  at streamwise position x                                 -> array  (ny, Nz)
The x (z) positions are rounded to the nearest point of the output grid.

Examples:
    python3 slices.py data/restart                     # HiCHAN (name of .bin/.hdr); y = -0.9, 0; z = 0; x = 0
    python3 slices.py data/restart.bin --grid tanh 4.5 # FDCHAN (the y-grid is not in the file)
    python3 slices.py data/samples/restart_t000105.0012.bin --grid tanh 4.5 --yplus 15 100 --retau 1000 \\
            --z 1.0 --x 2.0 --upsample 2 --fluct -o slices.mat

The y-grid: HiCHAN files contain it (grid_type 'file': --grid-file); for FDCHAN give --grid TYPE ALPHA (tanh,
uniform, cosine; FDCHAN's default is tanh 3.0) or --grid-file. Needs numpy; scipy for .mat files.
"""
import argparse
import os
import sys
import numpy as np


# ------------------------------------------------------------------------------------------- reading
def strip(name):
    return name[:-4] if name.endswith(('.bin', '.hdr')) else name


def is_hichan(name):
    return os.path.exists(strip(name) + '.hdr')


def read_hichan(name):
    base = strip(name)
    with open(base + '.hdr') as f:
        lines = [l for l in f.read().splitlines() if l.strip() and not l.lstrip().startswith('#')]
    nx, ny, nz, nscal = (int(t) for t in lines[0].split()[:4])
    lx, lz = (float(t.replace('D', 'E')) for t in lines[1].split()[:2])
    t = lines[2].split()
    h = dict(fmt='hichan', bin=base + '.bin', nx=nx, ny=ny, nz=nz, nscal=nscal, lx=lx, lz=lz,
             time=float(t[0].replace('D', 'E')), step=int(t[1]), dt=float(t[2].replace('D', 'E')), prec=8,
             grid_type=None, grid_alpha=None)
    if len(lines) > 3:
        g = lines[3].split()
        h['grid_type'], h['grid_alpha'] = g[0], float(g[1].replace('D', 'E'))
    need = (3 + nscal) * (nx // 2) * nz * ny * 16
    if os.path.getsize(h['bin']) < need:
        sys.exit('%s is smaller than its header says (%d < %d bytes)' % (h['bin'], os.path.getsize(h['bin']), need))
    return h


def read_fdchan(name):
    hd = np.fromfile(name, dtype='<i4', count=4)
    d = np.fromfile(name, dtype='<f8', count=3, offset=16)
    nx, ny, nz, prec = (int(v) for v in hd)
    if prec not in (4, 8) or nx <= 0 or ny <= 0 or nz <= 0 or nx % 2:
        sys.exit('%s does not look like an FDCHAN restart file (header %s)' % (name, hd))
    field = (nx // 2) * nz * ny * 2 * prec
    nf = (os.path.getsize(name) - 40) // field
    if nf < 3 or (os.path.getsize(name) - 40) % field:
        sys.exit('%s: the size is not a whole number (>= 3) of fields' % name)
    return dict(fmt='fdchan', bin=name, nx=nx, ny=ny, nz=nz, nscal=nf - 3, lx=d[0], lz=d[1], time=d[2], step=-1,
                dt=float('nan'), prec=prec, grid_type=None, grid_alpha=None)


def wall_normal_grid(ny, gtype, alpha, grid_file=None):
    """as both codes: y[0] = -1, y[ny-1] = 1"""
    if grid_file:
        y = np.loadtxt(grid_file).ravel()[:ny].astype(float)
        if y.size < ny:
            sys.exit('the grid file %s has fewer than ny = %d values' % (grid_file, ny))
        return y
    xi = np.arange(ny) / (ny - 1.0)
    if gtype == 'tanh':
        y = np.tanh(alpha * (xi - 0.5)) / np.tanh(0.5 * alpha)
    elif gtype == 'uniform':
        y = -1.0 + 2.0 * xi
    elif gtype == 'cosine':
        y = -np.cos(np.pi * xi)
    else:
        sys.exit("y-grid '%s': give --grid TYPE ALPHA (tanh, uniform, cosine) or --grid-file" % gtype)
    y[0], y[-1] = -1.0, 1.0
    return y


class Field:
    """memory-mapped access to one field; plane(j) gives the coefficients c[kx, m] (m: kz index in FFT order)"""

    def __init__(self, h, ifield):
        nxh, nz, ny = h['nx'] // 2, h['nz'], h['ny']
        self.h = h
        if h['fmt'] == 'hichan':
            self.data = np.memmap(h['bin'], dtype='<c16', mode='r', offset=ifield * nxh * nz * ny * 16,
                                  shape=(nxh, nz, ny), order='F')
        else:
            dt = '<c16' if h['prec'] == 8 else '<c8'
            self.data = np.memmap(h['bin'], dtype=dt, mode='r', offset=40 + ifield * nxh * nz * ny * 2 * h['prec'],
                                  shape=(ny, nz, nxh))

    def plane(self, j):
        if self.h['fmt'] == 'hichan':
            return np.array(self.data[:, :, j], dtype=complex)              # (nx/2, nz) coefficients
        p = np.array(self.data[j], dtype=complex).T                          # (nx/2, nz) values at z_k
        return np.fft.fft(p, axis=1) / self.h['nz']                          # -> coefficients


def zmodes(c, nz, Nz):
    """c[kx, m] (m in FFT order, any nz) -> the same modes in an Nz-point FFT array (zero padding; a Nyquist mode of
    an even nz is split between +nz/2 and -nz/2, which keeps the field real)"""
    if Nz == nz:
        return c
    g = np.zeros((c.shape[0], Nz), dtype=complex)
    npos = (nz - 1) // 2
    g[:, :npos + 1] = c[:, :npos + 1]
    g[:, Nz - npos:] = c[:, nz - npos:]
    if nz % 2 == 0:
        g[:, nz // 2] += 0.5 * c[:, nz // 2]
        g[:, Nz - nz // 2] += 0.5 * c[:, nz // 2]
    return g


def to_physical(c, nx, nz, Nx, Nz):
    """c[kx, m] (kx = 0..nx/2-1, m in FFT order, any nz) -> f(x, z) of shape (Nz, Nx) on Nx x Nz points"""
    g = np.fft.ifft(zmodes(c, nz, Nz), axis=1) * Nz                         # sum over kz -> f^(kx, z)
    G = np.zeros((Nx // 2 + 1, Nz), dtype=complex)
    G[:nx // 2, :] = g
    G[0, :] = G[0, :].real                                                   # kx = 0 is real in z
    return (np.fft.irfft(G, n=Nx, axis=0) * Nx).T                            # (Nz, Nx): rows z, columns x


def lines(c, nx, nz, Nx, Nz, ix, iz):
    """only f(x, z_iz) (length Nx) and f(x_ix, z) (length Nz): one sum over kz, one over kx and two 1D FFTs
    instead of the 2D transform of the whole plane"""
    g = zmodes(c, nz, Nz)
    h = g @ np.exp(2j * np.pi * np.fft.fftfreq(Nz, 1.0 / Nz) * iz / Nz)     # f^(kx, z_iz)
    H = np.zeros(Nx // 2 + 1, dtype=complex)
    H[:nx // 2] = h
    H[0] = H[0].real
    fx = np.fft.irfft(H, n=Nx) * Nx
    a = np.r_[1.0, np.full(nx // 2 - 1, 2.0)] * np.exp(2j * np.pi * np.arange(nx // 2) * ix / Nx)
    fz = np.real(np.fft.ifft(a @ g) * Nz)                                    # sum_kx a_i f^(kx, z) e^(i kx x_ix)
    return fx, fz


def even(n):
    n = int(round(n))
    return n + (n % 2)


# ------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('restart', help='HiCHAN restart (name of .bin/.hdr) or FDCHAN restart file (.bin)')
    ap.add_argument('-o', '--output', default=None,
                    help='output file, .mat (needs scipy) or .npz (default: <restart>_slices.mat, .npz without scipy)')
    ap.add_argument('--y', type=float, nargs='+', default=None, help='y of the xz-planes (default -0.9 0)')
    ap.add_argument('--yplus', type=float, nargs='+', default=None,
                    help='xz-planes at y+ from the bottom wall (needs --retau, or --re and --utau)')
    ap.add_argument('--retau', type=float, default=None, help='Re_tau (for --yplus)')
    ap.add_argument('--re', type=float, default=None, help='Reynolds number re (for --yplus with --utau)')
    ap.add_argument('--utau', type=float, default=None, help='friction velocity (for --yplus with --re)')
    ap.add_argument('--z', type=float, default=0.0, help='z of the xy-plane (default 0)')
    ap.add_argument('--x', type=float, default=0.0, help='x of the yz-plane (default 0)')
    ap.add_argument('--upsample', type=float, default=1.0,
                    help='output grid Nx = upsample*nx, Nz = upsample*nz (exact Fourier interpolation)')
    ap.add_argument('--fields', nargs='+', default=None, help='subset, e.g. u v s1 (default: all)')
    ap.add_argument('--fluct', action='store_true',
                    help="subtract the x-z mean (kx = kz = 0 mode): fluctuations u', v', ..., s'")
    ap.add_argument('--grid', nargs=2, metavar=('TYPE', 'ALPHA'), default=None,
                    help='y-grid (tanh, uniform, cosine) and stretching; needed for FDCHAN (default: HiCHAN header)')
    ap.add_argument('--grid-file', default=None, help='y-grid from a file (ny values)')
    a = ap.parse_args()

    h = read_hichan(a.restart) if is_hichan(a.restart) else read_fdchan(a.restart)
    nx, ny, nz = h['nx'], h['ny'], h['nz']
    if a.grid:
        gtype, alpha = a.grid[0], float(a.grid[1])
    elif h['grid_type'] is not None and h['grid_type'] != 'file':
        gtype, alpha = h['grid_type'], h['grid_alpha']
    elif h['grid_type'] == 'file' and not a.grid_file:
        sys.exit("grid_type = 'file': give the grid file with --grid-file")
    elif not a.grid_file:
        sys.exit('FDCHAN restart files do not contain the y-grid: give --grid TYPE ALPHA (as grid_type, grid_alpha '
                 'of the input; FDCHAN default tanh 3.0) or --grid-file')
    else:
        gtype, alpha = None, None
    y = wall_normal_grid(ny, gtype, alpha, a.grid_file)
    Nx, Nz = even(a.upsample * nx), int(round(a.upsample * nz))
    if a.upsample != 1:
        Nz = even(Nz)
    x = h['lx'] * np.arange(Nx) / Nx
    z = h['lz'] * np.arange(Nz) / Nz

    names = ['u', 'v', 'w'] + ['s%d' % (n + 1) for n in range(h['nscal'])]
    sel = names if a.fields is None else a.fields
    for f in sel:
        if f not in names:
            sys.exit('unknown field %s (available: %s)' % (f, ' '.join(names)))

    if a.yplus is not None:
        if a.retau is not None:
            ytarget = [-1.0 + yp / a.retau for yp in a.yplus]
        elif a.re is not None and a.utau is not None:
            ytarget = [-1.0 + yp / (a.utau * a.re) for yp in a.yplus]
        else:
            sys.exit('--yplus needs --retau, or --re and --utau')
    else:
        ytarget = a.y if a.y is not None else [-0.9, 0.0]
    jxz = [int(np.argmin(np.abs(y - yt))) for yt in ytarget]
    iz = int(round(a.z / h['lz'] * Nz)) % Nz
    ix = int(round(a.x / h['lx'] * Nx)) % Nx

    out = dict(x=x, y=y, z=z, time=h['time'], step=h['step'], lx=h['lx'], lz=h['lz'],
               nx=nx, ny=ny, nz=nz, Nx=Nx, Nz=Nz,
               y_xz=y[jxz], j_xz=np.array(jxz) + 1, z_xy=z[iz], x_yz=x[ix],
               fluct=int(a.fluct), fields=np.array(sel, dtype=object), source=h['fmt'])
    print('%s restart (%s precision): nx ny nz = %d %d %d, %d scalar(s), t = %.6g%s'
          % ('HiCHAN' if h['fmt'] == 'hichan' else 'FDCHAN', 'double' if h['prec'] == 8 else 'single', nx, ny, nz,
             h['nscal'], h['time'], ', step %d' % h['step'] if h['step'] >= 0 else ''))
    print('output grid Nx x Nz = %d x %d; xz-planes at y =' % (Nx, Nz), ' '.join('%.5f' % y[j] for j in jxz),
          '; xy-plane at z = %.5f; yz-plane at x = %.5f' % (z[iz], x[ix]))
    for name in sel:
        fld = Field(h, names.index(name))
        xz = np.zeros((Nz, Nx, len(jxz)))
        xy = np.zeros((ny, Nx))
        yz = np.zeros((ny, Nz))
        for j in range(ny):
            c = fld.plane(j)
            if a.fluct:
                c[0, 0] = 0.0
            if j in jxz:  # the whole plane
                p = to_physical(c, nx, nz, Nx, Nz)
                xy[j, :] = p[iz, :]
                yz[j, :] = p[:, ix]
                for m, jj in enumerate(jxz):
                    if jj == j:
                        xz[:, :, m] = p
            else:         # only the two lines
                xy[j, :], yz[j, :] = lines(c, nx, nz, Nx, Nz, ix, iz)
        out['xz_' + name] = xz
        out['xy_' + name] = xy
        out['yz_' + name] = yz
        print('  %-3s  min %12.5e  max %12.5e' % (name, min(xy.min(), yz.min(), xz.min()),
                                                  max(xy.max(), yz.max(), xz.max())))
    fout = a.output
    try:
        from scipy.io import savemat
    except ImportError:
        savemat = None
    if fout is None:
        fout = strip(a.restart) + ('_slices.mat' if savemat else '_slices.npz')
        if not savemat:
            print('(no scipy: writing a .npz file; plot with tools/plot_slices.py, or install scipy for .mat)')
    if fout.endswith('.npz'):
        npz = {k: np.asarray(v) for k, v in out.items()}
        npz['fields'] = np.array(sel, dtype=str)  # (no object arrays: np.load without allow_pickle)
        np.savez_compressed(fout, **npz)
    else:
        if savemat is None:
            sys.exit('scipy is required to write .mat files (pip install --user scipy), or use -o <name>.npz')
        savemat(fout, out, do_compression=True)
    print('wrote', fout)


if __name__ == '__main__':
    main()
