#!/usr/bin/env python3
"""Tests of tools/slices.py (python3 with numpy and scipy).

    python3 tests/test_slices.py [bin_dir]      (bin_dir with rsconv; default: the repository)

  * HiCHAN file: the slices against a direct evaluation of the Fourier series at the output points
    (also on a finer grid, --upsample 1.5), and with --fluct
  * FDCHAN file of the same field (rsconv h2f): the same slices; single precision: round-off
  * FDCHAN file with an odd nz and kz-Nyquist content (written directly): at --upsample 1 the stored z-point values
    are reproduced exactly; --upsample 2 at the even points gives the same values
  * --yplus/--retau, --y, --x, --z positions; .npz output; tools/plot_slices.py runs (if matplotlib is there)
"""
import os, subprocess, sys, tempfile
import numpy as np
try:
    from scipy.io import loadmat
except ImportError:
    print("slices: skipped (the test needs scipy: pip install --user scipy)")
    sys.exit(0)

B = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
T = os.path.join(B, "tools")
rng = np.random.default_rng(7)
ok = True


def check(name, err, tol):
    global ok
    good = err < tol
    ok &= good
    print(f"  {name:66s} {err:9.2e} (< {tol:.0e}) {'ok' if good else 'WRONG'}")


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout, r.stderr)
        raise SystemExit("failed: " + " ".join(cmd))
    return r.stdout


def spectrum(nkx, nz, ny):
    c = (rng.normal(size=(ny, nz, nkx)) + 1j * rng.normal(size=(ny, nz, nkx))) / (1 + np.arange(nkx))
    c0 = c[:, :, 0].copy()
    for m in range(nz):
        c[:, m, 0] = 0.5 * (c0[:, m] + np.conj(c0[:, (-m) % nz]))
    if nz % 2 == 0:
        c[:, nz // 2, :] = 0
    return c


def series(c, lx, lz, x, z):
    """direct evaluation: f(x, z) = sum over kx >= 0 (the conjugates of kx > 0 included) and all kz of c[m, i]"""
    nz, nkx = c.shape
    kx = 2 * np.pi * np.arange(nkx) / lx
    m = np.fft.fftfreq(nz, 1.0 / nz)
    kz = 2 * np.pi * m / lz
    f = np.zeros((z.size, x.size))
    for ii in range(nkx):
        w = 1.0 if ii == 0 else 2.0
        ez = np.exp(1j * np.outer(z, kz)) @ c[:, ii]          # sum over kz at each z
        f += w * np.real(np.outer(ez, np.exp(1j * kx[ii] * x)))
    return f


def write_hichan(prefix, fields, nx, nz, lx, lz, t):
    ny = fields[0].shape[0]
    with open(prefix + ".bin", "wb") as f:
        for c in fields:
            f.write(np.ascontiguousarray(c, dtype=np.complex128).tobytes())
    with open(prefix + ".hdr", "w") as f:
        f.write("# nx ny nz nscal lx lz time istep dt grid_type grid_alpha fd_order\n")
        f.write(f" {nx} {ny} {nz} {len(fields) - 3}\n {lx!r} {lz!r}\n {t!r} 77 0.001\n tanh 2.0 6\n")


def write_fdchan(fn, fields, nx, nz, lx, lz, t, prec=8):
    """fields: physical z-values [y, k, i] (complex, kx = 0 real)"""
    ny = fields[0].shape[0]
    with open(fn, "wb") as f:
        np.array([nx, ny, nz, prec], dtype="<i4").tofile(f)
        np.array([lx, lz, t], dtype="<f8").tofile(f)
        for v in fields:
            v.astype(np.complex128 if prec == 8 else np.complex64).tofile(f)


nx, ny, nz, lx, lz = 16, 11, 12, 2 * np.pi, np.pi
names = ["u", "v", "w", "s1"]
C = [spectrum(nx // 2, nz, ny) for _ in names]
yg = np.tanh(2.0 * (np.arange(ny) / (ny - 1) - 0.5)) / np.tanh(1.0)
print("slices tests")
with tempfile.TemporaryDirectory() as d:
    os.chdir(d)
    write_hichan("h", C, nx, nz, lx, lz, 4.5)
    for up in (1, 1.5):
        run([sys.executable, os.path.join(T, "slices.py"), "h", "--y", "-0.5", "0.2", "--x", "1.0", "--z", "0.7",
             "--upsample", str(up), "-o", "s.mat"])
        S = loadmat("s.mat", squeeze_me=True)
        x, z = S["x"], S["z"]
        e = 0
        for q, n in enumerate(names):
            jxz = [int(np.argmin(np.abs(yg - t))) for t in (-0.5, 0.2)]
            for k, j in enumerate(jxz):
                ref = series(C[q][j], lx, lz, x, z)
                e = max(e, np.max(np.abs(S["xz_" + n][:, :, k] - ref)) / np.max(np.abs(ref)))
            for j in range(ny):
                ref = series(C[q][j], lx, lz, x, np.array([float(S["z_xy"])]))[0]
                e = max(e, np.max(np.abs(S["xy_" + n][j] - ref)) / np.max(np.abs(ref)))
                ref = series(C[q][j], lx, lz, np.array([float(S["x_yz"])]), z)[:, 0]
                e = max(e, np.max(np.abs(S["yz_" + n][j] - ref)) / np.max(np.abs(ref)))
        check(f"HiCHAN, --upsample {up}: xz, xy, yz against the Fourier series", e, 1e-12)
    check("HiCHAN: x_yz, z_xy are output points nearest to --x 1.0, --z 0.7",
          float(abs(S["x_yz"] - 1.0) > lx / S["Nx"] / 2 + 1e-12 or abs(S["z_xy"] - 0.7) > lz / S["Nz"] / 2 + 1e-12), 0.5)
    run([sys.executable, os.path.join(T, "slices.py"), "h", "--y", "0.2", "--fluct", "-o", "sf.mat"])
    Sf = loadmat("sf.mat", squeeze_me=True)
    j = int(np.argmin(np.abs(yg - 0.2)))
    check("--fluct: the xz-plane has zero mean", np.max([abs(Sf["xz_" + n].mean()) for n in names]), 1e-12)
    ref = series(C[0][j], lx, lz, Sf["x"], Sf["z"]) - C[0][j][0, 0].real
    check("--fluct: the fluctuation", np.max(np.abs(Sf["xz_u"] - ref)), 1e-12)

    # FDCHAN file of the same field: the same slices
    run([os.path.join(B, "rsconv"), "h2f", "h", "f.bin"])
    run([os.path.join(B, "rsconv"), "f2f", "f.bin", "fs.bin", "--single"])
    for fn, tol, lab in (("f.bin", 1e-12, "double"), ("fs.bin", 1e-6, "single")):
        run([sys.executable, os.path.join(T, "slices.py"), fn, "--grid", "tanh", "2.0", "--y", "-0.5", "0.2",
             "--x", "1.0", "--z", "0.7", "--upsample", "1.5", "-o", "f.mat"])
        F = loadmat("f.mat", squeeze_me=True)
        e = max(np.max(np.abs(F[p + n] - S[p + n])) / np.max(np.abs(S[p + n])) for n in names for p in ("xz_", "xy_", "yz_"))
        check(f"FDCHAN file ({lab} precision) of the same field: the same slices", e, tol)
    r = subprocess.run([sys.executable, os.path.join(T, "slices.py"), "f.bin"], capture_output=True, text=True)
    check("FDCHAN file without --grid: refused (the grid is not in the file)", float(r.returncode == 0), 0.5)

    # FDCHAN, odd nz with Nyquist-like content (any z-values): upsample 1 reproduces the stored values
    nz2 = 9
    P = [rng.normal(size=(ny, nz2, nx // 2)) + 1j * rng.normal(size=(ny, nz2, nx // 2)) for _ in range(3)]
    for p in P:
        p[:, :, 0] = p[:, :, 0].real
    write_fdchan("odd.bin", P, nx, nz2, lx, lz, 1.0)
    run([sys.executable, os.path.join(T, "slices.py"), "odd.bin", "--grid", "tanh", "2.0", "--y", "0", "-o", "o1.mat"])
    O = loadmat("o1.mat", squeeze_me=True)
    j = int(np.argmin(np.abs(yg)))
    # f(x_i, z_k) = Re(c_0) + 2 Re(sum_{i>0} c_i exp(i kx x)) from the stored values c[k, i] at the z-points
    ref = np.array([[np.real(P[1][j][k, 0] + 2 * np.sum(P[1][j][k, 1:] * np.exp(2j * np.pi * np.arange(1, nx // 2) * i / nx)))
                     for i in range(nx)] for k in range(nz2)])
    check("FDCHAN odd nz = 9, --upsample 1: the stored z-point values", np.max(np.abs(O["xz_v"] - ref)) / np.max(np.abs(ref)), 1e-12)
    write_fdchan("even.bin", [np.concatenate([p, p[:, :1]], axis=1) for p in P], nx, nz2 + 1, lx, lz, 1.0)
    run([sys.executable, os.path.join(T, "slices.py"), "even.bin", "--grid", "tanh", "2.0", "--y", "0", "-o", "e1.mat"])
    run([sys.executable, os.path.join(T, "slices.py"), "even.bin", "--grid", "tanh", "2.0", "--y", "0", "--upsample", "2",
         "-o", "e2.mat"])
    E1, E2 = loadmat("e1.mat", squeeze_me=True), loadmat("e2.mat", squeeze_me=True)
    check("FDCHAN even nz with Nyquist content: --upsample 2 at the even points",
          np.max(np.abs(E2["xz_u"][::2, ::2] - E1["xz_u"])) / np.max(np.abs(E1["xz_u"])), 1e-12)

    # positions, npz, plot
    out = run([sys.executable, os.path.join(T, "slices.py"), "h", "--yplus", "10", "100", "--retau", "180", "-o", "p.npz"])
    Z = np.load("p.npz")
    want = sorted({int(np.argmin(np.abs(yg - (-1 + p / 180)))) for p in (10, 100)})
    check("--yplus 10 100 --retau 180 (.npz output): nearest grid points", float(sorted(Z["j_xz"] - 1) != want), 0.5)
    try:
        import matplotlib  # noqa
        run([sys.executable, os.path.join(T, "plot_slices.py"), "p.npz", "--fields", "u", "s1"])
        run([sys.executable, os.path.join(T, "plot_slices.py"), "s.mat", "--plane", "xz"])
        check("tools/plot_slices.py (.npz and .mat): the png files", float(not (os.path.exists("p_u.png") and
              os.path.exists("p_s1.png") and os.path.exists("s_w.png"))), 0.5)
    except ImportError:
        print("  (no matplotlib: tools/plot_slices.py not tested)")

print("slices: OK" if ok else "slices: FAIL")
sys.exit(0 if ok else 1)
