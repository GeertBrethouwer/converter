#!/usr/bin/env python3
"""Tests of rsspec (the output files have 9 significant digits: tolerances 1e-7) against an independent computation from the physical field (numpy).

    python3 tests/test_spectra.py [bin_dir]      (bin_dir with rsconv and rsspec; default: the repository)

A random HiCHAN restart (two scalars) is written; rsspec's output at a few y-positions is compared with
  * the variance <f'f'> of the physical field f(x, z) (inverse FFT of the coefficients);
  * Ex(kx) dkx = z-average of a_i |FFT_x(f)_i / nx|^2 (a_0 = 1, a_i = 2), the mean taken out;
  * Ez(kz) dkz = x-average of the folded |FFT_z(f)_m / nz|^2, the mean taken out;
  * E2: the folded |c|^2;
and also: the same spectra from the FDCHAN file (rsconv h2f), the average of two files, --sym, --yplus/--y.
"""
import os, subprocess, sys, tempfile
import numpy as np

B = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
rng = np.random.default_rng(5)
ok = True


def check(name, err, tol):
    global ok
    good = err < tol
    ok &= good
    print(f"  {name:66s} {err:9.2e} (< {tol:.0e}) {'ok' if good else 'WRONG'}")


def run(prog, *args):
    r = subprocess.run([os.path.join(B, prog), *map(str, args)], capture_output=True, text=True)
    if r.returncode:
        print(r.stdout, r.stderr)
        raise SystemExit(prog + " failed")
    return r.stdout


def spectrum(nkx, nz, ny):
    c = (rng.normal(size=(ny, nz, nkx)) + 1j * rng.normal(size=(ny, nz, nkx))) / (1 + np.arange(nkx))
    c0 = c[:, :, 0].copy()
    for m in range(nz):
        c[:, m, 0] = 0.5 * (c0[:, m] + np.conj(c0[:, (-m) % nz]))
    c[:, nz // 2, :] = 0
    c[:, 0, 0] += 3.0  # a mean
    return c


def write_hichan(prefix, fields, nx, nz, lx, lz, t):
    ny = fields[0].shape[0]
    with open(prefix + ".bin", "wb") as f:
        for c in fields:
            f.write(np.ascontiguousarray(c, dtype=np.complex128).tobytes())
    with open(prefix + ".hdr", "w") as f:
        f.write("# nx ny nz nscal lx lz time istep dt grid_type grid_alpha fd_order\n")
        f.write(f" {nx} {ny} {nz} {len(fields) - 3}\n {lx!r} {lz!r}\n {t!r} 0 0.001\n tanh 2.0 6\n")


def load(fn):
    return np.loadtxt(fn)


def physical(c, nx):
    """f(z, x) on the nz x nx grid from the coefficients c[m, i], i < nx/2"""
    nz, nkx = c.shape
    H = np.zeros((nz, nx // 2 + 1), dtype=complex)
    H[:, :nkx] = c
    return np.fft.irfft2(H, s=(nz, nx)) * nx * nz


nx, ny, nz, lx, lz = 24, 17, 20, 2 * np.pi, np.pi
dkx, dkz = 2 * np.pi / lx, 2 * np.pi / lz
names = ["u", "v", "w", "s1", "s2"]
A = [spectrum(nx // 2, nz, ny) for _ in names]
Bf = [spectrum(nx // 2, nz, ny) for _ in names]
print("rsspec tests")
with tempfile.TemporaryDirectory() as d:
    os.chdir(d)
    write_hichan("ha", A, nx, nz, lx, lz, 1.0)
    write_hichan("hb", Bf, nx, nz, lx, lz, 2.0)
    js = [0, 3, 8, 13]
    run("rsspec", "--j", ",".join(map(str, js)), "--out", "sa", "ha")
    ev, ex, ez, e2 = 0, 0, 0, 0
    for j in js:
        X, Z, T = load(f"sa/spec_x_j{j:04d}.dat"), load(f"sa/spec_z_j{j:04d}.dat"), load(f"sa/spec_2d_j{j:04d}.dat")
        for q, c in enumerate(A):
            f = physical(c[j], nx)
            fp = f - f.mean()
            var = np.mean(fp ** 2)
            ev = max(ev, abs(np.sum(X[:, 2 + q]) * dkx - var) / var)
            # streamwise: per z-line FFT in x
            Fx = np.fft.rfft(fp, axis=1) / nx
            sx = np.mean(np.abs(Fx[:, : nx // 2]) ** 2, axis=0) * np.where(np.arange(nx // 2) == 0, 1, 2)
            ex = max(ex, np.max(np.abs(X[:, 2 + q] * dkx - sx)) / var)
            # spanwise: per x-line FFT in z, folded
            Fz = np.abs(np.fft.fft(fp, axis=0) / nz) ** 2
            Fz = np.mean(Fz, axis=1)
            sz = np.array([Fz[0]] + [Fz[m] + Fz[nz - m] for m in range(1, nz // 2)] + [Fz[nz // 2]])
            ez = max(ez, np.max(np.abs(Z[:, 2 + q] * dkz - sz)) / var)
            # 2D
            E = np.abs(c[j]) ** 2
            E[0, 0] = 0
            fold = np.zeros((nx // 2, nz // 2 + 1))
            for m in range(nz):
                fold[:, min(m, nz - m)] += E[m, :] * np.where(np.arange(nx // 2) == 0, 1, 2)
            e2 = max(e2, np.max(np.abs(T[:, 4 + q].reshape(nx // 2, nz // 2 + 1) * dkx * dkz - fold)) / var)
    check("variance: sum Ex dkx = <f'f'> of the physical field", ev, 1e-7)
    check("Ex against FFTs in x of the physical field", ex, 1e-7)
    check("Ez against FFTs in z of the physical field", ez, 1e-7)
    check("E2 against the folded |c|^2", e2, 1e-7)
    # the same from the FDCHAN file
    run("rsconv", "h2f", "ha", "fa.bin")
    run("rsspec", "--j", ",".join(map(str, js)), "--out", "sf", "fa.bin")
    e = max(np.max(np.abs(load(f"sf/spec_2d_j{j:04d}.dat")[:, 4:] - load(f"sa/spec_2d_j{j:04d}.dat")[:, 4:]))
            / np.max(load(f"sa/spec_2d_j{j:04d}.dat")[:, 4:]) for j in js)
    check("FDCHAN file (rsconv h2f) gives the same spectra", e, 1e-7)
    run("rsconv", "f2f", "fa.bin", "fa_sp.bin", "--single")
    run("rsspec", "--j", ",".join(map(str, js)), "--out", "ss", "fa_sp.bin")
    e = max(np.max(np.abs(load(f"ss/spec_x_j{j:04d}.dat")[:, 2:] - load(f"sa/spec_x_j{j:04d}.dat")[:, 2:]))
            / np.max(load(f"sa/spec_x_j{j:04d}.dat")[:, 2:]) for j in js)
    check("single-precision FDCHAN file: single-precision round-off", e, 1e-6)
    # average of two files, also mixing the formats
    run("rsspec", "--j", ",".join(map(str, js)), "--out", "sb", "hb")
    run("rsspec", "--j", ",".join(map(str, js)), "--out", "sab", "fa.bin", "hb")
    e = max(np.max(np.abs(load(f"sab/spec_2d_j{j:04d}.dat")[:, 4:] - 0.5 * (load(f"sa/spec_2d_j{j:04d}.dat")[:, 4:]
            + load(f"sb/spec_2d_j{j:04d}.dat")[:, 4:]))) / np.max(load(f"sa/spec_2d_j{j:04d}.dat")[:, 4:]) for j in js)
    check("two files: the average of their spectra", e, 1e-7)
    # --sym: the average of j and ny-1-j
    run("rsspec", "--j", "3", "--sym", "--out", "sy", "ha")
    run("rsspec", "--j", str(ny - 1 - 3), "--out", "sm", "ha")
    e = np.max(np.abs(load("sy/spec_2d_j0003.dat")[:, 4:] - 0.5 * (load("sa/spec_2d_j0003.dat")[:, 4:] +
               load(f"sm/spec_2d_j{ny - 4:04d}.dat")[:, 4:]))) / np.max(load("sa/spec_2d_j0003.dat")[:, 4:])
    check("--sym: the average of j and ny-1-j", e, 1e-7)
    # positions: --y and --yplus pick the nearest grid points of the tanh grid from the HiCHAN header
    yg = np.tanh(2.0 * (np.arange(ny) / (ny - 1) - 0.5)) / np.tanh(1.0)
    out = run("rsspec", "--yplus", "10,100", "--retau", 180, "--out", "sp", "ha")
    want = sorted({int(np.argmin(np.abs(yg - (-1 + p / 180)))) for p in (10, 100)})
    got = sorted(int(f[8:12]) for f in os.listdir("sp") if f.startswith("spec_x"))
    check("--yplus 10,100 --retau 180: nearest grid points", float(want != got), 0.5)
    run("rsspec", "--y", "0", "--out", "sc", "ha")
    check("--y 0: the centre point", float(not os.path.exists(f"sc/spec_x_j{ny // 2:04d}.dat")), 0.5)
    sm = np.loadtxt("sa/summary.dat")
    means = np.array([[A[q][j, 0, 0].real for q in range(len(names))] for j in js])
    check("summary.dat: the means of the fields", np.max(np.abs(sm[:, 3 + len(names):] - means)), 1e-7)

print("rsspec: OK" if ok else "rsspec: FAIL")
sys.exit(0 if ok else 1)
