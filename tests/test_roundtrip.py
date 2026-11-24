#!/usr/bin/env python3
"""Round-trip and convention tests of rsconv with synthetic restart files (numpy).

    python3 tests/test_roundtrip.py [path/to/rsconv]

Checks, with random fields (Hermitian kz spectrum at kx = 0, zero Nyquist):
  * h2f: the FDCHAN values are f(z_k) = sum_m c_m exp(2 pi i m k / nz) (numpy's inverse FFT times nz)
  * h2f then f2h: the HiCHAN coefficients come back (round-off)
  * resolution: h2f --nz/--nx larger, then f2h back to the original: the coefficients come back
  * h2h --nz --nx: zero padding / truncation in spectral space
  * f2f --single, then --double: within single-precision round-off; f2f at the same precision: identical bytes
  * --nscal, the headers (time, step, dt, lx, lz) and odd nz for FDCHAN
"""
import os, subprocess, sys, tempfile
import numpy as np

RS = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "rsconv"))
rng = np.random.default_rng(3)
ok = True


def check(name, err, tol):
    global ok
    good = err < tol
    ok &= good
    print(f"  {name:62s} {err:9.2e} (< {tol:.0e}) {'ok' if good else 'WRONG'}")


def run(*args):
    r = subprocess.run([RS, *map(str, args)], capture_output=True, text=True)
    if r.returncode:
        print(r.stdout, r.stderr)
        raise SystemExit("rsconv failed: " + " ".join(map(str, args)))
    return r.stdout


def spectrum(nkx, nz, ny):
    """random coefficients c[y, m, i] (FFT order in m), Hermitian in m at i = 0, zero Nyquist"""
    c = rng.normal(size=(ny, nz, nkx)) + 1j * rng.normal(size=(ny, nz, nkx))
    c0 = c[:, :, 0].copy()
    for m in range(nz):
        c[:, m, 0] = 0.5 * (c0[:, m] + np.conj(c0[:, (-m) % nz]))
    if nz % 2 == 0:
        c[:, nz // 2, :] = 0
    return c


def write_hichan(prefix, fields, nx, nz, lx, lz, time, step, dt):
    ny = fields[0].shape[0]
    with open(prefix + ".bin", "wb") as f:
        for c in fields:  # Fortran (nkx, nz, ny): kx fastest, then m, then y = C order [y, m, i]
            f.write(np.ascontiguousarray(c, dtype=np.complex128).tobytes())
    with open(prefix + ".hdr", "w") as f:
        f.write("# nx ny nz nscal lx lz time istep dt grid_type grid_alpha fd_order\n")
        f.write(f" {nx} {ny} {nz} {len(fields) - 3}\n {lx!r} {lz!r}\n {time!r} {step} {dt!r}\n tanh 4.5 8\n")


def read_hichan(prefix):
    with open(prefix + ".hdr") as f:
        f.readline()
        nx, ny, nz, ns = map(int, f.readline().split())
        lx, lz = map(float, f.readline().split())
        t, st, dt = f.readline().split()
    a = np.fromfile(prefix + ".bin", dtype=np.complex128).reshape(3 + ns, ny, nz, nx // 2)
    return a, dict(nx=nx, ny=ny, nz=nz, lx=lx, lz=lz, time=float(t), step=int(st), dt=float(dt))


def read_fdchan(fname):
    h = np.fromfile(fname, dtype=np.int32, count=4)
    d = np.fromfile(fname, dtype=np.float64, count=3, offset=16)
    nx, ny, nz, prec = map(int, h)
    dtype = np.complex128 if prec == 8 else np.complex64
    a = np.fromfile(fname, dtype=dtype, offset=40)
    return a.reshape(-1, ny, nz, nx // 2), dict(nx=nx, ny=ny, nz=nz, prec=prec, lx=d[0], lz=d[1], time=d[2])


def resize(c, nkx_o, nz_o):
    """spectral resize as rsconv: kx < min, |m| <= (min(nz) - 1) // 2"""
    ny, nz, nkx = c.shape
    o = np.zeros((ny, nz_o, nkx_o), dtype=complex)
    mm = (min(nz, nz_o) - 1) // 2
    nk = min(nkx, nkx_o)
    for m in range(-mm, mm + 1):
        o[:, m % nz_o, :nk] = c[:, m % nz, :nk]
    return o


def rel(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


with tempfile.TemporaryDirectory() as d:
    os.chdir(d)
    nx, ny, nz, ns = 16, 9, 12, 2
    lx, lz, t, step, dt = 2 * np.pi, np.pi, 12.5, 4711, 1.6e-3
    C = [spectrum(nx // 2, nz, ny) for _ in range(3 + ns)]
    write_hichan("h0", C, nx, nz, lx, lz, t, step, dt)
    print("rsconv tests")
    run("info", "h0")

    # h2f: physical z-values from the coefficients
    run("h2f", "h0", "f1.bin")
    F, hf = read_fdchan("f1.bin")
    phys = [nz * np.fft.ifft(c, axis=1) for c in C]
    for p in phys:
        p[:, :, 0] = p[:, :, 0].real  # kx = 0 is real in z
    check("h2f: f(z_k) = sum_m c_m exp(2 pi i m k/nz)", max(rel(F[n], phys[n]) for n in range(3 + ns)), 1e-13)
    check("h2f: header (nx, ny, nz, prec, lx, lz, time)",
          float(not (hf["nx"] == nx and hf["ny"] == ny and hf["nz"] == nz and hf["prec"] == 8 and hf["lx"] == lx and
                     hf["lz"] == lz and hf["time"] == t)), 0.5)

    # f2h back
    run("f2h", "f1.bin", "h1", "--step", step, "--dt", dt)
    H, hh = read_hichan("h1")
    check("h2f + f2h: coefficients back", max(rel(H[n], C[n]) for n in range(3 + ns)), 1e-13)
    check("f2h: header (nx, ny, nz, lx, lz, time, step, dt)",
          float(not (hh["nx"] == nx and hh["nz"] == nz and hh["lx"] == lx and hh["lz"] == lz and hh["time"] == t and
                     hh["step"] == step and hh["dt"] == dt)), 0.5)

    # resolution: finer FDCHAN grid in z (odd), more x-modes, then back
    run("h2f", "h0", "f2.bin", "--nz", 27, "--nx", 24)
    F2, _ = read_fdchan("f2.bin")
    z27 = [27 * np.fft.ifft(resize(c, 12, 27), axis=1) for c in C]
    check("h2f --nz 27 --nx 24: spectral interpolation", max(rel(F2[n], z27[n]) for n in range(3 + ns)), 1e-13)
    run("f2h", "f2.bin", "h2", "--nz", nz, "--nx", nx)
    H2, _ = read_hichan("h2")
    check("... then f2h --nz 12 --nx 16: coefficients back", max(rel(H2[n], C[n]) for n in range(3 + ns)), 1e-13)

    # h2h: padding and truncation
    run("h2h", "h0", "h3", "--nz", 20, "--nx", 10)
    H3, hh3 = read_hichan("h3")
    check("h2h --nz 20 --nx 10 (pad z, truncate x)", max(rel(H3[n], resize(C[n], 5, 20)) for n in range(3 + ns)), 1e-15)
    check("h2h keeps step and dt", float(not (hh3["step"] == step and hh3["dt"] == dt)), 0.5)

    # precision
    run("f2f", "f1.bin", "f_sp.bin", "--single")
    run("f2f", "f_sp.bin", "f_dp.bin", "--double")
    Fs, hs = read_fdchan("f_sp.bin")
    Fd, hd = read_fdchan("f_dp.bin")
    check("f2f --single: precision in the header is 4", float(hs["prec"] != 4), 0.5)
    check("f2f --single then --double: single-precision round-off", max(rel(Fd[n], F[n]) for n in range(3 + ns)), 1e-7)
    check("f2f --double of single values: exact", max(rel(Fd[n], Fs[n].astype(complex)) for n in range(3 + ns)), 1e-300)
    run("f2f", "f1.bin", "f_copy.bin", "--double")
    check("f2f at the same precision: identical file", float(open("f1.bin", "rb").read() != open("f_copy.bin", "rb").read()), 0.5)

    # scalars
    run("h2f", "h0", "f_ns1.bin", "--nscal", 1)
    Fn, _ = read_fdchan("f_ns1.bin")
    check("--nscal 1: 4 fields, the same values", float(Fn.shape[0] != 4) + rel(Fn[3], F[3]), 1e-15)

    # FDCHAN single precision -> HiCHAN
    run("f2h", "f_sp.bin", "h4")
    H4, _ = read_hichan("h4")
    check("f2h from single precision: within single-precision round-off", max(rel(H4[n], C[n]) for n in range(3 + ns)), 1e-6)

    # refuses to overwrite
    r = subprocess.run([RS, "h2f", "h0", "f1.bin"], capture_output=True, text=True)
    check("refuses to overwrite an existing file without -f", float(r.returncode == 0), 0.5)

print("rsconv: OK" if ok else "rsconv: FAIL")
sys.exit(0 if ok else 1)
