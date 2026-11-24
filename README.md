# converter: restart files between HiCHAN and FDCHAN

`rsconv` converts restart files

* **HiCHAN → FDCHAN** (`h2f`) and **FDCHAN → HiCHAN** (`f2h`), velocities and scalars;
* **FDCHAN single ↔ double precision** (`f2f`);
* to a **different resolution in x and z** (`--nx`, `--nz`, any of the conversions, also HiCHAN → HiCHAN with
  `h2h`), by spectral interpolation (zero padding or truncation of the Fourier modes).

It is one C++17 file without libraries (its own FFT for any length), works on one y-plane at a time (a few
hundred MB of memory for production files) and uses threads over the kx-lines.

## Compile

Any C++17 compiler, on a login node:

```bash
cd $HOME/converter
make
```

This gives `./rsconv`. With another compiler: `make CXX=clang++` (or `CXX=nvc++` from the NVIDIA HPC SDK
module).

Tests (CPU, seconds; needs python3 with numpy):

```bash
make test
```

## Use

```bash
./rsconv info <file>                                  # header of a HiCHAN or FDCHAN restart file
./rsconv h2f <hichan_restart> <fdchan_out.bin> [options]
./rsconv f2h <fdchan_in.bin>  <hichan_restart> [options]
./rsconv f2f <fdchan_in.bin>  <fdchan_out.bin> [options]
./rsconv h2h <hichan_restart> <hichan_restart> [options]
```

A HiCHAN restart is a pair `<name>.bin` + `<name>.hdr` and is given by `<name>` (with or without the extension,
e.g. `data/restart`). An FDCHAN restart is one file (`data/restart.bin`).

| option | meaning |
|---|---|
| `--nx N` | nx of the output (nx/2 kx-modes); default: as the input |
| `--nz N` | nz of the output: FDCHAN z-points or HiCHAN Fourier modes; default: as the input |
| `--single`, `--double` | precision of an FDCHAN output; default: as an FDCHAN input, double for a HiCHAN input |
| `--nscal N` | keep only the first N scalars |
| `--time T` | time in the output header (default: as the input) |
| `--step N`, `--dt DT` | HiCHAN output: step number and time step in the header (default 0 and 1e-4 from FDCHAN, which does not store them; as the input for `h2h`) |
| `--grid TYPE ALPHA ORDER` | HiCHAN output: grid_type, grid_alpha, fd_order written to the header (for information; HiCHAN does not read them) |
| `--threads N` | threads (default: all cores) |
| `-f` | overwrite an existing output file |

The output is written to `<file>.tmp` and renamed at the end, so an interrupted conversion leaves no broken
restart file.

### Examples

HiCHAN → FDCHAN, same resolution, then restart FDCHAN from it:

```bash
./rsconv h2f /scratch/geert/hichan_run/data/restart restart_from_hichan.bin
```

In the FDCHAN input: `init_type = restart`, `restart_file = restart_from_hichan.bin`, and the same `nx`, `ny`,
`nz`, `lx`, `lz`, `grid_type`, `grid_alpha`, `nscal` as the file (`rsconv info` shows them).

HiCHAN → FDCHAN with more z-points (compact 8th order in z needs about 1.4× the Fourier modes for the same
resolution; nz must give at least 49 z-planes per GPU for `z_scheme = compact`, order 8):

```bash
./rsconv h2f data/restart restart_fd.bin --nz 2688
```

FDCHAN → HiCHAN (HiCHAN reads step and dt from the header; set them, e.g. to the last step of the FDCHAN run):

```bash
./rsconv f2h data/restart.bin restart_hc --step 0 --dt 1e-3 --grid tanh 4.5 8
```

Then in HiCHAN's `&init`: `init_type = 'restart'`, `restart_file = 'restart_hc'`.

FDCHAN single precision (`fdchan_sp`) → double precision, and back:

```bash
./rsconv f2f data/restart.bin restart_dp.bin --double
./rsconv f2f restart_dp.bin restart_sp.bin --single
```

Large files: `examples/arrhenius_convert.sh` runs a conversion as a batch job on one node
(`sbatch -A <your-account> $HOME/converter/examples/arrhenius_convert.sh`, after editing its last lines).
Speed: about 0.2 GB/s on 4 cores (compute bound); on a 72-core node the file system sets the speed, a few minutes
for a 2560 × 385 × 1920 restart with one scalar (76 GB).

## What the conversion does

Both codes represent a field by its Fourier coefficients in x, f(x, ·) = Σ f̂ exp(i kx x) over kx = 2πi/lx,
i = 0 … nx/2 − 1 (the conjugates implied), with the coefficient equal to the amplitude of the mode. They differ
in z:

| | HiCHAN | FDCHAN |
|---|---|---|
| z | `nz` Fourier modes in FFT order (0 … nz/2 − 1, the Nyquist mode nz/2 (zero), −nz/2 + 1 … −1) | `nz` points z_k = k lz/nz |
| layout of a field | Fortran (nx/2, nz, ny): kx fastest, then kz, then y | C (ny, nz, nx/2): kx fastest, then z, then y |
| fields | u, v, w, s₁ … s_nscal | u, v, w, s₁ … s_n |
| header | `<name>.hdr` (text): nx ny nz nscal / lx lz / time istep dt / grid_type grid_alpha fd_order | 40 bytes: int nx, ny, nz, sizeof(real); double lx, lz, time; n from the file size |
| precision | double | single or double |

`h2f` evaluates the Fourier series in z at the FDCHAN points (an inverse FFT per (kx, y)-line), `f2h` takes the
discrete Fourier transform of the z-points (divided by nz). Both are exact inverses of each other (to round-off)
when nz is the same and the FDCHAN field has no Nyquist content. y is not touched: **ny and the y-grid
(grid_type, grid_alpha) must be the same in both runs**; the files do not contain the grid, so check the inputs.
The y-points run from y = −1 (index 0) to y = 1 in both codes.

Things to know:

* **Nyquist mode.** HiCHAN has no kz Nyquist mode (it is zero); `f2h` drops it, and `rsconv` prints the fraction
  of the energy that is dropped per field. On a resolved grid it is tiny; it grows when the FDCHAN field has
  energy at the grid scale.
* **Resolution change.** `--nz`/`--nx` keep the modes |m| ≤ (min(nz_in, nz_out) − 1)/2 and kx < min(nx_in,
  nx_out)/2: zero padding for a finer grid (an exact interpolation of the Fourier series), truncation for a
  coarser one (the dropped energy fraction is printed).
* **Divergence.** A field that is divergence-free with FDCHAN's z-differences is not exactly divergence-free with
  HiCHAN's Fourier z-derivative, and the other way round (the derivatives differ at the high wavenumbers).
  FDCHAN projects the velocity after reading a restart (div ~ 1e-14 from the start). HiCHAN does not: its first
  `log.dat` line after the restart shows a larger `div_max` (on a coarse test grid 0.3), and from the first time
  step on it is at round-off again. The flow statistics are unchanged.
* **Step and Δt.** FDCHAN restart files contain only the time; `f2h` writes `--step` and `--dt` into the HiCHAN
  header (HiCHAN uses dt as the first time step if the time step is adaptive). FDCHAN starts with `dt` of its input
  file and grows it at `dt_growth` per step.
* **kx = 0.** FDCHAN keeps the kx = 0 values real in z; `h2f` sets their (round-off) imaginary parts to zero and
  prints the largest one, which shows that the HiCHAN data are Hermitian as they must be.

## Tests

| test | result |
|---|---|
| `tests/test_fft.cpp`: the FFT against a direct DFT, lengths 1–130, 1920, 2688, 2689 (prime, Bluestein), 2880 | 4e-15 |
| `tests/test_roundtrip.py`: h2f gives Σ c_m exp(2πi m k/nz) (numpy) | 6e-16 |
| h2f + f2h, also via a finer FDCHAN grid (`--nz 27 --nx 24`, odd nz) and back | coefficients back to 1e-15 |
| h2h `--nz 20 --nx 10` | exactly the padded/truncated spectrum |
| f2f `--single` then `--double` | single-precision round-off (4e-8); same precision: identical file |
| headers, `--nscal`, refusing to overwrite | ok |
| FDCHAN turbulent field (32 × 65 × 48, one scalar, 300 steps) → f2h → HiCHAN (CPU build) | HiCHAN's first log line: τ_w, Re_τ, dS/dy at both walls identical to FDCHAN's last; U_b 1.000007 (Nyquist energy dropped, coarse grid); div_max 0.34, then 1e-14 after one step |
| HiCHAN turbulent field (same grid, 300 steps) → h2f → FDCHAN | FDCHAN's first log line: U_b, τ_w, Re_τ, dS/dy identical to HiCHAN's last; div_max 7e-15; E_kin within 0.1 % |
