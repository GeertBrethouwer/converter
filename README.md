# converter_and_spectra: restart files between HiCHAN and FDCHAN, and their spectra

Two programs:

* `rsconv`: conversion of restart files (below);
* `rsspec`: streamwise, spanwise and two-dimensional spectra of the velocities and scalars from restart files of
  either code (section [Spectra](#spectra-rsspec)).
* `tools/slices.py`: 2D slices (xz-, xy-, yz-planes) of the velocities and scalars from restart files of either
  code, for visualisation (section [Slices](#slices-toolsslicespy)).

`rsconv` converts restart files

* **HiCHAN → FDCHAN** (`h2f`) and **FDCHAN → HiCHAN** (`f2h`), velocities and scalars;
* **FDCHAN single ↔ double precision** (`f2f`);
* to a **different resolution in x and z** (`--nx`, `--nz`, any of the conversions, also HiCHAN → HiCHAN with
  `h2h`), by spectral interpolation (zero padding or truncation of the Fourier modes).

Both are plain C++17 without libraries (their own FFT for any length, `src/rsio.hpp`), work on one y-plane at a
time (a few hundred MB of memory for production files) and use threads over the kx-lines.

## Get and compile

Any C++17 compiler, on a login node:

```bash
cd $HOME
git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
cd converter_and_spectra
make
```

The repository was called `converter` before. A clone under the old name keeps working (GitHub redirects), but
the job scripts look for `$HOME/converter_and_spectra`; rename the old clone and point it to the new name:

```bash
mv $HOME/converter $HOME/converter_and_spectra
cd $HOME/converter_and_spectra
git remote set-url origin https://github.com/GeertBrethouwer/converter_and_spectra.git
git pull
```

(or set `B=` / `RS=` to the directory when submitting a job script).

**Processor type.** The binaries run only on the type of processor they were built for: the login nodes of
Arrhenius are x86_64, the GH200 nodes ARM (aarch64). `make` remembers the type (`build/.arch`) and rebuilds when it
is run on the other one, and every job script runs `make` first on the node it runs on, so a job always has the
right binaries (a few seconds). `make` / `make test` on the login node gives x86 binaries for tests there.

This gives `./rsconv` and `./rsspec`. With another compiler: `make CXX=clang++` (or `CXX=nvc++` from the NVIDIA HPC SDK
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

Large files, as a batch job on one Arrhenius node:

* `examples/arrhenius_h2f.sh`: HiCHAN → FDCHAN. Checks the input and the free disk space, converts, and prints
  the lines for the FDCHAN input file (nx, ny, nz, lx, lz, grid, nscal, `init_type = restart`, `restart_file`).
  Settings as variables at the top or on the command line:

  ```bash
  sbatch -A <your-account> --export=ALL,HICHAN=/scratch/geert/run/data/restart,NZ=2688 $HOME/converter_and_spectra/examples/arrhenius_h2f.sh
  ```

  `HICHAN` the HiCHAN restart (without extension), `OUT` the FDCHAN file (default `restart_fdchan.bin` next to
  it), `NZ`, `NX` (default: as HiCHAN), `PREC=single` for fdchan_sp, `NSCAL`.
* `examples/arrhenius_f2h.sh`: FDCHAN (single or double precision) → HiCHAN. Checks the input, an even output nz
  and the free disk space, converts (printing the dropped Nyquist energy), and prints the `&grid`, `&init` and
  `&scalars` lines for the HiCHAN input file:

  ```bash
  sbatch -A <your-account> --export=ALL,FDCHAN=/scratch/geert/fd_run/data/restart.bin,NZ=1920,DT=1.5e-3 $HOME/converter_and_spectra/examples/arrhenius_f2h.sh
  ```

  `FDCHAN` the FDCHAN file, `OUT` the HiCHAN restart name (default `restart_hichan` next to it), `NZ` (even),
  `NX`, `NSCAL`, and for the HiCHAN header `STEP` (default 0), `DT` (default 1e-3; HiCHAN's first time step),
  `GRID` (default `"tanh 4.5 8"`: grid_type grid_alpha fd_order).
* `examples/arrhenius_f2f.sh`: FDCHAN single ↔ double precision (optionally another nx, nz):

  ```bash
  sbatch -A <your-account> --export=ALL,IN=/scratch/geert/fd_sp/data/restart.bin $HOME/converter_and_spectra/examples/arrhenius_f2f.sh
  ```

  `IN` the FDCHAN file, `PREC` the output precision (default `double`; `single` for fdchan_sp), `OUT` (default
  `restart_double.bin` / `restart_single.bin` next to it), `NZ`, `NX`, `NSCAL`.
* `examples/arrhenius_convert.sh`: any conversion (edit its last lines).
Speed: about 0.2 GB/s on 4 cores (compute bound); on a 72-core node the file system sets the speed, a few minutes
for a 2560 × 385 × 1920 restart with one scalar (76 GB).

## Spectra (`rsspec`)

```bash
./rsspec --yplus 5,15,30,100,200 --retau 1000 --grid tanh 4.5 --sym --out spectra  data/restart.bin  [more files ...]
```

reads the y-planes at the chosen positions of all fields (u, v, w, the scalars) from HiCHAN or FDCHAN restart files
(also mixed, single or double precision; same nx, ny, nz, lx, lz) and writes, per position, into `--out`:

| file | contents |
|---|---|
| `spec_x_jNNNN.dat` | kx, λx, then Ex(kx) of every field |
| `spec_z_jNNNN.dat` | kz, λz, then Ez(kz) of every field |
| `spec_2d_jNNNN.dat` | kx, kz, λx, λz, then E2(kx, kz) of every field; a block per kx (gnuplot `splot`; numpy: `reshape(nx/2, nz/2+1, -1)`) |
| `summary.dat` | j, y, y⁺, the variance ⟨f′f′⟩ and the mean of every field per position |

The header of each file gives j, y, y⁺, the times of the samples and the variance of every field.

Definitions (f′ = f minus its x-z mean at that y; one-sided spectra, kx, kz ≥ 0; densities):
Σ Ex(kx) Δkx = Σ Ez(kz) Δkz = Σ E2(kx, kz) Δkx Δkz = ⟨f′f′⟩, with Δkx = 2π/lx, Δkz = 2π/lz. From the Fourier
coefficients c: E2 Δkx Δkz = a_i (|c(i, m)|² + |c(i, −m)|²), a₀ = 1, a_i = 2 for kx > 0; the kz Nyquist mode of an
even nz counted once; the mean (kx = kz = 0) left out. The usual premultiplied spectra are kx Ex, kz Ez and
kx kz E2 against λ (`tools/plot_spectra.py`).

| option | meaning |
|---|---|
| `--yplus a,b,...` | positions in wall units from the bottom wall (needs `--retau`) |
| `--y a,b,...` | y-values (−1 bottom wall, 0 centre, 1 top wall) |
| `--j a,b,...` | grid indices (0 = bottom wall) |
| `--retau R` | Re_τ: for `--yplus` and y⁺ in the output |
| `--grid TYPE ALPHA` | y-grid of the run (`tanh`, `uniform`, `cosine`); default from a HiCHAN header, needed for FDCHAN files with `--y`/`--yplus` |
| `--yfile FILE` | y-grid from a file (ny values), e.g. HiCHAN's `grid_file` |
| `--sym` | average also the mirror points in the upper half (j and ny−1−j) |
| `--fields u,v,w,s1` | fields (default: all) |
| `--out DIR` | output directory (default `spectra`) |
| `--threads N` | threads (default: all cores) |

The nearest grid points are used; the actual y and y⁺ are printed. Several files are averaged (time samples): a
run overwrites its restart file every `iout_restart` steps, so keep copies for averages: both codes write them
with `tout_samples` (or `iout_samples`) into `data/samples/restart_t<time>.bin` (HiCHAN: `&output`, also `.hdr`). Only the needed planes are read:
seconds even for production files.

Plots (python3 with numpy and matplotlib):

```bash
python3 tools/plot_spectra.py spectra --retau 1000 --utau 0.0500
```

writes `spectra/spec_1d.png` (premultiplied kx Ex and kz Ez against λx⁺, λz⁺, a panel per field, a line per
position) and `spectra/spec_2d.png` (premultiplied kx kz E2, contours at 0.125, 0.25, 0.5, 0.75 of the maximum, a
panel per position and field). `--utau` (= Re_τ/Re for velocities in bulk units) scales the velocity spectra to
wall units.

On Arrhenius, `examples/arrhenius_spectra.sh` does both as a batch job:

```bash
sbatch -A <your-account> --export=ALL,FILES="/scratch/geert/run/samples/restart_*.bin",RETAU=1000,RE=20000,GRID="tanh 4.5" $HOME/converter_and_spectra/examples/arrhenius_spectra.sh
```

with `YPLUS` (default `5,15,30,100,200`; or `Y`, `J`), `SYM` (default 1), `FIELDS`, `OUT` (default `spectra` in
the submit directory). If python3 with matplotlib is not available on the node it prints the plot command.

For a HiCHAN run, `examples/arrhenius_spectra_hichan.sh` needs only the run directory:

```bash
sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run $HOME/converter_and_spectra/examples/arrhenius_spectra_hichan.sh
```

It takes the time samples `data/samples/restart_t*` (from `tout_samples`/`iout_samples` in `&output`), or
`data/restart` if there are none; the y-grid from the restart headers (a `grid_type = 'file'` grid from
`grid_file` of the input); Re from the input file (`*.nml` in `RUN`); and Re_τ as the mean of `log.dat` over the
times of the samples. Optional: `TMIN`, `TMAX` (only samples in this time range, e.g. after the transient),
`SAMPLES` (other files, globs work), `YPLUS` (default `5,15,30,100,200`; or `Y`, `J`), `RETAU`, `RE`, `SYM`,
`FIELDS`, `OUT` (default `RUN/spectra`), `INPUT`, `DATA`, `YFILE`.

For an FDCHAN run, `examples/arrhenius_spectra_fdchan.sh` does the same:

```bash
sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run $HOME/converter_and_spectra/examples/arrhenius_spectra_fdchan.sh
```

It takes the time samples `data/samples/restart_t*.bin` (from `tout_samples`/`iout_samples`), or
`data/restart.bin` if there are none (single or double precision); `data_dir`, the y-grid (`grid_type`,
`grid_alpha`; FDCHAN's defaults tanh 3.0 if not given) and Re from the input file (`RUN/input.txt`, since FDCHAN
restart files do not contain the grid); and Re_τ as the mean of `log.dat` over the times of the samples. The same
options as the HiCHAN script (`TMIN`, `TMAX`, `SAMPLES`, `YPLUS`/`Y`/`J`, `RETAU`, `RE`, `SYM`, `FIELDS`, `OUT`,
`INPUT`, `DATA`), and `GRID="tanh 4.5"` to set the grid.

## Slices (`tools/slices.py`)

`tools/slices.py` (python3 with numpy; scipy for `.mat` files) extracts 2D slices of u, v, w and the scalars from
a HiCHAN or FDCHAN restart file, single or double precision, as HiCHAN's `tools/restart_slices.py`:

* xz-planes at given wall-normal positions (`--y`, or `--yplus` with `--retau`), arrays (Nz, Nx, n_y);
* an xy-plane at `--z`, array (ny, Nx), and a yz-plane at `--x`, array (ny, Nz).

The fields are evaluated exactly from their Fourier series at the grid points y_j, optionally on a finer x-z grid
(`--upsample 2`: exact Fourier interpolation; for FDCHAN the z-points are transformed to their Fourier modes first,
and `--upsample 1` gives exactly the stored z-point values). `--fluct` removes the x-z mean. The output is a
MATLAB `.mat` file with the same contents as HiCHAN's (x, y, z, time, step, lx, lz, y_xz, z_xy, x_yz, fields,
xz_<f>, xy_<f>, yz_<f>; for a HiCHAN file identical to the output of `restart_slices.py`), or a `.npz` file
(`-o name.npz`, or by default when scipy is missing).

```bash
python3 tools/slices.py data/restart --y -0.9 0 --z 1.0 --x 2.0 --fluct                         # HiCHAN
python3 tools/slices.py data/restart.bin --grid tanh 4.5 --yplus 15 100 --retau 1000 --upsample 2 --fluct -o s.mat   # FDCHAN
```

FDCHAN restart files do not contain the y-grid: give `--grid TYPE ALPHA` (`grid_type`, `grid_alpha` of the
input; FDCHAN's default is `tanh 3.0`) or `--grid-file`. Other options: `--fields u v s1`, `--re` and `--utau`
(instead of `--retau`, as in HiCHAN's script). Only the xz-planes are transformed as a whole; of the other y-planes
only the two lines of the xy- and yz-planes are computed: at 2560 × 385 × 1920 about 0.2 s per plane and field
for an FDCHAN file (less for HiCHAN), a few minutes for all fields.

Pictures:

```bash
python3 tools/plot_slices.py s.mat            # (or .npz) s_<field>.png, one figure per field
```

`tools/plot_slices.py` (matplotlib) draws, per field, the xz-planes in the top row and the xy- and yz-plane below
(options `--fields`, `--clim lo hi`, `--plane xz|xy|yz`); `tools/plot_slices.m` (copied from HiCHAN) does the same
in MATLAB or Octave: `plot_slices('s.mat')`, `plot_slices('s.mat', {'u','s1'}, 'save', true, 'clim', [-0.3 0.3])`.

On Arrhenius, `examples/arrhenius_slices.sh` does both as a batch job:

```bash
sbatch -A <your-account> --export=ALL,FILE=/scratch/geert/fdchan_run/data/restart.bin,RUN=/scratch/geert/fdchan_run,YPLUS="15 100" $HOME/converter_and_spectra/examples/arrhenius_slices.sh
```

with `FILE` (FDCHAN `.bin` or HiCHAN restart name), `RUN` (optional: the FDCHAN grid from `RUN/input.txt` and Re_τ
from `RUN/data/log.dat`), `GRID`, `YPLUS` + `RETAU` or `Y` (default `"-0.9 0"`), `X`, `Z`, `UPSAMPLE`, `FLUCT`
(default 1), `FIELDS`, `OUT`, `PLOT` (default 1).

For a whole run, `examples/arrhenius_slices_hichan.sh` and `examples/arrhenius_slices_fdchan.sh` need only the run
directory, as the spectra scripts:

```bash
sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run $HOME/converter_and_spectra/examples/arrhenius_slices_hichan.sh
sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run,WHICH=all,TMIN=500,YPLUS="15 100",CLIM="-0.3 0.3" $HOME/converter_and_spectra/examples/arrhenius_slices_fdchan.sh
```

`WHICH=last` (default): the latest time sample (`data/samples/restart_t*`), else the restart file; `WHICH=all`:
every sample with `TMIN` ≤ t ≤ `TMAX` (frames, e.g. for a film: set `CLIM` for the same colour range in all
pictures); `WHICH=restart`: `data/restart`. They write `RUN/slices/<file>_slices.mat` (`FORMAT=npz`, or
automatically without scipy) and the png pictures (`<file>_slices_<field>.png`). The y-grid comes from the HiCHAN
headers (`grid_type = 'file'`: `grid_file` of the input), for FDCHAN from `RUN/input.txt` (with `data_dir`;
`GRID` to set it); Re_τ for `YPLUS` from `log.dat` at the time of each file (or `RETAU`). Other settings: `Y`
(default `"-0.9 0"`), `X`, `Z`, `UPSAMPLE`, `FLUCT` (default 1), `FIELDS`, `SAMPLES`, `OUTDIR` (default
`RUN/slices`), `PLOT`.

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

`make test` runs `tests/test_fft`, `tests/test_roundtrip.py` (rsconv) and `tests/test_spectra.py` (rsspec).

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
| `tests/test_spectra.py`: rsspec against numpy from the physical field (random HiCHAN file, two scalars) | variance, Ex (FFTs in x), Ez (FFTs in z), E2: to the 9 digits of the output; FDCHAN file of the same field: the same spectra; single precision: 3e-8; two files: their average; `--sym`, `--yplus`, `--y`, the means |
| `tests/test_slices.py`: slices.py on a random HiCHAN file against a direct evaluation of the Fourier series (`--upsample` 1 and 1.5, `--fluct`); the FDCHAN file of the same field (double, single precision); an FDCHAN file with odd nz and with Nyquist content (stored z-values reproduced, `--upsample 2`); positions, `.npz`, plot_slices.py | 4e-15 (single precision 4e-8); on a HiCHAN restart identical to HiCHAN's restart_slices.py |
