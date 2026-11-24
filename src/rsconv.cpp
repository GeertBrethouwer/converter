// rsconv: conversion of restart files between HiCHAN and FDCHAN, and between FDCHAN precisions.
//
//   rsconv info <file>
//   rsconv h2f  <hichan_restart> <fdchan_out.bin> [options]   HiCHAN  -> FDCHAN
//   rsconv f2h  <fdchan_in.bin>  <hichan_restart> [options]   FDCHAN  -> HiCHAN
//   rsconv f2f  <fdchan_in.bin>  <fdchan_out.bin> [options]   FDCHAN  -> FDCHAN (precision, resolution)
//   rsconv h2h  <hichan_restart> <hichan_restart> [options]   HiCHAN  -> HiCHAN (resolution)
//
// File formats (see README.md):
//   HiCHAN: <name>.hdr (text: nx ny nz nscal / lx lz / time istep dt / grid_type grid_alpha fd_order) and
//           <name>.bin: u, v, w, s_1..s_nscal, each a complex(8) array (nx/2, nz, ny) in Fortran order, i.e. kx
//           fastest, then the kz modes in FFT order (0..nz/2-1, Nyquist (zero), -nz/2+1..-1), then y.
//   FDCHAN: 40-byte header (int nx, ny, nz, sizeof(real); double lx, lz, time), then u, v, w, s_1..s_n, each a
//           complex array (ny, nz, nx/2) in C order (kx fastest, then the physical z-points, then y), in single or
//           double precision. The number of scalars follows from the file size.
// Both codes use f(x, z) = sum over the coefficients f^ exp(i (kx x + kz z)) with kx >= 0 (the conjugate of kx > 0
// implied), so a coefficient is the amplitude of its mode; the z-points of FDCHAN are z_k = k lz / nz. The y-grid
// (ny points from y = -1 to 1) is the same in both and is not changed by the conversion.
//
// The fields are converted one y-plane at a time (memory: a few planes), with threads over the kx-lines.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using cd = std::complex<double>;

[[noreturn]] static void fail(const std::string& m) {
  std::fprintf(stderr, "rsconv: %s\n", m.c_str());
  std::exit(1);
}

// ------------------------------------------------------------------------------------------- FFT (any length)
// Mixed-radix Cooley-Tukey (any factors, a direct DFT per factor) for lengths whose prime factors are <= 64, else
// Bluestein's algorithm with a power-of-two FFT. exec(x, sign): x_k <- sum_j x_j exp(sign 2 pi i j k / n),
// unnormalised.
class FFT {
 public:
  explicit FFT(int n) : n_(n) {
    int m = n;
    for (int p = 2; p * p <= m; ++p)
      while (m % p == 0) { fac_.push_back(p); m /= p; }
    if (m > 1) fac_.push_back(m);
    blue_ = !fac_.empty() && *std::max_element(fac_.begin(), fac_.end()) > 64;
    if (!blue_) {
      w_.resize(size_t(n));
      for (int k = 0; k < n; ++k) w_[size_t(k)] = std::polar(1.0, -2 * M_PI * k / n);
    } else {
      M_ = 1;
      while (M_ < 2 * n - 1) M_ *= 2;
      sub_.reset(new FFT(M_));
      chirp_.resize(size_t(n));
      for (int k = 0; k < n; ++k) {
        const long long kk = (long long)k * k % (2LL * n);  // exact phase for large k
        chirp_[size_t(k)] = std::polar(1.0, -M_PI * double(kk) / n);
      }
      bf_.assign(size_t(M_), cd(0, 0));
      bf_[0] = std::conj(chirp_[0]);
      for (int k = 1; k < n; ++k) bf_[size_t(k)] = bf_[size_t(M_ - k)] = std::conj(chirp_[size_t(k)]);
      sub_->exec(bf_.data(), -1);
    }
  }
  int size() const { return n_; }
  // work must hold at least work_size() values
  size_t work_size() const { return blue_ ? size_t(M_) + sub_->work_size() : size_t(n_); }
  void exec(cd* x, int sign, cd* work) const {
    if (n_ == 1) return;
    if (blue_) { bluestein(x, sign, work); return; }
    std::copy(x, x + n_, work);
    rec(x, work, n_, 1, 0, sign);
  }
  void exec(cd* x, int sign) const {
    std::vector<cd> w(work_size());
    exec(x, sign, w.data());
  }

 private:
  int n_, M_ = 0;
  bool blue_ = false;
  std::vector<int> fac_;
  std::vector<cd> w_, chirp_, bf_;
  std::unique_ptr<FFT> sub_;

  cd tw(long long t, int sign) const {  // exp(sign 2 pi i t / n)
    const cd z = w_[size_t(t % n_)];
    return sign < 0 ? z : std::conj(z);
  }
  // out[0..n) = DFT of in[0], in[s], ..., in[(n-1) s] (n divides n_)
  void rec(cd* out, const cd* in, int n, int s, int fi, int sign) const {
    if (n == 1) { out[0] = in[0]; return; }
    const int p = fac_[size_t(fi)], m = n / p;
    for (int r = 0; r < p; ++r) rec(out + size_t(r) * m, in + size_t(r) * s, m, s * p, fi + 1, sign);
    const int st = n_ / n;  // W_n^t = W_N^(t st)
    cd t[64];
    std::vector<cd> big;
    cd* tp = t;
    if (p > 64) { big.resize(size_t(p)); tp = big.data(); }
    for (int k = 0; k < m; ++k) {
      for (int r = 0; r < p; ++r) tp[r] = out[size_t(r) * m + k] * tw((long long)r * k * st, sign);
      for (int q = 0; q < p; ++q) {
        cd a = tp[0];
        for (int r = 1; r < p; ++r) a += tp[r] * tw((long long)r * q % p * (n_ / p), sign);
        out[size_t(q) * m + k] = a;
      }
    }
  }
  void bluestein(cd* x, int sign, cd* work) const {
    cd* a = work;
    cd* sw = work + M_;
    auto ch = [&](int k) { return sign < 0 ? chirp_[size_t(k)] : std::conj(chirp_[size_t(k)]); };
    std::fill(a, a + M_, cd(0, 0));
    for (int k = 0; k < n_; ++k) a[k] = x[k] * ch(k);
    sub_->exec(a, -1, sw);
    // the transform of conj(chirp) is bf_ for sign = -1; for sign = +1 use the conjugate chirp: its transform is
    // the reversed conjugate of bf_
    for (int k = 0; k < M_; ++k) a[k] *= (sign < 0) ? bf_[size_t(k)] : std::conj(bf_[size_t((M_ - k) % M_)]);
    sub_->exec(a, +1, sw);
    for (int k = 0; k < n_; ++k) x[k] = a[k] * ch(k) / double(M_);
  }
};

// ------------------------------------------------------------------------------------------- file descriptions
struct Desc {
  std::string fmt;          // "hichan" or "fdchan"
  std::string bin, hdr;     // data file, HiCHAN header file
  int nx = 0, ny = 0, nz = 0, nf = 0, prec = 8;  // nf: number of fields (3 + nscal); prec: bytes per real
  double lx = 0, lz = 0, time = 0, dt = 0;
  long step = 0;
  std::string grid_type = "unknown";
  double grid_alpha = 0;
  int fd_order = 0;
  int nkx() const { return nx / 2; }
  long long plane_bytes() const { return (long long)nkx() * nz * 2 * prec; }
  long long field_bytes() const { return plane_bytes() * ny; }
  long long offset() const { return fmt == "fdchan" ? 40 : 0; }
};

static long long file_size(const std::string& f) {
  std::FILE* p = std::fopen(f.c_str(), "rb");
  if (!p) fail("cannot open " + f);
  std::fseek(p, 0, SEEK_END);
  const long long s = (long long)ftello(p);
  std::fclose(p);
  return s;
}

static std::string strip(const std::string& name) {  // HiCHAN restart name without .bin / .hdr
  for (const char* e : {".bin", ".hdr"})
    if (name.size() > 4 && name.compare(name.size() - 4, 4, e) == 0) return name.substr(0, name.size() - 4);
  return name;
}

static Desc read_hichan(const std::string& name) {
  Desc d;
  d.fmt = "hichan";
  const std::string b = strip(name);
  d.bin = b + ".bin"; d.hdr = b + ".hdr";
  std::ifstream h(d.hdr);
  if (!h) fail("cannot open the HiCHAN header " + d.hdr);
  std::string line;
  std::getline(h, line);  // comment
  int ns = 0;
  if (!(h >> d.nx >> d.ny >> d.nz >> ns >> d.lx >> d.lz >> d.time >> d.step >> d.dt)) fail("cannot read " + d.hdr);
  if (!(h >> d.grid_type >> d.grid_alpha >> d.fd_order)) { d.grid_type = "unknown"; }
  d.nf = 3 + ns;
  d.prec = 8;
  const long long need = d.field_bytes() * d.nf, have = file_size(d.bin);
  if (have < need)
    fail(d.bin + " is too short: " + std::to_string(have) + " bytes, the header needs " + std::to_string(need));
  return d;
}

static Desc read_fdchan(const std::string& name) {
  Desc d;
  d.fmt = "fdchan";
  d.bin = name;
  std::FILE* f = std::fopen(name.c_str(), "rb");
  if (!f) fail("cannot open " + name);
  int h[4];
  double x[3];
  if (std::fread(h, sizeof(int), 4, f) != 4 || std::fread(x, sizeof(double), 3, f) != 3) fail("cannot read the header of " + name);
  std::fclose(f);
  d.nx = h[0]; d.ny = h[1]; d.nz = h[2]; d.prec = h[3];
  d.lx = x[0]; d.lz = x[1]; d.time = x[2];
  if ((d.prec != 4 && d.prec != 8) || d.nx <= 0 || d.ny <= 0 || d.nz <= 0 || d.nx % 2)
    fail(name + " does not look like an FDCHAN restart file (header " + std::to_string(h[0]) + " " +
         std::to_string(h[1]) + " " + std::to_string(h[2]) + " " + std::to_string(h[3]) + ")");
  const long long data = file_size(name) - 40;
  if (data % d.field_bytes()) fail(name + ": the size is not a whole number of fields");
  d.nf = int(data / d.field_bytes());
  if (d.nf < 3) fail(name + ": fewer than 3 fields");
  return d;
}

static bool is_hichan(const std::string& name) {
  std::ifstream h(strip(name) + ".hdr");
  return bool(h);
}

static void print_desc(const Desc& d, const char* what) {
  std::printf("%s: %s restart (%s)\n", what, d.fmt == "hichan" ? "HiCHAN" : "FDCHAN", d.bin.c_str());
  std::printf("  nx = %d (%d kx modes), ny = %d, nz = %d (%s), %d scalar(s), %s precision\n", d.nx, d.nkx(), d.ny, d.nz,
              d.fmt == "hichan" ? "Fourier modes" : "z-points", d.nf - 3, d.prec == 8 ? "double" : "single");
  std::printf("  lx = %.10g, lz = %.10g, time = %.10g", d.lx, d.lz, d.time);
  if (d.fmt == "hichan")
    std::printf(", step = %ld, dt = %.6g, grid %s alpha = %g, fd_order %d", d.step, d.dt, d.grid_type.c_str(),
                d.grid_alpha, d.fd_order);
  std::printf("\n  %.2f GB\n", double(d.field_bytes()) * d.nf / 1e9);
}

// ------------------------------------------------------------------------------------------- planes
// A plane in spectral form: S[mi * nkx + i], mi the kz index in FFT order (0..nz-1), i the kx index.
static int signed_mode(int mi, int n) { return mi <= n / 2 ? mi : mi - n; }

static void parallel_for(int n, int nthreads, const std::function<void(int, int)>& f) {  // f(item, thread)
  if (nthreads <= 1 || n < 2) {
    for (int i = 0; i < n; ++i) f(i, 0);
    return;
  }
  std::atomic<int> next(0);
  std::vector<std::thread> th;
  for (int t = 0; t < nthreads; ++t)
    th.emplace_back([&, t]() {
      for (int i; (i = next.fetch_add(1)) < n;) f(i, t);
    });
  for (auto& x : th) x.join();
}

struct Stats {
  double tot = 0, nyq = 0, cut = 0, imag0 = 0;  // energies: all, Nyquist (dropped), truncated; |Im| of kx = 0 in z
};

class Converter {
 public:
  Converter(const Desc& in, const Desc& out, int nthreads) : in_(in), out_(out), nth_(nthreads) {
    if (in.ny != out.ny) fail("ny cannot be changed");
    if (in_.fmt == "fdchan") fin_.reset(new FFT(in.nz));
    if (out_.fmt == "fdchan") fout_.reset(new FFT(out.nz));
    const size_t ws = std::max(fin_ ? fin_->work_size() : 0, fout_ ? fout_->work_size() : 0);
    work_.assign(size_t(nth_), std::vector<cd>(ws + size_t(std::max(in.nz, out.nz))));
  }
  // direct copy (FDCHAN to FDCHAN at the same nx, nz: only the precision changes, no transforms)
  bool direct() const { return in_.fmt == out_.fmt && in_.nx == out_.nx && in_.nz == out_.nz; }

  // in: raw plane of the input file; out: raw plane of the output file
  void plane(const std::vector<char>& rin, std::vector<char>& rout, Stats& st) {
    if (direct()) {
      const long long n = (long long)in_.nkx() * in_.nz * 2;
      for (long long k = 0; k < n; ++k) put(rout, out_.prec, k, get(rin, in_.prec, k));
      return;
    }
    const int nki = in_.nkx(), nko = out_.nkx(), nzi = in_.nz, nzo = out_.nz;
    std::vector<cd> S(size_t(nki) * nzi), T(size_t(nko) * nzo, cd(0, 0));
    // 1. input -> spectral S
    if (in_.fmt == "hichan") {
      for (long long k = 0; k < (long long)nki * nzi; ++k) S[size_t(k)] = cd(get(rin, 8, 2 * k), get(rin, 8, 2 * k + 1));
    } else {
      parallel_for(nki, nth_, [&](int i, int t) {
        cd* line = work_[size_t(t)].data();
        cd* w = line + nzi;
        for (int k = 0; k < nzi; ++k) {
          const long long e = 2 * ((long long)k * nki + i);
          line[k] = cd(get(rin, in_.prec, e), get(rin, in_.prec, e + 1));
        }
        fin_->exec(line, -1, w);
        for (int m = 0; m < nzi; ++m) S[size_t(m) * nki + i] = line[m] / double(nzi);
      });
    }
    // 2. resize: kx < min(nkx); |m| <= (min(nz) - 1) / 2 (the Nyquist mode of an even nz is not kept)
    const int nkm = std::min(nki, nko), mmax = (std::min(nzi, nzo) - 1) / 2;
    for (int mi = 0; mi < nzi; ++mi) {
      const int m = signed_mode(mi, nzi);
      const bool nyq = (nzi % 2 == 0 && mi == nzi / 2);
      for (int i = 0; i < nki; ++i) {
        const cd a = S[size_t(mi) * nki + i];
        const double e = std::norm(a) * (i == 0 ? 1.0 : 2.0);
        st.tot += e;
        if (nyq) { st.nyq += e; continue; }
        if (i >= nkm || std::abs(m) > mmax) { st.cut += e; continue; }
        const int mo = m >= 0 ? m : m + nzo;
        T[size_t(mo) * nko + i] = a;
      }
    }
    // 3. spectral T -> output
    if (out_.fmt == "hichan") {
      for (long long k = 0; k < (long long)nko * nzo; ++k) {
        put(rout, 8, 2 * k, T[size_t(k)].real());
        put(rout, 8, 2 * k + 1, T[size_t(k)].imag());
      }
    } else {
      std::vector<double> im0(size_t(nth_), 0.0);
      parallel_for(nko, nth_, [&](int i, int t) {
        cd* line = work_[size_t(t)].data();
        cd* w = line + nzo;
        for (int m = 0; m < nzo; ++m) line[m] = T[size_t(m) * nko + i];
        fout_->exec(line, +1, w);
        for (int k = 0; k < nzo; ++k) {
          cd a = line[k];
          if (i == 0) {  // real for kx = 0 (exactly, as FDCHAN keeps it)
            im0[size_t(t)] = std::max(im0[size_t(t)], std::abs(a.imag()));
            a = cd(a.real(), 0);
          }
          const long long e = 2 * ((long long)k * nko + i);
          put(rout, out_.prec, e, a.real());
          put(rout, out_.prec, e + 1, a.imag());
        }
      });
      for (double x : im0) st.imag0 = std::max(st.imag0, x);
    }
  }

 private:
  Desc in_, out_;
  int nth_;
  std::unique_ptr<FFT> fin_, fout_;
  std::vector<std::vector<cd>> work_;
  static double get(const std::vector<char>& r, int prec, long long k) {
    if (prec == 8) { double x; std::memcpy(&x, r.data() + k * 8, 8); return x; }
    float x; std::memcpy(&x, r.data() + k * 4, 4); return double(x);
  }
  static void put(std::vector<char>& r, int prec, long long k, double v) {
    if (prec == 8) { std::memcpy(r.data() + k * 8, &v, 8); return; }
    const float x = float(v);
    std::memcpy(r.data() + k * 4, &x, 4);
  }
};

// ------------------------------------------------------------------------------------------- main
static void usage() {
  std::puts(
      "usage:\n"
      "  rsconv info <file>                               print the header of a HiCHAN or FDCHAN restart file\n"
      "  rsconv h2f <hichan_restart> <fdchan_out.bin> [options]   HiCHAN -> FDCHAN\n"
      "  rsconv f2h <fdchan_in.bin> <hichan_restart> [options]    FDCHAN -> HiCHAN\n"
      "  rsconv f2f <fdchan_in.bin> <fdchan_out.bin> [options]    FDCHAN -> FDCHAN (precision, nx, nz)\n"
      "  rsconv h2h <hichan_restart> <hichan_restart> [options]   HiCHAN -> HiCHAN (nx, nz)\n"
      "HiCHAN restart names may be given with or without .bin/.hdr (e.g. data/restart).\n"
      "options:\n"
      "  --nx N        x-modes of the output (nx/2 kx modes; default: as the input)\n"
      "  --nz N        z of the output: FDCHAN z-points or HiCHAN Fourier modes (default: as the input)\n"
      "  --single      FDCHAN output in single precision (default: the precision of an FDCHAN input, else double)\n"
      "  --double      FDCHAN output in double precision\n"
      "  --nscal N     keep only the first N scalars (default: all)\n"
      "  --step N      HiCHAN output: time-step number in the header (default: 0)\n"
      "  --dt DT       HiCHAN output: time step in the header (default: 1e-4)\n"
      "  --time T      time in the output header (default: as the input)\n"
      "  --grid TYPE ALPHA ORDER   HiCHAN output: grid_type, grid_alpha, fd_order written to the header (for\n"
      "                information only; HiCHAN does not read them). Default: tanh 0 0, or as a HiCHAN input\n"
      "  --threads N   threads (default: all cores)\n"
      "  -f            overwrite the output file");
}

int main(int argc, char** argv) {
  if (argc < 3) { usage(); return 1; }
  const std::string cmd = argv[1];
  if (cmd == "info") {
    for (int a = 2; a < argc; ++a) {
      const std::string f = argv[a];
      print_desc(is_hichan(f) ? read_hichan(f) : read_fdchan(f), "file");
    }
    return 0;
  }
  if (cmd != "h2f" && cmd != "f2h" && cmd != "f2f" && cmd != "h2h") { usage(); return 1; }
  if (argc < 4) { usage(); return 1; }
  const std::string fin = argv[2], fout = argv[3];
  int nx = 0, nz = 0, nscal = -1, nth = int(std::max(1u, std::thread::hardware_concurrency())), prec = 0;
  long step = 0;
  double dt = 1e-4, time = NAN, alpha = NAN;
  std::string gtype;
  int order = -1;
  bool force = false;
  for (int a = 4; a < argc; ++a) {
    const std::string o = argv[a];
    auto val = [&](const char* what) -> std::string {
      if (a + 1 >= argc) fail(std::string("missing value for ") + what);
      return argv[++a];
    };
    if (o == "--nx") nx = std::atoi(val("--nx").c_str());
    else if (o == "--nz") nz = std::atoi(val("--nz").c_str());
    else if (o == "--nscal") nscal = std::atoi(val("--nscal").c_str());
    else if (o == "--single") prec = 4;
    else if (o == "--double") prec = 8;
    else if (o == "--step") step = std::atol(val("--step").c_str());
    else if (o == "--dt") dt = std::atof(val("--dt").c_str());
    else if (o == "--time") time = std::atof(val("--time").c_str());
    else if (o == "--threads") nth = std::max(1, std::atoi(val("--threads").c_str()));
    else if (o == "--grid") {
      gtype = val("--grid");
      alpha = std::atof(val("--grid").c_str());
      order = std::atoi(val("--grid").c_str());
    } else if (o == "-f") force = true;
    else fail("unknown option " + o);
  }
  const bool hin = (cmd[0] == 'h'), hout = (cmd[2] == 'h');
  const Desc in = hin ? read_hichan(fin) : read_fdchan(fin);
  Desc out = in;
  out.fmt = hout ? "hichan" : "fdchan";
  if (hout) { out.bin = strip(fout) + ".bin"; out.hdr = strip(fout) + ".hdr"; }
  else { out.bin = fout; out.hdr.clear(); }
  if (nx) {
    if (nx % 2 || nx < 2) fail("--nx must be even");
    out.nx = nx;
  }
  if (nz) out.nz = nz;
  if (nscal >= 0) {
    if (nscal > in.nf - 3) fail("--nscal: the input has only " + std::to_string(in.nf - 3) + " scalar(s)");
    out.nf = 3 + nscal;
  }
  out.prec = hout ? 8 : (prec ? prec : (hin ? 8 : in.prec));
  if (hout && prec == 4) fail("HiCHAN restart files are double precision (--single is for FDCHAN output)");
  if (!std::isnan(time)) out.time = time;
  if (hout) {
    if (!hin) { out.step = step; out.dt = dt; out.grid_type = "tanh"; out.grid_alpha = 0; out.fd_order = 0; }
    else {
      for (int a = 4; a < argc; ++a)
        if (std::string(argv[a]) == "--step") out.step = step;
        else if (std::string(argv[a]) == "--dt") out.dt = dt;
    }
    if (!gtype.empty()) { out.grid_type = gtype; out.grid_alpha = alpha; out.fd_order = order; }
  }
  if (out.bin == in.bin) fail("the output must be a different file than the input");
  if (!force) {
    std::FILE* t = std::fopen(out.bin.c_str(), "rb");
    if (t) { std::fclose(t); fail(out.bin + " exists (use -f to overwrite)"); }
  }
  print_desc(in, "input ");
  std::printf("output: %s restart (%s)\n", hout ? "HiCHAN" : "FDCHAN", out.bin.c_str());
  std::printf("  nx = %d, ny = %d, nz = %d (%s), %d scalar(s), %s precision, time = %.10g\n", out.nx, out.ny, out.nz,
              hout ? "Fourier modes" : "z-points", out.nf - 3, out.prec == 8 ? "double" : "single", out.time);
  if (!hin && !hout && in.nz == out.nz && in.nx == out.nx)
    std::printf("  (same nx and nz: the values are copied, only the precision changes)\n");
  std::fflush(stdout);

  Converter C(in, out, nth);
  std::FILE* fi = std::fopen(in.bin.c_str(), "rb");
  if (!fi) fail("cannot open " + in.bin);
  const std::string tmp = out.bin + ".tmp";
  std::FILE* fo = std::fopen(tmp.c_str(), "wb");
  if (!fo) fail("cannot write " + tmp);
  if (!hout) {
    const int h[4] = {out.nx, out.ny, out.nz, out.prec};
    const double x[3] = {out.lx, out.lz, out.time};
    std::fwrite(h, sizeof(int), 4, fo);
    std::fwrite(x, sizeof(double), 3, fo);
  }
  std::vector<char> rin(size_t(in.plane_bytes())), rout(size_t(out.plane_bytes()));
  static const char* fname[3] = {"u", "v", "w"};
  for (int f = 0; f < out.nf; ++f) {
    Stats st;
    if (fseeko(fi, off_t(in.offset() + in.field_bytes() * f), SEEK_SET)) fail("seek in " + in.bin);
    for (int j = 0; j < in.ny; ++j) {
      if (std::fread(rin.data(), 1, rin.size(), fi) != rin.size()) fail("read error in " + in.bin);
      C.plane(rin, rout, st);
      if (std::fwrite(rout.data(), 1, rout.size(), fo) != rout.size()) fail("write error in " + tmp + " (disk full?)");
    }
    const std::string nm = f < 3 ? fname[f] : "s" + std::to_string(f - 2);
    if (C.direct()) std::printf("  %-3s copied\n", nm.c_str());
    else {
      const double t = std::max(st.tot, 1e-300);
      std::printf("  %-3s energy fraction dropped: Nyquist kz %.2e, truncated modes %.2e", nm.c_str(), st.nyq / t, st.cut / t);
      if (!hout) std::printf("; max |Im| at kx = 0 (set to 0) %.1e", st.imag0);
      std::printf("\n");
    }
    std::fflush(stdout);
  }
  std::fclose(fi);
  if (std::fclose(fo)) fail("write error in " + tmp);
  if (std::rename(tmp.c_str(), out.bin.c_str())) fail("cannot rename " + tmp + " to " + out.bin);
  if (hout) {
    std::FILE* h = std::fopen(out.hdr.c_str(), "w");
    if (!h) fail("cannot write " + out.hdr);
    std::fprintf(h, "# nx ny nz nscal lx lz time istep dt grid_type grid_alpha fd_order\n");
    std::fprintf(h, " %d %d %d %d\n", out.nx, out.ny, out.nz, out.nf - 3);
    std::fprintf(h, " %.17E %.17E\n", out.lx, out.lz);
    std::fprintf(h, " %.17E %ld %.17E\n", out.time, out.step, out.dt);
    std::fprintf(h, " %s %.17E %d\n", out.grid_type.c_str(), out.grid_alpha, out.fd_order);
    std::fclose(h);
    std::printf("written %s and %s\n", out.bin.c_str(), out.hdr.c_str());
  } else {
    std::printf("written %s\n", out.bin.c_str());
  }
  return 0;
}
