// rsio.hpp: reading HiCHAN and FDCHAN restart files (headers, planes) and an FFT of any length; used by
// rsconv (conversion) and rsspec (spectra). File formats: see README.md.
#pragma once
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

[[noreturn]] inline void fail(const std::string& m) {
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

inline long long file_size(const std::string& f) {
  std::FILE* p = std::fopen(f.c_str(), "rb");
  if (!p) fail("cannot open " + f);
  std::fseek(p, 0, SEEK_END);
  const long long s = (long long)ftello(p);
  std::fclose(p);
  return s;
}

inline std::string strip(const std::string& name) {  // HiCHAN restart name without .bin / .hdr
  for (const char* e : {".bin", ".hdr"})
    if (name.size() > 4 && name.compare(name.size() - 4, 4, e) == 0) return name.substr(0, name.size() - 4);
  return name;
}

inline Desc read_hichan(const std::string& name) {
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

inline Desc read_fdchan(const std::string& name) {
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

inline bool is_hichan(const std::string& name) {
  std::ifstream h(strip(name) + ".hdr");
  return bool(h);
}

inline void print_desc(const Desc& d, const char* what) {
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
inline int signed_mode(int mi, int n) { return mi <= n / 2 ? mi : mi - n; }

inline void parallel_for(int n, int nthreads, const std::function<void(int, int)>& f) {  // f(item, thread)
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

