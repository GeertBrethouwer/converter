// rsspec: one-dimensional streamwise and spanwise spectra and two-dimensional (streamwise-spanwise) spectra of the
// velocity components and scalars at chosen wall-normal positions, from HiCHAN or FDCHAN restart files. Several
// files (time samples) are averaged; --sym also averages the mirror position in the upper half of the channel.
//
//   rsspec [options] <file> [<file> ...]          (see usage() or README.md)
//
// Definitions (f' = f - its x-z mean at the given y; one-sided in kx and kz, i.e. kx, kz >= 0):
//   E2(kx, kz) dkx dkz   the energy of f' in the band, sum over all bands = <f'f'>
//   Ex(kx) = sum_kz E2(kx, kz) dkz,   Ez(kz) = sum_kx E2(kx, kz) dkx,   sum Ex dkx = sum Ez dkz = <f'f'>
// with dkx = 2 pi/lx, dkz = 2 pi/lz. From the Fourier coefficients c(kx_i, kz_m) (amplitudes, both codes):
//   E2(kx_i, |kz_m|) dkx dkz = a_i sum_{+-m} |c(i, m)|^2, a_0 = 1, a_i = 2 (the conjugates, kx < 0),
// the kz Nyquist mode of an even nz counted once; the (0, 0) mode (the mean) is left out.
#include "rsio.hpp"

#include <map>
#include <sys/stat.h>

// y-grid as in both codes (y[0] = -1, y[ny-1] = 1)
static std::vector<double> ygrid(int ny, const std::string& type, double alpha, const std::string& file) {
  std::vector<double> y(static_cast<size_t>(ny));
  if (!file.empty()) {
    std::ifstream f(file);
    if (!f) fail("cannot open the y-grid file " + file);
    for (int j = 0; j < ny; ++j)
      if (!(f >> y[size_t(j)])) fail("the y-grid file " + file + " has fewer than ny = " + std::to_string(ny) + " values");
    return y;
  }
  for (int j = 0; j < ny; ++j) {
    const double xi = double(j) / (ny - 1);
    if (type == "tanh") y[size_t(j)] = std::tanh(alpha * (xi - 0.5)) / std::tanh(0.5 * alpha);
    else if (type == "uniform") y[size_t(j)] = -1 + 2 * xi;
    else if (type == "cosine") y[size_t(j)] = -std::cos(M_PI * xi);
    else fail("unknown grid type '" + type + "' (tanh, uniform, cosine; or --yfile)");
  }
  y[0] = -1;
  y[size_t(ny - 1)] = 1;
  return y;
}

static std::vector<double> numbers(const std::string& s) {
  std::vector<double> v;
  std::string t = s;
  for (char& c : t)
    if (c == ',') c = ' ';
  std::istringstream in(t);
  double x;
  while (in >> x) v.push_back(x);
  if (v.empty()) fail("empty list '" + s + "'");
  return v;
}

// one y-plane of one field as Fourier coefficients S[mi * nkx + i] (mi: kz index in FFT order)
class PlaneReader {
 public:
  PlaneReader(const Desc& d, int nth) : d_(d), nth_(nth) {
    f_ = std::fopen(d.bin.c_str(), "rb");
    if (!f_) fail("cannot open " + d.bin);
    if (d.fmt == "fdchan") {
      fft_.reset(new FFT(d.nz));
      work_.assign(size_t(nth), std::vector<cd>(fft_->work_size() + size_t(d.nz)));
    }
    raw_.resize(size_t(d.plane_bytes()));
  }
  ~PlaneReader() { std::fclose(f_); }
  void read(int field, int j, std::vector<cd>& S) {
    const long long off = d_.offset() + d_.field_bytes() * field + d_.plane_bytes() * j;
    if (fseeko(f_, off_t(off), SEEK_SET) || std::fread(raw_.data(), 1, raw_.size(), f_) != raw_.size())
      fail("read error in " + d_.bin);
    const int nkx = d_.nkx(), nz = d_.nz;
    S.assign(size_t(nkx) * nz, cd(0, 0));
    if (d_.fmt == "hichan") {
      std::memcpy(S.data(), raw_.data(), raw_.size());
      return;
    }
    parallel_for(nkx, nth_, [&](int i, int t) {
      cd* line = work_[size_t(t)].data();
      cd* w = line + nz;
      for (int k = 0; k < nz; ++k) {
        const long long e = 2 * ((long long)k * nkx + i);
        line[k] = cd(get(e), get(e + 1));
      }
      fft_->exec(line, -1, w);
      for (int m = 0; m < nz; ++m) S[size_t(m) * nkx + i] = line[m] / double(nz);
    });
  }

 private:
  Desc d_;
  int nth_;
  std::FILE* f_;
  std::unique_ptr<FFT> fft_;
  std::vector<std::vector<cd>> work_;
  std::vector<char> raw_;
  double get(long long k) const {
    if (d_.prec == 8) { double x; std::memcpy(&x, raw_.data() + k * 8, 8); return x; }
    float x; std::memcpy(&x, raw_.data() + k * 4, 4); return double(x);
  }
};

static void usage() {
  std::puts(
      "usage: rsspec [options] <restart> [<restart> ...]\n"
      "  Spectra of u, v, w and the scalars at chosen y-positions, averaged over the files (HiCHAN restarts as\n"
      "  <name> or <name>.bin/.hdr, FDCHAN restarts as <file>.bin; all with the same nx, ny, nz, lx, lz).\n"
      "positions (one of):\n"
      "  --yplus a,b,...     distances from the bottom wall in wall units (needs --retau)\n"
      "  --y a,b,...         y-values (-1 = bottom wall, 0 = centre, 1 = top wall)\n"
      "  --j a,b,...         grid indices (0 = bottom wall)\n"
      "  The nearest grid points are used.\n"
      "options:\n"
      "  --retau R           friction Reynolds number (for --yplus and the wall units in the output)\n"
      "  --grid TYPE ALPHA   y-grid: tanh, uniform or cosine and the stretching (default: from a HiCHAN header;\n"
      "                      FDCHAN files do not contain it)\n"
      "  --yfile FILE        y-grid from a file (ny values, e.g. HiCHAN's grid_file)\n"
      "  --sym               also average the mirror points in the upper half (j and ny-1-j)\n"
      "  --fields u,v,w,s1   fields (default: all)\n"
      "  --out DIR           output directory (default: spectra)\n"
      "  --threads N         threads (default: all cores)\n"
      "output (DIR/): per position j: spec_x_jNNNN.dat (kx, lambda_x, Ex of every field), spec_z_jNNNN.dat (kz,\n"
      "  lambda_z, Ez), spec_2d_jNNNN.dat (kx, kz, lambda_x, lambda_z, E2; gnuplot blocks per kx), and summary.dat\n"
      "  (y, y+, the variance of every field). E are densities: sum Ex dkx = sum Ez dkz = sum E2 dkx dkz = <f'f'>.");
}

int main(int argc, char** argv) {
  std::vector<std::string> files;
  std::vector<double> yplus, yv, jv;
  double retau = 0, alpha = NAN;
  std::string gtype, yfile, outdir = "spectra", fsel;
  bool sym = false;
  int nth = int(std::max(1u, std::thread::hardware_concurrency()));
  for (int a = 1; a < argc; ++a) {
    const std::string o = argv[a];
    auto val = [&]() -> std::string {
      if (a + 1 >= argc) fail("missing value for " + o);
      return argv[++a];
    };
    if (o == "--yplus") yplus = numbers(val());
    else if (o == "--y") yv = numbers(val());
    else if (o == "--j") jv = numbers(val());
    else if (o == "--retau") retau = std::atof(val().c_str());
    else if (o == "--grid") { gtype = val(); alpha = std::atof(val().c_str()); }
    else if (o == "--yfile") yfile = val();
    else if (o == "--sym") sym = true;
    else if (o == "--fields") fsel = val();
    else if (o == "--out") outdir = val();
    else if (o == "--threads") nth = std::max(1, std::atoi(val().c_str()));
    else if (o == "-h" || o == "--help") { usage(); return 0; }
    else if (!o.empty() && o[0] == '-') fail("unknown option " + o);
    else files.push_back(o);
  }
  if (files.empty()) { usage(); return 1; }
  if (int(!yplus.empty()) + int(!yv.empty()) + int(!jv.empty()) != 1) fail("give the positions with one of --yplus, --y, --j");
  if (!yplus.empty() && retau <= 0) fail("--yplus needs --retau");

  std::vector<Desc> D;
  for (const std::string& f : files) D.push_back(is_hichan(f) ? read_hichan(f) : read_fdchan(f));
  const Desc& d0 = D[0];
  for (const Desc& d : D)
    if (d.nx != d0.nx || d.ny != d0.ny || d.nz != d0.nz || std::abs(d.lx - d0.lx) > 1e-10 * d0.lx ||
        std::abs(d.lz - d0.lz) > 1e-10 * d0.lz)
      fail(d.bin + " has another grid or box than " + d0.bin);
  int nf = d0.nf;
  for (const Desc& d : D) nf = std::min(nf, d.nf);
  const int ny = d0.ny, nkx = d0.nkx(), nz = d0.nz, nzh = nz / 2;  // kz indices 0..nzh (nzh = Nyquist if even)
  const double dkx = 2 * M_PI / d0.lx, dkz = 2 * M_PI / d0.lz;

  // y-grid
  if (gtype.empty() && yfile.empty())
    for (const Desc& d : D)
      if (d.fmt == "hichan" && d.grid_type != "unknown") { gtype = d.grid_type; alpha = d.grid_alpha; break; }
  std::vector<double> y;
  if (!gtype.empty() || !yfile.empty()) y = ygrid(ny, gtype, alpha, yfile);
  else if (jv.empty()) fail("the y-grid is needed for --y/--yplus: give --grid TYPE ALPHA or --yfile");

  // positions -> grid indices (lower half for --yplus; --sym adds the mirror)
  std::vector<int> J;
  auto nearest = [&](double yy) {
    int b = 0;
    for (int j = 1; j < ny; ++j)
      if (std::abs(y[size_t(j)] - yy) < std::abs(y[size_t(b)] - yy)) b = j;
    return b;
  };
  for (double p : yplus) J.push_back(nearest(-1 + p / retau));
  for (double p : yv) J.push_back(nearest(p));
  for (double p : jv) {
    const int j = int(std::lround(p));
    if (j < 0 || j >= ny) fail("--j " + std::to_string(j) + " outside 0.." + std::to_string(ny - 1));
    J.push_back(j);
  }
  std::sort(J.begin(), J.end());
  J.erase(std::unique(J.begin(), J.end()), J.end());

  // fields
  std::vector<std::string> names = {"u", "v", "w"};
  for (int n = 1; n <= nf - 3; ++n) names.push_back("s" + std::to_string(n));
  std::vector<int> F;
  if (fsel.empty()) for (int f = 0; f < nf; ++f) F.push_back(f);
  else {
    std::string t = fsel;
    for (char& c : t) if (c == ',') c = ' ';
    std::istringstream in(t);
    std::string nm;
    while (in >> nm) {
      const auto it = std::find(names.begin(), names.end(), nm);
      if (it == names.end()) fail("unknown field " + nm + " (the files have u, v, w and " + std::to_string(nf - 3) + " scalar(s))");
      F.push_back(int(it - names.begin()));
    }
  }
  const int NF = int(F.size());

  std::printf("rsspec: %zu file(s), nx = %d, ny = %d, nz = %d, lx = %.6g, lz = %.6g, fields:", files.size(), d0.nx, ny,
              nz, d0.lx, d0.lz);
  for (int f : F) std::printf(" %s", names[size_t(f)].c_str());
  std::printf("\n");
  mkdir(outdir.c_str(), 0755);

  // E2[(p * NF + q) * (nkx * (nzh + 1)) + i * (nzh + 1) + m]: energy in band (i, m) (not yet a density)
  const size_t nb = size_t(nkx) * (nzh + 1);
  std::vector<double> E2(J.size() * size_t(NF) * nb, 0.0);
  std::vector<double> mean(J.size() * size_t(NF), 0.0);
  std::vector<cd> S;
  int nsamp = 0;
  std::vector<double> times;
  for (const Desc& d : D) {
    PlaneReader R(d, nth);
    times.push_back(d.time);
    for (size_t p = 0; p < J.size(); ++p) {
      std::vector<int> js = {J[p]};
      if (sym && ny - 1 - J[p] != J[p]) js.push_back(ny - 1 - J[p]);
      for (int q = 0; q < NF; ++q)
        for (int j : js) {
          R.read(F[size_t(q)], j, S);
          double* e = &E2[(p * NF + size_t(q)) * nb];
          for (int mi = 0; mi < nz; ++mi) {
            const int m = std::abs(signed_mode(mi, nz));
            for (int i = 0; i < nkx; ++i) {
              if (i == 0 && m == 0) continue;
              e[size_t(i) * (nzh + 1) + m] += (i == 0 ? 1.0 : 2.0) * std::norm(S[size_t(mi) * nkx + i]) / double(js.size());
            }
          }
          if (j == J[p]) mean[p * NF + size_t(q)] += S[0].real();  // the mean at the point itself (not the mirror)
        }
    }
    ++nsamp;
    std::printf("  %s (t = %.6g)\n", d.bin.c_str(), d.time);
    std::fflush(stdout);
  }
  for (double& x : E2) x /= nsamp;
  for (double& x : mean) x /= nsamp;

  // output
  auto hdr = [&](std::FILE* o, size_t p, const char* what) {
    const int j = J[p];
    std::fprintf(o, "# %s; rsspec, %d sample(s), t =", what, nsamp);
    for (size_t k = 0; k < times.size() && k < 20; ++k) std::fprintf(o, " %.6g", times[k]);
    if (times.size() > 20) std::fprintf(o, " ...");
    std::fprintf(o, "\n# j = %d%s", j, sym && ny - 1 - j != j ? (" (and its mirror " + std::to_string(ny - 1 - j) + ")").c_str() : "");
    if (!y.empty()) {
      std::fprintf(o, ", y = %.8g, 1 + y = %.6g", y[size_t(j)], 1 + y[size_t(j)]);
      if (retau > 0) std::fprintf(o, ", y+ = %.6g (Re_tau = %g)", (1 + y[size_t(j)]) * retau, retau);
    }
    std::fprintf(o, "\n# nx = %d, nz = %d, lx = %.10g, lz = %.10g; densities: sum E dk = <f'f'>\n# variance <f'f'>:", d0.nx,
                 nz, d0.lx, d0.lz);
    for (int q = 0; q < NF; ++q) {
      double v = 0;
      for (size_t b = 0; b < nb; ++b) v += E2[(p * NF + size_t(q)) * nb + b];
      std::fprintf(o, " %s %.8e", names[size_t(F[size_t(q)])].c_str(), v);
    }
    std::fprintf(o, "\n");
  };
  char fn[512];
  std::FILE* sm = std::fopen((outdir + "/summary.dat").c_str(), "w");
  if (!sm) fail("cannot write in " + outdir);
  std::fprintf(sm, "# rsspec: %d sample(s)%s; variance <f'f'> and mean (at j, not averaged with the mirror) of every field\n#   j            y           y+", nsamp,
               sym ? ", with the mirror points" : "");
  for (int q = 0; q < NF; ++q) std::fprintf(sm, "   var_%-10s", names[size_t(F[size_t(q)])].c_str());
  for (int q = 0; q < NF; ++q) std::fprintf(sm, "  mean_%-10s", names[size_t(F[size_t(q)])].c_str());
  std::fprintf(sm, "\n");
  for (size_t p = 0; p < J.size(); ++p) {
    const int j = J[p];
    std::vector<std::vector<double>> ex(size_t(NF), std::vector<double>(size_t(nkx), 0.0)),
        ez(size_t(NF), std::vector<double>(size_t(nzh + 1), 0.0));
    for (int q = 0; q < NF; ++q)
      for (int i = 0; i < nkx; ++i)
        for (int m = 0; m <= nzh; ++m) {
          const double e = E2[(p * NF + size_t(q)) * nb + size_t(i) * (nzh + 1) + m];
          ex[size_t(q)][size_t(i)] += e;
          ez[size_t(q)][size_t(m)] += e;
        }
    // streamwise
    std::snprintf(fn, sizeof(fn), "%s/spec_x_j%04d.dat", outdir.c_str(), j);
    std::FILE* o = std::fopen(fn, "w");
    hdr(o, p, "streamwise spectra Ex(kx)");
    std::fprintf(o, "#            kx        lambda_x");
    for (int q = 0; q < NF; ++q) std::fprintf(o, "       Ex_%-6s", names[size_t(F[size_t(q)])].c_str());
    std::fprintf(o, "\n");
    for (int i = 0; i < nkx; ++i) {
      std::fprintf(o, "%16.8e%16.8e", i * dkx, i ? d0.lx / i : INFINITY);
      for (int q = 0; q < NF; ++q) std::fprintf(o, "%16.8e", ex[size_t(q)][size_t(i)] / dkx);
      std::fprintf(o, "\n");
    }
    std::fclose(o);
    // spanwise
    std::snprintf(fn, sizeof(fn), "%s/spec_z_j%04d.dat", outdir.c_str(), j);
    o = std::fopen(fn, "w");
    hdr(o, p, "spanwise spectra Ez(kz)");
    std::fprintf(o, "#            kz        lambda_z");
    for (int q = 0; q < NF; ++q) std::fprintf(o, "       Ez_%-6s", names[size_t(F[size_t(q)])].c_str());
    std::fprintf(o, "\n");
    for (int m = 0; m <= nzh; ++m) {
      std::fprintf(o, "%16.8e%16.8e", m * dkz, m ? d0.lz / m : INFINITY);
      for (int q = 0; q < NF; ++q) std::fprintf(o, "%16.8e", ez[size_t(q)][size_t(m)] / dkz);
      std::fprintf(o, "\n");
    }
    std::fclose(o);
    // two-dimensional
    std::snprintf(fn, sizeof(fn), "%s/spec_2d_j%04d.dat", outdir.c_str(), j);
    o = std::fopen(fn, "w");
    hdr(o, p, "two-dimensional spectra E2(kx, kz); a block per kx (gnuplot splot / numpy reshape (nx/2, nz/2+1))");
    std::fprintf(o, "#            kx              kz        lambda_x        lambda_z");
    for (int q = 0; q < NF; ++q) std::fprintf(o, "       E2_%-6s", names[size_t(F[size_t(q)])].c_str());
    std::fprintf(o, "\n");
    for (int i = 0; i < nkx; ++i) {
      for (int m = 0; m <= nzh; ++m) {
        std::fprintf(o, "%16.8e%16.8e%16.8e%16.8e", i * dkx, m * dkz, i ? d0.lx / i : INFINITY, m ? d0.lz / m : INFINITY);
        for (int q = 0; q < NF; ++q)
          std::fprintf(o, "%16.8e", E2[(p * NF + size_t(q)) * nb + size_t(i) * (nzh + 1) + m] / (dkx * dkz));
        std::fprintf(o, "\n");
      }
      std::fprintf(o, "\n");
    }
    std::fclose(o);
    // summary
    std::fprintf(sm, "%5d %12.8f %12.5g", j, y.empty() ? NAN : y[size_t(j)],
                 (!y.empty() && retau > 0) ? (1 + y[size_t(j)]) * retau : NAN);
    for (int q = 0; q < NF; ++q) {
      double v = 0;
      for (int i = 0; i < nkx; ++i) v += ex[size_t(q)][size_t(i)];
      std::fprintf(sm, " %16.8e", v);
    }
    for (int q = 0; q < NF; ++q) std::fprintf(sm, " %16.8e", mean[p * NF + size_t(q)]);
    std::fprintf(sm, "\n");
    std::printf("  j = %4d", j);
    if (!y.empty()) std::printf(", y = %.6f", y[size_t(j)]);
    if (!y.empty() && retau > 0) std::printf(", y+ = %.2f", (1 + y[size_t(j)]) * retau);
    std::printf(": %s/spec_{x,z,2d}_j%04d.dat\n", outdir.c_str(), j);
  }
  std::fclose(sm);
  std::printf("written %s/summary.dat\n", outdir.c_str());
  return 0;
}
