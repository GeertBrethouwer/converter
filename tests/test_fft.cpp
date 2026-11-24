// FFT of rsconv against a direct DFT, lengths 1..130 and some production sizes (also Bluestein for large primes)
#define main rsconv_main
#include "../src/rsconv.cpp"  // (FFT from rsio.hpp)
#undef main
#include <random>
int main() {
  std::mt19937 g(1);
  std::normal_distribution<double> N(0, 1);
  std::vector<int> ns;
  for (int n = 1; n <= 130; ++n) ns.push_back(n);
  for (int n : {1920, 2688, 2689, 2880, 1031 * 2}) ns.push_back(n);
  double worst = 0;
  for (int n : ns) {
    FFT F(n);
    std::vector<cd> x(static_cast<size_t>(n)), y;
    for (auto& a : x) a = cd(N(g), N(g));
    for (int sign : {-1, 1}) {
      y = x;
      F.exec(y.data(), sign);
      double e = 0, m = 0;
      for (int k = 0; k < n; ++k) {
        cd s(0, 0);
        for (int j = 0; j < n; ++j) s += x[size_t(j)] * std::polar(1.0, sign * 2 * M_PI * double((long long)j * k % n) / n);
        e = std::max(e, std::abs(s - y[size_t(k)]));
        m = std::max(m, std::abs(s));
      }
      worst = std::max(worst, e / m);
      if (e / m > 1e-12) { std::printf("n = %d sign %d: error %.2e WRONG\n", n, sign, e / m); return 1; }
    }
  }
  std::printf("FFT vs direct DFT, n = 1..130, 1920, 2688, 2689 (prime), 2880, 2062: max relative error %.1e: OK\n", worst);
  return 0;
}
