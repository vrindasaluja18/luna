// SLAMCore-8G -- running statistics and deterministic RNG.
#pragma once

#include <cstdint>
#include <vector>

#include "slamcore/core/types.hpp"

namespace sc {

/// Welford accumulator: numerically stable mean/variance in one pass, O(1) RAM.
class RunningStats {
 public:
  void add(Scalar x) {
    ++n_;
    const Scalar d = x - mean_;
    mean_ += d / static_cast<Scalar>(n_);
    m2_ += d * (x - mean_);
    if (n_ == 1 || x < min_) min_ = x;
    if (n_ == 1 || x > max_) max_ = x;
  }
  Size count() const noexcept { return n_; }
  Scalar mean() const noexcept { return n_ ? mean_ : 0.0; }
  Scalar variance() const noexcept { return n_ > 1 ? m2_ / static_cast<Scalar>(n_ - 1) : 0.0; }
  Scalar stddev() const noexcept { return std::sqrt(variance()); }
  Scalar rms() const noexcept { return n_ ? std::sqrt((m2_ + static_cast<Scalar>(n_) * mean_ * mean_) / static_cast<Scalar>(n_)) : 0.0; }
  Scalar min() const noexcept { return n_ ? min_ : 0.0; }
  Scalar max() const noexcept { return n_ ? max_ : 0.0; }

  /// Percentile from a sorted copy of the samples (O(n log n), bench-only).
  static Scalar percentile(std::vector<Scalar> v, Scalar p);

 private:
  Size n_ = 0;
  Scalar mean_ = 0, m2_ = 0;
  Scalar min_ = 0, max_ = 0;
};

/// xorshift128+ : fast, tiny, and bit-reproducible across platforms, which is
/// what the deterministic replay mode depends on.
class Rng {
 public:
  explicit Rng(std::uint64_t s0 = 0x9E3779B97F4A7C15ull) { seed(s0); }

  void seed(std::uint64_t s) {
    // SplitMix64 expansion of the seed into the 128-bit state.
    s_ = splitmix(&s);
    s_ = s ? s : 0x9E3779B97F4A7C15ull;
    s1_ = splitmix(&s);
    if (!s_) s_ = 0xBF58476D1CE4E5B9ull;
    if (!s1_) s1_ = 0x94D049BB133111EBull;
  }

  std::uint64_t next() {
    std::uint64_t x = s_;
    const std::uint64_t y = s1_;
    s_ = y;
    x ^= x << 23;
    s1_ = x ^ y ^ (x >> 17) ^ (y >> 26);
    return s1_ + y;
  }

  /// Uniform in [0, 1).
  Scalar uniform() { return static_cast<Scalar>(next() >> 11) * (1.0 / 9007199254740992.0); }
  /// Uniform in [lo, hi).
  Scalar uniform(Scalar lo, Scalar hi) { return lo + (hi - lo) * uniform(); }
  /// Standard normal via Box-Muller (deterministic, no rejection loop).
  Scalar gaussian() {
    if (have_spare_) {
      have_spare_ = false;
      return spare_;
    }
    const Scalar u1 = std::max(uniform(), 1e-300);
    const Scalar u2 = uniform();
    const Scalar r = std::sqrt(-2.0 * std::log(u1));
    spare_ = r * std::sin(kTwoPi * u2);
    have_spare_ = true;
    return r * std::cos(kTwoPi * u2);
  }
  /// Gaussian with zero mean and given standard deviation.
  Scalar gaussian(Scalar sigma) { return sigma * gaussian(); }

 private:
  static std::uint64_t splitmix(std::uint64_t* x) {
    std::uint64_t z = (*x += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  std::uint64_t s_ = 0, s1_ = 0;
  Scalar spare_ = 0;
  bool have_spare_ = false;
};

}  // namespace sc
