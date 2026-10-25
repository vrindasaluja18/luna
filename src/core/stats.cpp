#include "luna/core/stats.hpp"

#include <algorithm>

namespace luna {

Scalar RunningStats::percentile(std::vector<Scalar> v, Scalar p) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const Scalar idx = clamp(p, 0.0, 1.0) * static_cast<Scalar>(v.size() - 1);
  const Size lo = static_cast<Size>(idx);
  const Size hi = std::min(lo + 1, v.size() - 1);
  const Scalar frac = idx - static_cast<Scalar>(lo);
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}

}  // namespace luna
