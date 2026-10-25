// SLAMCore-8G -- core scalar types, constants, small helpers.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace sc {

using Scalar = double;
using Index = std::int32_t;
using Size = std::size_t;

inline constexpr Scalar kPi = 3.14159265358979323846;
inline constexpr Scalar kTwoPi = 2.0 * kPi;
inline constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

/// Wrap an angle to (-pi, pi].
inline Scalar wrapPi(Scalar a) {
  a = std::fmod(a + kPi, kTwoPi);
  if (a < 0.0) a += kTwoPi;
  return a - kPi;
}

inline Scalar sq(Scalar x) { return x * x; }
inline Scalar clamp(Scalar x, Scalar lo, Scalar hi) { return x < lo ? lo : (x > hi ? hi : x); }
inline int sign(Scalar x) { return (x > 0.0) - (x < 0.0); }

}  // namespace sc
