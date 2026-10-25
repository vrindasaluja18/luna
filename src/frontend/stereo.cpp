// SLAMCore-8G -- stereo block matching for sparse feature depth.
//
// Dense depth is not needed: the frontend only needs a depth for the features
// it tracks, so this is a per-feature SAD search along the epipolar line. Three
// gates keep the output usable: a uniqueness ratio, left-right consistency, and
// a depth range. Without them a feature near a wall silently reports the wall's
// depth, which is a confident but wrong measurement.
#include "slamcore/frontend/features.hpp"

#include <limits>

namespace sc {
namespace frontend {

// ---------------------------------------------------------------------------
// Stereo block matching
// ---------------------------------------------------------------------------

void stereoDepth(const sim::Image& left, const sim::Image& right, const sim::CameraModel& cam,
                 const StereoOptions& opt, std::vector<Feature>* features) {
  const int h = opt.patch;

  // Sum of absolute differences for a window centred on `lx` in the left image
  // and `rx` in the right image. Returns INT_MAX when the window leaves either
  // image, which callers treat as "no such candidate".
  const auto sad = [&](int lx, int ly, int rx) {
    if (lx - h < 0 || ly - h < 0 || lx + h >= left.w || ly + h >= left.h) return INT_MAX;
    if (rx - h < 0 || rx + h >= right.w) return INT_MAX;
    int acc = 0;
    for (int dy = -h; dy <= h; ++dy)
      for (int dx = -h; dx <= h; ++dx)
        acc += std::abs(static_cast<int>(left.at(lx + dx, ly + dy)) -
                        static_cast<int>(right.at(rx + dx, ly + dy)));
    return acc;
  };

  // Best disparity for an anchor pixel, with the uniqueness gate and a
  // sub-pixel refinement by parabolic interpolation of the cost curve.
  //
  // `dir` selects which image holds the anchor: -1 anchors in the left image
  // and searches the right at `ax - d`; +1 anchors in the right image and
  // searches the left at `ax + d`. Both directions must be expressible because
  // the left-right consistency test below runs the matcher the other way round,
  // and a single-direction helper would silently compare the left image with
  // itself and "confirm" every match.
  //
  // Interpolation is not cosmetic: depth is inversely proportional to
  // disparity, so one whole pixel of quantization is already ~25% depth error
  // at 20 m, where the true disparity is only about 1.5 px.
  const auto match = [&](int ax, int ay, int dir, Scalar* out) {
    const auto cost_at = [&](int d) {
      return (dir < 0) ? sad(ax, ay, ax - d) : sad(ax + d, ay, ax);
    };
    Scalar best = std::numeric_limits<Scalar>::max();
    Scalar second = std::numeric_limits<Scalar>::max();
    int best_d = -1;
    for (int d = 1; d <= opt.max_disparity; ++d) {
      const int c = cost_at(d);
      if (c == INT_MAX) break;
      const Scalar s = static_cast<Scalar>(c);
      if (s < best) {
        second = best;
        best = s;
        best_d = d;
      } else if (s < second) {
        second = s;
      }
    }
    if (best_d <= 0) return false;
    // A flat cost curve means a textureless or repetitive match, which must be
    // rejected rather than fed to the filter.
    if (second > 0 && best > opt.sad_ratio * second) return false;
    Scalar d_sub = static_cast<Scalar>(best_d);
    if (best_d > 1 && best_d < opt.max_disparity) {
      const Scalar cm = static_cast<Scalar>(cost_at(best_d - 1));
      const Scalar c0 = static_cast<Scalar>(cost_at(best_d));
      const Scalar cp = static_cast<Scalar>(cost_at(best_d + 1));
      const Scalar denom = cm - 2.0 * c0 + cp;
      if (std::abs(denom) > 1e-6) {
        const Scalar delta = 0.5 * (cm - cp) / denom;
        if (std::abs(delta) < 1.0) d_sub += delta;
      }
    }
    *out = d_sub;
    return true;
  };

  for (Feature& f : *features) {
    const int lx = static_cast<int>(std::lround(f.x));
    const int ly = static_cast<int>(std::lround(f.y));
    Scalar disp = 0;
    if (!match(lx, ly, -1, &disp)) continue;
    if (disp <= 0.05) continue;

    // Left-right consistency: match back from the right image and require the
    // return trip to land on the original pixel. Without this, a feature near a
    // wall is matched to the *wall's* disparity, which looks like a perfectly
    // confident measurement but reports the wall's depth. The round-trip test
    // rejects exactly those ambiguous matches.
    const int rx = static_cast<int>(std::lround(static_cast<Scalar>(lx) - disp));
    if (rx < 0 || rx >= right.w) continue;
    Scalar back = 0;
    if (!match(rx, ly, +1, &back)) continue;
    if (std::abs((static_cast<Scalar>(rx) + back) - static_cast<Scalar>(lx)) > opt.lr_threshold)
      continue;

    const Scalar depth = cam.fx * cam.baseline / disp;
    if (depth < opt.min_depth || depth > opt.max_depth) continue;
    f.depth = depth;
  }
}

}  // namespace frontend
}  // namespace sc
