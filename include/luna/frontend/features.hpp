// LUNA -- BRIEF descriptors and KLT optical flow.
//
// BRIEF is implemented with a *fixed, seeded* sampling pattern so descriptors
// are reproducible across machines and runs. Descriptors are 256 bits packed
// into four 64-bit words, which is what makes the descriptor DB cheap: one
// XOR-popcount per candidate pair, and an LRU index that costs 32 bytes per
// keyframe.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "luna/core/stats.hpp"
#include "luna/core/types.hpp"
#include "luna/frontend/fast.hpp"
#include "luna/sim/sensors.hpp"

namespace luna {
namespace frontend {

/// 256-bit descriptor as four 64-bit words.
struct Descriptor {
  std::array<std::uint64_t, 4> w{{0, 0, 0, 0}};

  int hamming(const Descriptor& o) const {
#if defined(__GNUC__)
    return __builtin_popcountll(w[0] ^ o.w[0]) + __builtin_popcountll(w[1] ^ o.w[1]) +
           __builtin_popcountll(w[2] ^ o.w[2]) + __builtin_popcountll(w[3] ^ o.w[3]);
#else
    int d = 0;
    for (int i = 0; i < 4; ++i) d += popcount64(w[i] ^ o.w[i]);
    return d;
#endif
  }
};

/// A feature: pixel position, descriptor, and (when known) a depth in metres.
struct Feature {
  Scalar x = 0, y = 0;
  Descriptor desc;
  Scalar depth = -1;      ///< <0 when unknown
  Scalar response = 0;
  int track_id = -1;      ///< assigned by the tracker
  /// Ground-truth landmark id when the simulator provides one (evaluation only).
  std::uint32_t gt_id = 0;
  bool has_gt = false;
};

/// Precomputed BRIEF sampling pattern. Generating it once and sharing it
/// across every frame is what keeps descriptor extraction allocation-free.
class BriefPattern {
 public:
  explicit BriefPattern(std::uint64_t seed = 0xB4EEF1234ull, int patch = 31, int pairs = 256);

  /// Compute a descriptor for the pixel at (x, y). Returns false when the
  /// patch would leave the image.
  bool compute(const sim::Image& img, int x, int y, Descriptor* out) const;

  int patch() const { return patch_; }
  int halfPatch() const { return patch_ / 2; }
  Size bytes() const { return offsets_.size() * sizeof(std::int16_t); }

 private:
  int patch_ = 31;
  int pairs_ = 256;
  /// Packed (dx1, dy1, dx2, dy2) as int8 offsets, 4 bytes per pair.
  std::vector<std::array<std::int8_t, 4>> offsets_;
};

struct BriefOptions {
  int threshold = 10;     ///< reject features with no intensity contrast
  int max_features = 400;
};

/// Detect + describe in one pass. `gt_lookup` (optional) maps a pixel to the
/// ground-truth landmark id, which the evaluator uses to score the frontend.
void detectAndDescribe(const sim::Image& img, const BriefPattern& pattern, const FastOptions& fast,
                       const BriefOptions& brief, std::vector<Feature>* out);

/// KLT tracking options.
struct KltOptions {
  int window = 15;          ///< patch side length (odd)
  int max_level = 3;        ///< pyramid levels (4 levels = 8x scale range)
  int max_iter = 12;
  Scalar eps = 0.01;        ///< convergence threshold in pixels
  Scalar min_eigenvalue = 1e-4;
  Scalar max_residual = 80.0;
  Scalar fb_threshold = 2.0;  ///< forward-backward error rejection
  int max_track = 400;      ///< cap on features carried between frames
  /// Maximum Gauss-Newton step per pyramid level, in that level's pixels.
  /// This is the guard against a coarse level jumping a texture period.
  Scalar max_step_per_level = 2.0;
};

/// Pyramid-based Lucas-Kanade tracker with forward-backward outlier rejection.
///
/// Working in both directions is what makes tracking robust to the occlusion
/// and blackout faults: a track whose backward flow does not return to its
/// origin is discarded rather than poisoning the pose estimate.
class KltTracker {
 public:
  explicit KltTracker(const KltOptions& opt = KltOptions());

  /// Track features from `prev` to `cur`. On success, `out[i]` holds the new
  /// position and the input feature's descriptor is carried over. Features
  /// that fail are dropped, so `out` may be shorter than `in`.
  void track(const sim::Image& prev, const sim::Image& cur, const std::vector<Feature>& in,
             std::vector<Feature>* out);

  Size tracked() const { return tracked_; }
  Size rejected() const { return rejected_; }
  /// Rejection counters, in the order the gates run.
  Size rejected_coarse() const { return rejected_coarse_; }
  Size rejected_forward_backward() const { return rejected_fb_; }
  struct Stats { Size tracked = 0, coarse = 0, fb = 0; };
  Stats stats() const { return {tracked_, rejected_coarse_, rejected_fb_}; }

 private:
  void buildPyramid(const sim::Image& img);
  void buildCurrentPyramid(const sim::Image& img);
  /// Refine one point at one pyramid level.
  ///
  /// `x_ref` is where the reference patch lives in the *previous* image and
  /// must stay fixed for the whole track; `x_cur` is the running estimate of
  /// where that same point appears in the *current* image. Both are in this
  /// level's pixel coordinates. Keeping them separate is what makes
  /// coarse-to-fine work: each level solves for the *residual* displacement
  /// left over by the coarser level. Building the reference patch at the
  /// running estimate instead re-solves for the full displacement at every
  /// level, and the corrections then compound instead of shrinking.
  bool trackPoint(int level, const Vec2& x_ref, const Vec2& x_cur, Vec2* out);
  void downsample(const sim::Image& src, sim::Image* dst);
  Scalar sampleLevel(int level, Scalar x, Scalar y) const;
  Scalar samplePrev(int level, Scalar x, Scalar y) const;

  KltOptions opt_;
  std::vector<sim::Image> prev_pyr_, cur_pyr_;
  // Scratch buffers reused across every trackPoint() call so the tracker
  // performs zero heap allocation per feature.
  std::vector<Scalar> ref_, gx_, gy_;
  Size tracked_ = 0, rejected_ = 0;
  Size rejected_coarse_ = 0, rejected_fb_ = 0;
};

/// SIMD-friendly block matcher for stereo depth: for each left feature, search
/// a window in the right image along the epipolar line and return the best
/// disparity as a depth in metres.
struct StereoOptions {
  int max_disparity = 80;     ///< covers the 0.4-25 m depth band
  int patch = 4;              ///< half patch (9x9 window)
  Scalar sad_ratio = 0.6;     ///< best must beat second-best by this factor
  Scalar lr_threshold = 1.0;  ///< left-right consistency tolerance, pixels
  Scalar min_depth = 0.4;
  Scalar max_depth = 25.0;
};

void stereoDepth(const sim::Image& left, const sim::Image& right, const sim::CameraModel& cam,
                 const StereoOptions& opt, std::vector<Feature>* features);

}  // namespace frontend
}  // namespace luna
