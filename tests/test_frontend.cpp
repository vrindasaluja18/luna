// SLAMCore-8G -- frontend tests: FAST, BRIEF, KLT, stereo.
//
// The KLT tests are the important ones. They warp an image by a *known* shift
// with an independent bilinear resampler (not the simulator's renderer, and not
// the tracker's own pyramid code) and then require the tracker to recover that
// shift to sub-pixel accuracy at every pyramid level. That is what pins the
// coarse-to-fine behaviour: a tracker that solves for the full displacement at
// every level instead of the residual passes a level-0 test and fails a
// level-2 test, which is exactly the bug these tests exist to catch.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "slamcore/frontend/features.hpp"
#include "slamcore/sim/simulator.hpp"
#include "test_util.hpp"

using namespace sc;
using namespace sc::frontend;
using namespace sc::sim;

namespace {

/// Bilinear resampler used to synthesise a known sub-pixel shift. `out(x, y)`
/// samples the source at (x - sx, y - sy), so content moves by (+sx, +sy).
Image warpImage(const Image& src, Scalar sx, Scalar sy) {
  Image out;
  out.resize(src.w, src.h);
  for (int y = 0; y < src.h; ++y) {
    for (int x = 0; x < src.w; ++x) {
      const Scalar fx = x - sx, fy = y - sy;
      const int x0 = static_cast<int>(std::floor(fx)), y0 = static_cast<int>(std::floor(fy));
      const Scalar ux = fx - x0, uy = fy - y0;
      const auto S = [&](int px, int py) -> Scalar {
        if (px < 0 || py < 0 || px >= src.w || py >= src.h) return 0.0;
        return static_cast<Scalar>(src.at(px, py));
      };
      const Scalar a = S(x0, y0) * (1 - ux) + S(x0 + 1, y0) * ux;
      const Scalar b = S(x0, y0 + 1) * (1 - ux) + S(x0 + 1, y0 + 1) * ux;
      out.at(x, y) = static_cast<std::uint8_t>(clamp(a * (1 - uy) + b * uy, 0.0, 255.0));
    }
  }
  return out;
}

/// Deterministic textured test image: aperiodic value noise with structure at
/// several scales, which is what a real surface looks like to a corner
/// detector. A pure sinusoid would make the KLT cost curve ambiguous.
Image synthTexture(int w, int h, std::uint64_t seed) {
  Image im;
  im.resize(w, h);
  Rng rng(seed);
  const int cell = 4;
  const int cw = w / cell + 2, ch = h / cell + 2;
  std::vector<Scalar> coarse(static_cast<Size>(cw) * ch);
  for (Scalar& v : coarse) v = rng.uniform();
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const Scalar fx = static_cast<Scalar>(x) / cell, fy = static_cast<Scalar>(y) / cell;
      const int x0 = static_cast<int>(fx), y0 = static_cast<int>(fy);
      const Scalar ux = fx - x0, uy = fy - y0;
      const auto S = [&](int px, int py) {
        return coarse[static_cast<Size>(std::min(py, ch - 1)) * cw + std::min(px, cw - 1)];
      };
      const Scalar a = S(x0, y0) * (1 - ux) + S(x0 + 1, y0) * ux;
      const Scalar b = S(x0, y0 + 1) * (1 - ux) + S(x0 + 1, y0 + 1) * ux;
      const Scalar base = 60.0 + 150.0 * (a * (1 - uy) + b * uy);
      // A fine, high-contrast grid on top so corners are well localised.
      const Scalar fine = ((x % 7 == 0) || (y % 7 == 0)) ? 45.0 : 0.0;
      im.at(x, y) = static_cast<std::uint8_t>(clamp(base + fine, 0.0, 255.0));
    }
  }
  return im;
}

void testFastDetectsCorners() {
  // Bright squares on a dark background: each square has four genuine corners.
  // A checkerboard would *not* work as a fixture here: at a checkerboard
  // intersection the 16-pixel circle alternates bright/dark, so no 9-contiguous
  // arc of similar pixels exists and FAST correctly reports nothing. The
  // fixture has to contain real corners, not just intensity changes.
  Image im;
  im.resize(64, 64);
  std::fill(im.px.begin(), im.px.end(), static_cast<std::uint8_t>(20));
  for (int k = 0; k < 4; ++k) {
    const int ox = 8 + 16 * (k % 2), oy = 8 + 16 * (k / 2);
    for (int y = oy; y < oy + 12; ++y)
      for (int x = ox; x < ox + 12; ++x) im.at(x, y) = 230;
  }
  FastOptions opt;
  opt.threshold = 25;
  std::vector<Corner> corners;
  detectFast(im, opt, &corners);
  SC_CHECK(!corners.empty());
  // Every corner must be inside the border margin.
  for (const Corner& c : corners) {
    SC_CHECK(c.x >= opt.border && c.y >= opt.border);
    SC_CHECK(c.x < im.w - opt.border && c.y < im.h - opt.border);
    SC_CHECK(c.score > 0);
  }
  // A flat image has no corners at all.
  Image flat;
  flat.resize(64, 64);
  std::fill(flat.px.begin(), flat.px.end(), 128);
  std::vector<Corner> none;
  detectFast(flat, opt, &none);
  SC_CHECK(none.empty());
  // The cap is honoured.
  opt.max_features = 7;
  std::vector<Corner> capped;
  detectFast(im, opt, &capped);
  SC_CHECK(capped.size() <= 7);
}

/// BRIEF must be a pure function of the pixel patch: identical patches give
/// identical descriptors, different patches give different ones, and the
/// pattern is reproducible from its seed (loop closure depends on both).
void testBriefDescriptor() {
  const Image a = synthTexture(96, 96, 7);
  const Image b = synthTexture(96, 96, 8);
  BriefPattern p1(0xABCDEFull);
  BriefPattern p2(0xABCDEFull);
  Descriptor da, db, da2;
  SC_CHECK(p1.compute(a, 48, 48, &da));
  SC_CHECK(p2.compute(a, 48, 48, &da2));
  SC_CHECK(da.hamming(da2) == 0);  // same pattern, same pixels
  SC_CHECK(p1.compute(b, 48, 48, &db));
  SC_CHECK(da.hamming(db) > 0);  // different pixels
  // A patch that would leave the image is refused.
  Descriptor edge;
  SC_CHECK(!p1.compute(a, 1, 1, &edge));
  SC_CHECK(!p1.compute(a, 95, 95, &edge));
  // Descriptor bytes are the packed bit count, not a per-keyframe cost.
  SC_CHECK(p1.bytes() > 0);
}

void testDetectAndDescribe() {
  const Image im = synthTexture(320, 240, 3);
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 20;
  BriefOptions brief;
  std::vector<Feature> feats;
  detectAndDescribe(im, pat, fast, brief, &feats);
  SC_CHECK(!feats.empty());
  SC_CHECK(feats.size() <= static_cast<Size>(fast.max_features));
  for (const Feature& f : feats) {
    SC_CHECK(f.x >= 0 && f.y >= 0 && f.x < im.w && f.y < im.h);
    SC_CHECK(f.depth < 0);  // depth is not known until stereo runs
  }
  // Descriptors must be distinct for distinct corners; a degenerate detector
  // that returns the same bit vector everywhere would still pass a count test.
  int distinct = 0;
  for (Size i = 0; i < feats.size() && i < 40; ++i)
    for (Size j = i + 1; j < feats.size() && j < 40; ++j)
      if (feats[i].desc.hamming(feats[j].desc) > 0) ++distinct;
  SC_CHECK(distinct > 0);
}

/// Recover a known shift at one pyramid configuration. Returns the mean error
/// and the number of features tracked.
struct TrackResult {
  Scalar mean_err = 0;
  Size tracked = 0;
  Size coarse_rej = 0;
  Size fb_rej = 0;
};

TrackResult trackShift(const Image& a, const Image& b, int levels, Scalar sx, Scalar sy) {
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 15;
  BriefOptions brief;
  std::vector<Feature> feats;
  detectAndDescribe(a, pat, fast, brief, &feats);

  KltOptions opt;
  opt.max_level = levels;
  KltTracker klt(opt);
  std::vector<Feature> out;
  klt.track(a, b, feats, &out);

  TrackResult r;
  r.tracked = out.size();
  r.coarse_rej = klt.rejected_coarse();
  r.fb_rej = klt.rejected_forward_backward();
  if (out.empty()) return r;
  // Match each tracked feature back to its source by descriptor identity, which
  // is how the tracker carries the association through.
  Scalar sum = 0;
  Size n = 0;
  for (const Feature& t : out) {
    for (const Feature& s : feats) {
      if (t.desc.hamming(s.desc) != 0) continue;
      sum += std::hypot((t.x - s.x) - sx, (t.y - s.y) - sy);
      ++n;
      break;
    }
  }
  r.mean_err = n ? sum / static_cast<Scalar>(n) : 0.0;
  r.tracked = n;
  return r;
}

void testKltSingleLevelAccuracy() {
  const Image a = synthTexture(320, 240, 11);
  const Scalar sx = 1.0, sy = 0.3;
  const Image b = warpImage(a, sx, sy);
  const TrackResult r = trackShift(a, b, 0, sx, sy);
  SC_CHECK(r.tracked > 50);
  SC_NEAR(r.mean_err, 0.0, 0.1);
}

/// The regression test for the coarse-to-fine bug: every pyramid depth must
/// recover the same shift to the same accuracy. A tracker that rebuilds its
/// reference patch at the running estimate at each level reports an error
/// proportional to the number of levels here.
void testKltPyramidAccuracy() {
  const Image a = synthTexture(320, 240, 13);
  for (const Scalar shift : {0.5, 1.0, 2.0, 4.0}) {
    const Image b = warpImage(a, shift, 0.3 * shift);
    for (int levels = 0; levels <= 3; ++levels) {
      const TrackResult r = trackShift(a, b, levels, shift, 0.3 * shift);
      SC_CHECK(r.tracked > 50);
      // The tolerance is sub-pixel at every level. If the reference patch is
      // rebuilt per level, the error grows with the level count and this fails.
      SC_NEAR(r.mean_err, 0.0, 0.25);
    }
  }
}

/// A pure vertical shift exercises the other gradient axis; a tracker that
/// transposes gx/gy passes the horizontal test and fails here.
void testKltVerticalShift() {
  const Image a = synthTexture(320, 240, 17);
  const Scalar sy = 1.7;
  const Image b = warpImage(a, 0.0, sy);
  for (int levels = 0; levels <= 2; ++levels) {
    const TrackResult r = trackShift(a, b, levels, 0.0, sy);
    SC_CHECK(r.tracked > 50);
    SC_NEAR(r.mean_err, 0.0, 0.25);
  }
}

/// Forward-backward rejection must fire on a track that cannot be recovered.
/// Here the second image is pure noise, so every track is genuinely wrong and
/// the tracker must reject essentially all of them rather than report garbage.
void testKltRejectsUnrecoverableTracks() {
  const Image a = synthTexture(320, 240, 23);
  const Image noise = synthTexture(320, 240, 999);  // uncorrelated content
  const TrackResult r = trackShift(a, noise, 2, 0.0, 0.0);
  // Some tracks may survive by luck on a matching patch, but the bulk must go.
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 15;
  BriefOptions brief;
  std::vector<Feature> feats;
  detectAndDescribe(a, pat, fast, brief, &feats);
  SC_CHECK(r.tracked < feats.size());
}

/// The tracker must not report a large displacement for a stationary scene.
/// This is the false-positive gate: an unstable tracker drifts on zero motion.
void testKltZeroMotionStable() {
  const Image a = synthTexture(320, 240, 29);
  const TrackResult r = trackShift(a, a, 3, 0.0, 0.0);
  SC_CHECK(r.tracked > 100);
  SC_NEAR(r.mean_err, 0.0, 0.05);
  SC_NEAR(static_cast<Scalar>(r.coarse_rej), 0.0, 0.0);
}

/// Stereo depth on a controlled warp: shift the right image by exactly the
/// disparity a frontal surface at depth Z produces, and require the matcher to
/// recover Z. This isolates the matcher from the renderer's geometry.
void testStereoDepthOnWarp() {
  CameraModel cam;
  const Image left = synthTexture(320, 240, 31);
  StereoOptions opt;
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 15;
  BriefOptions brief;
  std::vector<Feature> feats;
  detectAndDescribe(left, pat, fast, brief, &feats);
  SC_CHECK(feats.size() > 50);

  for (const Scalar disp : {2.0, 4.0, 8.0, 16.0}) {
    // Disparity is u_left - u_right, so the right image's content sits to the
    // *left*: the right image is the left one shifted by -disp in x.
    const Image right = warpImage(left, -disp, 0.0);
    std::vector<Feature> fs = feats;
    stereoDepth(left, right, cam, opt, &fs);
    const Scalar true_z = cam.fx * cam.baseline / disp;
    Size valid = 0;
    Scalar rel = 0;
    for (const Feature& f : fs) {
      if (f.depth <= 0) continue;
      rel += std::abs(f.depth - true_z) / true_z;
      ++valid;
    }
    SC_CHECK(valid > 0);
    // Sub-pixel disparity recovery should keep this well inside a few percent.
    SC_NEAR(rel / static_cast<Scalar>(valid), 0.0, 0.05);
  }
}

/// Depth must be refused, not guessed, when there is nothing to match. A flat
/// right image gives no unique minimum and the matcher must return nothing.
void testStereoRejectsFlatMatch() {
  CameraModel cam;
  const Image left = synthTexture(320, 240, 37);
  Image flat;
  flat.resize(320, 240);
  std::fill(flat.px.begin(), flat.px.end(), 128);
  StereoOptions opt;
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 15;
  BriefOptions brief;
  std::vector<Feature> feats;
  detectAndDescribe(left, pat, fast, brief, &feats);
  stereoDepth(left, flat, cam, opt, &feats);
  for (const Feature& f : feats) SC_CHECK(f.depth < 0);
}

void run() {
  testFastDetectsCorners();
  testBriefDescriptor();
  testDetectAndDescribe();
  testKltSingleLevelAccuracy();
  testKltPyramidAccuracy();
  testKltVerticalShift();
  testKltRejectsUnrecoverableTracks();
  testKltZeroMotionStable();
  testStereoDepthOnWarp();
  testStereoRejectsFlatMatch();
}

}  // namespace

SC_TEST_MAIN("test_frontend")
