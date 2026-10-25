// LUNA -- BRIEF descriptors.
//
// BRIEF is implemented with a *fixed, seeded* sampling pattern so descriptors
// are reproducible across machines and runs. Descriptors are 256 bits packed
// into four 64-bit words, which is what makes the descriptor DB cheap: one
// XOR-popcount per candidate pair, and an LRU index that costs 32 bytes per
// keyframe.
#include "luna/frontend/features.hpp"

#include <cstdlib>

namespace luna {
namespace frontend {

// ---------------------------------------------------------------------------
// BRIEF
// ---------------------------------------------------------------------------

BriefPattern::BriefPattern(std::uint64_t seed, int patch, int pairs) : patch_(patch), pairs_(pairs) {
  Rng rng(seed);
  offsets_.resize(static_cast<Size>(pairs));
  const int h = patch / 2;
  // Isotropic Gaussian sampling around the patch centre, clipped to the patch.
  const Scalar sigma = 0.4 * patch;
  for (int i = 0; i < pairs; ++i) {
    for (int k = 0; k < 4; ++k) {
      Scalar v;
      do {
        v = rng.gaussian(sigma);
      } while (std::abs(v) > h);
      offsets_[static_cast<Size>(i)][k] = static_cast<std::int8_t>(std::lround(v));
    }
  }
}

bool BriefPattern::compute(const sim::Image& img, int x, int y, Descriptor* out) const {
  const int h = patch_ / 2;
  if (x - h < 0 || y - h < 0 || x + h >= img.w || y + h >= img.h) return false;
  out->w = {{0, 0, 0, 0}};
  for (int i = 0; i < pairs_; ++i) {
    const auto& o = offsets_[static_cast<Size>(i)];
    const int a = img.at(x + o[0], y + o[1]);
    const int b = img.at(x + o[2], y + o[3]);
    if (a < b) out->w[static_cast<Size>(i) >> 6] |= (1ull << (i & 63));
  }
  return true;
}

void detectAndDescribe(const sim::Image& img, const BriefPattern& pattern, const FastOptions& fast,
                       const BriefOptions& brief, std::vector<Feature>* out) {
  std::vector<Corner> corners;
  FastOptions fo = fast;
  fo.max_features = brief.max_features;
  fo.border = std::max(fo.border, pattern.halfPatch() + 1);
  detectFast(img, fo, &corners);

  out->clear();
  out->reserve(corners.size());
  const int h = pattern.halfPatch();
  for (const Corner& c : corners) {
    // Contrast gate: a corner whose local intensity range is flat carries no
    // descriptor information and only adds tracking noise.
    int mn = 255, mx = 0;
    for (int dy = -h; dy <= h; dy += 3)
      for (int dx = -h; dx <= h; dx += 3) {
        const int v = img.at(c.x + dx, c.y + dy);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      }
    if (mx - mn < brief.threshold) continue;

    Feature f;
    f.x = c.x;
    f.y = c.y;
    f.response = c.score;
    if (!pattern.compute(img, c.x, c.y, &f.desc)) continue;
    out->push_back(f);
  }
}
}  // namespace frontend
}  // namespace luna
