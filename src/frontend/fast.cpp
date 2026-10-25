#include "slamcore/frontend/fast.hpp"

#include <algorithm>
#include <array>

#if defined(__SSSE3__) || defined(__SSE2__)
#include <emmintrin.h>
#endif

namespace sc {
namespace frontend {

namespace {

/// Bresenham circle of radius 3, in the conventional FAST order starting at
/// the top pixel and going clockwise.
constexpr int kCircle[16][2] = {
    {0, -3},  {1, -3},  {2, -2},  {3, -1},  {3, 0},   {3, 1},   {2, 2},   {1, 3},
    {0, 3},   {-1, 3},  {-2, 2},  {-3, 1},  {-3, 0},  {-3, -1}, {-2, -2}, {-1, -3},
};

/// The 4 offsets that FAST first tests: 1, 5, 9, 13.
constexpr int kQuad[4] = {1, 5, 9, 13};

/// Returns true when at least 9 contiguous circle pixels are all brighter than
/// (p + t) or all darker than (p - t).
inline bool isCorner(const std::uint8_t* base, int stride, int p, int t, int* score) {
  const int hi = p + t;
  const int lo = p - t;

  // High-speed test on the 4 cardinal points. Any 9-contiguous arc must
  // include at least 3 of these 4.
  int brighter = 0, darker = 0;
  for (int i = 0; i < 4; ++i) {
    const int q = base[kCircle[kQuad[i]][1] * stride + kCircle[kQuad[i]][0]];
    if (q > hi) ++brighter;
    if (q < lo) ++darker;
  }
  if (brighter < 3 && darker < 3) return false;

  // Full test: count the longest contiguous run in both directions.
  std::array<int, 16> vals;
  for (int i = 0; i < 16; ++i) vals[i] = base[kCircle[i][1] * stride + kCircle[i][0]];

  int best_hi = 0, best_lo = 0;
  for (int start = 0; start < 16; ++start) {
    int run_hi = 0, run_lo = 0;
    for (int k = 0; k < 16; ++k) {
      const int v = vals[(start + k) & 15];
      if (v > hi) ++run_hi; else break;
    }
    for (int k = 0; k < 16; ++k) {
      const int v = vals[(start + k) & 15];
      if (v < lo) ++run_lo; else break;
    }
    best_hi = std::max(best_hi, run_hi);
    best_lo = std::max(best_lo, run_lo);
  }
  if (best_hi < 9 && best_lo < 9) return false;

  // Score: sum of absolute differences over the circle, saturated.
  int s = 0;
  for (int i = 0; i < 16; ++i) s += std::abs(vals[i] - p);
  *score = s;
  return true;
}

}  // namespace

void nonMaxSuppress(const std::vector<int>& scores, int w, int h, std::vector<Corner>* out,
                    int max_features) {
  out->clear();
  // Collect candidates above zero score, then keep only 3x3 maxima.
  struct Cand {
    int x, y, s;
  };
  std::vector<Cand> cands;
  cands.reserve(1024);
  for (int y = 1; y < h - 1; ++y) {
    for (int x = 1; x < w - 1; ++x) {
      const int s = scores[static_cast<Size>(y) * w + x];
      if (s <= 0) continue;
      bool is_max = true;
      for (int dy = -1; dy <= 1 && is_max; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0) continue;
          if (scores[static_cast<Size>(y + dy) * w + (x + dx)] > s) {
            is_max = false;
            break;
          }
        }
      if (is_max) cands.push_back({x, y, s});
    }
  }
  if (static_cast<int>(cands.size()) > max_features) {
    std::partial_sort(cands.begin(), cands.begin() + max_features, cands.end(),
                      [](const Cand& a, const Cand& b) { return a.s > b.s; });
    cands.resize(static_cast<Size>(max_features));
  }
  out->reserve(cands.size());
  for (const Cand& c : cands) out->push_back(Corner{c.x, c.y, c.s});
}

void detectFast(const sim::Image& img, const FastOptions& opt, std::vector<Corner>* out) {
  const int w = img.w, h = img.h;
  const int b = std::max(opt.border, 3);
  std::vector<int> scores(static_cast<Size>(w) * h, 0);
  const int stride = w;

  for (int y = b; y < h - b; ++y) {
    const std::uint8_t* row = img.px.data() + static_cast<Size>(y) * stride;
    for (int x = b; x < w - b; ++x) {
      const int p = row[x];
      int score = 0;
      if (isCorner(row + x, stride, p, opt.threshold, &score)) {
        scores[static_cast<Size>(y) * w + x] = score;
      }
    }
  }
  if (opt.nonmax) {
    nonMaxSuppress(scores, w, h, out, opt.max_features);
  } else {
    out->clear();
    for (int y = b; y < h - b; ++y)
      for (int x = b; x < w - b; ++x) {
        const int s = scores[static_cast<Size>(y) * w + x];
        if (s > 0) out->push_back(Corner{x, y, s});
      }
    if (static_cast<int>(out->size()) > opt.max_features) out->resize(opt.max_features);
  }
}

}  // namespace frontend
}  // namespace sc
