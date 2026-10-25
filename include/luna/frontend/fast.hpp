// LUNA -- FAST corner detector with SSE acceleration.
//
// The 16-pixel Bresenham circle test is branch-free up to the first three
// rejected offsets (the standard high-speed bail-out), and the surviving
// candidates are vectorised over four circle offsets at a time with SSE2/SSSE3
// integer compares. This keeps corner extraction well under 1 ms per 320x240
// frame on the reference CPU budget.
#pragma once

#include <vector>

#include "luna/core/types.hpp"
#include "luna/sim/sensors.hpp"

namespace luna {
namespace frontend {

struct FastOptions {
  int threshold = 25;      ///< intensity difference threshold
  bool nonmax = true;      ///< apply 3x3 non-maximum suppression
  int border = 4;          ///< pixels to skip at the image border
  int max_features = 400;  ///< hard cap (RAM discipline)
};

struct Corner {
  int x = 0, y = 0;
  int score = 0;           ///< FAST score (sum of absolute differences)
};

/// Detect corners in a grayscale image.
void detectFast(const sim::Image& img, const FastOptions& opt, std::vector<Corner>* out);

/// Non-maximum suppression over the score image; exposed for tests.
void nonMaxSuppress(const std::vector<int>& scores, int w, int h, std::vector<Corner>* out, int max_features);

}  // namespace frontend
}  // namespace luna
