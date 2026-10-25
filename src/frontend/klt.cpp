// LUNA -- pyramid Lucas-Kanade optical flow with forward-backward
// outlier rejection.
//
// The tracker works coarse-to-fine over a Gaussian pyramid. Two details carry
// the accuracy: the reference patch is built once, at the *previous* image's
// position, and each level then solves only for the residual displacement left
// by the coarser level; and the pyramid prefilter must be wide enough to band
// limit the decimated image, or the coarse levels match detail the fine level
// cannot reproduce and the estimate drifts instead of converging.
#include "luna/frontend/features.hpp"

#include <algorithm>

namespace luna {
namespace frontend {

// ---------------------------------------------------------------------------
// KLT
// ---------------------------------------------------------------------------

KltTracker::KltTracker(const KltOptions& opt) : opt_(opt) {}

void KltTracker::downsample(const sim::Image& src, sim::Image* dst) {
  const int w = std::max(1, src.w / 2), h = std::max(1, src.h / 2);
  dst->resize(w, h);
  // Separable Gaussian (sigma ~ 1.2 px) followed by 2x decimation. The
  // prefilter width is the load-bearing choice for the pyramid: it must cut
  // everything above the decimated Nyquist limit, or the coarse levels keep
  // high-frequency detail that has no counterpart in the fine levels. Then
  // coarse-to-fine descent converges against a pattern the fine level cannot
  // reproduce, and the estimate drifts by a pixel per level instead of
  // shrinking toward the true displacement. A 3-tap binomial is too weak here.
  static const Scalar k[5] = {0.0706, 0.2448, 0.3691, 0.2448, 0.0706};
  const auto cl = [&](int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); };
  for (int y = 0; y < h; ++y) {
    const int sy = 2 * y;
    for (int x = 0; x < w; ++x) {
      const int sx = 2 * x;
      Scalar acc = 0;
      for (int dy = -2; dy <= 2; ++dy) {
        const int py = cl(sy + dy, 0, src.h - 1);
        Scalar row = 0;
        for (int dx = -2; dx <= 2; ++dx) {
          row += k[dx + 2] * src.at(cl(sx + dx, 0, src.w - 1), py);
        }
        acc += k[dy + 2] * row;
      }
      dst->at(x, y) = static_cast<std::uint8_t>(clamp(acc, 0.0, 255.0));
    }
  }
}

void KltTracker::buildPyramid(const sim::Image& img) {
  prev_pyr_.resize(static_cast<Size>(opt_.max_level + 1));
  prev_pyr_[0] = img;
  for (int l = 1; l <= opt_.max_level; ++l)
    downsample(prev_pyr_[static_cast<Size>(l - 1)], &prev_pyr_[static_cast<Size>(l)]);
}

void KltTracker::buildCurrentPyramid(const sim::Image& img) {
  cur_pyr_.resize(static_cast<Size>(opt_.max_level + 1));
  cur_pyr_[0] = img;
  for (int l = 1; l <= opt_.max_level; ++l)
    downsample(cur_pyr_[static_cast<Size>(l - 1)], &cur_pyr_[static_cast<Size>(l)]);
}

Scalar KltTracker::samplePrev(int level, Scalar x, Scalar y) const {
  const sim::Image& img = prev_pyr_[static_cast<Size>(level)];
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const Scalar fx = x - x0, fy = y - y0;
  const auto s = [&](int px, int py) -> Scalar {
    return static_cast<Scalar>(img.at(clamp(px, 0, img.w - 1), clamp(py, 0, img.h - 1)));
  };
  const Scalar a = s(x0, y0) * (1 - fx) + s(x0 + 1, y0) * fx;
  const Scalar b = s(x0, y0 + 1) * (1 - fx) + s(x0 + 1, y0 + 1) * fx;
  return a * (1 - fy) + b * fy;
}

Scalar KltTracker::sampleLevel(int level, Scalar x, Scalar y) const {
  const sim::Image& img = cur_pyr_[static_cast<Size>(level)];
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const Scalar fx = x - x0, fy = y - y0;
  const auto s = [&](int px, int py) -> Scalar {
    return static_cast<Scalar>(img.at(clamp(px, 0, img.w - 1), clamp(py, 0, img.h - 1)));
  };
  const Scalar a = s(x0, y0) * (1 - fx) + s(x0 + 1, y0) * fx;
  const Scalar b = s(x0, y0 + 1) * (1 - fx) + s(x0 + 1, y0 + 1) * fx;
  return a * (1 - fy) + b * fy;
}

bool KltTracker::trackPoint(int level, const Vec2& x_ref, const Vec2& x_cur, Vec2* out) {
  const sim::Image& prev = prev_pyr_[static_cast<Size>(level)];
  const sim::Image& cur = cur_pyr_[static_cast<Size>(level)];
  // Window size is in pixels *at this level*, so a coarse level sees the same
  // physical support as the fine level. A constant pixel window would make
  // coarse levels cover 2^L times more of the scene and average the signal away.
  const int half = opt_.window / 2;
  const int rx = static_cast<int>(std::lround(x_ref.x()));
  const int ry = static_cast<int>(std::lround(x_ref.y()));
  if (rx - half - 2 < 0 || ry - half - 2 < 0 || rx + half + 2 >= prev.w ||
      ry + half + 2 >= prev.h)
    return false;
  if (x_cur.x() - half - 2 < 0 || x_cur.y() - half - 2 < 0 ||
      x_cur.x() + half + 2 >= cur.w || x_cur.y() + half + 2 >= cur.h)
    return false;

  // Reference patch and its spatial gradients, sampled bilinearly at the exact
  // sub-pixel reference position. Sampling at a rounded integer position would
  // bake that rounding offset into the answer.
  const int n = opt_.window * opt_.window;
  ref_.resize(static_cast<Size>(n));
  gx_.resize(static_cast<Size>(n));
  gy_.resize(static_cast<Size>(n));
  Scalar gxx = 0, gxy = 0, gyy = 0;
  for (int dy = -half; dy <= half; ++dy) {
    for (int dx = -half; dx <= half; ++dx) {
      const int i = (dy + half) * opt_.window + (dx + half);
      const Scalar cx = x_ref.x() + dx, cy = x_ref.y() + dy;
      ref_[static_cast<Size>(i)] = samplePrev(level, cx, cy);
      const Scalar gxv = 0.5 * (samplePrev(level, cx + 1, cy) - samplePrev(level, cx - 1, cy));
      const Scalar gyv = 0.5 * (samplePrev(level, cx, cy + 1) - samplePrev(level, cx, cy - 1));
      gx_[static_cast<Size>(i)] = gxv;
      gy_[static_cast<Size>(i)] = gyv;
      gxx += gxv * gxv;
      gxy += gxv * gyv;
      gyy += gyv * gyv;
    }
  }
  // Shrinkage keeps the normal matrix well conditioned on edges, so the
  // eigenvalue gate is a real texture test rather than a numerical accident.
  gxx += 1e-3 * n;
  gyy += 1e-3 * n;
  const Scalar det = gxx * gyy - gxy * gxy;
  if (det <= 0) return false;
  const Scalar min_eig = 0.5 * (gxx + gyy) - std::sqrt(0.25 * sq(gxx - gyy) + sq(gxy));
  if (min_eig < opt_.min_eigenvalue) return false;
  const Scalar inv_det = 1.0 / det;

  Vec2 p = x_cur;
  for (int it = 0; it < opt_.max_iter; ++it) {
    Scalar bx = 0, by = 0;
    for (int dy = -half; dy <= half; ++dy)
      for (int dx = -half; dx <= half; ++dx) {
        const int i = (dy + half) * opt_.window + (dx + half);
        const Scalar d = sampleLevel(level, p.x() + dx, p.y() + dy) - ref_[static_cast<Size>(i)];
        bx += gx_[static_cast<Size>(i)] * d;
        by += gy_[static_cast<Size>(i)] * d;
      }
    // Gauss-Newton step on the photometric residual: solve J'J delta = J'r.
    // The step is clamped per level: unclamped, a coarse level can jump a whole
    // texture period and land in the neighbouring minimum, a silent large
    // error rather than a visible failure.
    Scalar ddx = -(gyy * bx - gxy * by) * inv_det;
    Scalar ddy = -(gxx * by - gxy * bx) * inv_det;
    const Scalar step = std::hypot(ddx, ddy);
    if (step > opt_.max_step_per_level) {
      const Scalar sc = opt_.max_step_per_level / step;
      ddx *= sc;
      ddy *= sc;
    }
    p.x() += ddx;
    p.y() += ddy;
    if (std::abs(ddx) + std::abs(ddy) < opt_.eps) break;
  }

  Scalar sad = 0;
  for (int dy = -half; dy <= half; ++dy)
    for (int dx = -half; dx <= half; ++dx) {
      const int i = (dy + half) * opt_.window + (dx + half);
      sad += std::abs(sampleLevel(level, p.x() + dx, p.y() + dy) - ref_[static_cast<Size>(i)]);
    }
  if (sad / n > opt_.max_residual) return false;
  *out = p;
  return true;
}

void KltTracker::track(const sim::Image& prev, const sim::Image& cur, const std::vector<Feature>& in,
                       std::vector<Feature>* out) {
  tracked_ = 0;
  rejected_ = 0;
  rejected_coarse_ = 0;
  rejected_fb_ = 0;
  out->clear();
  if (in.empty()) return;

  buildPyramid(prev);
  buildCurrentPyramid(cur);

  out->reserve(in.size());
  for (const Feature& f : in) {
    if (static_cast<int>(out->size()) >= opt_.max_track) break;

    // Coarse-to-fine. `x_ref` is the fixed reference position in the previous
    // image; only the estimate in the current image is refined. Each level
    // therefore solves for the residual left by the level above it.
    const Vec2 ref_full(f.x, f.y);
    Vec2 p = ref_full;
    bool ok = true;
    for (int l = opt_.max_level; l >= 0; --l) {
      const Scalar scale = static_cast<Scalar>(1 << l);
      Vec2 q;
      if (!trackPoint(l, Vec2(ref_full.x() / scale, ref_full.y() / scale),
                      Vec2(p.x() / scale, p.y() / scale), &q)) {
        // A coarse-level failure is recoverable; the finest level is not.
        if (l == 0) {
          ok = false;
          break;
        }
        continue;
      }
      p = Vec2(q.x() * scale, q.y() * scale);
    }
    if (!ok) {
      ++rejected_;
      ++rejected_coarse_;
      continue;
    }

    // Forward-backward consistency: track the estimate back into the previous
    // image and require the round trip to return to the start. This is the gate
    // that makes occlusion and blackout faults survivable instead of corrupting
    // the pose.
    //
    // The images are swapped, so two positions are in play and they live in
    // different frames: `b_ref` stays fixed at `p` in the *current* image (that
    // is the patch being matched), while `b_est` walks through the *previous*
    // image and must end up at `f`. Rebuilding the reference at the running
    // estimate instead would compare the estimate against itself, and the check
    // would pass no matter how badly the track failed.
    std::swap(prev_pyr_, cur_pyr_);
    const Vec2 b_ref = p;
    Vec2 b_est = p;
    bool fb_ok = true;
    for (int l = opt_.max_level; l >= 0; --l) {
      const Scalar scale = static_cast<Scalar>(1 << l);
      Vec2 q;
      if (!trackPoint(l, Vec2(b_ref.x() / scale, b_ref.y() / scale),
                      Vec2(b_est.x() / scale, b_est.y() / scale), &q)) {
        if (l == 0) {
          fb_ok = false;
          break;
        }
        continue;
      }
      b_est = Vec2(q.x() * scale, q.y() * scale);
    }
    std::swap(prev_pyr_, cur_pyr_);
    if (!fb_ok || (b_est - ref_full).norm() > opt_.fb_threshold) {
      ++rejected_;
      ++rejected_fb_;
      continue;
    }

    Feature nf = f;
    nf.x = p.x();
    nf.y = p.y();
    out->push_back(nf);
    ++tracked_;
  }
}

}  // namespace frontend
}  // namespace luna
