#include "luna/backend/loop_close.hpp"

#include <algorithm>
#include <cmath>

#include <Eigen/SVD>

#include "luna/core/stats.hpp"

namespace luna {
namespace backend {

namespace {

/// Build a rotation from the correlation matrix by SVD, with the reflection
/// guard.
///
/// The naive R = V U^T can come out with determinant -1, which is a reflection
/// and not a rotation. On noisy data this happens often enough to matter, and
/// a reflected "alignment" of two point sets can still fit many inliers, so it
/// would be accepted as a loop closure and then mirror the map. Flipping the
/// smallest singular vector's sign is the standard fix.
Mat3 rotationFromCorrelation(const Mat3& H) {
  Eigen::JacobiSVD<Mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Mat3 U = svd.matrixU();
  Mat3 V = svd.matrixV();
  Mat3 R = V * U.transpose();
  if (R.determinant() < 0) {
    U.col(2) *= -1.0;
    R = V * U.transpose();
  }
  return R;
}

}  // namespace

bool alignPointSets(const std::vector<Vec3>& src, const std::vector<Vec3>& dst, Mat3* R, Vec3* t) {
  if (src.size() != dst.size() || src.size() < 3) return false;
  const Size n = src.size();
  Vec3 cs = Vec3::Zero(), cd = Vec3::Zero();
  for (Size i = 0; i < n; ++i) {
    cs += src[i];
    cd += dst[i];
  }
  cs /= static_cast<Scalar>(n);
  cd /= static_cast<Scalar>(n);

  Mat3 H = Mat3::Zero();
  for (Size i = 0; i < n; ++i) H += (src[i] - cs) * (dst[i] - cd).transpose();

  *R = rotationFromCorrelation(H);
  *t = cd - (*R) * cs;
  return true;
}

bool verifyRansac(const std::vector<Vec3>& a, const std::vector<Vec3>& b, const LoopCloseOptions& opt,
                  Mat4* T, std::vector<int>* inliers, Size* attempts) {
  if (a.size() != b.size() || a.size() < 3) return false;
  const Size n = a.size();
  if (attempts) *attempts = 0;

  Rng rng(0xC0FFEEull ^ static_cast<std::uint64_t>(n));
  const Scalar thresh2 = opt.inlier_dist * opt.inlier_dist;

  std::vector<int> best;
  Mat3 best_R = Mat3::Identity();
  Vec3 best_t = Vec3::Zero();

  std::vector<Vec3> sa(3), sb(3);
  for (int iter = 0; iter < opt.ransac_iterations; ++iter) {
    // Sample three distinct indices.
    int i0 = static_cast<int>(rng.uniform() * n) % n;
    int i1 = static_cast<int>(rng.uniform() * n) % n;
    int i2 = static_cast<int>(rng.uniform() * n) % n;
    if (i0 == i1 || i1 == i2 || i0 == i2) continue;
    sa[0] = a[i0]; sa[1] = a[i1]; sa[2] = a[i2];
    sb[0] = b[i0]; sb[1] = b[i1]; sb[2] = b[i2];

    Mat3 R;
    Vec3 t;
    if (!alignPointSets(sa, sb, &R, &t)) continue;

    std::vector<int> inl;
    for (Size k = 0; k < n; ++k) {
      const Vec3 pred = R * a[k] + t;
      if ((pred - b[k]).squaredNorm() <= thresh2) inl.push_back(static_cast<int>(k));
    }
    if (inl.size() > best.size()) {
      best = inl;
      best_R = R;
      best_t = t;
    }
    if (attempts) ++(*attempts);
  }

  if (best.size() < 3) return false;

  // Refit on the inlier set, then re-count. The three-point model is noisy;
  // refitting on all inliers is what makes the final inlier count meaningful
  // rather than a product of a lucky sample.
  std::vector<Vec3> ra, rb;
  ra.reserve(best.size());
  rb.reserve(best.size());
  for (int idx : best) {
    ra.push_back(a[idx]);
    rb.push_back(b[idx]);
  }
  Mat3 R;
  Vec3 t;
  if (!alignPointSets(ra, rb, &R, &t)) return false;

  std::vector<int> final_inl;
  for (Size k = 0; k < n; ++k) {
    const Vec3 pred = R * a[k] + t;
    if ((pred - b[k]).squaredNorm() <= thresh2) final_inl.push_back(static_cast<int>(k));
  }

  *T = Mat4::Identity();
  T->topLeftCorner<3, 3>() = R;
  T->topRightCorner<3, 1>() = t;
  if (inliers) *inliers = final_inl;
  return true;
}

LoopCloser::LoopCloser(const LoopCloseOptions& opt) : opt_(opt) {}

void LoopCloser::clear() {
  db_.clear();
  stats_ = LoopCloseStats();
}

void LoopCloser::addKeyframe(const Keyframe& kf) {
  db_.push_back(kf);
  ++stats_.keyframes;
  // Bounded database: evict the oldest once the cap is exceeded. A long
  // sequence would otherwise grow this without limit.
  while (db_.size() > opt_.max_keyframes) {
    db_.erase(db_.begin());
    ++stats_.db_evictions;
  }
}

Size LoopCloser::matchDescriptors(const std::vector<MapPoint>& a, const std::vector<MapPoint>& b,
                                  std::vector<std::pair<int, int>>* matches) const {
  if (matches) matches->clear();
  if (a.empty() || b.empty()) return 0;

  Size n_match = 0;
  for (Size i = 0; i < a.size(); ++i) {
    int best_d = opt_.max_hamming + 1, second_d = opt_.max_hamming + 1;
    int best_j = -1;
    for (Size j = 0; j < b.size(); ++j) {
      const int d = a[i].desc.hamming(b[j].desc);
      if (d < best_d) {
        second_d = best_d;
        best_d = d;
        best_j = static_cast<int>(j);
      } else if (d < second_d) {
        second_d = d;
      }
    }
    // Lowe's ratio test. The guard on second_d keeps a lone descriptor in the
    // database from being rejected simply because there is no runner-up.
    if (best_j < 0 || best_d > opt_.max_hamming) continue;
    if (second_d <= opt_.max_hamming && best_d > opt_.ratio * second_d) continue;
    if (matches) matches->push_back({static_cast<int>(i), best_j});
    ++n_match;
  }
  return n_match;
}

bool LoopCloser::query(const Keyframe& query, LoopClosure* out) {
  if (db_.empty() || query.points.size() < 6) return false;

  for (const Keyframe& kf : db_) {
    if (std::abs(kf.id - query.id) < opt_.min_index_gap) continue;
    if (kf.points.size() < 6) continue;
    ++stats_.candidates_tested;

    std::vector<std::pair<int, int>> matches;
    const Size n = matchDescriptors(query.points, kf.points, &matches);
    if (n < static_cast<Size>(opt_.min_inliers)) {
      ++stats_.rejected_inliers;
      continue;
    }

    std::vector<Vec3> a, b;
    a.reserve(matches.size());
    b.reserve(matches.size());
    for (const auto& m : matches) {
      a.push_back(query.points[m.first].p);
      b.push_back(kf.points[m.second].p);
    }

    ++stats_.verification_attempts;
    Mat4 T;
    std::vector<int> inliers;
    Size attempts = 0;
    if (!verifyRansac(a, b, opt_, &T, &inliers, &attempts)) {
      ++stats_.rejected_inliers;
      continue;
    }
    const Scalar ratio =
        a.empty() ? 0.0 : static_cast<Scalar>(inliers.size()) / static_cast<Scalar>(a.size());
    if (static_cast<int>(inliers.size()) < opt_.min_inliers || ratio < opt_.min_inlier_ratio) {
      ++stats_.rejected_ratio;
      continue;
    }

    if (out) {
      out->from = kf.id;
      out->to = query.id;
      // The measurement is the transform RANSAC fitted to the matched 3D
      // points, NOT the relative pose of the current estimates. The two are
      // very different things: the estimates already contain the drift the
      // closure exists to correct, so using them makes the constraint's
      // residual zero and the loop edge adds no information at all -- it looks
      // like a working loop closure while doing nothing, and every accepted
      // closure scores as a false positive because the graph never moves.
      //
      // `verifyRansac` aligns the query's points (a) onto the stored
      // keyframe's points (b), so the fitted transform maps the query frame
      // into the stored frame, which is exactly x_from^{-1} x_to.
      const Mat3 R = T.topLeftCorner<3, 3>();
      const Vec3 t = T.topRightCorner<3, 1>();
      // RANSAC fits a rigid transform T = [R | t] acting on raw 3D points, but
      // a graph edge's `z` is an SE(2) *tangent* vector: se2Exp(z) is what the
      // residual compares against. The two differ by the left Jacobian of
      // SE(2), which scales the translation by (sin th / th) and rotates it by
      // half the angle. Storing the raw translation instead leaves an error
      // proportional to the drift being corrected, and it grows with the
      // relative rotation -- the opposite of what a closure is for.
      Mat3 Tse2 = Mat3::Identity();
      Tse2.topLeftCorner<2, 2>() = R.topLeftCorner<2, 2>();
      Tse2.topRightCorner<2, 1>() = t.head<2>();
      out->z = se2Log(Tse2);
      out->inliers = static_cast<int>(inliers.size());
      out->correspondences = static_cast<int>(a.size());
      out->inlier_ratio = ratio;

      // The measurement noise is the spread of the inlier residuals, not the
      // RANSAC inlier threshold.
      //
      // The threshold only says the points agree to within `inlier_dist`; it
      // says nothing about how accurate the fitted transform is. Deriving sigma
      // from it makes the closure look an order of magnitude more precise than
      // the 3D reconstruction it came from, so the graph over-trusts it and
      // *worsens* a trajectory whose odometry was already good. Measuring the
      // residuals directly makes the weight self-tuning: a clean fit earns a
      // stiff constraint, a marginal one earns a soft one.
      const Mat3 Ri = T.topLeftCorner<3, 3>();
      const Vec3 ti = T.topRightCorner<3, 1>();
      Scalar sq = 0;
      for (int k : inliers) {
        sq += (Ri * a[static_cast<Size>(k)] + ti - b[static_cast<Size>(k)]).squaredNorm();
      }
      const Scalar rms = inliers.empty()
                             ? opt_.inlier_dist
                             : std::sqrt(sq / static_cast<Scalar>(inliers.size()));
      // The residual RMS, used as-is. It would be tempting to divide by
      // sqrt(n) on the grounds that the pose is fitted from n points, but the
      // residuals here are dominated by a *systematic* error -- the depth bias
      // of the stereo-lite reconstruction, shared by every point in the
      // keyframe -- so averaging over more points does not reduce it. Measured
      // against ground truth the fitted transform is wrong by about the RMS,
      // not by RMS/sqrt(n), and claiming the smaller figure makes the graph
      // over-trust the closure and drag a good trajectory off course.
      const Scalar sigma = std::max<Scalar>(opt_.min_sigma, rms);
      out->info = Mat3::Identity() / (sigma * sigma);
    }
    ++stats_.accepted;
    return true;
  }
  return false;
}

}  // namespace backend
}  // namespace luna
