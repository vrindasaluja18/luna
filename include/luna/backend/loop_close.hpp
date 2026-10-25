// LUNA -- loop closure: appearance retrieval plus geometric verification.
//
// Two stages, and both are load-bearing:
//
//   1. Retrieval. A bounded descriptor database returns candidate keyframes
//      that look like the current one. This is recall-oriented: it is allowed
//      to be wrong, because the next stage rejects.
//   2. Verification. RANSAC over 3D-3D correspondences fits a rigid transform
//      and requires a high inlier ratio before the constraint is accepted.
//      This is precision-oriented: a wrong loop closure corrupts every pose
//      between the two ends of the loop, which is far worse than missing one.
//
// The matcher deliberately does *not* use the simulator's ground-truth
// landmark ids. Those ids exist in the frontend features (`gt_id`/`has_gt`) and
// using them here would make the loop-closure precision figure meaningless:
// the pipeline would be matching on the answer key. Ground truth is used only
// afterwards, by the evaluator, to score what the descriptor matcher found.
#pragma once

#include <cstdint>
#include <vector>

#include "luna/backend/pose_graph.hpp"
#include "luna/core/se2.hpp"
#include "luna/frontend/features.hpp"

namespace luna {
namespace backend {

/// One observed 3D point in a keyframe's body frame.
struct MapPoint {
  Vec3 p = Vec3::Zero();
  frontend::Descriptor desc;
};

/// A keyframe: a pose plus the 3D structure observed from it.
struct Keyframe {
  Index id = -1;
  Pose2 pose;
  Scalar t = 0;
  std::vector<MapPoint> points;
};

struct LoopCloseOptions {
  /// Minimum Hamming distance improvement of the best match over the second
  /// best (Lowe's ratio test). Without this, a repetitive corridor matches
  /// every door to every other door.
  Scalar ratio = 0.8;
  /// Descriptors this far apart are not considered matches at all.
  int max_hamming = 80;
  /// RANSAC iterations for the verification stage.
  int ransac_iterations = 200;
  /// Inlier threshold, metres, for 3D-3D verification.
  Scalar inlier_dist = 0.15;
  /// Minimum inliers for a closure to be accepted.
  int min_inliers = 12;
  /// Minimum inlier fraction of the candidate correspondences.
  Scalar min_inlier_ratio = 0.5;
  /// Floor on the closure's position/rotation sigma, metres and radians. Keeps a
  /// perfect fit from claiming infinite precision, which would let one closure
  /// weld two poses together and override the odometry backbone. Raise it to
  /// distrust closures relative to odometry; lower it only if the measurements
  /// are genuinely better than the odometry they correct.
  Scalar min_sigma = 0.10;
  /// Keyframes closer than this in *index* are never loop candidates: a
  /// nearby keyframe is odometry, not a loop, and treating it as one adds a
  /// redundant constraint that hides real drift.
  int min_index_gap = 20;
  /// Cap on the number of keyframes retained in the descriptor database. This
  /// is the RAM discipline: without a cap the database grows without bound on
  /// a long sequence, which is exactly what the 8 GB budget forbids.
  Size max_keyframes = 400;
  /// Keyframe spacing, metres travelled between keyframes.
  Scalar keyframe_distance = 0.5;
};

/// A verified loop closure.
struct LoopClosure {
  Index from = -1;      ///< keyframe id
  Index to = -1;        ///< keyframe id
  Vec3 z = Vec3::Zero();  ///< relative transform in tangent form, x_from^{-1} x_to
  Mat3 info = Mat3::Identity();
  int inliers = 0;
  int correspondences = 0;
  Scalar inlier_ratio = 0;
  Scalar chi2 = 0;      ///< residual after the graph is optimised
};

struct LoopCloseStats {
  Size candidates_tested = 0;
  Size verification_attempts = 0;
  Size accepted = 0;
  Size rejected_inliers = 0;
  Size rejected_ratio = 0;
  Size db_evictions = 0;
  Size keyframes = 0;
};

/// Bounded keyframe database with descriptor-based retrieval and RANSAC
/// verification.
class LoopCloser {
 public:
  explicit LoopCloser(const LoopCloseOptions& opt = LoopCloseOptions());

  /// Add a keyframe. Evicts the oldest entry when the cap is reached, so the
  /// database is bounded by construction rather than by hoping the sequence is
  /// short.
  void addKeyframe(const Keyframe& kf);

  /// Look for a verified closure between `query` and the stored keyframes.
  /// Returns true and fills `out` when one is accepted.
  bool query(const Keyframe& query, LoopClosure* out);

  const LoopCloseStats& stats() const { return stats_; }
  Size size() const { return db_.size(); }
  void clear();

  /// Correspondence count between two descriptor sets, applying the ratio
  /// test. Exposed for testing.
  Size matchDescriptors(const std::vector<MapPoint>& a, const std::vector<MapPoint>& b,
                        std::vector<std::pair<int, int>>* matches) const;

 private:
  LoopCloseOptions opt_;
  std::vector<Keyframe> db_;
  LoopCloseStats stats_;
};

/// Rigid alignment of two 3D point sets (Horn's method, no scale). Returns
/// false when either set has fewer than three points or is degenerate.
bool alignPointSets(const std::vector<Vec3>& src, const std::vector<Vec3>& dst, Mat3* R, Vec3* t);

/// Verification: RANSAC over 3D-3D correspondences. Fills `inliers` with the
/// indices of the agreeing correspondences and returns the fitted transform.
bool verifyRansac(const std::vector<Vec3>& a, const std::vector<Vec3>& b,
                  const LoopCloseOptions& opt, Mat4* T, std::vector<int>* inliers,
                  Size* attempts);

}  // namespace backend
}  // namespace luna
