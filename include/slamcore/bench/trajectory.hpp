// SLAMCore-8G -- trajectory accuracy metrics (ATE, RPE).
//
// Absolute Trajectory Error is the standard SLAM yardstick: align the estimate
// to ground truth, then measure the position error at matched timestamps.
// Relative Pose Error measures local drift over a fixed interval, which is the
// metric that exposes scale and yaw drift even when the estimate happens to
// return to the right place at the end.
//
// Both are reported in metres *and* as a fraction of the ground-truth path
// length, because "0.3 % of distance" is the contract and a raw metre figure
// means nothing without the distance it was accumulated over.
#pragma once

#include <string>
#include <vector>

#include "slamcore/core/se3.hpp"

namespace sc {
namespace bench {

/// A timestamped planar pose. Timestamps are what make an estimate comparable
/// to ground truth sampled at a different rate.
struct TrajectoryPoint {
  Scalar t = 0;
  Pose2 pose;
};

using Trajectory = std::vector<TrajectoryPoint>;

struct AteResult {
  Size matched = 0;          ///< number of estimate poses that found ground truth
  Size unmatched = 0;        ///< estimate poses with no ground truth in range
  Scalar rmse = 0;           ///< metres
  Scalar mean = 0;
  Scalar median = 0;
  Scalar max = 0;
  Scalar path_length = 0;    ///< ground-truth path length, metres
  Scalar percent_of_distance = 0;  ///< rmse / path_length * 100
  bool ok = false;
};

struct RpeResult {
  Size pairs = 0;
  Scalar rmse_trans = 0;     ///< metres over the interval
  Scalar mean_trans = 0;
  Scalar rmse_rot = 0;       ///< degrees over the interval
  Scalar mean_rot = 0;
  Scalar interval = 0;       ///< the fixed interval used, seconds
  bool ok = false;
};

struct EvalOptions {
  /// Ground truth is interpolated at each estimate timestamp. An estimate
  /// outside the ground-truth time span by more than this is dropped rather
  /// than extrapolated, because extrapolated ground truth is a guess and would
  /// flatter or penalise the estimate arbitrarily.
  Scalar max_time_gap = 0.05;
  /// Fixed interval for RPE, seconds.
  Scalar rpe_interval = 1.0;
};

/// Position-only ATE. No alignment is performed: the simulator's estimate and
/// its ground truth share a frame by construction, so the metric measures
/// drift, not an unknown rigid transform. Applying a best-fit transform here
/// would hide exactly the error the contract is about.
AteResult absoluteTrajectoryError(const Trajectory& estimate, const Trajectory& truth,
                                  const EvalOptions& opt = EvalOptions());

/// Relative pose error over a fixed time interval.
RpeResult relativePoseError(const Trajectory& estimate, const Trajectory& truth,
                            const EvalOptions& opt = EvalOptions());

/// Total path length of a trajectory, in metres.
Scalar pathLength(const Trajectory& traj);

/// Read/write TUM trajectory format: "t x y z qx qy qz qw" per line. The
/// simulator writes z = 0 and a yaw-only quaternion, which is exactly what a
/// planar estimate has, so the files are directly comparable with the public
/// TUM tooling.
bool writeTum(const std::string& path, const Trajectory& traj);
bool readTum(const std::string& path, Trajectory* out);

}  // namespace bench
}  // namespace sc
