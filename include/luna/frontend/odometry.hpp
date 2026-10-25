// LUNA -- wheel-odometry / IMU dead reckoning.
//
// This is the baseline the pose graph improves on, and it is kept separate from
// the graph so the benchmark can report both numbers. Folding dead reckoning
// into the optimiser would leave no way to show that the optimiser does
// anything.
//
// The sensor model matters. Wheel odometry reports *forward distance and yaw
// only*: the lateral component of the relative transform is not measured at all
// (see OdomEmulator::measure, which drops it). So the measurement observes two
// of a transform's three degrees of freedom, and the information matrix must say
// so -- it is very *soft* across-track. Modelling the unobserved direction as
// tight is what makes a dead-reckoning run look far better than it is on any
// route with lateral motion.
//
// Two clocks feed this: odometry ticks at 20 Hz and the IMU at 200 Hz, while
// the caller may ask for an increment at any point. Both sensors are therefore
// accumulated raw over the whole interval and blended once, in take(). Blending
// per call instead is a real error: a call that contains no odometry tick has
// d_yaw = 0 but a full interval of gyro yaw, and averaging those two halves the
// true rotation on every other frame.
#pragma once

#include <vector>

#include "luna/core/se2.hpp"
#include "luna/sim/sensors.hpp"

namespace luna {
namespace frontend {

struct OdomOptions {
  /// Standard deviation of the forward-distance *random* error, metres per
  /// metre travelled.
  Scalar dist_sigma = 0.01;
  /// Standard deviation of the forward-distance *scale* error, metres per metre.
  ///
  /// Separate from `dist_sigma` because it is a different failure mode: a wheel
  /// that reads 1% short is wrong by 1% of every metre, every time. Over one
  /// keyframe interval that looks like random noise and the estimator cannot
  /// tell them apart, but it does not average away over a trajectory -- it
  /// accumulates as a proportional error, which is exactly what the reported
  /// drift is made of. A covariance that omits it is overconfident by
  /// construction, and a chance constraint built on it would be too.
  Scalar dist_scale_sigma = 0.01;
  /// Standard deviation of the yaw error, radians per radian turned.
  Scalar angle_sigma = 0.005;
  /// Yaw *scale* error, radians per radian turned. Same argument as
  /// `dist_scale_sigma`: a systematic mis-calibration that the least-squares
  /// problem treats as noise and that therefore biases the covariance low.
  Scalar angle_scale_sigma = 0.01;
  /// Gyro-bias random-walk allowance, radians per second.
  Scalar gyro_bias_sigma = 1e-4;
  /// Floor on every standard deviation, so a stationary step does not produce
  /// an infinitely stiff constraint that welds two poses together.
  Scalar min_sigma = 1e-3;
  /// Weight of the IMU yaw rate against odometry yaw in the fused heading.
  /// Odometry yaw is accurate but drifts with the wheel scale error; the gyro
  /// is drift-free over short windows but biased. A complementary blend uses
  /// each where it is strong.
  Scalar imu_yaw_weight = 0.5;
};

/// A relative-pose measurement with its information matrix.
struct OdomConstraint {
  Pose2 delta;      ///< relative transform, x_i^{-1} x_j
  Mat3 info;        ///< inverse covariance
  Scalar dt = 0;    ///< time spanned
  Scalar distance = 0;
};

/// Integrates odometry and IMU yaw rate into relative-pose constraints while
/// also maintaining a running dead-reckoned pose.
class OdomIntegrator {
 public:
  explicit OdomIntegrator(const OdomOptions& opt = OdomOptions());

  /// Fold one odometry sample and the IMU batch that arrived with it into the
  /// pending increment, and advance the running pose.
  void add(const sim::OdomSample& odom, const std::vector<sim::ImuSample>& imu);

  /// Take the accumulated increment as a constraint and reset the accumulator.
  /// The running pose is *not* reset: the caller wants both a per-frame
  /// trajectory and per-keyframe increments, and resetting the pose here would
  /// force the caller to integrate a second time.
  OdomConstraint take();

  /// True when an increment is pending.
  bool pending() const { return have_pending_; }

  /// The dead-reckoned pose, i.e. the composition of every increment taken so
  /// far. This is what the odometry-only baseline trajectory is built from.
  const Pose2& pose() const { return pose_; }

  /// Total distance integrated, for diagnostics.
  Scalar totalDistance() const { return total_distance_; }

 private:
  OdomOptions opt_;
  /// Accumulated transform for the pending interval, composed from the same
  /// per-window increments that advance the running pose.
  Pose2 pending_{0, 0, 0};
  /// Bookkeeping for the pending interval, used to size the noise.
  Scalar acc_yaw_ = 0, acc_dt_ = 0, acc_distance_ = 0;
  // Running dead-reckoned pose.
  Pose2 pose_{0, 0, 0};
  Scalar total_distance_ = 0;
  bool have_pending_ = false;
};

}  // namespace frontend
}  // namespace luna
