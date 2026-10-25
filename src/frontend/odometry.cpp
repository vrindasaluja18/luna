#include "luna/frontend/odometry.hpp"

#include <algorithm>
#include <cmath>

namespace luna {
namespace frontend {

OdomIntegrator::OdomIntegrator(const OdomOptions& opt) : opt_(opt) {}

void OdomIntegrator::add(const sim::OdomSample& odom, const std::vector<sim::ImuSample>& imu) {
  const Scalar d_forward = odom.valid ? odom.d_forward : 0.0;
  const Scalar d_yaw = odom.valid ? odom.d_yaw : 0.0;

  // Integrate the gyro over exactly the window this odometry sample covers.
  Scalar gyro_yaw = 0;
  Scalar imu_dt = 0;
  for (Size i = 1; i < imu.size(); ++i) {
    const Scalar dt = imu[i].t - imu[i - 1].t;
    if (dt <= 0 || dt > 1.0) continue;
    gyro_yaw += 0.5 * (imu[i].gyro.z() + imu[i - 1].gyro.z()) * dt;  // trapezoidal
    imu_dt += dt;
  }

  // Blend the two yaw sources once per window. Odometry yaw is accurate but
  // drifts with wheel scale error; the gyro is drift-free over short windows
  // but biased. Both describe the same interval here, which is what makes the
  // blend valid -- blending a windowed aggregate against a single sub-tick
  // gyro reading is what halves the rotation.
  const bool have_gyro = imu_dt > 1e-6;
  const Scalar yaw =
      have_gyro ? (1.0 - opt_.imu_yaw_weight) * d_yaw + opt_.imu_yaw_weight * gyro_yaw : d_yaw;

  const Pose2 inc = matToPose2(se2Exp(Vec3(d_forward, 0.0, yaw)));
  pending_ = pose2Compose(pending_, inc);
  // The running pose advances by the same increments the constraints are built
  // from, so the dead-reckoning baseline and the graph disagree only where the
  // optimiser moves things -- never because they integrated differently.
  pose_ = pose2Compose(pose_, inc);
  total_distance_ += std::abs(d_forward);
  acc_yaw_ += yaw;
  acc_dt_ += have_gyro ? imu_dt : 0.05;
  acc_distance_ += std::abs(d_forward);
  have_pending_ = true;
}

OdomConstraint OdomIntegrator::take() {
  OdomConstraint c;
  c.delta = pending_;
  c.dt = acc_dt_;
  c.distance = acc_distance_;

  // Noise scales with the distance travelled, which is the right model for
  // wheel slip and for the gyro's short-term drift. A constant sigma would make
  // a long interval look as trustworthy as a short one.
  const Scalar travel = std::max(c.distance, 0.05);
  const Scalar sig_fwd =
      std::max(std::hypot(opt_.dist_sigma, opt_.dist_scale_sigma) * travel, opt_.min_sigma);
  const Scalar sig_yaw =
      std::max(std::hypot(opt_.angle_sigma, opt_.angle_scale_sigma) * std::abs(acc_yaw_) +
                   opt_.gyro_bias_sigma * acc_dt_,
               opt_.min_sigma);

  // Cross-track uncertainty is *derived*, not assumed. A wheel odometry sensor
  // measures forward distance and yaw; it does not measure lateral motion, so
  // there is no lateral noise to model. The lateral error of the resulting
  // increment comes from the yaw error rotating the forward motion: an error of
  // sig_yaw in heading displaces the endpoint sideways by `travel * sig_yaw`.
  //
  // Modelling the unobserved axis as its own large noise term is the usual
  // shortcut and it is wrong in a way that matters. It makes the cross-track
  // variance independent of how far the vehicle actually went and of how
  // uncertain its heading was, so a stationary step and a 10 m step get the
  // same lateral uncertainty, and a heading estimate with a milliradian of
  // error is treated as being as bad as one with a tenth of a radian. The
  // covariance then over-reports by orders of magnitude, and every consumer --
  // here, the controller's chance constraint -- either stops the vehicle or
  // ignores the number entirely.
  const Scalar sig_lat = std::max(travel * sig_yaw, opt_.min_sigma);

  Mat3 cov = Mat3::Zero();
  cov(0, 0) = sig_fwd * sig_fwd;
  cov(1, 1) = sig_lat * sig_lat;
  cov(2, 2) = sig_yaw * sig_yaw;
  c.info = cov.inverse();

  pending_ = Pose2{0, 0, 0};
  acc_yaw_ = acc_dt_ = acc_distance_ = 0;
  have_pending_ = false;
  return c;
}

}  // namespace frontend
}  // namespace luna
