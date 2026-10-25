// SLAMCore-8G -- kinematic vehicle model and trajectory ground truth.
//
// The plant is a bicycle model with velocity and steering-rate inputs, which is
// the same model the MPC assumes. Keeping plant and controller identical makes
// the closed-loop tracking metric measure the *controller*, not a plant
// mismatch; plant mismatch is injected deliberately by the fault module.
#pragma once

#include <functional>
#include <vector>

#include "slamcore/core/se3.hpp"
#include "slamcore/core/stats.hpp"
#include "slamcore/sim/world.hpp"

namespace sc {
namespace sim {

struct VehicleParams {
  // Geometry is chosen so the body's *diagonal* fits the generator's default
  // 1.2 m doorway with margin to spare. A 0.8 x 0.45 m body has a 0.92 m
  // diagonal, which leaves only 0.14 m each side of a doorway and turns every
  // doorway into a collision even with a perfect path. The body is sized for
  // the world it is meant to drive in, not for the world to be re-sized around
  // an arbitrarily large body.
  Scalar wheelbase = 0.5;         ///< metres, front-to-rear axle
  Scalar max_speed = 2.5;         ///< m/s
  Scalar max_accel = 2.0;         ///< m/s^2
  Scalar max_steer = 0.75;        ///< rad
  Scalar max_steer_rate = 2.0;    ///< rad/s
  Scalar max_yaw_rate = 1.5;      ///< rad/s
  Scalar track_width = 0.34;      ///< body width used for clearance checks
  Scalar length = 0.62;           ///< body length
};

/// Commanded control input.
struct Control {
  Scalar v = 0;      ///< longitudinal speed (m/s)
  Scalar steer = 0;  ///< steering angle (rad)
  Scalar accel = 0;  ///< commanded acceleration (m/s^2), for the MPC plant
};

/// Full ground-truth state, including the derivatives the IMU emulator needs.
struct VehicleState {
  Pose2 pose;
  Scalar v = 0;
  Scalar steer = 0;
  Scalar accel = 0;
  Scalar t = 0;
};

/// Integrates the bicycle model with fixed-step RK4.
class Vehicle {
 public:
  Vehicle(const VehicleParams& params, const Pose2& initial, std::uint64_t seed);

  /// Advance by `dt` under a constant command. Returns the new state.
  const VehicleState& step(const Control& cmd, Scalar dt);

  /// Advance along a reference path (used to synthesise dataset-like runs).
  const VehicleState& follow(const std::vector<Vec2>& path, Scalar speed, Scalar dt);

  const VehicleState& state() const { return state_; }
  const VehicleParams& params() const { return params_; }
  void setState(const VehicleState& s) { state_ = s; }

  /// True when the body footprint overlaps an occupied cell.
  bool inCollision(const GridMap& grid) const;

  /// Accumulated path length in metres.
  Scalar distanceTravelled() const { return travelled_; }

 private:
  void integrate(const Control& cmd, Scalar dt);

  VehicleParams params_;
  VehicleState state_;
  Scalar travelled_ = 0;
  Size path_idx_ = 0;
};

/// Ground-truth logger: records the full state history and, optionally, the
/// per-frame landmark observations that the evaluator uses to score SLAM.
struct GroundTruthFrame {
  Scalar t = 0;
  Pose2 pose;
  Scalar v = 0;
  std::vector<std::uint32_t> landmark_ids;
  std::vector<Vec2> landmark_px;
};

class GroundTruthLog {
 public:
  void add(const GroundTruthFrame& f) { frames_.push_back(f); }
  const std::vector<GroundTruthFrame>& frames() const { return frames_; }
  void clear() { frames_.clear(); }
  Size size() const { return frames_.size(); }

  /// Total travelled distance of the logged trajectory.
  Scalar pathLength() const;

  /// Save in TUM trajectory format: `t x y z qx qy qz qw`.
  bool saveTum(const std::string& path) const;
  /// Save the landmark id/pixel stream for the evaluator.
  bool saveLandmarks(const std::string& path) const;

 private:
  std::vector<GroundTruthFrame> frames_;
};

}  // namespace sim
}  // namespace sc
