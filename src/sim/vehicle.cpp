#include "slamcore/sim/vehicle.hpp"

#include <cmath>
#include <fstream>

namespace sc {
namespace sim {

Vehicle::Vehicle(const VehicleParams& params, const Pose2& initial, std::uint64_t seed)
    : params_(params) {
  state_.pose = initial;
  (void)seed;
}

void Vehicle::integrate(const Control& cmd, Scalar dt) {
  // Bicycle model with the rear axle as reference point:
  //   x' = v cos(th),  y' = v sin(th),  th' = v tan(steer) / L
  const Scalar L = params_.wheelbase;
  const auto deriv = [&](const VehicleState& s, Scalar steer, Scalar accel) {
    const Scalar yaw_rate = clamp(s.v * std::tan(steer) / L, -params_.max_yaw_rate, params_.max_yaw_rate);
    return std::array<Scalar, 4>{s.v * std::cos(s.pose.th), s.v * std::sin(s.pose.th), yaw_rate, accel};
  };

  // Clamp the command to the actuator envelope.
  Control c = cmd;
  c.v = clamp(c.v, -params_.max_speed, params_.max_speed);
  c.steer = clamp(c.steer, -params_.max_steer, params_.max_steer);
  c.accel = clamp(c.accel, -params_.max_accel, params_.max_accel);

  // Two input conventions are supported, because both callers are legitimate:
  // a velocity command `v` (the planner's natural output) is tracked at the
  // acceleration limit, while a nonzero `accel` is taken as an explicit
  // acceleration input (what the MPC produces). Treating `v` as inert, as an
  // accel-only plant does, makes a speed command silently do nothing, which is
  // exactly the kind of failure that looks like a controller bug downstream.
  Scalar accel = c.accel;
  if (accel == 0.0) {
    const Scalar dv = c.v - state_.v;
    accel = clamp(dv / std::max(dt, 1e-9), -params_.max_accel, params_.max_accel);
  }

  // RK4 on (x, y, th, v); steering slews toward its command at the rate limit.
  VehicleState s = state_;
  const Scalar steer_step = clamp(c.steer - s.steer, -params_.max_steer_rate * dt, params_.max_steer_rate * dt);
  const Scalar steer_new = s.steer + steer_step;

  const auto k1 = deriv(s, s.steer, accel);
  VehicleState s2 = s;
  s2.pose.x += 0.5 * dt * k1[0];
  s2.pose.y += 0.5 * dt * k1[1];
  s2.pose.th += 0.5 * dt * k1[2];
  s2.v += 0.5 * dt * accel;
  const auto k2 = deriv(s2, s.steer + 0.5 * steer_step, accel);

  VehicleState s3 = s;
  s3.pose.x += 0.5 * dt * k2[0];
  s3.pose.y += 0.5 * dt * k2[1];
  s3.pose.th += 0.5 * dt * k2[2];
  s3.v += 0.5 * dt * accel;
  const auto k3 = deriv(s3, s.steer + 0.5 * steer_step, accel);

  VehicleState s4 = s;
  s4.pose.x += dt * k3[0];
  s4.pose.y += dt * k3[1];
  s4.pose.th += dt * k3[2];
  s4.v += dt * accel;
  const auto k4 = deriv(s4, steer_new, accel);

  state_.pose.x += dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  state_.pose.y += dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  state_.pose.th = wrapPi(state_.pose.th + dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]));
  state_.v += dt / 6.0 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
  state_.v = clamp(state_.v, -params_.max_speed, params_.max_speed);
  state_.steer = steer_new;
  state_.accel = accel;
  state_.t += dt;
  travelled_ += std::abs(state_.v) * dt;
}

const VehicleState& Vehicle::step(const Control& cmd, Scalar dt) {
  integrate(cmd, dt);
  return state_;
}

const VehicleState& Vehicle::follow(const std::vector<Vec2>& path, Scalar speed, Scalar dt) {
  if (path.size() < 2) return state_;
  // Pure-pursuit style waypoint chaser: pick the first waypoint at least
  // `lookahead` ahead and steer toward it.
  const Scalar lookahead = clamp(0.6 + 0.5 * std::abs(speed), 0.8, 3.0);
  Vec2 target = path.back();
  for (Size i = path_idx_; i < path.size(); ++i) {
    if ((path[i] - Vec2(state_.pose.x, state_.pose.y)).norm() >= lookahead) {
      target = path[i];
      path_idx_ = i;
      break;
    }
  }
  const Vec2 rel = target - Vec2(state_.pose.x, state_.pose.y);
  const Scalar c = std::cos(state_.pose.th), s = std::sin(state_.pose.th);
  const Scalar dx = c * rel.x() + s * rel.y();
  const Scalar dy = -s * rel.x() + c * rel.y();
  const Scalar dist = std::max(rel.norm(), 1e-3);
  const Scalar alpha = std::atan2(dy, dx);
  Scalar steer = std::atan2(2.0 * params_.wheelbase * std::sin(alpha), dist);
  Control cmd;
  cmd.steer = clamp(steer, -params_.max_steer, params_.max_steer);
  const Scalar dv = clamp(speed - state_.v, -params_.max_accel * dt, params_.max_accel * dt);
  cmd.v = state_.v + dv;
  cmd.accel = dv / dt;
  integrate(cmd, dt);
  return state_;
}

bool Vehicle::inCollision(const GridMap& grid) const {
  // Test the four corners plus the centre of the body rectangle.
  const Scalar hl = params_.length * 0.5;
  const Scalar hw = params_.track_width * 0.5;
  const Vec2 corners[5] = {
      Vec2(hl, hw), Vec2(hl, -hw), Vec2(-hl, hw), Vec2(-hl, -hw), Vec2(0, 0),
  };
  for (const Vec2& c : corners) {
    if (grid.atWorld(pose2Apply(state_.pose, c)) == Cell::kOccupied) return true;
  }
  return false;
}

Scalar GroundTruthLog::pathLength() const {
  Scalar d = 0;
  for (Size i = 1; i < frames_.size(); ++i) {
    const Vec2 a(frames_[i - 1].pose.x, frames_[i - 1].pose.y);
    const Vec2 b(frames_[i].pose.x, frames_[i].pose.y);
    d += (b - a).norm();
  }
  return d;
}

bool GroundTruthLog::saveTum(const std::string& path) const {
  std::ofstream f(path);
  if (!f) return false;
  for (const auto& fr : frames_) {
    const Quat q(Eigen::AngleAxis<Scalar>(fr.pose.th, Vec3::UnitZ()));
    f << fr.t << ' ' << fr.pose.x << ' ' << fr.pose.y << " 0 "
      << q.x() << ' ' << q.y() << ' ' << q.z() << ' ' << q.w() << '\n';
  }
  return static_cast<bool>(f);
}

bool GroundTruthLog::saveLandmarks(const std::string& path) const {
  std::ofstream f(path);
  if (!f) return false;
  for (const auto& fr : frames_) {
    f << fr.t;
    for (Size i = 0; i < fr.landmark_ids.size(); ++i) {
      f << ' ' << fr.landmark_ids[i] << ' ' << fr.landmark_px[i].x() << ' ' << fr.landmark_px[i].y();
    }
    f << '\n';
  }
  return static_cast<bool>(f);
}

}  // namespace sim
}  // namespace sc
