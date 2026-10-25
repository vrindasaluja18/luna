#include "luna/bench/trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "luna/core/se2.hpp"
#include "luna/core/stats.hpp"
#include <sstream>

namespace luna {
namespace bench {

namespace {

/// Ground truth sampled at an arbitrary time, by linear interpolation of the
/// position and shortest-arc interpolation of the yaw.
///
/// Yaw must be interpolated on the circle: between a pose at +179 deg and one
/// at -179 deg the true heading passes through 180, whereas interpolating the
/// raw numbers passes through 0 and reports a 180 deg error for a pose that is
/// nearly correct.
bool sampleAt(const Trajectory& traj, Scalar t, Pose2* out) {
  if (traj.empty()) return false;
  if (t < traj.front().t || t > traj.back().t) return false;
  // Upper bound on the first element with time > t, then step back one.
  Size hi = static_cast<Size>(
      std::upper_bound(traj.begin(), traj.end(), t,
                       [](Scalar v, const TrajectoryPoint& p) { return v < p.t; }) -
      traj.begin());
  if (hi == 0) {
    *out = traj.front().pose;
    return true;
  }
  if (hi >= traj.size()) {
    *out = traj.back().pose;
    return true;
  }
  const TrajectoryPoint& a = traj[hi - 1];
  const TrajectoryPoint& b = traj[hi];
  const Scalar dt = b.t - a.t;
  const Scalar u = dt > 1e-12 ? (t - a.t) / dt : 0.0;
  out->x = a.pose.x + u * (b.pose.x - a.pose.x);
  out->y = a.pose.y + u * (b.pose.y - a.pose.y);
  out->th = wrapPi(a.pose.th + u * wrapPi(b.pose.th - a.pose.th));
  return true;
}

}  // namespace

Scalar pathLength(const Trajectory& traj) {
  Scalar len = 0;
  for (Size i = 1; i < traj.size(); ++i) {
    len += std::hypot(traj[i].pose.x - traj[i - 1].pose.x, traj[i].pose.y - traj[i - 1].pose.y);
  }
  return len;
}

AteResult absoluteTrajectoryError(const Trajectory& estimate, const Trajectory& truth,
                                  const EvalOptions& opt) {
  AteResult res;
  if (estimate.empty() || truth.size() < 2) return res;

  res.path_length = pathLength(truth);

  std::vector<Scalar> errs;
  errs.reserve(estimate.size());
  for (const TrajectoryPoint& e : estimate) {
    Pose2 gt;
    if (!sampleAt(truth, e.t, &gt)) {
      ++res.unmatched;
      continue;
    }
    errs.push_back(std::hypot(e.pose.x - gt.x, e.pose.y - gt.y));
  }
  if (errs.empty()) return res;

  RunningStats stats;
  for (Scalar v : errs) stats.add(v);
  std::sort(errs.begin(), errs.end());

  res.matched = errs.size();
  res.rmse = std::sqrt(stats.mean() * stats.mean() + stats.variance());
  res.mean = stats.mean();
  res.median = errs[errs.size() / 2];
  res.max = errs.back();
  res.percent_of_distance = res.path_length > 1e-9 ? 100.0 * res.rmse / res.path_length : 0.0;
  res.ok = true;
  return res;
}

RpeResult relativePoseError(const Trajectory& estimate, const Trajectory& truth,
                            const EvalOptions& opt) {
  RpeResult res;
  res.interval = opt.rpe_interval;
  if (estimate.size() < 2 || truth.size() < 2) return res;

  RunningStats trans, rot;
  for (const TrajectoryPoint& e : estimate) {
    Pose2 gt_a, gt_b;
    const Scalar t_b = e.t + opt.rpe_interval;
    if (!sampleAt(truth, e.t, &gt_a) || !sampleAt(truth, t_b, &gt_b)) continue;
    // Find the estimate pose closest to t_b, within half a sample period.
    Pose2 est_b;
    if (!sampleAt(estimate, t_b, &est_b)) continue;

    const Pose2 rel_gt = pose2Relative(gt_a, gt_b);
    const Pose2 rel_est = pose2Relative(e.pose, est_b);
    trans.add(std::hypot(rel_est.x - rel_gt.x, rel_est.y - rel_gt.y));
    rot.add(std::abs(wrapPi(rel_est.th - rel_gt.th)) * 180.0 / kPi);
    ++res.pairs;
  }
  if (res.pairs == 0) return res;
  res.rmse_trans = std::sqrt(trans.mean() * trans.mean() + trans.variance());
  res.mean_trans = trans.mean();
  res.rmse_rot = std::sqrt(rot.mean() * rot.mean() + rot.variance());
  res.mean_rot = rot.mean();
  res.ok = true;
  return res;
}

bool writeTum(const std::string& path, const Trajectory& traj) {
  std::ofstream f(path);
  if (!f) return false;
  // Full double precision. The default 6 significant digits is enough to plot a
  // trajectory but not enough to round-trip one: at 20 m the quantisation is
  // ~1e-4 m, which is the same order as the RPE the evaluator is asked to
  // resolve. A reproducibility kit whose own output cannot be read back
  // exactly is not reproducible.
  f << std::setprecision(17);
  for (const TrajectoryPoint& p : traj) {
    const Scalar half = 0.5 * p.pose.th;
    // Yaw-only quaternion: q = (0, 0, sin(theta/2), cos(theta/2)), with the
    // rotation about +z.
    f << p.t << ' ' << p.pose.x << ' ' << p.pose.y << " 0 0 0 " << std::sin(half) << ' '
      << std::cos(half) << '\n';
  }
  return f.good();
}

bool readTum(const std::string& path, Trajectory* out) {
  if (!out) return false;
  std::ifstream f(path);
  if (!f) return false;
  out->clear();
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    TrajectoryPoint p;
    Scalar qx = 0, qy = 0, qz = 0, qw = 1;
    if (!(ss >> p.t >> p.pose.x >> p.pose.y)) continue;
    // z is ignored (planar), but must be consumed if present.
    Scalar z = 0;
    if (ss >> z) {
      if (ss >> qx >> qy >> qz >> qw) {
        // Recover yaw from the quaternion's z/w components.
        p.pose.th = wrapPi(2.0 * std::atan2(qz, qw));
      }
    }
    out->push_back(p);
  }
  return !out->empty();
}

}  // namespace bench
}  // namespace luna
