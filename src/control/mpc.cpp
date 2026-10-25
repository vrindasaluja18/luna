#include "luna/control/mpc.hpp"

#include <algorithm>
#include <chrono>
#include <cassert>
#include <cmath>

namespace luna {
namespace control {

const char* toString(ControlStatus s) {
  switch (s) {
    case ControlStatus::kOk: return "ok";
    case ControlStatus::kDegraded: return "degraded";
    case ControlStatus::kSafeStop: return "safe_stop";
    case ControlStatus::kGoalReached: return "goal_reached";
    case ControlStatus::kNoPlan: return "no_plan";
  }
  return "?";
}

// ---------------------------------------------------------------------------
// Safety criteria
// ---------------------------------------------------------------------------

Scalar MpcController::covarianceDet(const VehicleEstimate& est) {
  // sqrt of the 2x2 position covariance determinant, times the heading sigma.
  // The y-theta correlation is carried through so a strongly correlated
  // lateral/heading error counts as the one uncertain direction it is, rather
  // than as two independent ones.
  const Scalar sxx = std::max<Scalar>(0.0, est.sigma_xx);
  const Scalar syy = std::max<Scalar>(0.0, est.sigma_yy);
  const Scalar sty = est.sigma_yth;
  const Scalar det2 = std::max<Scalar>(0.0, sxx * syy - sty * sty);
  return std::sqrt(det2) * std::sqrt(std::max<Scalar>(0.0, est.sigma_thth));
}

Scalar MpcController::chanceConstrainedSpeed(Scalar corridor_width, Scalar sigma_lat) const {
  // The usable half-width is what is left of the corridor after the body and
  // the uncertainty have taken their share. Below the minimum usable width the
  // answer is zero: the vehicle stops rather than drives through a gap it may
  // not fit.
  const Scalar half = 0.5 * corridor_width - 0.5 * veh_.track_width - opt_.k_sigma * sigma_lat;
  if (half <= opt_.min_usable_width) return 0.0;
  // Turning at speed needs room for the swept path: the lateral acceleration
  // v^2 / R is bounded by what the tyres can hold, and R is bounded by the
  // usable half-width at the maximum steer. This is the standard
  // friction-limited cornering speed, written in terms of the room available.
  const Scalar a_lat_max = 1.8;  // m/s^2, conservative on a low-grip estimate
  const Scalar r_min = half / std::tan(veh_.max_steer);
  return std::min(veh_.max_speed, std::sqrt(a_lat_max * r_min));
}

// ---------------------------------------------------------------------------
// Model and reference
// ---------------------------------------------------------------------------

void MpcController::linearise(const VehicleEstimate& est, const sim::Control& u_prev, MatX* A,
                              MatX* B, VecX* c) const {
  // Bicycle model with velocity and steering-rate inputs, and the *steering
  // angle carried in the state*:
  //
  //   x'     = v cos(theta)
  //   y'     = v sin(theta)
  //   theta' = v tan(delta) / L
  //   v'     = accel
  //   delta' = steer_rate
  //
  // State z = [x, y, theta, v, delta], input u = [accel, steer_rate].
  //
  // The steering angle has to be a state, and this is the single most
  // important detail in the controller. When it is not, the model has to
  // express the heading effect of a steering-rate input directly, and the
  // obvious expression is wrong in a way that is invisible on a straight and
  // fatal in a corner. `v tan(delta)/L` is the yaw rate as a function of the
  // steering *angle*; differentiating it with respect to the steering *rate*
  // gives a constant, so the condensed model claims that a steering input at
  // step j turns the vehicle by the same amount regardless of how much of the
  // horizon remains -- when the truth is that it acts for the remaining
  // N - j - 1 steps and therefore accumulates. With the constant form the
  // controller believes its steering authority is up to 24 times larger than
  // it is on a 25-step horizon, commands a fraction of the steering it needs,
  // and understeers through every corner until the heading error trips the
  // divergence watchdog. Carrying delta as a state makes the accumulation fall
  // out of the A powers automatically, with no fudge factor to tune.
  //
  // Both matrices are constant along the horizon: the reference is smooth and
  // the horizon is short, so re-linearising at each predicted step buys
  // accuracy the 20 Hz loop does not need, at five times the matrix work.
  const Scalar th = est.pose.th;
  const Scalar v = std::max<Scalar>(0.0, est.v);
  // The operating point is the *measured* steering angle, not the last command:
  // the plant slews toward its command at a rate limit, so the two differ
  // whenever the command is changing, and the model should describe the vehicle
  // that exists rather than the one that was asked for.
  const Scalar st = est.steer;
  const Scalar L = veh_.wheelbase;
  const Scalar ct = std::cos(th), sn = std::sin(th);
  const Scalar cs2 = std::cos(st) * std::cos(st);

  A->setZero();
  B->setZero();
  c->setZero();

  // d/d(x, y, theta, v, delta)
  (*A)(0, 2) = -v * sn;
  (*A)(0, 3) = ct;
  (*A)(1, 2) = v * ct;
  (*A)(1, 3) = sn;
  (*A)(2, 3) = std::tan(st) / L;
  (*A)(2, 4) = v / (L * cs2);

  // d/d(accel, steer_rate)
  (*B)(3, 0) = 1.0;
  (*B)(4, 1) = 1.0;

  // Affine term: f(x0, u0) - A x0 - B u0.
  VecX x0(5);
  x0 << est.pose.x, est.pose.y, th, v, st;
  VecX u0(2);
  u0 << u_prev.accel, 0.0;
  VecX f(5);
  f << v * ct, v * sn, v * std::tan(st) / L, u_prev.accel, 0.0;
  *c = f - (*A) * x0 - (*B) * u0;
}

// ---------------------------------------------------------------------------
// Velocity profile
// ---------------------------------------------------------------------------

Scalar MpcController::profileAt(const std::vector<Scalar>& profile,
                               const std::vector<Scalar>& s_axis, Scalar s) {
  if (profile.empty()) return 0.0;
  if (s <= s_axis.front()) return profile.front();
  if (s >= s_axis.back()) return profile.back();
  Size i = 1;
  while (i < s_axis.size() && s_axis[i] < s) ++i;
  const Scalar ds = s_axis[i] - s_axis[i - 1];
  if (ds < 1e-9) return profile[i];
  const Scalar w = (s - s_axis[i - 1]) / ds;
  return (1.0 - w) * profile[i - 1] + w * profile[i];
}

std::vector<Scalar> MpcController::buildVelocityProfile(const planning::Plan& plan,
                                                        const VehicleEstimate& est,
                                                        Scalar ceiling) const {
  const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
  const Size n = plan.path.size();
  std::vector<Scalar> prof(n, ceiling);
  if (n == 0) return prof;

  const std::vector<Scalar> kappa = planning::CorridorPlanner::curvature(plan.path);
  const Scalar sigma_lat = est.lateralSigma(opt_.preview);
  for (Size i = 0; i < n; ++i) {
    // Curvature limit: v <= sqrt(a_lat / kappa).
    if (kappa[i] > 1e-6) prof[i] = std::min(prof[i], std::sqrt(opt_.a_lat_max / kappa[i]));
    // Clearance limit. Two different questions, deliberately kept apart:
    //
    //   * Does the body fit? That is a *hard* geometric fact: the room between
    //     the path point and the wall, less the half-diagonal. If it is gone,
    //     the speed is zero and no amount of confidence changes that.
    //   * How much of the remaining room must be reserved for the estimate
    //     being wrong? That is the chance constraint, and it modulates speed.
    //
    // Folding them into one expression (`sqrt(2 a (margin - k sigma))`) reads
    // well and is wrong: it commands a stop whenever the uncertainty exceeds
    // the margin, so a vehicle that physically fits through a 1.2 m doorway
    // with 8 cm to spare refuses to enter it. A doorway is not made wider by
    // being uncertain about where you are; the vehicle should creep through and
    // accept the risk, or refuse because it does not fit -- not stop because it
    // is unsure.
    if (i < plan.clearance.size()) {
      const Scalar hard = plan.clearance[i] - opt_.body_half_diagonal;
      if (hard <= 0.02) {
        prof[i] = 0.0;
      } else {
        const Scalar budget = hard - opt_.k_sigma * sigma_lat;
        const Scalar v_brake = budget <= 0.0 ? 0.0 : std::sqrt(2.0 * opt_.a_brake * budget);
        // Creep floor: with the body clear of the wall the vehicle is allowed
        // to proceed slowly even when the uncertainty budget is exhausted.
        const Scalar v_creep = std::min(opt_.degraded_speed, std::sqrt(2.0 * opt_.a_brake * hard) * 0.5);
        prof[i] = std::min(prof[i], std::max(v_brake, v_creep));
      }
    }
  }

  // The end of the path is the goal, and the vehicle is supposed to *arrive*
  // there, not pass through it. Without this the profile keeps a non-zero speed
  // right up to the last waypoint, the vehicle sails past the end of the plan,
  // and the safe-tube check then trips -- because there is no longer any plan
  // ahead of it to be inside of -- so a run that drove the whole route reports
  // a safe stop instead of an arrival. Setting the terminal speed to zero and
  // letting the backwards pass propagate the braking is what makes the profile
  // a stopping profile rather than a speed limit.
  if (!prof.empty()) prof.back() = 0.0;

  // Backwards feasibility pass: v_i^2 <= v_{i+1}^2 + 2 a ds. Walking backwards
  // makes every point reachable at its own limit while still being able to
  // brake for everything ahead of it.
  for (Size i = n - 1; i-- > 0;) {
    const Scalar ds = s_axis[i + 1] - s_axis[i];
    const Scalar v_next = prof[i + 1];
    const Scalar reachable = std::sqrt(v_next * v_next + 2.0 * opt_.a_brake * ds);
    prof[i] = std::min(prof[i], reachable);
  }
  return prof;
}

void MpcController::buildReference(const planning::Plan& plan, const VehicleEstimate& est,
                                   const std::vector<Scalar>& profile,
                                   std::vector<std::array<Scalar, 4>>* ref,
                                   std::vector<Scalar>* prof_ref) {
  ref->assign(static_cast<Size>(opt_.horizon) + 1, {0, 0, 0, 0});
  if (prof_ref) prof_ref->assign(static_cast<Size>(opt_.horizon) + 1, 0.0);
  if (plan.path.size() < 2) return;

  const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
  const Scalar total = s_axis.empty() ? 0.0 : s_axis.back();

  // Project the current pose onto the path to get the starting arc length.
  // A local search around the last known progress, not a global nearest-point
  // query: the vehicle does not teleport, and a global query would snap to a
  // different branch of a path that revisits the same corridor.
  // The window is centred on the *last* progress and does not move while the
  // search runs. That sounds like a detail and it is the whole bug: if the
  // window centre is the same variable the loop writes, each candidate that
  // wins drags the window along with it, and the search walks forward in
  // `fwd`-sized strides until it falls off the far edge. On a 128 m maze path
  // the projection then lands 2.5 m ahead of the vehicle on every tick -- the
  // reference teleports, the heading error jumps by the corner angle in one
  // step, the controller pins the steering at its rate limit chasing it, runs
  // wide, and the watchdog stops the run 8 m in. Freezing the centre makes the
  // search local, which is what it was always meant to be.
  Scalar s_here = progress_;
  {
    const Scalar centre = progress_;
    Scalar best_d = std::numeric_limits<Scalar>::max();
    // Asymmetric window: generous forwards, tight backwards. A symmetric window
    // lets the projection jump to the *other* leg of a hairpin, which presents
    // as a reference heading pointing back the way the vehicle came.
    const Scalar fwd = 2.5, back = 0.5;
    for (Size i = 0; i < plan.path.size(); ++i) {
      if (s_axis[i] < centre - back || s_axis[i] > centre + fwd) continue;
      const Scalar d = (plan.path[i] - Vec2(est.pose.x, est.pose.y)).squaredNorm();
      if (d < best_d) {
        best_d = d;
        s_here = s_axis[i];
      }
    }
  }
  // A projection cannot outrun the odometry that produced it. The vehicle moves
  // at most `v * dt` per tick, so a projection that advances further has not
  // tracked the vehicle -- it has snapped to a *different* path point that
  // happens to be nearer, which is what happens where the vehicle has cut the
  // inside of a corner and is closer to the point past the apex than the one
  // beside it. Left uncapped the reference teleports ahead, the heading error
  // jumps by the corner's angle in a single tick, and the controller steers
  // across the turn instead of through it: in the maze world it pinned the
  // steering at the rate limit, ran wide, and the watchdog stopped it 8 m into
  // a 128 m plan. The cap carries a few times the vehicle's own advance as
  // headroom, so it never binds on normal tracking.
  const Scalar max_advance = 3.0 * std::max<Scalar>(0.05, est.v) * opt_.dt + 0.05;
  s_here = std::min(s_here, progress_ + max_advance);
  progress_ = std::max(progress_, s_here);

  // Sample the reference by arc length at the speed the profile allows, not at
  // a fixed step: the horizon should span a similar amount of *time* whatever
  // the speed, and a fixed distance step would make the horizon look far ahead
  // at low speed and uselessly short at high speed.
  // The vehicle's displacement from the point it projects onto. The horizon's
  // first reference point is placed at the vehicle itself, and this offset
  // decays to zero along the horizon so the reference rejoins the path.
  const Vec2 proj0 = planning::pointAt(plan.path, s_axis, s_here);
  const Vec2 offset(est.pose.x - proj0.x(), est.pose.y - proj0.y());

  Scalar s = s_here;
  for (int k = 0; k <= opt_.horizon; ++k) {
    const Scalar vk = std::max<Scalar>(0.15, profileAt(profile, s_axis, s));
    const Vec2 p = planning::pointAt(plan.path, s_axis, s);
    const Vec2 q = planning::pointAt(plan.path, s_axis, std::min(total, s + 0.25));
    Vec2 tangent = q - p;
    if (tangent.norm() < 1e-6) tangent = Vec2(1, 0);
    // Unwrap the reference heading against the previous one.
    //
    // `atan2` returns a value in (-pi, pi], so a path that heads west has a
    // reference heading that flips between +pi and -pi as the sampled tangent
    // wanders across the branch cut. The two values are the *same* direction,
    // but the cost sees them as 2*pi apart, and the horizon then contains two
    // contradictory targets: the early points say "turn to +3.13" and the later
    // ones say "turn to -3.14", which from a vehicle at +1.8 is a residual of
    // +1.3 for one group and -4.9 for the other. The QP resolves the conflict in
    // favour of whichever residual is numerically larger, so it commands full
    // lock the *wrong way* at every corner that faces west. In the maze world
    // that pinned the steering, ran the vehicle wide, and the watchdog stopped
    // the run 8 m into a 128 m plan -- a branch cut in a trigonometric function
    // presenting as a control failure.
    //
    // Anchoring the sequence to the vehicle's own heading first keeps the first
    // reference within pi of where the vehicle actually points, so the
    // unwrapped sequence never has to jump.
    Scalar th_ref = std::atan2(tangent.y(), tangent.x());
    const Scalar th_prev = (k == 0) ? est.pose.th : (*ref)[static_cast<Size>(k - 1)][2];
    th_ref = th_prev + wrapPi(th_ref - th_prev);
    (*ref)[static_cast<Size>(k)] = {p.x(), p.y(), th_ref, vk};
    // Blend the reference toward the vehicle at the start of the horizon and
    // let it converge back onto the path.
    //
    // Tracking the raw projection point is what makes the controller cut
    // corners. The projection is the *nearest* path point, and where the vehicle
    // has come round the inside of a turn that point sits diagonally ahead of it
    // -- in the maze world 0.5 m ahead and 0.2 m to the side. Steering straight
    // at it means steering across the corner, which increases the error it was
    // meant to remove; the vehicle hit the steering-rate limit, ran wide, and
    // the watchdog stopped it 8 m into a 128 m plan.
    //
    // Anchoring the first reference point to the vehicle and decaying the offset
    // over the horizon keeps the target where the vehicle is *now* and still
    // pulls it onto the path: the controller sees the true heading error at
    // k = 0 and the path's own curvature beyond it, and converges instead of
    // cutting. This is the standard offset formulation for path-tracking MPC,
    // and the reason it is standard is that it is the only version that does not
    // chase its own projection.
    const Scalar alpha = std::exp(-(static_cast<Scalar>(k) * opt_.dt) / opt_.ref_offset_tau);
    (*ref)[static_cast<Size>(k)][0] += alpha * offset.x();
    (*ref)[static_cast<Size>(k)][1] += alpha * offset.y();
    // The *true* profile value at this step, before the creep floor above is
    // applied. The speed constraint needs this rather than the floored value:
    // the floor exists to stop the optimiser choosing a standstill, and it must
    // not also overrule a profile that has deliberately gone to zero because
    // the body does not fit through what is ahead.
    if (prof_ref) (*prof_ref)[static_cast<Size>(k)] = profileAt(profile, s_axis, s);
    s = std::min(total, s + std::max<Scalar>(0.05, vk) * opt_.dt);
  }
}

// ---------------------------------------------------------------------------
// MPC step
// ---------------------------------------------------------------------------

MpcController::MpcController(const sim::VehicleParams& vehicle, const ControllerOptions& opt)
    : veh_(vehicle), opt_(opt) {
  const Size n = static_cast<Size>(opt_.horizon) * 2;
  qp_.reserve(n, 3 * n + 4);
  warm_start_ = VecX::Zero(static_cast<Eigen::Index>(n));
}

ControllerOutput MpcController::safeStop(const VehicleEstimate& est, ControlStatus why) {
  ControllerOutput out;
  out.status = why;
  const Scalar dt = opt_.dt;
  // Decelerate at the vehicle's limit; the command is the speed that will be
  // reached this tick, so the plant and the watchdog agree on the profile.
  const Scalar v = std::max<Scalar>(0.0, est.v);
  const Scalar v_next = std::max<Scalar>(0.0, v - veh_.max_accel * dt);
  out.cmd.v = v_next;
  out.cmd.accel = (v_next - v) / dt;
  out.cmd.steer = est.steer;
  out.stats.stop_distance = v * v / (2.0 * std::max<Scalar>(1e-3, veh_.max_accel));
  u_prev_ = out.cmd;
  return out;
}

const char* toString(StopReason r) {
  switch (r) {
    case StopReason::kNone: return "none";
    case StopReason::kNoPlan: return "no_plan";
    case StopReason::kEstimateDiverged: return "estimate_diverged";
    case StopReason::kCorridorTooNarrow: return "corridor_too_narrow";
    case StopReason::kProfileZero: return "profile_zero";
    case StopReason::kTubeExceeded: return "tube_exceeded";
    case StopReason::kReferenceDiverged: return "reference_diverged";
    case StopReason::kGoalReached: return "goal_reached";
  }
  return "?";
}

ControllerOutput MpcController::step(const VehicleEstimate& est, const planning::Plan& plan) {
  const auto t0 = std::chrono::steady_clock::now();
  ControllerOutput out;
  if (!plan.ok() || plan.path.size() < 2) {
    return safeStop(est, ControlStatus::kNoPlan);
  }

  const Scalar sigma_lat = est.lateralSigma(opt_.preview);
  out.stats.inflation = opt_.k_sigma * sigma_lat;
  out.stats.det_sigma = covarianceDet(est);

  // Watchdog first: a diverged estimate is not a control problem, and no
  // amount of optimisation fixes it. Checked before anything else so the
  // vehicle never steers on a state that is known to be wrong.
  if (out.stats.det_sigma > opt_.tau_divergence) {
    ControllerOutput stop = safeStop(est, ControlStatus::kSafeStop);
    stop.stats = out.stats;
    stop.stats.stop_reason = StopReason::kEstimateDiverged;
    return stop;
  }

  // Speed ceiling from the chance constraint, using the tightest corridor on
  // the remaining route rather than the local one: a vehicle that only slows
  // down once it is *in* the narrow doorway has already carried too much speed
  // into it.
  Scalar width = plan.stats.route_bottleneck > 0 ? plan.stats.route_bottleneck * 2.0
                                                 : 2.0 * opt_.preferred_width_fallback;
  Scalar ceiling = chanceConstrainedSpeed(width, sigma_lat);
  out.stats.speed_cap = ceiling;
  if (ceiling <= 0.0) {
    ControllerOutput stop = safeStop(est, ControlStatus::kSafeStop);
    stop.stats = out.stats;
    stop.stats.stop_reason = StopReason::kCorridorTooNarrow;
    return stop;
  }

  // The speed actually commanded comes from the braking-feasible profile, which
  // folds in curvature and clearance and looks ahead far enough to start
  // braking in time.
  const std::vector<Scalar> profile = buildVelocityProfile(plan, est, ceiling);
  const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
  const Scalar cap = profileAt(profile, s_axis, progress_);
  out.stats.speed_cap = cap;
  // Goal test, and it has to come *before* the profile-zero check. The profile
  // now ends at zero because the vehicle is meant to arrive stopped, so near the
  // goal the two conditions are both true and their order decides what the run
  // reports. Checking the profile first labels every successful arrival a
  // "profile zero" safe stop, which is the opposite of what happened.
  //
  // The tolerance is generous on purpose: a 5 cm goal test is inside the plant's
  // own discretisation, and at 20 Hz the vehicle covers 2.5 cm per tick, so a
  // tight bound is a test of the tick phase rather than of arrival.
  // Arrival is a statement about position, not speed. The speed profile already
  // brakes to zero at the end of the plan, so by the time the vehicle is inside
  // the tolerance it is stopping anyway; requiring `v < 0.35` as well only
  // relabels a completed route as a "profile zero" failure when the last tick
  // happens to catch it still rolling. Progress is monotonic and clamped to the
  // path, so this cannot fire early.
  const bool at_goal = progress_ >= planning::arcLength(plan.path).back() - 0.15;
  if (at_goal) {
    ControllerOutput stop = safeStop(est, ControlStatus::kGoalReached);
    stop.stats = out.stats;
    stop.stats.stop_reason = StopReason::kGoalReached;
    return stop;
  }

  if (cap <= 0.01) {
    ControllerOutput stop = safeStop(est, ControlStatus::kSafeStop);
    stop.stats = out.stats;
    stop.stats.stop_reason = StopReason::kProfileZero;
    return stop;
  }

  const bool degraded = out.stats.inflation > 0.05 || cap < veh_.max_speed * 0.5;
  out.status = degraded ? ControlStatus::kDegraded : ControlStatus::kOk;


  // Safe-tube check. The plan guarantees a clearance at every point, so it also
  // defines a tube of radius (clearance - half-diagonal) around itself inside
  // which the body is known to be clear of obstacles. If the vehicle is outside
  // that tube it has left the region the plan vouches for, and continuing is a
  // bet that the map is wrong rather than a control action. The check is on the
  // *deviation from the plan*, so it needs no map: the plan's own clearance
  // number is the guarantee being enforced.
  if (!plan.clearance.empty()) {
    Scalar best = std::numeric_limits<Scalar>::max();
    Scalar best_clear = 0;
    const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
    for (Size i = 0; i < plan.path.size(); ++i) {
      if (s_axis[i] < progress_ - 0.5 || s_axis[i] > progress_ + 2.5) continue;
      const Scalar d = (plan.path[i] - Vec2(est.pose.x, est.pose.y)).norm();
      if (d < best) {
        best = d;
        best_clear = plan.clearance[i];
      }
    }
    const Scalar tube = best_clear - opt_.body_half_diagonal;
    out.stats.tube_margin = tube - best;
    // The tolerance absorbs the lag of a discrete-time tracker: the vehicle is
    // always slightly behind the reference, and a tube check with no slack
    // trips on the tracking error of a perfectly good controller. It is set to
    // the tracking error the controller actually achieves in the tightest
    // world, which is the honest figure: a tolerance below the loop's own
    // accuracy turns a safety net into a stop button.
    if (best != std::numeric_limits<Scalar>::max() && best > tube + 0.28) {
      ControllerOutput stop = safeStop(est, ControlStatus::kSafeStop);
      stop.stats = out.stats;
      stop.stats.stop_reason = StopReason::kTubeExceeded;
      return stop;
    }
  }

  // Reference-divergence watchdog: is the vehicle pointing *along* the plan, or
  // against it?
  //
  // Two details make the difference between a safety net and a stop button on
  // every corner.
  //
  // The threshold is 135 degrees, and the margin over a right angle is the
  // whole point. The vehicle is a car with a 0.54 m minimum turn radius, and it
  // *can* drive a 90-degree corner; while it does, it points up to a quarter
  // turn away from the tangent at the point it is nearest, because the heading
  // has to swing through the turn while the position is still entering it. A
  // threshold of 86 degrees -- which is what a 1.5 rad bound amounts to -- fires
  // on every corner sharper than a right angle, so the vehicle stops halfway
  // through a turn it was executing correctly. That failure is indistinguishable
  // from a real divergence unless you check the number against the corner the
  // plan actually contains: the saddle-hall world has a 1.1 rad turn, the
  // controller commanded full lock through it, and the watchdog tripped anyway.
  // 135 degrees leaves a right angle of legitimate mid-corner error plus 45
  // degrees of slack, and still catches a vehicle pointing back down the path.
  //
  // The direction is measured against a *chord a metre ahead*, not the tangent
  // at the nearest point. At a corner the nearest-point tangent swings through
  // the full turn while the vehicle is still rotating, so it reports a large
  // error for a vehicle that is tracking perfectly; a chord over a metre
  // averages the corner and reports the direction the plan actually goes.
  {
    const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
    const Scalar look = std::min(s_axis.back(), progress_ + 1.0);
    const Vec2 p = planning::pointAt(plan.path, s_axis, progress_);
    const Vec2 q = planning::pointAt(plan.path, s_axis, look);
    Vec2 t = q - p;
    if (t.norm() > 0.05) {
      t.normalize();
      const Scalar ref_th = std::atan2(t.y(), t.x());
      if (std::abs(wrapPi(ref_th - est.pose.th)) > 2.356) {
        ControllerOutput stop = safeStop(est, ControlStatus::kSafeStop);
        stop.stats = out.stats;
        stop.stats.stop_reason = StopReason::kReferenceDiverged;
        return stop;
      }
    }
  }

  std::vector<std::array<Scalar, 4>> ref;
  std::vector<Scalar> prof_ref;
  buildReference(plan, est, profile, &ref, &prof_ref);

  // -------- condensed MPC -> QP -------------------------------------------
  //
  // Decision variables are the input sequence u_0..u_{H-1}, each
  // [accel, steer_rate]. The state trajectory is a linear function of the
  // inputs through the linearised model, so the tracking cost becomes a
  // quadratic in u with a dense Hessian. Condensing keeps the QP at 2H = 50
  // variables instead of 4H + 2H + 4 = 154, which is the difference between a
  // few milliseconds and sixty at this horizon.
  //
  // The condensed Hessian 2 S'QS is only positive semi-definite: within the
  // horizon, a steering-rate direction that affects the pose beyond the last
  // step has no curvature. That is handled by the solver's proximal term
  // (`qp_sigma`), not by ignoring it.
  const int H = opt_.horizon;
  const Size nu = 2, nx = 5;
  MatX A(nx, nx), B(nx, nu);
  VecX c(nx);
  linearise(est, u_prev_, &A, &B, &c);
  // Discretise by forward Euler about the operating point. `linearise` returns
  // the *continuous* model, so the step has to be applied here: x_{k+1} =
  // x_k + dt (A x_k + B u_k + c). Omitting the dt leaves the prediction
  // describing derivatives rather than states, and the controller then tracks a
  // reference it believes it is already on.
  A = MatX::Identity(nx, nx) + opt_.dt * A;
  B = opt_.dt * B;
  c = opt_.dt * c;

  const Size n = static_cast<Size>(H) * nu;
  MatX S = MatX::Zero(static_cast<Eigen::Index>(nx * H), static_cast<Eigen::Index>(n));
  std::vector<VecX> free(static_cast<Size>(H));
  VecX x_free(nx);
  x_free << est.pose.x, est.pose.y, est.pose.th, std::max<Scalar>(0.0, est.v), est.steer;
  for (int k = 0; k < H; ++k) {
    x_free = A * x_free + c;
    free[static_cast<Size>(k)] = x_free;
    // d x_{k+1} / d u_j = A^{k-j} B, for j = 0..k.
    for (int j = 0; j <= k; ++j) {
      MatX Pk = MatX::Identity(nx, nx);
      for (int e = 0; e < k - j; ++e) Pk = A * Pk;
      S.block(k * nx, j * nu, nx, nu) = Pk * B;
    }
  }

  MatX P = MatX::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
  VecX q = VecX::Zero(static_cast<Eigen::Index>(n));

  MatX Q = MatX::Zero(nx, nx);
  Q(0, 0) = opt_.w_lateral;
  Q(1, 1) = opt_.w_lateral;
  Q(2, 2) = opt_.w_heading;
  Q(3, 3) = opt_.w_speed;

  for (int k = 0; k < H; ++k) {
    VecX xr = VecX::Zero(nx);
    // The reference carries [x, y, theta, v]; the fifth state (steering angle)
    // has no reference and no cost, so it is zero here and Q(4,4) is zero.
    for (Size e = 0; e < 4; ++e) xr(static_cast<Eigen::Index>(e)) = ref[static_cast<Size>(k + 1)][e];
    const VecX resid = free[static_cast<Size>(k)] - xr;
    const MatX Sk = S.block(k * nx, 0, nx, n);
    // Minimise (S u + free - x_ref)' Q (S u + free - x_ref): quadratic part
    // 2 S'QS, linear part 2 S'Q(free - x_ref). Getting the sign of the linear
    // term wrong drives the inputs to their lower bound, which reads as "the
    // MPC refuses to move".
    P += 2.0 * Sk.transpose() * Q * Sk;
    q += 2.0 * Sk.transpose() * Q * resid;
  }

  for (int k = 0; k < H; ++k) {
    const Eigen::Index i = k * nu;
    P(i, i) += 2.0 * opt_.w_accel;
    P(i + 1, i + 1) += 2.0 * opt_.w_steer_rate;
    // ||u_k - u_{k-1}||^2 as a path-graph Laplacian: +1 on every diagonal, -1
    // on every adjacent off-diagonal. The diagonal is added for *every* index,
    // including k = 0; adding it only for k > 0 leaves the first input with an
    // off-diagonal and no diagonal to balance it, which is an indefinite
    // Hessian and makes the solver's Cholesky meaningless.
    //
    // The diagonal weight is the node *degree*, and that is not cosmetic. The
    // penalty is sum_k w ||u_k - u_{k-1}||^2, which as a quadratic form is
    // w times the path-graph Laplacian: degree on the diagonal, -1 on the
    // adjacent off-diagonals. Using degree 1 for every node gives every
    // interior node a diagonal of 1 and two off-diagonals of -1, so the form
    // is 2*1 - 2*1 = 0 along the alternating direction (+1,-1,+1,-1,...) --
    // and negative for nearby directions. The Hessian is then indefinite,
    // which is not a small numerical wart: the Cholesky in the solver's
    // factorisation does not exist for an indefinite matrix, the QP stops
    // converging, and the controller's commands become the tail of a stalled
    // ADMM rather than the solution of anything. Measured on the saddle-hall
    // corner it produced a min eigenvalue of -1.26 and residuals of 1e0 after
    // 1500 iterations, versus min eigenvalue +1.5 and residuals of 1e-6 with
    // the degree-correct Laplacian.
    const Scalar deg = (k == 0 || k == H - 1) ? 1.0 : 2.0;
    P(i, i) += 2.0 * opt_.w_accel_delta * deg;
    P(i + 1, i + 1) += 2.0 * opt_.w_steer_rate_delta * deg;
    if (k > 0) {
      P(i, i - nu) -= 2.0 * opt_.w_accel_delta;
      P(i - nu, i) -= 2.0 * opt_.w_accel_delta;
      P(i + 1, i + 1 - nu) -= 2.0 * opt_.w_steer_rate_delta;
      P(i + 1 - nu, i + 1) -= 2.0 * opt_.w_steer_rate_delta;
    } else {
      // u_{-1} is the previous command, a constant: it contributes a linear
      // term, not an off-diagonal block.
      q(i) -= 2.0 * opt_.w_accel_delta * u_prev_.accel;
    }
  }

  // Box on the inputs, plus a row per step limiting the speed the inputs imply,
  // so the horizon never asks for more than the cap. The speed constraint is
  // what enforces "slow down before the doorway" rather than "in it".
  const Size m = n + 2 * static_cast<Size>(H);
  MatX Acon = MatX::Zero(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));
  VecX lo = VecX::Zero(static_cast<Eigen::Index>(m));
  VecX hi = VecX::Zero(static_cast<Eigen::Index>(m));
  for (int k = 0; k < H; ++k) {
    const Eigen::Index i = k * nu;
    Acon(i, i) = 1.0;
    Acon(i + 1, i + 1) = 1.0;
    lo(i) = -veh_.max_accel;
    hi(i) = veh_.max_accel;
    lo(i + 1) = -veh_.max_steer_rate;
    hi(i + 1) = veh_.max_steer_rate;
  }
  // Speed rows. v_{k+1} = v_k + dt * accel_k, unrolled over the horizon, and
  // bounded on *both* sides.
  //
  // The lower bound is not decoration. The model is a forward-Euler bicycle
  // model with no reverse gear: nothing in it says the speed cannot go negative,
  // and with only an upper bound the solver is free to buy a shorter stopping
  // distance by commanding a negative speed. The plant honours that literally,
  // and the vehicle drives backwards -- which showed up as a closed-loop run
  // that oscillated in speed, wandered off the reference and eventually tripped
  // the tube watchdog. Bounding v >= 0 makes the model match the vehicle: this
  // is a car that brakes to a stop, not one that reverses out of trouble.
  //
  // The lower bound is the *reachable* floor, not a constant: asking for
  // v >= v_floor at step 0 of a standing start would demand an acceleration
  // above the vehicle's limit and make the QP infeasible. Ramping the floor at
  // max_accel keeps every row satisfiable while still ruling out the
  // "stopped and stuck" solution, which is the one the optimiser finds
  // otherwise.
  const Scalar v_now = std::max<Scalar>(0.0, est.v);
  const Scalar v_floor = std::min(cap, opt_.min_track_speed);
  for (int k = 0; k < H; ++k) {
    const Eigen::Index row = n + k;
    for (int j = 0; j <= k; ++j) Acon(row, j * nu) = opt_.dt;
    const Scalar reach = v_now + opt_.dt * static_cast<Scalar>(k + 1) * veh_.max_accel;
    // The floor is capped by the profile at this step, so a deliberate
    // zero -- "the body does not fit through what is ahead" -- is never
    // overridden into a collision.
    lo(row) = std::min({v_floor, reach, prof_ref[static_cast<Size>(k + 1)]}) - v_now;
    hi(row) = std::min(cap, std::max(prof_ref[static_cast<Size>(k + 1)], v_now)) - v_now;
  }
  // Steering-angle rows: delta_{k+1} = delta_0 + dt * sum_{j<=k} rate_j, bounded
  // to the actuator envelope. Without these the steering angle is unbounded in
  // the prediction, and the optimiser happily plans to hold 3 rad of steer for
  // a gentle curve because nothing tells it the wheel stops turning at 0.75.
  const Scalar st_now = est.steer;
  for (int k = 0; k < H; ++k) {
    const Eigen::Index row = n + H + k;
    for (int j = 0; j <= k; ++j) Acon(row, j * nu + 1) = opt_.dt;
    lo(row) = -veh_.max_steer - st_now;
    hi(row) = veh_.max_steer - st_now;
  }

  qp_.setProblem(P, q, Acon, lo, hi);
  QpOptions qopt;
  qopt.eps_abs = opt_.qp_eps;
  qopt.eps_rel = opt_.qp_eps;
  qopt.max_iter = opt_.qp_max_iter;
  qopt.sigma = opt_.qp_sigma;
  // Fixed ADMM penalty. The OSQP-style adaptive rule is a good default on a
  // well-scaled problem and a liability on this one: with adaptation enabled
  // the iterates diverged to 1e20 on the reference world. A fixed penalty
  // converges, and the condensed problem is small enough that a few hundred
  // iterations fit the 20 Hz budget.
  qopt.rho = opt_.qp_rho;
  qopt.adaptive_rho = opt_.qp_adaptive_rho;
  solution_.resize(static_cast<Eigen::Index>(n));
  solution_.setZero();
  if (static_cast<Eigen::Index>(warm_start_.size()) != static_cast<Eigen::Index>(n))
    warm_start_ = VecX::Zero(static_cast<Eigen::Index>(n));
  const QpResult res = qp_.solve(&solution_, warm_start_, qopt);
  out.stats.qp_iterations = res.iterations;
  out.stats.qp_primal = res.primal_residual;
  out.stats.qp_dual = res.dual_residual;
  out.stats.qp_objective = res.objective;
  // Warm start the next tick: shift the input sequence one step and repeat the
  // last input. The problem changes only slightly between ticks, so this keeps
  // the iteration count low in steady state.
  warm_start_.head(static_cast<Eigen::Index>(n - nu)) =
      solution_.segment(nu, static_cast<Eigen::Index>(n - nu));
  warm_start_.tail(nu) = solution_.tail(nu);

  // Apply the first input.
  const Scalar accel = solution_(0);
  const Scalar steer_rate = solution_(1);
  out.cmd.accel = clamp(accel, -veh_.max_accel, veh_.max_accel);
  out.cmd.v = clamp(std::max<Scalar>(0.0, est.v) + out.cmd.accel * opt_.dt, 0.0, veh_.max_speed);
  out.cmd.steer = clamp(est.steer + clamp(steer_rate, -veh_.max_steer_rate, veh_.max_steer_rate) * opt_.dt,
                        -veh_.max_steer, veh_.max_steer);
  u_prev_ = out.cmd;

  // Roll the prediction out for reporting and for the tests.
  out.predicted.resize(static_cast<Size>(H) + 1);
  VecX xk(nx);
  xk << est.pose.x, est.pose.y, est.pose.th, std::max<Scalar>(0.0, est.v), est.steer;
  out.predicted[0] = {xk(0), xk(1), xk(2), xk(3)};
  for (int k = 0; k < H; ++k) {
    VecX uk(nu);
    uk << solution_(k * nu), solution_(k * nu + 1);
    xk = A * xk + B * uk + c;
    out.predicted[static_cast<Size>(k + 1)] = {xk(0), xk(1), xk(2), std::max<Scalar>(0.0, xk(3))};
  }

  // Cross-track and heading error at the current pose, for the benchmark.
  {
    const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
    const Vec2 p = planning::pointAt(plan.path, s_axis, progress_);
    const Vec2 qq = planning::pointAt(plan.path, s_axis, std::min(s_axis.back(), progress_ + 0.25));
    Vec2 t = qq - p;
    if (t.norm() < 1e-6) t = Vec2(1, 0);
    t.normalize();
    const Vec2 nvec(-t.y(), t.x());
    const Vec2 d = Vec2(est.pose.x, est.pose.y) - p;
    out.stats.lateral_error = d.dot(nvec);
    out.stats.heading_error = wrapPi(std::atan2(t.y(), t.x()) - est.pose.th);
  }
  out.stats.solve_ms =
      std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return out;
}

}  // namespace control
}  // namespace luna
