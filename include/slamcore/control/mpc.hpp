// SLAMCore-8G -- bicycle-model MPC, chance constraints, safety watchdogs.
//
// The controller closes the loop around the plan. Three things make it more
// than a path follower, and all three are consequences of running on a SLAM
// estimate rather than on ground truth:
//
//   1. Chance constraints. The pose the controller steers by is an estimate
//      with a covariance. A corridor that is 1.2 m wide is wide only relative
//      to a *known* position; at 3 sigma of 0.4 m of lateral uncertainty it is
//      not wide at all. The controller therefore inflates the corridor by
//      `k_sigma * sigma_lateral` and, when the inflated corridor no longer
//      fits, slows down -- the vehicle needs less room when it stops sooner.
//
//   2. Graceful degradation. As covariance grows the speed cap tightens
//      continuously; the vehicle does not fall off a cliff from "driving" to
//      "stopped". That is what lets a fault be survivable rather than fatal.
//
//   3. A divergence watchdog. `sqrt(det Sigma)` above a threshold means the
//      estimate is no longer trustworthy, and no amount of control fixes a
//      broken state estimate. The watchdog commands a safe stop instead of
//      steering confidently in a direction that may be wrong.
//
// The optimisation itself is a linear time-varying MPC: the bicycle model is
// linearised about the current state each tick and condensed into a dense QP
// over the input sequence, solved by the ADMM solver in core/qp.hpp. Horizon
// 25 at 20 Hz on four CPU cores.
#pragma once

#include <array>
#include <vector>

#include "slamcore/core/qp.hpp"
#include "slamcore/core/se3.hpp"
#include "slamcore/core/types.hpp"
#include "slamcore/planning/planner.hpp"
#include "slamcore/sim/vehicle.hpp"

namespace sc {
namespace control {

/// The controller's view of where it is and how sure it is. Decoupled from the
/// SLAM types on purpose: the controller is testable with a hand-written
/// covariance, and the estimator is testable without the controller.
struct VehicleEstimate {
  Pose2 pose;
  Scalar v = 0;        ///< forward speed, m/s
  Scalar steer = 0;    ///< current steering angle, rad
  /// 3x3 position/heading covariance in the local frame (x forward, y left,
  /// theta). Only the diagonal and the y-theta block are used.
  Scalar sigma_xx = 0, sigma_yy = 0, sigma_thth = 0, sigma_yth = 0;
  /// Worst-case lateral standard deviation, metres: what the chance constraint
  /// budgets for. sqrt(sigma_yy) alone understates the error, because a heading
  /// error also displaces the vehicle laterally over the preview distance.
  Scalar lateralSigma(Scalar preview_m) const {
    return std::sqrt(std::max<Scalar>(0.0, sigma_yy)) +
           std::abs(preview_m) * std::sqrt(std::max<Scalar>(0.0, sigma_thth));
  }
};

struct ControllerOptions {
  Scalar dt = 0.05;          ///< 20 Hz
  int horizon = 25;
  /// Tracking weights: cross-track, heading, speed, then input effort.
  Scalar w_lateral = 20.0;
  Scalar w_heading = 25.0;
  Scalar w_speed = 4.0;
  Scalar w_accel = 0.6;
  Scalar w_steer_rate = 1.0;
  /// Input regularisation: penalise changes between consecutive inputs so the
  /// commands are smooth. A jerk to the steering is a jerk to the passengers
  /// and to the odometry.
  Scalar w_accel_delta = 1.5;
  Scalar w_steer_rate_delta = 1.5;
  /// Preview distance for reference sampling, metres.
  /// Lowest speed the optimiser is allowed to plan, m/s, when the profile
  /// permits motion at all.
  ///
  /// This exists to break a deadlock that is intrinsic to linear MPC on a
  /// nonholonomic vehicle. The heading response to a steering input is
  /// proportional to speed (`theta' = v tan(delta) / L`), so at v = 0 the
  /// linearised model says steering does nothing, and a vehicle that has
  /// stopped with a heading error has no predicted way to fix it. Stopping is
  /// then a genuine local optimum of the QP -- not a solver failure; the
  /// residuals converge cleanly -- and the vehicle sits still forever, driving
  /// 12 m of a 43 m plan. Keeping the speed off zero keeps the linearisation
  /// valid and the vehicle always able to turn its way out.
  Scalar min_track_speed = 0.2;
  /// Time constant, seconds, over which the reference rejoins the path from the
  /// vehicle's current offset. Short enough to correct promptly, long enough
  /// that the correction is a turn rather than a step.
  Scalar ref_offset_tau = 0.4;
  Scalar preview = 2.0;
  /// Chance-constraint confidence: the corridor is inflated by this many
  /// standard deviations of lateral uncertainty.
  Scalar k_sigma = 2.0;
  /// Lateral acceleration the controller is willing to use, m/s^2. Below the
  /// tyre limit on purpose: this is a planning figure for a low-grip estimate,
  /// not a performance target.
  Scalar a_lat_max = 1.0;
  /// Braking deceleration assumed when building the speed profile, m/s^2.
  Scalar a_brake = 1.2;
  /// Watchdog: safe-stop when sqrt(det Sigma) exceeds this (m^2 * rad).
  Scalar tau_divergence = 0.35;
  /// Hard cap on speed while degraded, m/s.
  Scalar degraded_speed = 0.35;
  /// Residual half-width the vehicle insists on keeping after the body and the
  /// uncertainty have taken their share, metres.
  ///
  /// This is a *margin on what is left*, not a width, and the distinction is
  /// the whole point of the number. The usable half-width is
  ///
  ///     half = W/2 - body_half_width - k_sigma * sigma_lat
  ///
  /// and the vehicle drives only if `half` exceeds this margin. Setting the
  /// margin to a body-scale figure -- a third of a metre, say -- silently
  /// demands a corridor of `2 * (body + k sigma + 0.35)` metres, which on a
  /// 1.8 m doorway with a well-estimated pose is more than the doorway has, so
  /// the vehicle refuses to move at all. The failure is quiet: the plan is
  /// valid, the geometry fits, and the controller stops with no indication that
  /// its own margin, not the world, was the obstacle.
  ///
  /// Zero is the right default: the chance constraint has already reserved
  /// `k_sigma` standard deviations, so a strictly positive residual means the
  /// gap is wider than the body plus the uncertainty budget, which is exactly
  /// the condition for proceeding. Asking for more than that is asking for a
  /// margin on a margin. The knob exists for callers who want to be
  /// deliberately more timid than their own covariance justifies.
  Scalar min_usable_width = 0.0;
  /// Corridor width assumed when the plan does not report a bottleneck. A plan
  /// that crossed no doorway has no bottleneck to report, and treating that as
  /// "zero width" would stop the vehicle in an open room.
  Scalar preferred_width_fallback = 1.5;
  /// Body half-diagonal, metres: the distance from the centre to the furthest
  /// corner. Clearance is spent against this, not against the half-width,
  /// because the corner leads when the vehicle is not aligned with the wall.
  /// The default matches the simulator's 0.34 x 0.62 m body.
  Scalar body_half_diagonal = 0.353;
  Scalar qp_eps = 1e-5;
  /// Proximal regularisation handed to the QP solver. See the comment at the
  /// solve site for why a condensed MPC needs it.
  Scalar qp_sigma = 1e-3;
  /// ADMM penalty. Fixed, not adapted: the adaptive rule diverges on this
  /// problem (measured residual 1e20 with adaptation on, 1e-6 with a fixed
  /// rho), and a fixed penalty at this magnitude converges in well under a
  /// hundred iterations. See the comment at the solve site.
  Scalar qp_rho = 10.0;
  bool qp_adaptive_rho = false;
  Size qp_max_iter = 1500;
};

enum class ControlStatus {
  kOk,          ///< tracking normally
  kDegraded,    ///< covariance inflated the corridor; slowed down
  kSafeStop,    ///< watchdog tripped or no usable corridor: stopping
  kGoalReached,
  kNoPlan,
};

const char* toString(ControlStatus s);

/// Which safety check produced a stop. A safe stop has several distinct causes
/// -- a diverged estimate, a corridor too narrow for the body, a reference the
/// vehicle cannot follow -- and they call for different fixes. Reporting only
/// "safe stop" makes a controller that is working exactly as designed look
/// identical to one that is broken.
enum class StopReason {
  kNone,
  kNoPlan,
  kEstimateDiverged,
  kCorridorTooNarrow,
  kProfileZero,
  kTubeExceeded,
  kReferenceDiverged,
  kGoalReached,
};

const char* toString(StopReason r);

struct ControllerStats {
  Size qp_iterations = 0;
  Scalar qp_primal = 0;
  Scalar qp_dual = 0;
  Scalar qp_objective = 0;
  StopReason stop_reason = StopReason::kNone;
  Scalar solve_ms = 0;
  Scalar lateral_error = 0;      ///< cross-track error at the current step
  Scalar heading_error = 0;
  Scalar speed_cap = 0;          ///< the cap actually applied
  Scalar inflation = 0;          ///< corridor inflation from covariance, metres
  Scalar det_sigma = 0;
  Scalar stop_distance = 0;
  Scalar clearance_margin = 0;   ///< plan clearance minus the body half-diagonal
  Scalar tube_margin = 0;        ///< safe-tube radius minus the deviation from the plan
};

struct ControllerOutput {
  sim::Control cmd;
  ControlStatus status = ControlStatus::kNoPlan;
  /// Predicted state trajectory over the horizon: [x, y, theta, v] per step.
  std::vector<std::array<Scalar, 4>> predicted;
  ControllerStats stats;
  bool stopping() const { return status == ControlStatus::kSafeStop; }
};

/// Linear time-varying MPC over a planned path.
class MpcController {
 public:
  MpcController(const sim::VehicleParams& vehicle, const ControllerOptions& opt = ControllerOptions());

  /// Compute the command for this tick. `plan` may be replaced by the caller
  /// between ticks; the controller holds no state about it beyond the current
  /// progress along it.
  ControllerOutput step(const VehicleEstimate& est, const planning::Plan& plan);

  /// Safe stop with no plan: decelerate along the current heading.
  ControllerOutput safeStop(const VehicleEstimate& est, ControlStatus why);

  /// Progress along the plan: arc length of the closest point, metres.
  Scalar progress() const { return progress_; }
  void resetProgress() { progress_ = 0; }
  const ControllerOptions& options() const { return opt_; }

  /// Speed cap implied by the chance constraint for a given corridor width and
  /// lateral uncertainty. Exposed so the benchmark can report it and the tests
  /// can assert the monotonicity that the whole scheme depends on.
  Scalar chanceConstrainedSpeed(Scalar corridor_width, Scalar sigma_lat) const;

  /// sqrt(det Sigma) from the estimate's covariance.
  static Scalar covarianceDet(const VehicleEstimate& est);

  /// Per-point speed limit along the plan: the tightest of the chance
  /// constraint, the curvature limit, and the clearance limit, then made
  /// feasible by a backwards pass that respects the braking limit. Point i is
  /// reachable at profile[i] and stoppable before any point j > i.
  ///
  /// The backwards pass is what makes the difference between "slows down in
  /// the corner" and "slows down before it". A forward-only limit is only
  /// correct if the vehicle is already at or below the limit when it arrives;
  /// a vehicle that is above it has to brake, and the brake has to start
  /// earlier. This is the standard speed-profile construction, and it is much
  /// easier to get right than tuning a lookahead distance.
  std::vector<Scalar> buildVelocityProfile(const planning::Plan& plan, const VehicleEstimate& est,
                                           Scalar ceiling) const;

  /// Speed limit at arc length `s`, linearly interpolated between profile
  /// points.
  static Scalar profileAt(const std::vector<Scalar>& profile, const std::vector<Scalar>& s_axis,
                          Scalar s);


 private:
  void linearise(const VehicleEstimate& est, const sim::Control& u_prev, MatX* A, MatX* B,
                 VecX* c) const;
  /// Reference trajectory over the horizon: position, heading, speed.
  /// Non-const: it advances `progress_` along the plan as a side effect, which
  /// is the controller's only piece of cross-tick state.
  void buildReference(const planning::Plan& plan, const VehicleEstimate& est,
                      const std::vector<Scalar>& profile, std::vector<std::array<Scalar, 4>>* ref,
                      std::vector<Scalar>* prof_ref = nullptr);

  sim::VehicleParams veh_;
  ControllerOptions opt_;
  QpSolver qp_;
  VecX solution_;
  VecX warm_start_;
  sim::Control u_prev_;
  Scalar progress_ = 0;
};

}  // namespace control
}  // namespace sc
