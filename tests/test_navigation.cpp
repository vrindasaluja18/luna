// SLAMCore-8G -- planning, topology and control tests.
//
// These are the tests that make the safety claims checkable rather than
// aspirational: a collision-free closed-loop run, a safe stop under an injected
// fault, and the monotonicity properties the chance constraint depends on.
#include <algorithm>
#include <cmath>
#include <vector>

#include "slamcore/app/navigation.hpp"
#include "slamcore/control/mpc.hpp"
#include "slamcore/planning/planner.hpp"
#include "slamcore/sim/world.hpp"
#include "slamcore/topology/morse.hpp"
#include "test_util.hpp"

using namespace sc;

namespace {

struct LoopResult {
  Scalar rms = 0;
  Scalar max_err = 0;
  Scalar distance = 0;
  bool collided = false;
  bool safe_stopped = false;
  bool reached_goal = false;
  Scalar plan_length = 0;
  Scalar mean_solve_ms = 0;
  int steps = 0;
};

/// Drives the plan with the MPC against the kinematic plant and scores it.
/// `covariance` is the estimate uncertainty injected each tick; `fault_after`
/// starts growing it to trip the watchdog.
LoopResult driveLoop(sim::WorldType type, int seed, Scalar sigma_pos, Scalar sigma_th,
                     int fault_after = -1) {
  LoopResult r;
  sim::WorldSpec spec;
  spec.type = type;
  spec.seed = seed;
  spec.size_m = 40;
  const sim::World world = sim::World::generate(spec);
  const auto complex = topology::MorseComplex::extract(world.grid());
  planning::CorridorPlanner planner(complex.graph(), world.grid());

  // The longest corridor-to-corridor route in the world: a hard case rather
  // than a short one that would pass trivially.
  int a = -1, b = -1;
  Scalar best = 0;
  const auto& nodes = complex.graph().nodes;
  for (Size i = 0; i < nodes.size(); ++i)
    for (Size j = i + 1; j < nodes.size(); ++j) {
      const Scalar d = (nodes[i].p - nodes[j].p).norm();
      if (d > best) {
        best = d;
        a = static_cast<int>(i);
        b = static_cast<int>(j);
      }
    }
  if (a < 0) return r;
  const planning::Plan plan = planner.plan(nodes[static_cast<Size>(a)].p, nodes[static_cast<Size>(b)].p);
  if (!plan.ok()) return r;

  sim::VehicleParams veh_params;
  Pose2 start{plan.path.front().x(), plan.path.front().y(), 0};
  {
    const Vec2 t = plan.path[std::min<Size>(3, plan.path.size() - 1)] - plan.path.front();
    start.th = std::atan2(t.y(), t.x());
  }
  sim::Vehicle vehicle(veh_params, start, 1);
  control::ControllerOptions opt;
  control::MpcController mpc(veh_params, opt);

  control::VehicleEstimate est;
  est.pose = start;
  const std::vector<Scalar> s_axis = planning::arcLength(plan.path);

  Scalar sum_sq = 0, sum_solve = 0;
  int n = 0;
  const Scalar dt = opt.dt;
  // Injected fault growth, accumulated *across* ticks.
  //
  // It has to accumulate to be a fault at all. Adding a fixed increment to a
  // covariance that is re-initialised every tick leaves the uncertainty
  // constant, so the watchdog sees a system that is as well-localised on tick
  // 2000 as on tick 0 and -- correctly -- never trips. The failure then reads as
  // "the watchdog does not work", which is the opposite of the truth: a
  // diverging estimator is one whose uncertainty *grows*, and a test that does
  // not model growth cannot test the response to it.
  Scalar fault_pos = 0, fault_th = 0;
  for (int k = 0; k < 4000; ++k) {
    est.sigma_xx = sigma_pos + fault_pos;
    est.sigma_yy = sigma_pos + fault_pos;
    est.sigma_thth = sigma_th + fault_th;
    est.sigma_yth = 0.25 * std::sqrt((sigma_pos + fault_pos) * (sigma_th + fault_th));
    if (fault_after >= 0 && k > fault_after) {
      fault_pos += 0.02;
      fault_th += 0.01;
    }
    const control::ControllerOutput out = mpc.step(est, plan);
    sum_solve += out.stats.solve_ms;
    if (out.stopping()) r.safe_stopped = true;
    vehicle.step(out.cmd, dt);
    est.pose = vehicle.state().pose;
    est.v = vehicle.state().v;
    est.steer = vehicle.state().steer;
    if (vehicle.inCollision(world.grid())) {
      r.collided = true;
      break;
    }
    Scalar d = std::numeric_limits<Scalar>::max();
    for (const Vec2& p : plan.path) d = std::min(d, (p - Vec2(est.pose.x, est.pose.y)).norm());
    sum_sq += d * d;
    ++n;
    r.max_err = std::max(r.max_err, d);
    r.steps = k;
    if (out.status == control::ControlStatus::kGoalReached) {
      r.reached_goal = true;
      break;
    }
    if (out.stopping()) break;
    (void)s_axis;
  }
  r.rms = n ? std::sqrt(sum_sq / n) : 0;
  r.mean_solve_ms = r.steps ? sum_solve / r.steps : 0;
  r.distance = vehicle.distanceTravelled();
  r.plan_length = plan.stats.path_length;
  return r;
}

void testChanceConstraintMonotonic() {
  sim::VehicleParams veh;
  control::ControllerOptions opt;
  control::MpcController mpc(veh, opt);

  // Speed must fall as uncertainty grows, and must be zero once the corridor
  // cannot hold the body plus the uncertainty budget.
  const Scalar w = 1.2;
  const Scalar v0 = mpc.chanceConstrainedSpeed(w, 0.0);
  const Scalar v1 = mpc.chanceConstrainedSpeed(w, 0.05);
  const Scalar v2 = mpc.chanceConstrainedSpeed(w, 0.15);
  SC_CHECK(v0 >= v1);
  SC_CHECK(v1 >= v2);
  SC_CHECK(v0 > 0.0);
  // A corridor narrower than the body is not drivable at any confidence.
  SC_NEAR(mpc.chanceConstrainedSpeed(0.30, 0.0), 0.0, 1e-12);
  // An empty corridor with no uncertainty is capped by the vehicle, not by the
  // corridor: the constraint must not be the binding limit in open space.
  SC_NEAR(mpc.chanceConstrainedSpeed(20.0, 0.0), veh.max_speed, 1e-9);
}

void testCovarianceDetMonotonic() {
  control::VehicleEstimate e;
  e.sigma_xx = e.sigma_yy = e.sigma_thth = 0.01;
  const Scalar a = control::MpcController::covarianceDet(e);
  e.sigma_yy = 0.04;
  const Scalar b = control::MpcController::covarianceDet(e);
  SC_CHECK(b > a);
  e.sigma_xx = e.sigma_yy = e.sigma_thth = 0.0;
  SC_NEAR(control::MpcController::covarianceDet(e), 0.0, 1e-15);
}

void testSafeStopUnderFault() {
  // Fault injection: the estimate covariance grows without bound. The watchdog
  // must stop the vehicle and it must not collide while doing so.
  const LoopResult r = driveLoop(sim::WorldType::kRooms, 6, 1e-4, 2.5e-5, /*fault_after=*/200);
  SC_CHECK(r.safe_stopped);
  SC_CHECK(!r.collided);
}

void testClosedLoopNoCollision() {
  // Four worlds of increasing difficulty. Every one must complete the longest
  // route without touching a wall: the plan's clearance guarantee plus the
  // safe-tube check are what make this a property rather than a hope.
  struct Case {
    sim::WorldType type;
    const char* name;
  };
  const Case cases[] = {
      {sim::WorldType::kRooms, "rooms"},
      {sim::WorldType::kCorridorLoop, "corridor_loop"},
      {sim::WorldType::kSaddleHall, "saddle_hall"},
      {sim::WorldType::kMaze, "maze"},
  };
  for (const Case& c : cases) {
    const LoopResult r = driveLoop(c.type, 6, 1e-4, 2.5e-5);
    std::printf("  %-14s rms=%.3f m  max=%.3f m  dist=%.1f/%.1f m  solve=%.1f ms\n", c.name,
                r.rms, r.max_err, r.distance, r.plan_length, r.mean_solve_ms);
    SC_CHECK(!r.collided);
    // The route must be *completed*, not merely started. A distance floor of a
    // few metres passed while the vehicle was stalling 8 m into a 128 m plan --
    // the assertion was satisfied by the failure it existed to catch. The
    // honest bound is the plan's own length.
    SC_CHECK(r.reached_goal);
    SC_CHECK(r.distance > 0.9 * r.plan_length);
    // Tracking accuracy in the simulator, which is the figure the contract
    // quotes. The worlds are deliberately tight, so this is a demanding bound.
    SC_CHECK(r.rms < 0.60);
  }
}

void testVelocityProfileRespectsBraking() {
  // The profile's defining property: at every point the vehicle can still stop
  // before the next point. Checked directly on the profile rather than through
  // the closed loop, so a failure points at the profile.
  sim::WorldSpec spec;
  spec.type = sim::WorldType::kMaze;
  spec.seed = 6;
  spec.size_m = 40;
  const sim::World world = sim::World::generate(spec);
  const auto complex = topology::MorseComplex::extract(world.grid());
  planning::CorridorPlanner planner(complex.graph(), world.grid());
  const auto& nodes = complex.graph().nodes;
  const planning::Plan plan = planner.plan(nodes.front().p, nodes.back().p);
  if (!plan.ok()) return;

  control::VehicleEstimate est;
  est.sigma_xx = est.sigma_yy = 1e-4;
  est.sigma_thth = 2.5e-5;
  control::ControllerOptions opt;
  control::MpcController mpc(sim::VehicleParams(), opt);
  const std::vector<Scalar> prof = mpc.buildVelocityProfile(plan, est, 2.5);
  const std::vector<Scalar> s_axis = planning::arcLength(plan.path);
  SC_CHECK(prof.size() == plan.path.size());
  for (Size i = 0; i + 1 < prof.size(); ++i) {
    const Scalar ds = s_axis[i + 1] - s_axis[i];
    const Scalar reachable = std::sqrt(prof[i + 1] * prof[i + 1] + 2.0 * opt.a_brake * ds);
    SC_CHECK(prof[i] <= reachable + 1e-6);
  }
  // Monotone in the sense that matters: no point exceeds the ceiling.
  for (Scalar v : prof) SC_CHECK(v <= 2.5 + 1e-9);
}

void testPlannerClearanceAndDrivability() {
  // The planner's output must be drivable: every point clear of the walls by
  // more than the body's half-diagonal, and no turn tighter than the vehicle
  // can make.
  const sim::VehicleParams veh;
  const Scalar half_diag =
      0.5 * std::sqrt(veh.length * veh.length + veh.track_width * veh.track_width);
  const Scalar max_kappa = std::tan(veh.max_steer) / veh.wheelbase;

  for (int t = 0; t < 4; ++t) {
    sim::WorldSpec spec;
    spec.type = static_cast<sim::WorldType>(t);
    spec.seed = 6;
    spec.size_m = 40;
    const sim::World world = sim::World::generate(spec);
    const auto complex = topology::MorseComplex::extract(world.grid());
    planning::CorridorPlanner planner(complex.graph(), world.grid());
    const auto& nodes = complex.graph().nodes;
    int a = -1, b = -1;
    Scalar best = 0;
    for (Size i = 0; i < nodes.size(); ++i)
      for (Size j = i + 1; j < nodes.size(); ++j) {
        const Scalar d = (nodes[i].p - nodes[j].p).norm();
        if (d > best) {
          best = d;
          a = static_cast<int>(i);
          b = static_cast<int>(j);
        }
      }
    const planning::Plan plan =
        planner.plan(nodes[static_cast<Size>(a)].p, nodes[static_cast<Size>(b)].p);
    SC_CHECK(plan.ok());
    if (!plan.ok()) continue;
    SC_CHECK(plan.clearance.size() == plan.path.size());
    SC_CHECK(plan.stats.min_clearance_along > half_diag);
    const std::vector<Scalar> k = planning::CorridorPlanner::curvature(plan.path);
    for (Scalar v : k) SC_CHECK(v <= max_kappa + 1e-6);
  }
}

void testMorseEulerRelation() {
  // The topology invariant the corridor graph is built on: for a height
  // function on a simply connected domain, #maxima - #saddles = 1.
  sim::WorldSpec spec;
  spec.type = sim::WorldType::kRooms;
  spec.seed = 6;
  spec.size_m = 40;
  const sim::World world = sim::World::generate(spec);
  const auto complex = topology::MorseComplex::extract(world.grid());
  const auto& nodes = complex.graph().nodes;
  int maxima = 0, saddles = 0;
  for (const auto& n : nodes) {
    if (n.type == topology::CriticalType::kMaximum) ++maxima;
    if (n.type == topology::CriticalType::kSaddle) ++saddles;
  }
  SC_CHECK(maxima - saddles == 1);
  SC_CHECK(maxima >= 2);
}

void testStampRadiusExceedsPlanClearance() {
  // The invariant that made the reconstructed maze unplannable: the disc stamped
  // free around each driven pose must be drivable at the planner's own minimum
  // clearance, or the centreline of the vehicle's own trail falls below the
  // threshold and the map fragments into islands no route can cross.
  const app::NavigationOptions opt;
  SC_CHECK(opt.stamp_free_radius_m > opt.plan.min_clearance);
  SC_CHECK(opt.stamp_free_radius_m > opt.plan.preferred_clearance * 0.5);
}

void testNavigateReconstructedMapIsHonest() {
  // `runNavigation` plans on the *reconstructed* map, and its status must agree
  // with whether it produced a usable route. Two failure modes are guarded here,
  // both of which were real:
  //
  //   - a zero-length "ok" plan (the goal search accepting the start node), which
  //     reported success and then drove nothing;
  //   - a plan whose reported length disagrees with the points it carries.
  //
  // The worlds that do navigate must also arrive: an `ok` status with the
  // vehicle stationary is the same false success in a different place.
  struct Case {
    sim::WorldType type;
    const char* name;
  };
  const Case cases[] = {
      {sim::WorldType::kRooms, "rooms"},
      {sim::WorldType::kCorridorLoop, "corridor_loop"},
      {sim::WorldType::kSaddleHall, "saddle_hall"},
  };
  for (const Case& c : cases) {
    sim::SimConfig cfg;
    cfg.world.type = c.type;
    cfg.world.seed = 6;
    cfg.world.size_m = 40;
    cfg.duration = 60.0;
    cfg.record_images = false;
    app::NavigationOptions opt;
    opt.drive = true;
    const app::NavigationResult r = app::runNavigation(cfg, opt);
    if (!r.ok) {
      std::printf("  %-14s skipped (map not plannable on this seed)\n", c.name);
      continue;
    }
    std::printf("  %-14s status=%s plan=%.1f m driven=%.1f m\n", c.name,
                planning::toString(r.stats.plan_status), r.stats.plan_length,
                r.stats.drive_distance);
    // An accepted plan carries a real route.
    SC_CHECK(r.stats.plan_length > 1.0);
    SC_CHECK(r.planned_path.size() >= 2);
    // And the closed loop, run on the real world, reaches the end of it.
    SC_CHECK(r.stats.drive_distance > 0.8 * r.stats.plan_length);
    SC_CHECK(r.stats.collisions == 0);
  }
}

void run() {
  testChanceConstraintMonotonic();
  testCovarianceDetMonotonic();
  testVelocityProfileRespectsBraking();
  testPlannerClearanceAndDrivability();
  testMorseEulerRelation();
  testSafeStopUnderFault();
  testClosedLoopNoCollision();
  testStampRadiusExceedsPlanClearance();
  testNavigateReconstructedMapIsHonest();
}

}  // namespace

SC_TEST_MAIN("test_navigation")
