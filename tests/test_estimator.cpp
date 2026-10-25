// SLAMCore-8G -- end-to-end estimator tests.
//
// These run the real pipeline over the simulator and assert on the numbers it
// produces. The important property is the *comparison* between dead reckoning
// and the pose graph: a test that only checked "the optimised ATE is small"
// would pass even if the graph did nothing, and a test that only checked
// "optimisation succeeds" would pass even if it made things worse. The
// assertions below are therefore mostly of the form "the graph does not
// degrade the odometry it was given", which is what the backend is for.
#include <cmath>
#include <cstdio>
#include <string>

#include "slamcore/app/estimator.hpp"
#include "slamcore/bench/trajectory.hpp"
#include "slamcore/core/stats.hpp"
#include "slamcore/sim/simulator.hpp"
#include "test_util.hpp"

using namespace sc;

namespace {

sim::SimConfig worldConfig(sim::WorldType type, unsigned seed, Scalar duration) {
  sim::SimConfig cfg;
  cfg.world.type = type;
  cfg.world.seed = seed;
  cfg.duration = duration;
  cfg.record_images = false;
  return cfg;
}

}  // namespace

/// The pipeline must run to completion and produce trajectories with the same
/// time base as the ground truth. Timestamp alignment is the quiet failure
/// mode: an optimised trajectory stamped with keyframe indices instead of times
/// still evaluates, and reports a huge error for a perfect estimate.
void testPipelineRunsAndAligns() {
  const app::EstimatorResult r =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 1, 120));
  SC_CHECK(r.ok);
  SC_CHECK(r.stats.keyframes > 20);
  SC_CHECK(r.optimised_trajectory.size() == r.stats.keyframes);
  SC_CHECK(!r.odom_trajectory.empty());
  SC_CHECK(!r.truth.empty());

  // Every trajectory must span the same time window, and be monotone in time.
  auto spans = [](const bench::Trajectory& t) {
    return std::make_pair(t.front().t, t.back().t);
  };
  const auto [ot0, ot1] = spans(r.odom_trajectory);
  const auto [kt0, kt1] = spans(r.optimised_trajectory);
  SC_NEAR(ot0, kt0, 1e-6);
  SC_CHECK(kt1 > kt0);
  SC_CHECK(kt1 <= r.truth.back().t + 1e-6);
  for (Size i = 1; i < r.optimised_trajectory.size(); ++i) {
    SC_CHECK(r.optimised_trajectory[i].t > r.optimised_trajectory[i - 1].t);
  }
}

/// The graph is given odometry increments that are already consistent, so with
/// no loop closures it must reproduce dead reckoning rather than improve or
/// degrade it. Both outcomes would be a bug: better means the graph invented
/// information, worse means it is not solving the problem it was handed.
void testGraphReproducesOdometryWithoutLoops() {
  app::EstimatorOptions opt;
  opt.loop.min_index_gap = 1 << 30;  // no closures
  const app::EstimatorResult r =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 1, 120), opt);
  SC_CHECK(r.ok);
  SC_CHECK(r.stats.loop_closures == 0);
  SC_CHECK(r.stats.graph_failures == 0);

  // Dead reckoning must already be accurate on this short, well-textured run,
  // otherwise the comparison below is meaningless.
  SC_CHECK(r.odom_ate.percent_of_distance < 2.0);
  // The optimised result must stay within a small margin of it.
  SC_CHECK(r.optimised_ate.rmse < r.odom_ate.rmse + 0.5);
}

/// A genuine revisit must be detected, and the closures it produces must be
/// correct and must not corrupt the trajectory.
///
/// The precision assertions are the load-bearing ones. An earlier version of
/// this code reported closures while storing the *estimated* relative pose as
/// the measurement, which made every edge residual zero: the graph did nothing,
/// and every accepted closure scored as false against ground truth.
///
/// On the well-textured reference world the odometry is already accurate to
/// 0.09 % of distance, so the total drift the closures could remove is about
/// 0.2 m -- comparable to the error of a single stereo-lite loop measurement. At
/// that ratio closures cannot move the global metric, and their noise can move it
/// slightly the wrong way. The assertion here is that they do not *corrupt* the
/// trajectory; the test below covers the case where drift is large enough for
/// closures to pay off.
void testLoopClosureHelpsOrHolds() {
  const app::EstimatorResult r =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 6, 300));
  SC_CHECK(r.ok);
  SC_CHECK(r.stats.loop_closures > 0);
  SC_CHECK(r.stats.graph_failures == 0);
  // Every accepted closure must agree with the ground-truth relationship
  // between the two keyframes it joined.
  SC_CHECK(r.stats.loop_true_positives == r.stats.loop_closures);
  SC_CHECK(r.stats.loop_false_positives == 0);
  // With closures the estimate must not be materially worse than the odometry.
  SC_CHECK(r.optimised_ate.rmse <= r.odom_ate.rmse + 0.2);
  // And it must be accurate in absolute terms.
  SC_CHECK(r.optimised_ate.percent_of_distance < 1.0);
}

/// When the odometry drifts enough for the closures to have something to
/// correct, they must measurably improve the trajectory. This is the property
/// the module exists for, and it is invisible on a run whose odometry is already
/// good -- so it needs its own world with a deliberate wheel-scale error.
void testLoopClosureCorrectsRealDrift() {
  sim::SimConfig cfg = worldConfig(sim::WorldType::kRooms, 6, 300);
  cfg.odom_noise.dist_scale = 0.05;
  cfg.odom_noise.angle_scale = 0.05;
  cfg.odom_noise.dist_sigma = 0.05;
  cfg.odom_noise.angle_sigma = 0.05;

  app::EstimatorOptions no_loop;
  no_loop.loop.min_index_gap = 1 << 30;

  const app::EstimatorResult with_loops = app::runEstimator(cfg);
  const app::EstimatorResult without = app::runEstimator(cfg, no_loop);

  SC_CHECK(with_loops.ok && without.ok);
  SC_CHECK(with_loops.stats.loop_closures > 0);
  SC_CHECK(with_loops.stats.loop_false_positives == 0);
  // The drift is what makes the comparison meaningful.
  SC_CHECK(without.optimised_ate.percent_of_distance > 0.3);
  // The closures must reduce it.
  SC_CHECK(with_loops.optimised_ate.rmse < without.optimised_ate.rmse);
}

/// The contract's latency budget is measured per optimisation. If a single
/// optimisation overruns, the real-time claim is false regardless of the mean.
void testOptimiseLatencyBudget() {
  const app::EstimatorResult r =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 6, 300));
  SC_CHECK(r.ok);
  SC_CHECK(r.stats.optimise_ms_p99 < 450.0);
  SC_CHECK(r.stats.optimise_ms_max < 450.0);
}

/// Determinism: the same seed must produce bit-identical results. Without this
/// no benchmark number is reproducible, which is the contract's central claim.
void testDeterminism() {
  const app::EstimatorResult a =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 4, 120));
  const app::EstimatorResult b =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 4, 120));
  SC_CHECK(a.optimised_trajectory.size() == b.optimised_trajectory.size());
  SC_CHECK(a.stats.loop_closures == b.stats.loop_closures);
  for (Size i = 0; i < a.optimised_trajectory.size(); ++i) {
    SC_CHECK(a.optimised_trajectory[i].t == b.optimised_trajectory[i].t);
    SC_CHECK(a.optimised_trajectory[i].pose.x == b.optimised_trajectory[i].pose.x);
    SC_CHECK(a.optimised_trajectory[i].pose.y == b.optimised_trajectory[i].pose.y);
    SC_CHECK(a.optimised_trajectory[i].pose.th == b.optimised_trajectory[i].pose.th);
  }
}

/// The trajectory writer must round-trip through the TUM format the evaluator
/// and the external tooling both read.
void testTumRoundTrip() {
  const app::EstimatorResult r =
      app::runEstimator(worldConfig(sim::WorldType::kRooms, 1, 60));
  const std::string path = "estimator_tum_roundtrip.tum";
  SC_CHECK(bench::writeTum(path, r.optimised_trajectory));
  bench::Trajectory back;
  SC_CHECK(bench::readTum(path, &back));
  std::remove(path.c_str());
  SC_CHECK(back.size() == r.optimised_trajectory.size());
  for (Size i = 0; i < back.size(); ++i) {
    SC_NEAR(back[i].pose.x, r.optimised_trajectory[i].pose.x, 1e-6);
    SC_NEAR(back[i].pose.y, r.optimised_trajectory[i].pose.y, 1e-6);
  }
}

void run() {
  testPipelineRunsAndAligns();
  testGraphReproducesOdometryWithoutLoops();
  testLoopClosureHelpsOrHolds();
  testLoopClosureCorrectsRealDrift();
  testOptimiseLatencyBudget();
  testDeterminism();
  testTumRoundTrip();
}

SC_TEST_MAIN("test_estimator")
