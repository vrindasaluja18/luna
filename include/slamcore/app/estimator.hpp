// SLAMCore-8G -- the estimator that ties the pieces together.
//
//   odometry  ->  pose graph  <-  loop closures
//                    |
//                    +-> optimised trajectory (TUM), ATE/RPE
//
// This is deliberately a thin orchestrator: every part of the pipeline is
// usable and testable on its own, and this file is the only place that knows
// the order they run in. A pipeline that is only reachable through its
// top-level entry point cannot be benchmarked piecewise, and the contract
// needs the dead-reckoning number and the optimised number side by side.
#pragma once

#include <string>
#include <vector>

#include "slamcore/backend/loop_close.hpp"
#include "slamcore/backend/pose_graph.hpp"
#include "slamcore/bench/trajectory.hpp"
#include "slamcore/frontend/features.hpp"
#include "slamcore/frontend/odometry.hpp"
#include "slamcore/sim/simulator.hpp"

namespace sc {
namespace app {

struct EstimatorOptions {
  frontend::OdomOptions odom;
  backend::LoopCloseOptions loop;
  /// Consistency gate on an accepted closure. The measured transform must agree
  /// with what the graph already believes about the relative pose of the two
  /// keyframes, within these bounds. The graph's belief is essentially the
  /// odometry chain, whose drift is small over a single loop, while a
  /// perceptually aliased match (the same door, matched to a different door)
  /// claims an arbitrary transform. Without this gate those aliased matches pass
  /// RANSAC with a dozen inliers, because the correspondence set really is
  /// geometrically consistent -- it is just consistently wrong.
  Scalar loop_gate_translation = 0.6;  ///< metres
  Scalar loop_gate_rotation = 0.3;    ///< radians
  backend::GraphOptions graph;
  /// FAST corner threshold. The simulator renders band-limited procedural
  /// texture, which has far less high-frequency energy than a real camera, so
  /// the default threshold tuned for real imagery finds only a handful of
  /// corners per frame. A handful of descriptors cannot support loop-closure
  /// retrieval, and the failure looks like a broken matcher rather than a
  /// threshold that is simply too high for the image statistics.
  int fast_threshold = 10;
  /// Optimise after every keyframe rather than only at the end. This is what
  /// makes the reported latency a real per-keyframe figure instead of a
  /// one-off cost, and it is the mode the real-time contract applies to.
  bool incremental = true;
  /// Run a global optimisation at the end over the whole graph.
  bool final_optimise = true;
};

struct EstimatorStats {
  Size keyframes = 0;
  Size frames = 0;
  Size graph_edges = 0;
  Size loop_closures = 0;
  Size loop_rejected = 0;
  /// Accepted closures that agree with ground truth, and those that do not.
  /// Precision is the contract's figure; a closure is counted as correct when
  /// the two keyframes really were at nearly the same place, judged from the
  /// simulator's ground truth and never from anything the matcher produced.
  Size loop_true_positives = 0;
  Size loop_false_positives = 0;
  Size dropped_frames = 0;
  Size graph_failures = 0;
  Scalar odom_distance = 0;
  Scalar keyframe_distance = 0;
  Scalar optimise_ms_total = 0;
  Scalar optimise_ms_max = 0;
  Scalar optimise_ms_p99 = 0;
  Scalar estimate_ms_total = 0;
  int graph_iterations = 0;
};

struct EstimatorResult {
  /// Ground truth, captured from the frames the estimator already consumes.
  /// Collected here so a benchmark run needs one simulation rather than two,
  /// and so the estimate and its reference are guaranteed to come from the
  /// same run.
  bench::Trajectory truth;
  /// Dead reckoning only: the baseline the optimiser must beat.
  bench::Trajectory odom_trajectory;
  /// After pose-graph optimisation.
  bench::Trajectory optimised_trajectory;
  bench::AteResult odom_ate;
  bench::AteResult optimised_ate;
  bench::RpeResult odom_rpe;
  bench::RpeResult optimised_rpe;
  /// Marginal covariance of the final pose, from the pose graph's information
  /// matrix (Laplace approximation). This is the estimator's own statement of
  /// how well it knows where it is, and it is what the controller's chance
  /// constraint consumes. Using the ATE instead would be circular: ATE is
  /// computed *against ground truth*, so a controller driven by it would be
  /// using the answer key.
  Mat3 final_covariance = Mat3::Identity();
  bool final_covariance_ok = false;
  /// Covariance of the *recent* pose relative to a keyframe a few metres back,
  /// which is the bounded, locally-relevant uncertainty. See
  /// `PoseGraph::relativeCovariance` for why a controller must not use the
  /// unbounded marginal.
  Mat3 local_covariance = Mat3::Identity();
  bool local_covariance_ok = false;
  EstimatorStats stats;
  bool ok = false;
};

/// Run the full pipeline over one simulated session.
///
/// No ground-truth landmark id or pose ever feeds a constraint: the frames'
/// `ground_truth` field is written to the result for reporting only, so the
/// accuracy figures measure the pipeline rather than its access to the answer
/// key.
EstimatorResult runEstimator(const sim::SimConfig& cfg,
                             const EstimatorOptions& opt = EstimatorOptions());

}  // namespace app
}  // namespace sc
