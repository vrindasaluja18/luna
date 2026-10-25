#include "slamcore/app/estimator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace sc {
namespace app {

namespace {

/// Unproject a pixel at known depth into the camera's body frame.
///
/// The camera looks along body +x with +y to the left and +z up, which is the
/// convention sim::project uses. The two signs below are mirrored from that
/// function deliberately: getting them backwards places every map point behind
/// the vehicle, so every loop closure fails and the failure looks like a
/// descriptor problem rather than a frame-convention problem.
Vec3 unproject(const sim::CameraModel& cam, Scalar u, Scalar v, Scalar depth) {
  const Scalar y = -(u - cam.cx) * depth / cam.fx;  // image +u is body -y
  const Scalar z = -(v - cam.cy) * depth / cam.fy;  // image +v is body -z (down)
  return Vec3(depth, y, z + cam.body_height);
}

/// Build the 3D structure a keyframe observed from the features that have a
/// depth. Features without depth carry no position and are skipped.
Size buildMapPoints(const sim::CameraModel& cam, const std::vector<frontend::Feature>& feats,
                    std::vector<backend::MapPoint>* out) {
  out->clear();
  out->reserve(feats.size());
  for (const frontend::Feature& f : feats) {
    if (f.depth <= 0) continue;
    backend::MapPoint mp;
    mp.p = unproject(cam, f.x, f.y, f.depth);
    mp.desc = f.desc;
    out->push_back(mp);
  }
  return out->size();
}

}  // namespace

EstimatorResult runEstimator(const sim::SimConfig& cfg, const EstimatorOptions& opt) {
  EstimatorResult res;

  backend::PoseGraph graph;
  graph.addNode(Pose2{0, 0, 0});
  frontend::OdomIntegrator odom(opt.odom);
  backend::LoopCloser closer(opt.loop);

  const frontend::BriefPattern pattern;
  frontend::FastOptions fast;
  fast.threshold = opt.fast_threshold;
  frontend::BriefOptions brief;
  brief.max_features = fast.max_features;

  sim::Simulator sim(cfg);
  Size kf_id = 0;
  // Dead-reckoned pose at the last keyframe. The odometry increment must be the
  // relative transform between two *dead-reckoned* poses. Using the graph's
  // optimised pose on one side and the raw pose on the other mixes two frames
  // and injects the optimiser's own correction as a bogus measurement, which
  // makes the graph strictly worse than the odometry it started from.
  Pose2 last_kf_dr{0, 0, 0};
  std::vector<Scalar> optimise_times;
  // Real keyframe timestamps, one per graph node and in node order. The
  // optimised trajectory must carry these and not keyframe indices: ATE matches
  // the estimate against ground truth *by time*, so index-as-time would compare
  // the tenth keyframe against the tenth second of the run and report a large
  // error for a perfect estimate.
  std::vector<Scalar> kf_times;
  // Ground-truth pose of each keyframe, used only to score accepted closures.
  // Kept out of every constraint: the matcher must not see it.
  std::vector<Pose2> kf_truth;
  // The simulator's ground truth is in world coordinates, but odometry has no
  // absolute reference and necessarily starts at the origin. Ground truth is
  // therefore re-expressed relative to its own first pose, so the two
  // trajectories share a frame. This is the standard odometry-evaluation
  // convention: the metric measures accumulated drift, not the unknown constant
  // offset between the world frame and the estimator's start frame.
  bool have_truth0 = false;
  Pose2 truth0{0, 0, 0};

  const auto t_start = std::chrono::steady_clock::now();

  sim.run([&](const sim::SensorFrame& f) {
    ++res.stats.frames;
    // Reporting-only ground truth. It is never read by any constraint.
    if (!have_truth0) {
      truth0 = f.ground_truth;
      have_truth0 = true;
    }
    res.truth.push_back(bench::TrajectoryPoint{f.t, pose2Relative(truth0, f.ground_truth)});

    odom.add(f.odom, f.imu);
    res.stats.odom_distance = odom.totalDistance();
    if (!f.camera_ok) {
      ++res.stats.dropped_frames;
      return;
    }

    // The dead-reckoned pose after this frame. This is the odometry-only
    // baseline, but it is *not* recorded yet: the graph only has a node per
    // keyframe, so the optimised trajectory is sampled at keyframe times. ATE
    // averages over whatever samples it is given, so comparing a full-rate
    // odometry trajectory against a keyframe-rate optimised one compares two
    // different sample sets, and the odometry number is the one that changes.
    // Both trajectories are therefore built from the same keyframe instants,
    // below.
    const Pose2 dr = odom.pose();

    // Keyframes are chosen on distance travelled rather than frame count, so
    // map density does not depend on the camera rate.
    const Scalar moved = std::hypot(dr.x - last_kf_dr.x, dr.y - last_kf_dr.y);
    if (moved < opt.loop.keyframe_distance && kf_id > 0) return;

    // -------- new keyframe ------------------------------------------------
    std::vector<frontend::Feature> feats;
    frontend::detectAndDescribe(f.image, pattern, fast, brief, &feats);
    if (cfg.stereo && !f.image_right.px.empty()) {
      frontend::StereoOptions so;
      frontend::stereoDepth(f.image, f.image_right, cfg.camera, so, &feats);
    }

    const frontend::OdomConstraint inc = odom.take();

    // Extend the graph *before* building the keyframe, so the keyframe's own
    // node exists and its pose is well defined.
    if (kf_id > 0) {
      graph.addNode(dr);
      backend::PoseEdge e;
      e.i = static_cast<Index>(kf_id - 1);
      e.j = static_cast<Index>(kf_id);
      // The increment since the previous keyframe, expressed in that
      // keyframe's frame. Both sides are dead reckoning, so the measurement
      // contains only sensor information.
      e.z = se2Log(se2Inverse(pose2ToSE2(last_kf_dr)) * pose2ToSE2(dr));
      e.info = inc.info;
      // Odometry is the graph's backbone: a bad wheel sample should show up in
      // the cost rather than be silently down-weighted, so no robust kernel.
      e.robust = false;
      graph.addEdge(e);
    }

    // -------- loop closure ------------------------------------------------
    backend::Keyframe kf;
    kf.id = static_cast<Index>(kf_id);
    kf.pose = graph.pose(static_cast<Index>(kf_id));
    kf.t = f.t;
    buildMapPoints(cfg.camera, feats, &kf.points);

    backend::LoopClosure lc;
    // Reject closures that disagree with the graph's current belief about the
    // relative pose. Checked before the edge is added so a bad match never
    // enters the problem.
    bool loop_accepted = false;
    if (closer.query(kf, &lc)) {
      const Vec3 predicted =
          se2Log(se2Inverse(pose2ToSE2(graph.pose(lc.from))) * pose2ToSE2(graph.pose(lc.to)));
      const Vec3 disagree = se2Log(se2Inverse(se2Exp(predicted)) * se2Exp(lc.z));
      loop_accepted = std::hypot(disagree.x(), disagree.y()) <= opt.loop_gate_translation &&
                      std::abs(disagree.z()) <= opt.loop_gate_rotation;
    }
    if (loop_accepted) {
      backend::PoseEdge e;
      e.i = lc.from;
      e.j = lc.to;
      e.z = lc.z;
      e.info = lc.info;
      // The kernel stays on even though the closure passed verification. It is
      // not the primary defence -- the geometric check and the gate above are --
      // but it bounds the damage if a wrong match slips through, and a single
      // wrong closure corrupts every pose between its two ends. The kernel is
      // annealed inside the optimiser (see `PoseGraph::optimize`) so it does not
      // suppress the large, honest correction a real closure carries.
      e.robust = true;
      if (graph.addEdge(e)) {
        ++res.stats.loop_closures;
        // Score the closure against ground truth. "Correct" means the two
        // keyframes really were at nearly the same pose, so the constraint the
        // matcher added is one the world actually supports.
        if (lc.from >= 0 && lc.from < static_cast<Index>(kf_truth.size())) {
          const Pose2 rel = pose2Relative(kf_truth[static_cast<Size>(lc.from)], f.ground_truth);
          // Same convention as the graph's edge residual: r = Log(Z^-1 *
          // x_from^-1 * x_to). Using Z^-1 * rel rather than rel^-1 * Z matters
          // for large rotations: the two differ by the adjoint, which is
          // enormous when the relative rotation is near pi, and a correct
          // closure would then be scored as false.
          const Vec3 err = se2Log(se2Inverse(se2Exp(lc.z)) * pose2ToSE2(rel));
          if (std::hypot(err.x(), err.y()) < 0.5 && std::abs(err.z()) < 0.3) {
            ++res.stats.loop_true_positives;
          } else {
            ++res.stats.loop_false_positives;
          }
        }
      } else {
        ++res.stats.loop_rejected;
      }
    }

    closer.addKeyframe(kf);
    last_kf_dr = dr;
    ++kf_id;
    kf_times.push_back(f.t);
    kf_truth.push_back(f.ground_truth);
    // Same instant, same sample set as the optimised trajectory.
    res.odom_trajectory.push_back(bench::TrajectoryPoint{f.t, dr});
    ++res.stats.keyframes;

    if (opt.incremental && kf_id > 1) {
      const auto t0 = std::chrono::steady_clock::now();
      const backend::GraphResult gr = graph.optimize(opt.graph);
      const Scalar ms =
          std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0)
              .count();
      optimise_times.push_back(ms);
      res.stats.optimise_ms_total += ms;
      res.stats.optimise_ms_max = std::max(res.stats.optimise_ms_max, ms);
      res.stats.graph_iterations += gr.iterations;
      if (!gr.ok) ++res.stats.graph_failures;
    }
  });

  res.stats.estimate_ms_total =
      std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t_start)
          .count();

  if (opt.final_optimise) {
    const auto t0 = std::chrono::steady_clock::now();
    const backend::GraphResult gr = graph.optimize(opt.graph);
    const Scalar ms =
        std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0).count();
    optimise_times.push_back(ms);
    res.stats.optimise_ms_total += ms;
    res.stats.optimise_ms_max = std::max(res.stats.optimise_ms_max, ms);
    if (!gr.ok) ++res.stats.graph_failures;
  }

  if (!optimise_times.empty()) {
    std::vector<Scalar> s = optimise_times;
    std::sort(s.begin(), s.end());
    const Size idx = std::min(s.size() - 1, static_cast<Size>(std::floor(0.99 * s.size())));
    res.stats.optimise_ms_p99 = s[idx];
  }
  res.stats.graph_edges = graph.edgeCount();

  // One optimised pose per keyframe, stamped with the keyframe's real time.
  for (Size i = 0; i < graph.nodes() && i < kf_times.size(); ++i) {
    res.optimised_trajectory.push_back(
        bench::TrajectoryPoint{kf_times[i], graph.pose(static_cast<Index>(i))});
  }

  // Marginal covariance of the last pose. Taken from the graph the pipeline
  // actually built, not from a model of it, so the number the controller sees
  // is the number the optimiser believes.
  if (graph.nodes() > 0) {
    const Index last = static_cast<Index>(graph.nodes() - 1);
    res.final_covariance_ok = graph.poseCovariance(last, &res.final_covariance);
    // Local window: the relative pose over the last `kLocalWindow` keyframes.
    // Large enough that the relative estimate is genuinely constrained by
    // independent observations, small enough that it does not accumulate the
    // drift the marginal does.
    const int kLocalWindow = 5;
    if (last >= kLocalWindow) {
      res.local_covariance_ok = graph.relativeCovariance(last - kLocalWindow, last,
                                                        &res.local_covariance);
    }
  }

  const bench::EvalOptions ev;
  const bench::Trajectory& truth = res.truth;
  res.odom_ate = bench::absoluteTrajectoryError(res.odom_trajectory, truth, ev);
  res.odom_rpe = bench::relativePoseError(res.odom_trajectory, truth, ev);
  res.optimised_ate = bench::absoluteTrajectoryError(res.optimised_trajectory, truth, ev);
  res.optimised_rpe = bench::relativePoseError(res.optimised_trajectory, truth, ev);
  res.ok = res.optimised_ate.ok;
  return res;
}

}  // namespace app
}  // namespace sc
