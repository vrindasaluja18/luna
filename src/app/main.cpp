// LUNA -- command line interface.
//
//   luna sim   -- generate a world, run the scripted route, render frames
//   luna fe    -- run the frontend (FAST+BRIEF+KLT+stereo) over a sim run
//   luna world -- print a world as ASCII and write a PGM occupancy map
//
// Everything here is CPU-only and self-contained: no datasets, no ROS, no GPU.
// The `fe` subcommand is the end-to-end vertical slice the contract asks for,
// and it reports the same quantities the benchmark tables use.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "luna/app/estimator.hpp"
#include "luna/app/navigation.hpp"
#include "luna/bench/trajectory.hpp"
#include "luna/frontend/features.hpp"
#include "luna/sim/faults.hpp"
#include "luna/sim/simulator.hpp"

namespace {

using luna::Scalar;
using luna::Size;
using luna::Vec2;
using luna::frontend::BriefOptions;
using luna::frontend::BriefPattern;
using luna::frontend::FastOptions;
using luna::frontend::Feature;
using luna::frontend::KltOptions;
using luna::frontend::KltTracker;
using luna::frontend::StereoOptions;
using luna::sim::SimConfig;
using luna::sim::Simulator;
using luna::sim::SimResult;
using luna::sim::SensorFrame;
using luna::sim::World;
using luna::sim::WorldSpec;
using luna::sim::WorldType;

int usage() {
  std::fprintf(stderr,
               "usage: luna <command> [options]\n"
               "\n"
               "commands:\n"
               "  world [--type T] [--size M] [--seed N] [--rooms N] [--pgm FILE]\n"
               "        generate a world and print it as ASCII\n"
               "  sim   [--type T] [--seed N] [--duration S] [--speed MPS] [--faults]\n"
               "        [--out PREFIX] [--pgm FILE]\n"
               "        run the scripted route and write ground truth\n"
               "  fe    [--type T] [--seed N] [--duration S] [--speed MPS] [--faults]\n"
               "        [--stereo 0|1] [--fast N] [--klt-levels N]\n"
               "        run the frontend over a sim run and report its metrics\n"
               "  estimate [--type T] [--seed N] [--duration S] [--speed MPS] [--faults]\n"
               "        [--stereo 0|1] [--fast N] [--kf-dist M] [--out PREFIX]\n"
               "        [--no-loop] [--no-incremental]\n"
               "        run SLAM end to end (odometry -> graph -> loop closure) and\n"
               "        report ATE/RPE against ground truth, plus dead reckoning for\n"
               "        comparison\n"
               "\n"
               "world types: rooms corridor_loop saddle_hall maze cluttered textureless\n");
  return 2;
}

bool parseWorld(const std::string& s, WorldType* out) {
  if (s == "rooms") return *out = WorldType::kRooms, true;
  if (s == "corridor_loop" || s == "corridor") return *out = WorldType::kCorridorLoop, true;
  if (s == "saddle_hall" || s == "saddle") return *out = WorldType::kSaddleHall, true;
  if (s == "maze") return *out = WorldType::kMaze, true;
  if (s == "cluttered") return *out = WorldType::kCluttered, true;
  if (s == "textureless") return *out = WorldType::kTextureless, true;
  return false;
}

/// Pull "--name value" out of argv. Returns nullptr when absent.
const char* opt(int argc, char** argv, const char* name) {
  for (int i = 0; i < argc; ++i)
    if (std::strcmp(argv[i], name) == 0 && i + 1 < argc) return argv[i + 1];
  return nullptr;
}

bool flag(int argc, char** argv, const char* name) {
  for (int i = 0; i < argc; ++i)
    if (std::strcmp(argv[i], name) == 0) return true;
  return false;
}

SimConfig baseConfig(int argc, char** argv) {
  SimConfig cfg;
  cfg.world.type = WorldType::kRooms;
  if (const char* t = opt(argc, argv, "--type")) {
    WorldType wt;
    if (!parseWorld(t, &wt)) {
      std::fprintf(stderr, "unknown world type '%s'\n", t);
      std::exit(2);
    }
    cfg.world.type = wt;
  }
  if (const char* v = opt(argc, argv, "--seed")) cfg.world.seed = std::strtoull(v, nullptr, 10);
  if (const char* v = opt(argc, argv, "--size")) cfg.world.size_m = std::atof(v);
  if (const char* v = opt(argc, argv, "--rooms")) cfg.world.rooms = std::atoi(v);
  if (const char* v = opt(argc, argv, "--duration")) cfg.duration = std::atof(v);
  if (const char* v = opt(argc, argv, "--speed")) cfg.target_speed = std::atof(v);
  if (const char* v = opt(argc, argv, "--stereo")) cfg.stereo = std::atoi(v) != 0;
  if (flag(argc, argv, "--faults")) {
    // The standard suite starts once the vehicle is moving and is spaced so
    // each fault is individually observable in the metrics.
    cfg.faults = luna::sim::FaultSchedule::standardSuite(3.0, std::max(1.0, cfg.duration / 12.0));
  }
  return cfg;
}

int cmdWorld(int argc, char** argv) {
  const SimConfig cfg = baseConfig(argc, argv);
  const World world = World::generate(cfg.world);
  const auto& st = world.stats();
  std::printf("world type=%s seed=%llu size=%.1fm\n", luna::sim::toString(cfg.world.type),
              static_cast<unsigned long long>(cfg.world.seed), cfg.world.size_m);
  std::printf("  grid      : %dx%d cells @ %.2f m\n", world.grid().w, world.grid().h, world.grid().res);
  std::printf("  free area : %.1f m^2 (%zu cells)\n", st.freeArea, st.freeCells);
  std::printf("  landmarks : %zu\n", st.landmarks);
  std::printf("  corridors : min %.2f m, median %.2f m\n", st.corridorWidthMin,
              st.corridorWidthMedian);
  std::printf("  grid RAM  : %.1f KB\n", world.grid().bytes() / 1024.0);
  std::printf("\n%s\n", world.toAscii(96).c_str());
  if (const char* pgm = opt(argc, argv, "--pgm")) {
    std::printf("wrote %s: %s\n", pgm, luna::sim::writePgm(pgm, world.grid()) ? "ok" : "FAILED");
  }
  return 0;
}

int cmdSim(int argc, char** argv) {
  const SimConfig cfg = baseConfig(argc, argv);
  Simulator sim(cfg);
  const SimResult r = sim.run(nullptr);
  std::printf("sim type=%s seed=%llu duration=%.1fs\n", luna::sim::toString(cfg.world.type),
              static_cast<unsigned long long>(cfg.world.seed), cfg.duration);
  std::printf("  frames    : %zu (%zu dropped)\n", r.frames, r.dropped_frames);
  std::printf("  imu       : %zu samples\n", r.imu_samples);
  std::printf("  path      : %.1f m\n", r.path_length);
  std::printf("  collisions: %zu\n", r.collisions);
  std::printf("  faults    : %zu injected\n", r.faults_injected);
  std::printf("  wall time : %.2f s (%.1fx real time)\n", r.wall_time_s,
              r.wall_time_s > 0 ? cfg.duration / r.wall_time_s : 0.0);
  if (const char* prefix = opt(argc, argv, "--out")) {
    const std::string gt = std::string(prefix) + "_gt.tum";
    const std::string lm = std::string(prefix) + "_landmarks.txt";
    std::printf("  wrote %s: %s\n", gt.c_str(), r.truth.saveTum(gt) ? "ok" : "FAILED");
    std::printf("  wrote %s: %s\n", lm.c_str(), r.truth.saveLandmarks(lm) ? "ok" : "FAILED");
  }
  if (const char* pgm = opt(argc, argv, "--pgm")) {
    std::printf("  wrote %s: %s\n", pgm, luna::sim::writePgm(pgm, sim.world().grid()) ? "ok" : "FAILED");
  }
  return 0;
}

/// Runs the feature frontend over a simulated run and reports its metrics.
///
/// The interesting numbers are: how many features survive tracking frame to
/// frame (a low rate means the tracker, not the scene, is the bottleneck), how
/// accurate the tracked positions are against the simulator's exact landmark
/// projections, and how accurate the stereo depth is against the true range.
/// Those three together are what the backend will later consume.
int cmdFrontend(int argc, char** argv) {
  SimConfig cfg = baseConfig(argc, argv);
  cfg.record_images = false;

  FastOptions fast;
  fast.threshold = opt(argc, argv, "--fast") ? std::atoi(opt(argc, argv, "--fast")) : 12;
  BriefOptions brief;
  brief.max_features = fast.max_features;
  KltOptions klt_opt;
  if (const char* v = opt(argc, argv, "--klt-levels")) klt_opt.max_level = std::atoi(v);
  StereoOptions stereo_opt;

  const BriefPattern pattern;
  KltTracker klt(klt_opt);

  Simulator sim(cfg);
  luna::sim::Image prev;
  std::vector<Feature> prev_feats;
  bool have_prev = false;

  Size frames = 0, detected = 0, tracked = 0, depth_valid = 0, stereo_frames = 0;
  Size coarse_rej = 0, fb_rej = 0;
  Size track_err_n = 0;
  Scalar track_err_sum = 0, track_err_sq = 0;
  Scalar detect_ms = 0, track_ms = 0, depth_ms = 0;
  std::vector<Scalar> track_errs;

  const auto t0 = std::chrono::steady_clock::now();
  const SimResult r = sim.run([&](const SensorFrame& f) {
    ++frames;

    std::vector<Feature> feats;
    auto a = std::chrono::steady_clock::now();
    luna::frontend::detectAndDescribe(f.image, pattern, fast, brief, &feats);
    auto b = std::chrono::steady_clock::now();

    if (cfg.stereo && !f.image_right.px.empty()) {
      luna::frontend::stereoDepth(f.image, f.image_right, cfg.camera, stereo_opt, &feats);
      ++stereo_frames;
    }
    auto c = std::chrono::steady_clock::now();

    if (have_prev) {
      std::vector<Feature> out;
      klt.track(prev, f.image, prev_feats, &out);
      const auto st = klt.stats();
      tracked += st.tracked;
      coarse_rej += st.coarse;
      fb_rej += st.fb;

      // Position accuracy against the simulator's exact landmark projections.
      // This is the ground truth a dataset cannot give: the true pixel a
      // landmark projects to, independent of any detection or matching.
      for (const Feature& tf : out) {
        Scalar best = 1e9;
        for (Size i = 0; i < f.visible_pixels.size(); ++i) {
          const Scalar d2 = (f.visible_pixels[i] - Vec2(tf.x, tf.y)).squaredNorm();
          if (d2 < best) best = d2;
        }
        if (best > 9.0) continue;  // no landmark projection nearby
        const Scalar e = std::sqrt(best);
        track_err_sum += e;
        track_err_sq += e * e;
        track_errs.push_back(e);
        ++track_err_n;
      }
    }
    auto d = std::chrono::steady_clock::now();

    if (have_prev) {
      detect_ms += std::chrono::duration<double, std::milli>(b - a).count();
      depth_ms += std::chrono::duration<double, std::milli>(c - b).count();
      track_ms += std::chrono::duration<double, std::milli>(d - c).count();
    }

    detected += feats.size();
    for (const Feature& x : feats)
      if (x.depth > 0) ++depth_valid;

    prev = f.image;
    prev_feats = feats;
    have_prev = true;
  });
  const auto t1 = std::chrono::steady_clock::now();
  const Scalar wall = std::chrono::duration<double>(t1 - t0).count();

  const Scalar n_track_frames = frames > 0 ? static_cast<Scalar>(frames - 1) : 1.0;
  std::printf("frontend type=%s seed=%llu frames=%zu faults=%s\n", luna::sim::toString(cfg.world.type),
              static_cast<unsigned long long>(cfg.world.seed), frames,
              flag(argc, argv, "--faults") ? "on" : "off");
  std::printf("  features  : %.1f detected/frame, %.1f tracked/frame (%.1f%% of detected)\n",
              detected / n_track_frames, tracked / n_track_frames,
              100.0 * static_cast<Scalar>(tracked) / std::max<Scalar>(1.0, static_cast<Scalar>(detected)));
  std::printf("  rejects   : %.1f coarse/frame, %.1f forward-backward/frame\n",
              coarse_rej / n_track_frames, fb_rej / n_track_frames);
  std::printf("  track err : mean %.4f px, rms %.4f px, p95 %.4f px (n=%zu)\n",
              track_err_n ? track_err_sum / track_err_n : 0.0,
              track_err_n ? std::sqrt(track_err_sq / track_err_n) : 0.0,
              luna::RunningStats::percentile(track_errs, 0.95), track_err_n);
  std::printf("  stereo    : %.1f valid depths/frame over %zu stereo frames\n",
              depth_valid / std::max<Scalar>(1.0, static_cast<Scalar>(stereo_frames)), stereo_frames);
  std::printf("  timing    : detect %.2f ms, track %.2f ms, stereo %.2f ms, total %.2f ms/frame\n",
              detect_ms / n_track_frames, track_ms / n_track_frames, depth_ms / n_track_frames,
              wall * 1000.0 / n_track_frames);
  std::printf("  sim       : path %.1f m, collisions %zu, dropped frames %zu\n", r.path_length,
              r.collisions, r.dropped_frames);
  return 0;
}

}  // namespace

/// Runs the whole pipeline and reports accuracy against the simulator's own
/// ground truth.
///
/// Two numbers matter and they are printed together on purpose: dead reckoning
/// alone, and the pose graph. Reporting only the optimised figure would make it
/// impossible to tell whether the backend is contributing anything, which is
/// the entire point of having a backend.
int cmdEstimate(int argc, char** argv) {
  SimConfig cfg = baseConfig(argc, argv);
  cfg.record_images = false;
  luna::app::EstimatorOptions eopt;
  if (flag(argc, argv, "--no-incremental")) eopt.incremental = false;
  if (flag(argc, argv, "--no-loop")) {
    // Push the candidate gap beyond any possible index, so no pair is ever
    // considered. This isolates the effect of loop closure rather than merely
    // weakening it.
    eopt.loop.min_index_gap = 1 << 30;
  }
  if (const char* v = opt(argc, argv, "--kf-dist")) eopt.loop.keyframe_distance = std::atof(v);
  if (const char* v = opt(argc, argv, "--fast")) eopt.fast_threshold = std::atoi(v);

  const luna::app::EstimatorResult er = luna::app::runEstimator(cfg, eopt);

  std::printf("estimate type=%s seed=%llu duration=%.1fs\n", luna::sim::toString(cfg.world.type),
              static_cast<unsigned long long>(cfg.world.seed), cfg.duration);
  std::printf("  keyframes : %zu (of %zu frames, %zu dropped)\n", er.stats.keyframes,
              er.stats.frames, er.stats.dropped_frames);
  std::printf("  graph     : %zu nodes, %zu edges, %zu loop closures (%zu rejected)\n",
              er.stats.keyframes, er.stats.graph_edges, er.stats.loop_closures,
              er.stats.loop_rejected);
  if (er.stats.loop_closures > 0) {
    const Size tp = er.stats.loop_true_positives;
    const Size fp = er.stats.loop_false_positives;
    std::printf("  closures  : %zu verified correct, %zu false, precision %.1f%%\n", tp, fp,
                tp + fp > 0 ? 100.0 * static_cast<double>(tp) / static_cast<double>(tp + fp) : 0.0);
  }
  std::printf("  distance  : %.1f m dead reckoning\n", er.stats.odom_distance);
  std::printf("\n");
  std::printf("  %-14s %10s %10s %10s %10s\n", "estimate", "ATE(m)", "ATE(%)", "RPE_t(m)",
              "RPE_r(deg)");
  std::printf("  %-14s %10.4f %10.4f %10.4f %10.4f\n", "dead reckoning", er.odom_ate.rmse,
              er.odom_ate.percent_of_distance, er.odom_rpe.mean_trans, er.odom_rpe.mean_rot);
  std::printf("  %-14s %10.4f %10.4f %10.4f %10.4f\n", "pose graph", er.optimised_ate.rmse,
              er.optimised_ate.percent_of_distance, er.optimised_rpe.mean_trans,
              er.optimised_rpe.mean_rot);
  std::printf("\n");
  std::printf("  timing    : optimise p99 %.1f ms, max %.1f ms, total %.1f ms\n",
              er.stats.optimise_ms_p99, er.stats.optimise_ms_max, er.stats.optimise_ms_total);
  std::printf("  pipeline  : %.1f ms wall (%.1fx real time)\n", er.stats.estimate_ms_total,
              er.stats.estimate_ms_total > 0 ? 1000.0 * cfg.duration / er.stats.estimate_ms_total
                                             : 0.0);
  if (er.stats.graph_failures) {
    std::printf("  WARNING   : %zu graph optimisation failures\n", er.stats.graph_failures);
  }

  if (const char* prefix = opt(argc, argv, "--out")) {
    const std::string od = std::string(prefix) + "_odom.tum";
    const std::string op = std::string(prefix) + "_slam.tum";
    std::printf("\n  wrote %s: %s\n", od.c_str(),
                luna::bench::writeTum(od, er.odom_trajectory) ? "ok" : "FAILED");
    std::printf("  wrote %s: %s\n", op.c_str(),
                luna::bench::writeTum(op, er.optimised_trajectory) ? "ok" : "FAILED");
  }
  return er.ok ? 0 : 1;
}

int cmdNavigate(int argc, char** argv) {
  SimConfig cfg = baseConfig(argc, argv);
  luna::app::NavigationOptions nopt;
  nopt.drive = !flag(argc, argv, "--no-drive");
  if (const char* v = opt(argc, argv, "--inflation")) nopt.map_inflation_m = std::atof(v);
  if (const char* v = opt(argc, argv, "--kf-dist")) nopt.estimator.loop.keyframe_distance = std::atof(v);
  if (const char* v = opt(argc, argv, "--drive-duration")) nopt.drive_duration = std::atof(v);

  const luna::app::NavigationResult nr = luna::app::runNavigation(cfg, nopt);
  const auto& s = nr.stats;
  std::printf("navigate type=%s seed=%llu duration=%.1fs stereo=%d\n",
              luna::sim::toString(cfg.world.type), static_cast<unsigned long long>(cfg.world.seed),
              cfg.duration, cfg.stereo ? 1 : 0);
  std::printf("\n  [SLAM]\n");
  std::printf("    keyframes      : %zu\n", s.keyframes);
  std::printf("    ATE            : %.4f m (%.4f%% of distance)\n", s.slam_ate, s.slam_ate_percent);
  std::printf("    loop closures  : %zu, precision %.1f%%\n", s.loop_closures,
              100.0 * s.loop_precision);
  if (s.slam_covariance_ok) {
    std::printf("    final sigma    : %.4f m position, %.4f rad heading (Laplace)\n",
                s.slam_sigma_xy, s.slam_sigma_th);
  } else {
    std::printf("    final sigma    : UNAVAILABLE (graph singular)\n");
  }
  std::printf("\n  [MAP]\n");
  std::printf("    occupied cells : %zu, known %.1f%% of grid\n", s.map_occupied_cells,
              100.0 * s.map_known_fraction);
  std::printf("    tiles resident : %zu (%.1f MB)\n", s.map_tiles_resident,
              static_cast<double>(s.map_bytes_resident) / 1048576.0);
  std::printf("\n  [TOPOLOGY]\n");
  std::printf("    Morse complex  : %d nodes, %d edges (%d maxima, %d saddles), %.1f ms\n",
              s.morse_nodes, s.morse_edges, s.maxima, s.saddles, s.topology_ms);
  std::printf("\n  [PLAN]\n");
  std::printf("    status         : %s\n", luna::planning::toString(s.plan_status));
  std::printf("    length         : %.1f m, min clearance %.3f m, %.1f ms\n", s.plan_length,
              s.plan_clearance, s.plan_ms);
  std::printf("\n  [CONTROL]\n");
  std::printf("    ticks          : %zu, driven %.1f m\n", s.control_ticks, s.drive_distance);
  std::printf("    tracking       : rms %.3f m, max %.3f m\n", s.tracking_rms, s.tracking_max);
  std::printf("    solve latency  : p99 %.1f ms, max %.1f ms\n", s.control_ms_p99, s.control_ms_max);
  std::printf("    goal error     : %.3f m\n", s.goal_error);
  std::printf("    safe stop      : %s\n", s.safe_stopped ? "yes" : "no");
  std::printf("    COLLISIONS     : %zu\n", s.collisions);

  if (const char* prefix = opt(argc, argv, "--out")) {
    const std::string pl = std::string(prefix) + "_plan.tum";
    const std::string dr = std::string(prefix) + "_drive.tum";
    luna::bench::Trajectory pt, dt;
    for (Size i = 0; i < nr.planned_path.size(); ++i)
      pt.push_back(luna::bench::TrajectoryPoint{static_cast<Scalar>(i) * 0.1,
                                              luna::Pose2{nr.planned_path[i].x(), nr.planned_path[i].y(), 0}});
    for (Size i = 0; i < nr.driven_path.size(); ++i)
      dt.push_back(luna::bench::TrajectoryPoint{static_cast<Scalar>(i) * 0.05,
                                              luna::Pose2{nr.driven_path[i].x(), nr.driven_path[i].y(), 0}});
    std::printf("\n  wrote %s: %s\n", pl.c_str(), luna::bench::writeTum(pl, pt) ? "ok" : "FAILED");
    std::printf("  wrote %s: %s\n", dr.c_str(), luna::bench::writeTum(dr, dt) ? "ok" : "FAILED");
  }
  return nr.ok ? 0 : 1;
}

int main(int argc, char** argv) {
  if (argc < 2) return usage();
  const std::string cmd = argv[1];
  if (cmd == "world") return cmdWorld(argc - 1, argv + 1);
  if (cmd == "sim") return cmdSim(argc - 1, argv + 1);
  if (cmd == "fe" || cmd == "frontend") return cmdFrontend(argc - 1, argv + 1);
  if (cmd == "estimate" || cmd == "slam") return cmdEstimate(argc - 1, argv + 1);
  if (cmd == "navigate" || cmd == "nav") return cmdNavigate(argc - 1, argv + 1);
  if (cmd == "-h" || cmd == "--help" || cmd == "help") return usage();
  std::fprintf(stderr, "unknown command '%s'\n", cmd.c_str());
  return usage();
}
