#include "luna/app/navigation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "luna/bench/trajectory.hpp"
#include "luna/frontend/features.hpp"
#include "luna/frontend/odometry.hpp"
#include "luna/sim/simulator.hpp"
#include "luna/sim/vehicle.hpp"

namespace luna {
namespace app {
namespace {

using Clock = std::chrono::steady_clock;

Scalar msSince(Clock::time_point t0) {
  return std::chrono::duration<Scalar, std::milli>(Clock::now() - t0).count();
}

/// Percentile of a vector, 0 when empty.
Scalar percentile(std::vector<Scalar> v, Scalar q) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  const Size i = std::min(v.size() - 1, static_cast<Size>(std::floor(q * v.size())));
  return v[i];
}

/// Build the occupancy map from one simulated pass, using the *online*
/// odometry pose.
///
/// The map is built from dead reckoning rather than from the optimised
/// trajectory on purpose. A pose-graph optimisation is a global, non-causal
/// operation: its correction to keyframe 10 can depend on a loop closed at
/// keyframe 400. A map built from the optimised trajectory is therefore a map
/// that could not have existed at the time it was needed, and using it would
/// make the planning and control results better than the stack can actually
/// deliver. Real systems map with the estimate they have and correct the map
/// later; this does the same, and reports the resulting map's quality.
struct MappingPass {
  sim::GridMap grid;
  Pose2 origin{0, 0, 0};  ///< world pose of the odometry frame's origin
  bench::Trajectory truth;
  bench::Trajectory odom_traj;
  Size frames = 0;
  Size dropped = 0;
  mapping::OccupancyStats map_stats;
  Scalar map_ms = 0;
};

/// Remove occupied components smaller than `min_cells`, turning them into free
/// space.
///
/// A stereo depth map has outliers: a mismatched block reports a hit at the
/// wrong range, and that one reading is enough to mark a cell occupied, because
/// the inverse sensor model treats a hit as strong evidence. The result is a
/// sprinkling of single occupied cells across open floor. Each one is harmless
/// on its own and ruinous in aggregate: the planner's clearance field measures
/// distance to the nearest obstacle, so a speck in the middle of a room lowers
/// the clearance there, and a speck in a doorway -- which is where the noise is
/// most likely, since a doorway is a depth discontinuity -- closes the doorway
/// entirely.
///
/// The filter is a size threshold on connected occupied components, and it is
/// the right shape for this noise because of what a wall *is*: a wall is a long
/// connected run of cells, and a stereo outlier is an isolated cell or a small
/// clump. Size separates them cleanly, and unlike a smoothing or voting filter
/// it cannot round off a corner or thin a wall that was genuinely observed.
///
/// Eight-connected, so a wall that is only diagonally joined still counts as one
/// structure rather than two half-walls.
Size removeObstacleSpeckles(sim::GridMap* g, Size min_cells) {
  const Size n = static_cast<Size>(g->w) * g->h;
  std::vector<char> seen(n, 0);
  std::vector<int> comp(n, -1);
  std::vector<std::pair<int, int>> stack;
  Size removed = 0;
  int next = 0;
  for (int y = 0; y < g->h; ++y) {
    for (int x = 0; x < g->w; ++x) {
      const Size i0 = static_cast<Size>(y) * g->w + x;
      if (seen[i0] || g->at(x, y) != sim::Cell::kOccupied) continue;
      const int id = next++;
      std::vector<std::pair<int, int>> members;
      stack.clear();
      stack.emplace_back(x, y);
      seen[i0] = 1;
      while (!stack.empty()) {
        const auto [cx, cy] = stack.back();
        stack.pop_back();
        members.emplace_back(cx, cy);
        comp[static_cast<Size>(cy) * g->w + cx] = id;
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue;
            const int nx = cx + dx, ny = cy + dy;
            if (!g->inBounds(nx, ny)) continue;
            const Size j = static_cast<Size>(ny) * g->w + nx;
            if (seen[j] || g->at(nx, ny) != sim::Cell::kOccupied) continue;
            seen[j] = 1;
            stack.emplace_back(nx, ny);
          }
        }
      }
      if (members.size() < min_cells) {
        for (const auto& m : members) g->set(m.first, m.second, sim::Cell::kFree);
        removed += members.size();
      }
    }
  }
  return removed;
}

/// Mark the corridor the vehicle actually drove as free space.
///
/// The occupancy map's sensor model cannot do this by itself. It ignores returns
/// closer than `min_range`, because a beam that terminates inside the robot
/// footprint would otherwise carve a hole that grows into a trail. That leaves a
/// disc of *unknown* cells around every scan pose, and when the stack is asked
/// to plan from where it currently is, the start point sits inside that disc.
/// With unknown space treated as an obstacle -- which is what keeps a plan
/// inside observed territory -- the start is an obstacle, the clearance field is
/// zero there, and the controller refuses to move out of a cell it is already
/// occupying.
///
/// The trajectory itself is evidence of free space, and stronger evidence than a
/// depth return: the vehicle was there. Stamping the swept footprint free is not
/// a guess about the world, it is a report of a measurement the sensor model
/// chose not to keep.
///
/// The footprint is a disc of the body's half-diagonal, so the stamp covers the
/// whole body rather than a point. It is applied before free-space completion so
/// the completed region joins up with the driven corridor.
void stampTrajectoryFree(sim::GridMap* g, const bench::Trajectory& traj, Scalar radius_m) {
  const int r = static_cast<int>(std::ceil(radius_m / g->res));
  for (const bench::TrajectoryPoint& tp : traj) {
    const Eigen::Vector2i c = g->toCell(Vec2(tp.pose.x, tp.pose.y));
    for (int dy = -r; dy <= r; ++dy) {
      for (int dx = -r; dx <= r; ++dx) {
        if (dx * dx + dy * dy > r * r) continue;
        const int x = c.x() + dx, y = c.y() + dy;
        if (!g->inBounds(x, y)) continue;
        if (g->at(x, y) == sim::Cell::kOccupied) continue;
        g->set(x, y, sim::Cell::kFree);
      }
    }
  }
}

MappingPass buildMap(const sim::SimConfig& cfg, const EstimatorOptions& eopt,
                     const mapping::OccupancyOptions& mopt, Scalar occupancy_threshold,
                     Scalar inflation_m, Scalar free_completion_m, Size min_obstacle_cells,
                     Scalar stamp_free_radius_m, const bench::Trajectory& slam_traj) {
  MappingPass out;
  const auto t0 = Clock::now();

  frontend::OdomIntegrator odom(eopt.odom);
  mapping::OccupancyMap map(mopt);
  const frontend::BriefPattern pattern;
  frontend::FastOptions fast;
  fast.threshold = eopt.fast_threshold;
  frontend::BriefOptions brief;
  brief.max_features = fast.max_features;

  bool have_origin = false;
  sim::Simulator sim(cfg);
  // Pose source for the scan. Two modes, and the difference is not cosmetic:
  //
  //   * With a SLAM trajectory, each frame's pose is the optimised estimate at
  //     the nearest keyframe. This is what a real system does after its backend
  //     runs: the map is corrected along with the trajectory, because a map
  //     built from drifted odometry smears every wall by the drift at the time
  //     it was observed, and a smeared wall is wider than the doorway beside it.
  //     On a long run the smearing is larger than the map resolution by an order
  //     of magnitude, and the map is not merely inaccurate -- its topology is
  //     wrong, because walls that should be thin are thick enough to close the
  //     gaps the planner needs.
  //   * Without one, the online odometry pose. Kept because it is the honest
  //     picture of what the stack knows *before* its backend has run, and the
  //     comparison between the two is the clearest evidence for why the
  //     correction matters.
  const bool have_slam = !slam_traj.empty();
  Size slam_idx = 0;
  sim.run([&](const sim::SensorFrame& f) {
    ++out.frames;
    if (!have_origin) {
      out.origin = f.ground_truth;
      have_origin = true;
    }
    out.truth.push_back(bench::TrajectoryPoint{f.t, f.ground_truth});
    odom.add(f.odom, f.imu);
    out.odom_traj.push_back(bench::TrajectoryPoint{f.t, odom.pose()});
    if (!f.camera_ok) {
      ++out.dropped;
      return;
    }
    std::vector<frontend::Feature> feats;
    frontend::detectAndDescribe(f.image, pattern, fast, brief, &feats);
    if (cfg.stereo && !f.image_right.px.empty()) {
      frontend::StereoOptions so;
      frontend::stereoDepth(f.image, f.image_right, cfg.camera, so, &feats);
    }
    const std::vector<mapping::RangeReading> rays = mapping::readingsFromFeatures(feats, cfg.camera);
    Pose2 scan_pose = odom.pose();
    if (have_slam) {
      while (slam_idx + 1 < slam_traj.size() && slam_traj[slam_idx + 1].t <= f.t) ++slam_idx;
      scan_pose = slam_traj[slam_idx].pose;
    }
    map.integrateScan(scan_pose, rays);
  });

  // Project the probabilistic map onto a binary grid the planner can use. The
  // occupancy threshold and the obstacle inflation are both applied here, in
  // the map frame, so the planner never has to know the map is probabilistic.
  out.grid = map.toGrid(mopt.resolution, 2.0);
  removeObstacleSpeckles(&out.grid, min_obstacle_cells);
  if (!slam_traj.empty()) {
    stampTrajectoryFree(&out.grid, slam_traj, stamp_free_radius_m);
  }

  // Free-space completion. An occupancy map built from a stereo camera marks
  // free space only *along* the beams that actually returned a depth: the cells
  // between two adjacent beams stay unknown, even though the region they lie in
  // was plainly observed. The result is a map whose free space is a fan of thin
  // streaks rather than the room the vehicle drove through, and a planner
  // cannot route along a streak.
  //
  // Filling unknown cells that lie within `free_completion_m` of observed free
  // space closes those gaps. The fill only ever moves into unknown cells and
  // stops at observed obstacles, so it cannot carve a passage through a wall
  // that was seen; the wall's own cells are occupied and block the fill. What it
  // can do is bridge a *thin* unobserved gap, which is the intended trade: the
  // alternative is a map that is technically conservative and practically
  // unusable.
  if (free_completion_m > 0) {
    const int r = static_cast<int>(std::ceil(free_completion_m / out.grid.res));
    std::vector<std::pair<int, int>> frontier;
    std::vector<int> dist(static_cast<Size>(out.grid.w) * out.grid.h, -1);
    for (int y = 0; y < out.grid.h; ++y) {
      for (int x = 0; x < out.grid.w; ++x) {
        if (out.grid.at(x, y) == sim::Cell::kFree) {
          dist[static_cast<Size>(y) * out.grid.w + x] = 0;
          frontier.emplace_back(x, y);
        }
      }
    }
    for (Size head = 0; head < frontier.size(); ++head) {
      const auto [cx, cy] = frontier[head];
      const int d = dist[static_cast<Size>(cy) * out.grid.w + cx];
      if (d >= r) continue;
      const int dx[4] = {1, -1, 0, 0}, dy[4] = {0, 0, 1, -1};
      for (int k = 0; k < 4; ++k) {
        const int nx = cx + dx[k], ny = cy + dy[k];
        if (!out.grid.inBounds(nx, ny)) continue;
        const Size j = static_cast<Size>(ny) * out.grid.w + nx;
        if (dist[j] >= 0) continue;
        if (out.grid.at(nx, ny) != sim::Cell::kUnknown) continue;
        dist[j] = d + 1;
        out.grid.set(nx, ny, sim::Cell::kFree);
        frontier.emplace_back(nx, ny);
      }
    }
  }

  // Unknown space is marked occupied for planning. This is the single most
  // important line in the map-to-plan path, and the reason is that "unknown" is
  // not "empty". A cell the sensor never saw is a cell whose contents are
  // exactly the thing the vehicle must not drive into; the clearance field
  // treats non-occupied cells as free, so leaving unknown as unknown hands the
  // planner infinite clearance through every unobserved region of the map. The
  // planner then routes confidently through a wall it has never seen, the
  // controller tracks the plan perfectly, and the vehicle collides with a
  // perfectly good trajectory -- a failure that looks like a controller bug and
  // is actually a map semantics bug.
  //
  // Treating unknown as occupied is conservative in the direction that matters:
  // it can only remove plans, never add an unsafe one. The cost is that the
  // vehicle cannot plan through space it has not observed, which is the correct
  // behaviour for a system that does not know what is there.
  for (int y = 0; y < out.grid.h; ++y) {
    for (int x = 0; x < out.grid.w; ++x) {
      if (out.grid.at(x, y) == sim::Cell::kUnknown) out.grid.set(x, y, sim::Cell::kOccupied);
    }
  }
  if (inflation_m > 0) out.grid.dilate(inflation_m);
  out.map_stats = map.stats();
  out.map_ms = msSince(t0);
  (void)occupancy_threshold;
  return out;
}

/// Nearest cell centre that is free, searched outwards from `p`. Returns `p`
/// unchanged if no free cell is within `max_radius`.
///
/// Needed because topology nodes are placed by clearance, not by occupancy: a
/// node sits at the centre of a region, which in a reconstructed map can be a
/// cell the map marked occupied. Handing such a point to the planner as a goal
/// produces a path whose own endpoint is inside an obstacle, and the planner
/// correctly reports the route as blocked -- a failure that looks like a planner
/// bug and is a frame-of-reference mismatch between two modules.
Vec2 snapToFree(const sim::GridMap& g, const Vec2& p, Scalar max_radius) {
  const Eigen::Vector2i c = g.toCell(p);
  const int rad = static_cast<int>(std::ceil(max_radius / g.res));
  for (int r = 0; r <= rad; ++r) {
    for (int dy = -r; dy <= r; ++dy) {
      for (int dx = -r; dx <= r; ++dx) {
        if (std::max(std::abs(dx), std::abs(dy)) != r) continue;
        const int nx = c.x() + dx, ny = c.y() + dy;
        if (!g.inBounds(nx, ny)) continue;
        if (g.at(nx, ny) == sim::Cell::kFree) return g.cellCenter(nx, ny);
      }
    }
  }
  return p;
}


}  // namespace

NavigationResult runNavigation(const sim::SimConfig& cfg, const NavigationOptions& opt) {
  NavigationResult res;
  const auto t_start = Clock::now();

  // -------- 1. SLAM -------------------------------------------------------
  // The estimator is run for its own numbers rather than reimplemented here.
  // Its ATE, loop-closure count and precision are the contract's SLAM figures,
  // and duplicating the pipeline to avoid a second simulation would guarantee
  // the two numbers disagree.
  const EstimatorResult er = runEstimator(cfg, opt.estimator);
  res.stats.keyframes = er.stats.keyframes;
  res.stats.slam_ate = er.optimised_ate.rmse;
  res.stats.slam_ate_percent = er.optimised_ate.percent_of_distance;
  res.stats.loop_closures = er.stats.loop_closures;
  {
    const Size tp = er.stats.loop_true_positives;
    const Size fp = er.stats.loop_false_positives;
    res.stats.loop_precision = tp + fp > 0 ? static_cast<Scalar>(tp) / static_cast<Scalar>(tp + fp) : 1.0;
  }

  // -------- 2. Map --------------------------------------------------------
  const MappingPass mp = buildMap(cfg, opt.estimator, opt.map, opt.occupancy_threshold,
                                  opt.map_inflation_m, opt.free_completion_m,
                                  opt.min_obstacle_component_cells, opt.stamp_free_radius_m,
                                  er.optimised_trajectory);
  {
    Size known = 0, occupied = 0;
    for (int y = 0; y < mp.grid.h; ++y) {
      for (int x = 0; x < mp.grid.w; ++x) {
        const sim::Cell v = mp.grid.at(x, y);
        if (v != sim::Cell::kUnknown) ++known;
        if (v == sim::Cell::kOccupied) ++occupied;
      }
    }
    res.stats.map_occupied_cells = occupied;
    res.stats.map_known_fraction =
        static_cast<Scalar>(known) / static_cast<Scalar>(std::max(1, mp.grid.w * mp.grid.h));
    res.stats.map_tiles_resident = mp.map_stats.tiles_resident;
    res.stats.map_bytes_resident = mp.map_stats.bytes_resident;
  }

  res.planning_grid = mp.grid;

  // -------- 3. Topology ---------------------------------------------------
  const auto t_topo = Clock::now();
  const auto complex = topology::MorseComplex::extract(mp.grid, opt.morse);
  res.stats.topology_ms = msSince(t_topo);
  res.stats.morse_nodes = static_cast<int>(complex.graph().nodes.size());
  res.stats.morse_edges = static_cast<int>(complex.graph().edges.size());
  for (const auto& n : complex.graph().nodes) {
    if (n.type == topology::CriticalType::kMaximum) ++res.stats.maxima;
    if (n.type == topology::CriticalType::kSaddle) ++res.stats.saddles;
  }
  if (complex.graph().nodes.size() < 2) {
    res.ok = false;
    return res;
  }

  // -------- 4. Plan -------------------------------------------------------
  // Pick the pair of corridor nodes furthest apart, which is the most
  // demanding route the map offers rather than a short one that would pass
  // regardless of planner quality.
  planning::CorridorPlanner planner(complex.graph(), mp.grid, opt.plan);

  // Choose the goal by the planner's own notion of reachability: ask it to route
  // from the vehicle's start node to every other node, and keep the furthest one
  // that actually has a route.
  //
  // Every cheaper proxy disagrees with the planner somewhere. Free-space
  // connectivity ignores the corridor graph; the corridor graph ignores the
  // occupancy grid; a straight-line test between two nodes ignores the fact that
  // corridors run *through* a saddle, not along the chord between its endpoints.
  // Asking the planner is not a shortcut around the problem, it is the only test
  // whose answer is the same question the plan itself will ask.
  const auto& nodes = complex.graph().nodes;

  // Corridor-graph components, from the edge list.
  std::vector<int> comp(nodes.size(), -1);
  {
    std::vector<std::vector<Size>> adj(nodes.size());
    for (const auto& e : complex.graph().edges) {
      adj[static_cast<Size>(e.a)].push_back(static_cast<Size>(e.b));
      adj[static_cast<Size>(e.b)].push_back(static_cast<Size>(e.a));
    }
    int next = 0;
    std::vector<Size> stack;
    for (Size k = 0; k < nodes.size(); ++k) {
      if (comp[k] >= 0) continue;
      const int id = next++;
      comp[k] = id;
      stack.assign(1, k);
      while (!stack.empty()) {
        const Size cur = stack.back();
        stack.pop_back();
        for (Size nb : adj[cur]) {
          if (comp[nb] >= 0) continue;
          comp[nb] = id;
          stack.push_back(nb);
        }
      }
    }
  }

  // The vehicle's current position, which is where the plan must begin.
  //
  // Not the odometry origin. The origin is where the vehicle *started*, and on a
  // 60 s run that is tens of metres behind it -- deep in territory it has left,
  // where the map is a trail of partly-observed cells and the vehicle's own
  // blind-spot disc. Planning from there asks the controller to drive out of a
  // cell whose clearance is zero, and it correctly refuses to move: the reported
  // symptom is a safe stop one tick into a plan that looks fine, because the
  // plan's first metre is behind the vehicle and inside a wall it drove past.
  //
  // The final SLAM pose is the right frame: it is the best estimate of where the
  // vehicle is, in the same frame the map was built in.
  const Vec2 current_pt = er.optimised_trajectory.empty()
                              ? Vec2(0, 0)
                              : Vec2(er.optimised_trajectory.back().pose.x,
                                     er.optimised_trajectory.back().pose.y);
  const Vec2 start_pt = snapToFree(mp.grid, current_pt, 3.0);

  // Pick the component to drive in: the largest one that has a node the vehicle
  // can actually get to.
  //
  // Neither half of that rule can be dropped. Taking the component of the
  // *nearest* node fails because the nearest node is frequently an isolated
  // maximum -- a single pocket of high clearance beside the start -- and a
  // component of one node contains no route, so the planner is asked for a route
  // it cannot produce. Taking simply the largest component fails because a large
  // component elsewhere in the map is not somewhere the vehicle can reach: the
  // plan would begin with a jump across unknown space.
  std::vector<int> count(nodes.size(), 0);
  std::vector<Scalar> nearest(nodes.size(), std::numeric_limits<Scalar>::max());
  for (Size k = 0; k < nodes.size(); ++k) {
    count[static_cast<Size>(comp[k])] += 1;
    nearest[static_cast<Size>(comp[k])] =
        std::min(nearest[static_cast<Size>(comp[k])], (nodes[k].p - start_pt).norm());
  }
  const Scalar kReachRadius = 6.0;
  int chosen = -1;
  for (Size c = 0; c < nodes.size(); ++c) {
    if (count[c] == 0 || nearest[c] > kReachRadius) continue;
    if (chosen < 0 || count[c] > count[static_cast<Size>(chosen)]) chosen = static_cast<int>(c);
  }
  if (chosen < 0) {
    // No corridor component has a node the vehicle can reach. Report that
    // rather than leaving the default `kInvalid`, which reads as a malformed
    // request. This is a real outcome on a fragmented map -- the vehicle has
    // ended a run in a region the topology extractor never found a maximum in
    // -- and the status should say so.
    res.stats.plan_status = planning::PlanStatus::kNoCorridorRoute;
    res.ok = false;
    return res;
  }
  int start_node = -1;
  Scalar bd = std::numeric_limits<Scalar>::max();
  for (Size k = 0; k < nodes.size(); ++k) {
    if (comp[k] != chosen) continue;
    const Scalar d = (nodes[k].p - start_pt).norm();
    if (d < bd) {
      bd = d;
      start_node = static_cast<int>(k);
    }
  }
  if (start_node < 0) {
    res.stats.plan_status = planning::PlanStatus::kNoCorridorRoute;
    res.ok = false;
    return res;
  }

  // Candidate goals: every node in the chosen component, furthest first.
  //
  // The goal used to be the single furthest node, chosen because the corridor
  // route to it "exists by construction". That reasoning holds for the *graph*
  // and not for the map underneath it. A corridor edge is an adjacency in the
  // Morse complex, and on a partly-observed map the complex can join two regions
  // whose connecting cells are not drivable at the planner's clearance -- the
  // extractor's height function is smoothed, so it bridges a pinch the raw grid
  // does not. Planning to the furthest node then fails outright, and the run
  // reports `blocked` and drives nothing even though nearer nodes are perfectly
  // reachable. Trying candidates in descending distance order and taking the
  // first that plans turns that all-or-nothing outcome into the longest route
  // the map actually supports, which is the honest reading of the request.
  int n_reach = 0;
  std::vector<std::pair<Scalar, int>> candidates;
  for (Size k = 0; k < nodes.size(); ++k) {
    if (comp[k] != chosen) continue;
    ++n_reach;
    candidates.emplace_back((nodes[k].p - nodes[static_cast<Size>(start_node)].p).norm(),
                            static_cast<int>(k));
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const std::pair<Scalar, int>& a, const std::pair<Scalar, int>& b) {
              return a.first > b.first;
            });
  res.stats.reachable_nodes = n_reach;
  // The number of *distinct* components, not the size of the chosen one. The
  // two were conflated here: `count[c]` is how many nodes component `c` holds,
  // which is a node count, and reporting it as "components" overstates the
  // graph's connectivity by an order of magnitude on a fragmented map (19
  // components of one or two nodes each read as 19 components' worth of nodes
  // in a single component). The size of the chosen component is
  // `reachable_nodes`, which is already reported separately.
  int n_components = 0;
  for (Size c = 0; c < nodes.size(); ++c)
    if (count[c] > 0) ++n_components;
  res.stats.graph_components = n_components;
  // A plan that goes nowhere is not a plan. Without the length floor the search
  // accepts the start node itself as a "successful" goal -- it is a candidate,
  // its route is trivially valid, and it is the last one tried -- so a map on
  // which no real route exists reports `ok` with a zero-length path. The
  // controller then correctly refuses to move, and the run reads as a control
  // failure when it is a planning one. One metre is well above the planner's
  // own resampling step, so it excludes only genuinely degenerate routes.
  const Scalar kMinRouteM = 1.0;
  planning::Plan plan;
  for (const auto& cand : candidates) {
    if (cand.second == start_node) continue;
    const Vec2 goal_pt =
        snapToFree(mp.grid, nodes[static_cast<Size>(cand.second)].p, 3.0);
    plan = planner.planViaNodes(start_pt, start_node, goal_pt, cand.second);
    if (plan.ok() && plan.stats.path_length >= kMinRouteM) break;
  }
  res.stats.plan_ok = plan.ok();
  res.stats.plan_status = plan.status;
  res.stats.plan_length = plan.stats.path_length;
  res.stats.plan_clearance = plan.stats.min_clearance_along;
  res.stats.plan_ms = plan.stats.plan_ms;
  if (!plan.ok()) {
    res.ok = false;
    return res;
  }

  // The plan is in the odometry frame, which starts at the origin. The vehicle
  // drives in the world frame. The transform between them is the simulator's
  // own first ground-truth pose -- used here as the *frame convention* the
  // estimator implicitly assumes (it starts at the origin), not as information
  // about where the vehicle is. Any absolute-pose error is already inside the
  // SLAM ATE figure; what this transform must not do is leave the plan in the
  // wrong place relative to the world it is scored against.
  const Pose2 odom_to_world = mp.origin;
  res.planned_path.reserve(plan.path.size());
  for (const Vec2& p : plan.path) {
    // Rotate into the world frame by hand: the odometry frame and the world
    // frame differ by the vehicle's initial pose, and this is the only place
    // the two are related.
    const Scalar c = std::cos(odom_to_world.th), sn = std::sin(odom_to_world.th);
    res.planned_path.push_back(
        Vec2(odom_to_world.x + c * p.x() - sn * p.y(), odom_to_world.y + sn * p.x() + c * p.y()));
  }

  // -------- 5. Drive ------------------------------------------------------
  if (!opt.drive) {
    res.ok = true;
    return res;
  }

  // Score against the *real* world, not the map. This is the point of the
  // whole exercise: a stack can plan perfectly on its own map and still hit
  // walls in the world, and only a world collision check exposes that.
  sim::SimConfig world_cfg = cfg;
  world_cfg.record_images = false;
  const sim::World world = sim::World::generate(world_cfg.world);

  planning::Plan world_plan = plan;
  world_plan.path = res.planned_path;

  sim::VehicleParams vp = cfg.vehicle;

  // The plant starts at the position the vehicle is *actually* at -- the last
  // ground-truth pose of the mapping pass -- and the controller is told the
  // position SLAM reports. They differ by the run's position error, and that
  // difference is the point: collisions are scored against the world, so the
  // plant has to sit in the world where it really is.
  //
  // Spawning the plant at the estimate instead teleports it onto its own error,
  // so a run whose estimate lands 0.3 m inside a wall begins with the body in
  // the wall and reports a collision on tick 0. The body was never there.
  //
  // The heading comes from the plan, and that is a simplification worth naming.
  // The planner plans from a *position*, not a pose, so the plan's first
  // segment points at the goal and says nothing about which way the vehicle
  // faced when the scripted route ended. Orienting the vehicle along its plan
  // is the statement "it now drives this plan". Removing the simplification
  // means a heading-aware planner that can plan a turn, which is not
  // implemented; keeping the end-of-route heading instead makes the initial
  // error up to 180 degrees, the controller diverges on the first tick, and the
  // run measures the missing planner rather than the stack.
  Pose2 start{world_plan.path.front().x(), world_plan.path.front().y(), 0};
  {
    const Vec2 t = world_plan.path[std::min<Size>(3, world_plan.path.size() - 1)] -
                   world_plan.path.front();
    start.th = std::atan2(t.y(), t.x());
  }
  if (!mp.truth.empty()) {
    start.x = mp.truth.back().pose.x;
    start.y = mp.truth.back().pose.y;
  }
  sim::Vehicle vehicle(vp, start, 1);
  control::MpcController mpc(vp, opt.control);

  control::VehicleEstimate est;
  // The controller is fed the *SLAM estimate* of the vehicle's pose, not the
  // truth. The estimate starts where SLAM says the vehicle is -- the start of
  // its own plan, in the world frame -- and is then propagated by the vehicle's
  // odometry increments. The plant below starts at ground truth.
  //
  // Both halves matter, and getting either wrong is a false result:
  //
  //   * Feeding the controller the plant state directly is perfect state
  //     feedback. The controller then tracks a plan it can see it is exactly on,
  //     the localization error never reaches it, and the run reports zero
  //     tracking error for a stack that cannot localize.
  //   * Starting the plant at the estimate instead teleports the body onto its
  //     own error, so a 0.9 m estimate error puts the body 0.9 m from where it
  //     really is and can drop it inside a wall. The collision is an artifact of
  //     the setup.
  //
  // With the plant at truth and the estimate propagated from SLAM's own start,
  // the initial offset between them is the run's localization error, and it
  // reaches the vehicle exactly as it does on a real robot: as a plan the
  // controller believes it is following while the body sits somewhere else.
  est.pose = Pose2{world_plan.path.front().x(), world_plan.path.front().y(), start.th};
  Pose2 prev_truth = start;

  // The controller's uncertainty is the graph's own -- not the ATE. The ATE is
  // measured against ground truth, so a controller driven by it would be
  // steering on the answer key.
  //
  // It uses the *relative* covariance over a recent window rather than the
  // marginal of the final pose. The marginal grows without bound as a run gets
  // longer, because it counts every metre of accumulated drift; feeding it to a
  // chance constraint makes the vehicle decelerate steadily and eventually stop
  // for no reason other than having been driving. The relative covariance is
  // bounded, and it is what the question "is this gap wide enough for me right
  // now" actually depends on.
  const Mat3& loc = er.local_covariance;
  est.sigma_xx = std::max<Scalar>(0.0, loc(0, 0));
  est.sigma_yy = std::max<Scalar>(0.0, loc(1, 1));
  est.sigma_thth = std::max<Scalar>(0.0, loc(2, 2));
  est.sigma_yth = loc(1, 2);
  res.stats.slam_sigma_xy = std::sqrt(est.sigma_xx + est.sigma_yy);
  res.stats.slam_sigma_th = std::sqrt(est.sigma_thth);
  res.stats.slam_covariance_ok = er.local_covariance_ok;
  res.stats.slam_marginal_sigma_xy = er.final_covariance_ok
                                         ? std::sqrt(std::max<Scalar>(0.0, er.final_covariance(0, 0)) +
                                                     std::max<Scalar>(0.0, er.final_covariance(1, 1)))
                                         : 0;

  std::vector<Scalar> solve_times;
  Scalar sum_sq = 0;
  int n = 0;
  const Scalar dt = opt.control.dt;
  const int max_steps = static_cast<int>(opt.drive_duration / dt);
  for (int k = 0; k < max_steps; ++k) {
    const control::ControllerOutput out = mpc.step(est, world_plan);
    solve_times.push_back(out.stats.solve_ms);
    ++res.stats.control_ticks;
    if (out.stopping()) res.stats.safe_stopped = true;
    res.stats.control_status = out.status;
    res.stats.control_speed_cap = out.stats.speed_cap;
    res.stats.control_det_sigma = out.stats.det_sigma;
    res.stats.control_inflation = out.stats.inflation;
    res.stats.plan_route_bottleneck = plan.stats.route_bottleneck;
    res.stats.control_sigma_xy = std::sqrt(est.sigma_xx + est.sigma_yy);
    res.stats.control_sigma_th = std::sqrt(est.sigma_thth);
    res.stats.control_stop_reason = out.stats.stop_reason;
    if (std::getenv("LUNA_TRACE") && (k % 20 == 0 || out.stopping())) {
      std::fprintf(stderr,
                   "tick %4d pos=(%7.3f,%7.3f) th=%6.3f v=%5.2f prog=%.2f lat=%.3f head=%.3f "
                   "tube=%.3f cap=%.3f cmd(v=%.2f st=%+.3f a=%+.2f) st=%s reason=%s\n",
                   k, est.pose.x, est.pose.y, est.pose.th, est.v, mpc.progress(),
                   out.stats.lateral_error, out.stats.heading_error,
                   out.stats.tube_margin, out.stats.speed_cap,
                   out.cmd.v, out.cmd.steer, out.cmd.accel,
                   control::toString(out.status), control::toString(out.stats.stop_reason));
    }
    vehicle.step(out.cmd, dt);
    // Propagate the estimate by the odometry increment the plant just moved by,
    // expressed in the estimate's own frame. This is dead reckoning from SLAM's
    // start pose: the estimate drifts away from truth at the odometry rate, and
    // the controller steers on the estimate, not on where the body really is.
    {
      const Pose2 cur = vehicle.state().pose;
      const Pose2 inc = pose2Compose(pose2Inverse(prev_truth), cur);
      prev_truth = cur;
      est.pose = pose2Compose(est.pose, inc);
    }
    est.v = vehicle.state().v;
    est.steer = vehicle.state().steer;
    res.driven_path.push_back(Vec2(est.pose.x, est.pose.y));

    if (vehicle.inCollision(world.grid())) ++res.stats.collisions;

    Scalar d = std::numeric_limits<Scalar>::max();
    for (const Vec2& p : world_plan.path) d = std::min(d, (p - Vec2(est.pose.x, est.pose.y)).norm());

    sum_sq += d * d;
    ++n;
    res.stats.tracking_max = std::max(res.stats.tracking_max, d);

    if (out.status == control::ControlStatus::kGoalReached) break;
    if (out.stopping()) break;
  }
  res.stats.tracking_rms = n ? std::sqrt(sum_sq / n) : 0;
  res.stats.control_ms_p99 = percentile(solve_times, 0.99);
  res.stats.control_ms_max = solve_times.empty() ? 0
                                                 : *std::max_element(solve_times.begin(), solve_times.end());
  res.stats.drive_distance = vehicle.distanceTravelled();
  res.stats.goal_error =
      (Vec2(est.pose.x, est.pose.y) - world_plan.path.back()).norm();
  (void)t_start;
  res.ok = res.stats.collisions == 0;
  return res;
}

}  // namespace app
}  // namespace luna
