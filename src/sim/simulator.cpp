#include "luna/sim/simulator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <deque>
#include <fstream>

namespace luna {
namespace sim {

namespace {

/// Order a set of waypoints into a short tour by nearest-neighbour, then
/// close the loop. Good enough to script a repeatable route; the point is
/// coverage and determinism, not optimality.

/// Greedy nearest-neighbour ordering of the sampled anchors into a closed tour.
std::vector<Vec2> orderTour(const std::vector<Vec2>& pts, const Vec2& start) {
  if (pts.empty()) return {};
  std::vector<Vec2> remaining = pts;
  std::vector<Vec2> tour;
  Vec2 cur = start;
  while (!remaining.empty()) {
    Size best = 0;
    Scalar bestd = std::numeric_limits<Scalar>::max();
    for (Size i = 0; i < remaining.size(); ++i) {
      const Scalar d = (remaining[i] - cur).squaredNorm();
      if (d < bestd) {
        bestd = d;
        best = i;
      }
    }
    cur = remaining[best];
    tour.push_back(cur);
    remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(best));
  }
  tour.push_back(start);

  return tour;
}

}  // namespace

std::vector<Vec2> buildRoute(const World& world, Rng* rng, Scalar clearance_m,
                             Scalar min_clearance_m, Scalar min_turn_radius) {
  // Sample free space broadly, then chain the samples into a closed tour and
  // densify each leg into a collision-free grid path.
  auto anchors =
      sampleFreePoints(world.grid(), 10, 4.0, rng, clearance_m * 0.75, &world.clearance());
  // Drop anchors near the world border. The tour has to reverse direction at
  // the end of a dead end, and a dead end against the border gives the vehicle
  // nowhere to complete the turn: it swings into the wall. Anchors further in
  // leave room for the turn.
  const Scalar inset = std::max<Scalar>(1.5, 3.0 * clearance_m);
  const Scalar lim = world.spec().size_m - inset;
  anchors.erase(std::remove_if(anchors.begin(), anchors.end(),
                               [&](const Vec2& p) {
                                 return p.x() < inset || p.y() < inset || p.x() > lim || p.y() > lim;
                               }),
                anchors.end());
  if (anchors.size() < 4) {
    // Fall back to the unfiltered set rather than collapsing the route.
    anchors = sampleFreePoints(world.grid(), 10, 4.0, rng, clearance_m * 0.75,
                               &world.clearance());
  }
  if (anchors.empty()) return {};
  const Vec2 start = anchors.front();
  const auto tour = orderTour(anchors, start);

  // Plan the tour at a margin, returning an empty route when too few legs can
  // be planned to make a usable tour.
  const Scalar min_turn_radius_ = min_turn_radius;
  const auto attempt = [&](Scalar clearance) {
    std::vector<Vec2> route;
    Vec2 cur = start;
    int planned = 0;
    for (const Vec2& next : tour) {
      const auto leg = world.gridPath(cur, next, clearance);
      if (leg.empty()) continue;
      ++planned;
      // Smooth the staircase corners before shortening. A raw 4-connected grid
      // path turns in 90-degree steps, which no finite-turn-radius vehicle can
      // follow, so the tracker cuts every corner and collides even though the
      // path itself is clear. Smoothing first makes the polyline drivable, then
      // shortcutting removes the points smoothing left collinear.
      const auto smooth = world.smoothPath(leg, clearance, min_turn_radius_);
      const auto tight = world.shortcut(smooth, clearance);
      for (const Vec2& p : tight) {
        // Skip points that repeat the previous one. A duplicate is not
        // harmless: it creates a zero-length segment, and the initial heading
        // is taken from the first segment, so a duplicated start makes the
        // heading atan2(0, 0) = 0 and the vehicle leaves in an arbitrary
        // direction it can never correct from.
        if (!route.empty() && (route.back() - p).squaredNorm() < 1e-12) continue;
        route.push_back(p);
      }
      cur = next;
    }
    return std::make_pair(std::move(route), planned);
  };

  // Relax the margin from the preferred value down to the floor, and keep the
  // *first* margin that connects enough of the tour. Planning at a single
  // margin is not robust: one narrow pinch fails every leg through it and the
  // tour collapses to its start point, which the vehicle then leaves in a
  // straight line into a wall. But relaxing too eagerly is just as wrong: the
  // narrowest route is the one most likely to be undrivable, so the search
  // stops at the first margin that works rather than at the last.
  const int want = static_cast<int>(std::ceil(0.6 * static_cast<Scalar>(tour.size())));
  const Scalar floor_m = std::min(clearance_m, std::max(min_clearance_m, clearance_m * 0.3));
  std::vector<Vec2> best_route;
  int best_planned = -1;
  for (Scalar c = clearance_m; c >= floor_m - 1e-9; c = std::max(floor_m, c * 0.8)) {
    auto [route, planned] = attempt(c);
    if (planned > best_planned) {
      best_planned = planned;
      best_route = std::move(route);
    }
    if (planned >= want) break;
    if (c <= floor_m + 1e-9) break;
  }
  return best_route;
}

Simulator::Simulator(const SimConfig& config) : config_(config) {
  world_ = World::generate(config.world);
  renderer_ = std::make_unique<Renderer>(world_, config.camera, config.camera_noise,
                                        config.world.seed ^ 0xA5A5A5A5ull);
  imu_ = std::make_unique<ImuEmulator>(config.imu_noise, config.imu_rate, config.world.seed ^ 0x1234ull);
  odom_ = std::make_unique<OdomEmulator>(config.odom_noise, config.world.seed ^ 0x5678ull);
  faults_ = std::make_unique<FaultInjector>(config.faults);
}

SimResult Simulator::run(const FrameSink& frame_sink) {
  const auto wall0 = std::chrono::steady_clock::now();
  SimResult res;
  res.world_stats = world_.stats();

  Rng rng(config_.world.seed ^ 0xDEADBEEFull);
  // The route margin must exceed the vehicle's half-diagonal footprint, or a
  // clean run collides at wall corners. The preferred margin is deliberately
  // wide: planning at a wide margin keeps the route in the middle of the open
  // space and away from the thin slits the generator leaves between walls,
  // which are free but far too narrow to drive. Capping the margin by the
  // world's narrowest clearance looks prudent but does the opposite: that
  // minimum is set by those same irrelevant slits, so capping it lets the
  // planner route straight through them.
  const Scalar half_diag = 0.5 * std::sqrt(config_.vehicle.length * config_.vehicle.length +
                                          config_.vehicle.track_width * config_.vehicle.track_width);
  const Scalar drive_floor = half_diag + 0.12;
  // Smooth to slightly more than the vehicle's physical minimum turn radius.
  // Smoothing to much more than that is counter-productive: a doorway forces a
  // turn, and if the target radius exceeds what the doorway can hold, smoothing
  // stops at the clearance limit and leaves the corner as sharp as it found it.
  const Scalar min_turn_radius =
      1.15 * config_.vehicle.wheelbase / std::tan(config_.vehicle.max_steer);
  res.route = buildRoute(world_, &rng, std::max<Scalar>(0.75, drive_floor), drive_floor,
                         min_turn_radius);


  // Start at the route origin, or a free cell if no route exists.
  Pose2 start{2.0, 2.0, 0.0};
  if (res.route.size() >= 2) {
    start.x = res.route.front().x();
    start.y = res.route.front().y();
    // Use the first segment long enough to define a direction. A degenerate
    // leading segment would otherwise fix the heading to atan2(0, 0) = 0.
    for (Size i = 1; i < res.route.size(); ++i) {
      const Vec2 dir = res.route[i] - res.route[0];
      if (dir.squaredNorm() > 1e-9) {
        start.th = std::atan2(dir.y(), dir.x());
        break;
      }
    }
  } else if (!res.route.empty()) {
    start.x = res.route.front().x();
    start.y = res.route.front().y();
  } else {
    const auto pts = sampleFreePoints(world_.grid(), 1, 0.1, &rng);
    if (!pts.empty()) {
      start.x = pts[0].x();
      start.y = pts[0].y();
    }
  }

  Vehicle vehicle(config_.vehicle, start, config_.world.seed);
  const Scalar ctrl_dt = 1.0 / config_.control_rate;
  const Scalar cam_dt = 1.0 / config_.camera_rate;
  const Scalar odom_dt = 1.0 / config_.odom_rate;
  Scalar next_cam = 0.0;
  Scalar next_ctrl = 0.0;
  Scalar next_odom = 0.0;
  Scalar t = 0.0;
  Control active_cmd;
  Pose2 last_odom_pose = start;
  // Odometry ticks faster than the camera. The deltas from every tick inside a
  // frame interval are accumulated here and handed over as one increment: a
  // frame that only received the most recent tick would see roughly half the
  // distance travelled, which shows up as a dead-reckoning scale error of
  // almost exactly 2x and looks like an estimator bug rather than a harness
  // bug.
  Scalar odom_acc_forward = 0, odom_acc_yaw = 0;
  bool odom_acc_valid = true;
  const Scalar sim_dt = 1.0 / 200.0;

  SensorFrame frame;
  frame.image.resize(config_.camera.width, config_.camera.height);
  if (config_.stereo) frame.image_right.resize(config_.camera.width, config_.camera.height);

  // Actuator lag is implemented as a short command delay queue.
  std::deque<Control> cmd_queue;
  bool was_colliding = false;

  // Track progress along the route by arc length rather than by waypoint index.
  //
  // Index-based pursuit breaks at any reversal in the path: the nearest
  // waypoint at least `lookahead` away can sit *behind* the vehicle, so the
  // controller steers backwards into the wall it just rounded. Projecting the
  // vehicle onto the polyline, advancing a monotone arc-length cursor, and
  // aiming a fixed distance further along makes the target advance with the
  // vehicle by construction, so a reversal cannot invert the steering.
  std::vector<Scalar> route_s(res.route.size(), 0.0);
  for (Size i = 1; i < res.route.size(); ++i)
    route_s[i] = route_s[i - 1] + (res.route[i] - res.route[i - 1]).norm();
  const Scalar route_len = route_s.empty() ? 0.0 : route_s.back();
  Scalar s_cursor = 0.0;

  // Point at arc length `s`, interpolated within its segment.
  const auto point_at = [&](Scalar s) {
    s = clamp(s, 0.0, route_len);
    Size lo = 0, hi = route_s.size() - 1;
    while (lo + 1 < hi) {
      const Size mid = (lo + hi) / 2;
      if (route_s[mid] <= s) lo = mid; else hi = mid;
    }
    const Scalar seg = route_s[hi] - route_s[lo];
    const Scalar t = seg > 1e-9 ? (s - route_s[lo]) / seg : 0.0;
    return Vec2(res.route[lo] + (res.route[hi] - res.route[lo]) * t);
  };

  while (t < config_.duration) {
    const VehicleState& vs = vehicle.state();

    // -------- control at the control rate --------------------------------
    // The plant integrates at the inner rate below; the command is only
    // refreshed here, which is what makes camera motion between frames equal
    // to the physical distance actually travelled.
    if (t + 1e-9 >= next_ctrl) {
      Control cmd;
      Scalar speed_cap = config_.target_speed;
      bool route_done = false;
      const std::vector<Vec2>& path = res.route;
      if (path.size() >= 2) {
        // Pure pursuit cuts the inside of every corner by roughly
        // lookahead^2 / (2 * turn_radius). The vehicle's body half-width must
        // stay inside that cut, so the lookahead has to shrink when the local
        // clearance is small. A lookahead tuned for open rooms cuts 0.5 m off a
        // corner, which is wider than a 0.9 m corridor's usable half-width and
        // puts the body in the wall.
        const Scalar clr = world_.clearanceAtWorld(Vec2(vs.pose.x, vs.pose.y));
        const Scalar lookahead = clamp(0.20 + 0.50 * clr, 0.30, 0.90);
        const Vec2 here(vs.pose.x, vs.pose.y);

        // Advance the arc-length cursor to the vehicle's projection on the
        // path, searching forward only so progress can never go backwards.
        Scalar best_s = s_cursor, best_d = kInf;
        for (Size i = 0; i + 1 < res.route.size(); ++i) {
          if (route_s[i + 1] < s_cursor - 1e-9) continue;
          const Vec2 a = res.route[i], ab = res.route[i + 1] - a;
          const Scalar l2 = ab.squaredNorm();
          if (l2 < 1e-12) continue;
          const Scalar t = clamp((here - a).dot(ab) / l2, 0.0, 1.0);
          const Scalar d = (a + ab * t - here).norm();
          const Scalar s = route_s[i] + t * (route_s[i + 1] - route_s[i]);
          if (s < s_cursor - 1e-9) continue;
          if (d < best_d) {
            best_d = d;
            best_s = s;
          }
          // Only look a short distance ahead. A route that doubles back on
          // itself can otherwise project the vehicle onto a *later* pass that
          // happens to be nearer, and the controller then steers toward a point
          // behind the vehicle. The window is generous compared with the
          // per-step travel of v / control_rate.
          if (route_s[i] > s_cursor + 1.5) break;
        }
        s_cursor = std::max(s_cursor, best_s);
        // The route is a closed tour, so finishing it means arriving back where
        // the vehicle started. Hold station there. Without this the projection
        // clamps to the final waypoint, the cross-track error keeps growing as
        // the vehicle drives past it, and the steering saturates into a circle
        // that carries the vehicle off the map -- which silently destroys the
        // very revisit the loop-closure stage is supposed to detect.
        route_done = s_cursor >= route_len - 0.05;

        // Curvature-based speed scheduling, previewed far enough ahead to slow
        // down *before* the corner rather than in it. Sampling the curvature
        // only over the lookahead means the cap reacts when the corner is
        // already underneath the vehicle, which is too late: the vehicle enters
        // at full speed, overshoots, and leaves the route. The preview distance
        // is the braking distance from the current speed, so the speed cap is
        // applied early enough to be achievable.
        const auto curvature_over = [&](Scalar span) {
          const Vec2 a = point_at(s_cursor);
          const Vec2 b = point_at(s_cursor + span);
          const Vec2 c = point_at(s_cursor + 2.0 * span);
          const Scalar cross = std::abs((b - a).x() * (c - b).y() - (b - a).y() * (c - b).x());
          const Scalar denom = (b - a).norm() * (b - c).norm();
          return denom > 1e-9 ? 2.0 * cross / denom : 0.0;
        };
        const Scalar a_lat_max = 0.7;
        const Scalar brake_dist =
            vs.v * vs.v / std::max<Scalar>(2.0 * config_.vehicle.max_accel, 1e-3);
        // The preview must cover the distance actually travelled while the
        // speed command is being reduced, which is dominated by the tracker's
        // own lag rather than by the brakes. A preview of one lookahead only
        // sees a corner once the vehicle is already in it.
        const Scalar preview = std::max<Scalar>(3.0, brake_dist + 3.0 * lookahead);
        const Scalar kappa = std::max({curvature_over(0.4 * preview),
                                       curvature_over(0.7 * preview),
                                       curvature_over(preview)});
        if (kappa > 1e-6) speed_cap = std::min(speed_cap, std::sqrt(a_lat_max / kappa));

        // Clearance-based speed scaling, the same idea a safety-gated planner
        // uses: the tighter the space, the less speed the tracker may carry.
        // This is not cosmetic. The command is only refreshed at the control
        // rate, so the vehicle travels v / control_rate between updates; at
        // 1.2 m/s and 20 Hz that is 6 cm of uncommanded motion inside a
        // corridor whose usable half-width is under 10 cm. Slowing down is what
        // makes the remaining margin survive a discretely-updated controller.
        // Scale speed by the room left *over the body*, not by the raw
        // clearance: below half_diag the body already overlaps the wall, so the
        // margin that matters is clr - half_diag.
        const Scalar v_clear = 0.25 + 1.2 * std::max<Scalar>(0.0, clr - half_diag);
        speed_cap = std::min(speed_cap, std::max<Scalar>(0.25, v_clear));
        // Stanley steering: heading error plus a cross-track correction,
        // plus a curvature feed-forward.
        //
        // Pure pursuit was replaced here because it converges to the path only
        // asymptotically: it always aims at a point ahead, so it lags on entry
        // to a turn and cuts the inside of it. In a 0.9 m corridor that cut is
        // the wall. Stanley drives the cross-track error to zero directly, and
        // the feed-forward removes the steady-state lag on a constant curve, so
        // the body stays on the centreline where the clearance actually is.
        const Scalar s_here = s_cursor;
        Size seg = 0;
        while (seg + 2 < route_s.size() && route_s[seg + 1] < s_here) ++seg;
        Vec2 tangent = res.route[seg + 1] - res.route[seg];
        if (tangent.norm() > 1e-9) tangent.normalize();
        const Vec2 on_path = point_at(s_here);
        const Vec2 offset = here - on_path;
        // Positive when the vehicle sits to the left of the path tangent.
        const Scalar e_ct = tangent.x() * offset.y() - tangent.y() * offset.x();
        const Scalar psi_e = wrapPi(std::atan2(tangent.y(), tangent.x()) - vs.pose.th);
        const Scalar k_ct = 2.2;
        const Scalar v_ref = std::max(vs.v, 0.4);
        cmd.steer = clamp(psi_e - std::atan2(k_ct * e_ct, v_ref) +
                              std::atan(config_.vehicle.wheelbase * kappa),
                          -config_.vehicle.max_steer, config_.vehicle.max_steer);
      }
      if (route_done) {
        speed_cap = 0;
        cmd.steer = 0;
      }
      cmd.v = speed_cap;
      cmd.accel = clamp(cmd.v - vs.v, -config_.vehicle.max_accel, config_.vehicle.max_accel);
      active_cmd = cmd;
      const Scalar lag = faults_->actuatorLag(t);
      if (lag > 0) {
        cmd_queue.push_back(cmd);
        const Size delay_steps = std::max<Size>(1, static_cast<Size>(lag / ctrl_dt));
        if (cmd_queue.size() > delay_steps) {
          active_cmd = cmd_queue.front();
          cmd_queue.pop_front();
        } else {
          active_cmd = Control{};
        }
      }
      next_ctrl += ctrl_dt;
    }

    // -------- plant integration at the inner rate -------------------------
    vehicle.step(active_cmd, sim_dt);
    const VehicleState& vs2 = vehicle.state();
    // Count collision *events* (rising edge), not per-step samples: a single
    // glancing contact would otherwise dominate the metric.
    const bool colliding = vehicle.inCollision(world_.grid());
    if (colliding && !was_colliding) ++res.collisions;
    was_colliding = colliding;

    // -------- IMU at 200 Hz ----------------------------------------------
    std::vector<ImuSample> imu_batch;
    imu_->step(t, vs2.v, vs2.v * std::tan(vs2.steer) / config_.vehicle.wheelbase, vs2.accel, &imu_batch);
    faults_->applyImu(t, &imu_batch);
    res.imu_samples += imu_batch.size();

    // -------- odometry at 20 Hz ------------------------------------------
    if (t + 1e-9 >= next_odom) {
      const Pose2 d = pose2Compose(pose2Inverse(last_odom_pose), vs2.pose);
      OdomSample tick = odom_->measure(t, d.x, d.th);
      faults_->applyOdom(t, &tick);
      last_odom_pose = vs2.pose;
      next_odom += odom_dt;
      odom_acc_forward += tick.d_forward;
      odom_acc_yaw += tick.d_yaw;
      odom_acc_valid = odom_acc_valid && tick.valid;
    }
    OdomSample os;
    os.t = t;

    // -------- camera at the camera rate -----------------------------------
    if (t + 1e-9 >= next_cam) {
      frame.t = t;
      frame.ground_truth = vs2.pose;
      frame.imu = imu_batch;
      frame.odom.d_forward = odom_acc_forward;
      frame.odom.d_yaw = odom_acc_yaw;
      frame.odom.valid = odom_acc_valid;
      frame.odom.t = t;
      odom_acc_forward = 0;
      odom_acc_yaw = 0;
      odom_acc_valid = true;

      if (config_.stereo) {
        renderer_->renderStereo(vs2.pose, &frame.image, &frame.image_right);
      } else {
        renderer_->render(vs2.pose, &frame.image);
      }
      frame.camera_ok = faults_->applyCamera(t, &frame.image);

      frame.visible_landmarks.clear();
      frame.visible_pixels.clear();
      if (!faults_->suppressLandmarks(t)) {
        const auto& ids = renderer_->visibleLandmarks();
        const auto& px = renderer_->visiblePixels();
        for (Size i = 0; i < ids.size(); ++i) {
          if (config_.camera_noise.occlusion_prob > 0 && rng.uniform() < config_.camera_noise.occlusion_prob)
            continue;
          frame.visible_landmarks.push_back(ids[i]);
          frame.visible_pixels.push_back(px[i]);
        }
      }
      if (!frame.camera_ok) ++res.dropped_frames;

      GroundTruthFrame gt;
      gt.t = t;
      gt.pose = vs2.pose;
      gt.v = vs2.v;
      gt.landmark_ids = frame.visible_landmarks;
      gt.landmark_px = frame.visible_pixels;
      res.truth.add(gt);

      if (frame_sink) frame_sink(frame);
      ++res.frames;
      next_cam += cam_dt;
    }

    // Advance simulated time with a fixed step so the fault schedule fires at
    // identical times regardless of host speed.
    t += sim_dt;
  }

  res.path_length = res.truth.pathLength();
  res.faults_injected = faults_->injectedCount();
  const auto wall1 = std::chrono::steady_clock::now();
  res.wall_time_s = std::chrono::duration<Scalar>(wall1 - wall0).count();
  return res;
}

bool writePgm(const std::string& path, const GridMap& grid, Scalar scale) {
  const int step = std::max(1, static_cast<int>(std::round(1.0 / scale)));
  const int w = (grid.w + step - 1) / step;
  const int h = (grid.h + step - 1) / step;
  std::ofstream f(path, std::ios::binary);
  if (!f) return false;
  f << "P5\n" << w << ' ' << h << "\n255\n";
  for (int y = h - 1; y >= 0; --y) {
    for (int x = 0; x < w; ++x) {
      const Cell c = grid.at(x * step, y * step);
      const std::uint8_t v = (c == Cell::kOccupied) ? 30 : ((c == Cell::kFree) ? 245 : 128);
      f.write(reinterpret_cast<const char*>(&v), 1);
    }
  }
  return static_cast<bool>(f);
}

}  // namespace sim
}  // namespace luna
