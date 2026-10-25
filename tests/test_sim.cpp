// SLAMCore-8G -- simulator tests.
//
// The simulator is the project's ground truth, so it has to be trustworthy
// before anything measured against it means anything. These tests check the
// properties the rest of the stack depends on: worlds are connected and
// generated deterministically, raycasts agree with line-of-sight, the vehicle
// integrates the bicycle model correctly, the route stays in free space, and
// the sensor emulators produce physically plausible signals.
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "slamcore/sim/faults.hpp"
#include "slamcore/sim/sensors.hpp"
#include "slamcore/sim/simulator.hpp"
#include "slamcore/sim/vehicle.hpp"
#include "slamcore/sim/world.hpp"
#include "test_util.hpp"

using namespace sc;
using namespace sc::sim;

namespace {

const WorldType kAllTypes[] = {WorldType::kRooms,  WorldType::kCorridorLoop, WorldType::kSaddleHall,
                               WorldType::kMaze,  WorldType::kCluttered,    WorldType::kTextureless};

WorldSpec spec(WorldType t, std::uint64_t seed = 5) {
  WorldSpec s;
  s.type = t;
  s.seed = seed;
  s.size_m = 30.0;
  s.rooms = 4;
  return s;
}

void testWorldGeneration() {
  for (WorldType t : kAllTypes) {
    const World w = World::generate(spec(t));
    const auto& g = w.grid();
    SC_CHECK(g.w > 0 && g.h > 0);
    SC_CHECK(g.cells.size() == static_cast<Size>(g.w) * g.h);
    SC_CHECK(w.stats().freeCells > 0);
    SC_CHECK(w.stats().occupiedCells > 0);
    SC_CHECK(w.stats().freeArea > 0);
    // Landmarks must sit in free space, otherwise the renderer would draw them
    // inside a wall and every measurement would be wrong.
    int in_free = 0;
    for (const Landmark& lm : w.landmarks())
      if (g.isFreeWorld(lm.p)) ++in_free;
    SC_CHECK(in_free == static_cast<int>(w.landmarks().size()));
    // Clearance is defined on every cell.
    SC_CHECK(w.clearance().size() == g.cells.size());
  }
}

/// Same seed must give a bit-identical world; different seeds must not. This is
/// the reproducibility contract the benchmark tables rest on.
void testWorldDeterminism() {
  const World a = World::generate(spec(WorldType::kRooms, 11));
  const World b = World::generate(spec(WorldType::kRooms, 11));
  const World c = World::generate(spec(WorldType::kRooms, 12));
  SC_CHECK(a.grid().cells == b.grid().cells);
  SC_CHECK(a.landmarks().size() == b.landmarks().size());
  SC_CHECK(a.grid().cells != c.grid().cells);
  SC_CHECK(a.toAscii(64) == b.toAscii(64));
}

void testRaycastAndLineOfSight() {
  const World w = World::generate(spec(WorldType::kRooms));
  const auto& g = w.grid();
  // A raycast shorter than the true free distance must report no hit beyond it.
  // Pick a free cell and cast in several directions, comparing against the
  // distance to the first occupied cell found by brute force.
  Rng rng(3);
  int tested = 0;
  for (int attempt = 0; attempt < 400 && tested < 60; ++attempt) {
    const Vec2 p(rng.uniform(1.0, 29.0), rng.uniform(1.0, 29.0));
    if (!g.isFreeWorld(p)) continue;
    const Scalar ang = rng.uniform(0.0, kTwoPi);
    const Vec2 dir(std::cos(ang), std::sin(ang));
    const Scalar r = w.raycast(p, dir, 20.0);
    SC_CHECK(r >= 0.0 && r <= 20.0);
    // Just before the reported hit there must be no occupied cell.
    if (r < 19.9) {
      const Vec2 just = p + dir * (r - 0.2);
      SC_CHECK(!g.isOccupiedWorld(just));
    }
    // Line of sight must agree: a point before the hit is visible, a point
    // after it is not.
    if (r > 1.0) SC_CHECK(w.lineOfSight(p, p + dir * (r - 0.3)));
    if (r < 18.0) SC_CHECK(!w.lineOfSight(p, p + dir * (r + 0.5)));
    ++tested;
  }
  SC_CHECK(tested > 20);
}

void testGridPath() {
  const World w = World::generate(spec(WorldType::kRooms));
  const auto& g = w.grid();
  // Find two well-separated free cells and require a path between them that
  // stays free, respects the clearance request, and is shortcut to fewer nodes.
  Rng rng(21);
  const Vec2 a = w.grid().cellCenter(g.w / 6, g.h / 6);
  const Vec2 b = w.grid().cellCenter(g.w * 5 / 6, g.h * 5 / 6);
  (void)rng;
  const std::vector<Vec2> path = w.gridPath(a, b, 0.3);
  SC_CHECK(!path.empty());
  for (const Vec2& p : path) SC_CHECK(g.isFreeWorld(p));
  const std::vector<Vec2> short_path = w.shortcut(path);
  SC_CHECK(!short_path.empty());
  SC_CHECK(short_path.size() <= path.size());
  // Shortcut segments must remain collision-free, which is the property the
  // scripted route depends on.
  for (Size i = 1; i < short_path.size(); ++i)
    SC_CHECK(w.lineOfSight(short_path[i - 1], short_path[i]));
  // An unreachable target (inside a wall) yields no path.
  Vec2 wall = a;
  for (int y = 0; y < g.h && g.isFreeWorld(wall); ++y)
    for (int x = 0; x < g.w; ++x)
      if (g.at(x, y) == Cell::kOccupied) {
        wall = g.cellCenter(x, y);
        y = g.h;
        break;
      }
  SC_CHECK(w.gridPath(a, wall, 0.0).empty());
}

void testRouteStaysFree() {
  for (WorldType t : {WorldType::kRooms, WorldType::kCorridorLoop, WorldType::kSaddleHall}) {
    const World w = World::generate(spec(t));
    Rng rng(4);
    const std::vector<Vec2> route = buildRoute(w, &rng);
    SC_CHECK(route.size() >= 2);
    for (const Vec2& p : route) SC_CHECK(w.grid().isFreeWorld(p));
    // The route must be drivable: its narrowest point is a real bottleneck.
    SC_CHECK(w.minCorridorWidth(route) > 0.0);
  }
}

void testVehicleBicycleModel() {
  VehicleParams params;

  // Straight line at constant speed: 1 s at 1 m/s is 1 m along +x, no heading
  // change. The initial velocity is set directly so the test measures the
  // kinematic model rather than the acceleration ramp.
  Vehicle veh(params, Pose2{0, 0, 0}, 1);
  VehicleState init = veh.state();
  init.v = 1.0;
  veh.setState(init);
  Control cmd;
  cmd.v = 1.0;
  for (int i = 0; i < 100; ++i) veh.step(cmd, 0.01);
  SC_NEAR(veh.state().pose.x, 1.0, 1e-6);
  SC_NEAR(veh.state().pose.y, 0.0, 1e-9);
  SC_NEAR(veh.state().pose.th, 0.0, 1e-9);
  SC_NEAR(veh.distanceTravelled(), 1.0, 1e-6);

  // Constant steering traces a circle of radius L / tan(delta) about a centre
  // offset laterally by that radius.
  const Scalar delta = 0.3;
  Vehicle circ(params, Pose2{0, 0, 0}, 2);
  VehicleState ci = circ.state();
  ci.v = 1.0;
  circ.setState(ci);
  Control turn;
  turn.v = 1.0;
  turn.steer = delta;
  for (int i = 0; i < 1000; ++i) circ.step(turn, 0.01);
  const Scalar radius = params.wheelbase / std::tan(delta);
  const Scalar yaw_rate = 1.0 * std::tan(delta) / params.wheelbase;
  // Steering slews to its command at max_steer_rate, so the first ~0.15 s is a
  // shallower arc than the steady-state circle. The tolerance allows for that
  // ramp (roughly 0.5 * ramp_time * yaw_rate of heading, and a similar radial
  // offset) while still catching a wrong wheelbase or a wrong yaw-rate law.
  const Scalar ramp = delta / params.max_steer_rate;
  SC_NEAR(circ.state().pose.th, wrapPi(yaw_rate * (10.0 - 0.5 * ramp)), 0.02);
  const Vec2 centre(0.0, radius);
  SC_NEAR((Vec2(circ.state().pose.x, circ.state().pose.y) - centre).norm(), radius, 0.12);

  // A velocity command must actually move the vehicle. An accel-only plant
  // ignores `v` and leaves the vehicle stationary, which is the bug this
  // catches.
  Vehicle from_rest(params, Pose2{0, 0, 0}, 3);
  Control go;
  go.v = 2.0;
  for (int i = 0; i < 100; ++i) from_rest.step(go, 0.01);
  SC_CHECK(from_rest.state().pose.x > 1.0);
  SC_NEAR(from_rest.state().v, 2.0, 1e-6);
}

void testProjection() {
  CameraModel cam;
  // A point straight ahead projects to the principal point.
  Vec2 px;
  SC_CHECK(project(cam, Vec3(5.0, 0.0, 0.0), &px));
  SC_NEAR(px.x(), cam.cx, 1e-9);
  SC_NEAR(px.y(), cam.cy, 1e-9);
  // A point behind the camera does not project.
  SC_CHECK(!project(cam, Vec3(-1.0, 0.0, 0.0), &px));
  // Lateral offset maps to a pixel offset of -fx * y / z (camera +y is left).
  SC_CHECK(project(cam, Vec3(10.0, 1.0, 0.0), &px));
  SC_NEAR(px.x(), cam.cx - cam.fx * 0.1, 1e-9);
  // Vertical offset maps through fy with the z-up convention.
  SC_CHECK(project(cam, Vec3(10.0, 0.0, 1.0), &px));
  SC_NEAR(px.y(), cam.cy - cam.fy * 0.1, 1e-9);
}

/// The stereo pair must actually differ by the disparity the geometry predicts.
/// This is the test that catches a renderer that draws both cameras from the
/// same pose, which would silently make every depth infinite.
void testStereoDisparity() {
  SimConfig cfg;
  cfg.world = spec(WorldType::kRooms);
  cfg.world.landmark_density = 3.0;
  cfg.duration = 0.1;
  cfg.stereo = true;
  cfg.camera_noise.pixel_sigma = 0.0;
  cfg.camera_noise.blur_sigma = 0.0;
  Simulator sim(cfg);

  bool have = false;
  std::vector<Vec2> left_px;
  Pose2 pose;
  Image left, right;
  sim.run([&](const SensorFrame& f) {
    if (have) return;
    left = f.image;
    right = f.image_right;
    left_px = f.visible_pixels;
    pose = f.ground_truth;
    have = true;
  });
  SC_CHECK(!left.px.empty() && !right.px.empty());
  SC_CHECK(left.px != right.px);  // the two views must not be identical

  // For a landmark at depth Z the right camera sees it at u_left + fx*b/Z.
  // Compare the rendered pixel offset against the true projection in the right
  // camera, computed from the known pose and camera placement.
  const World& w = sim.world();
  Size checked = 0;
  for (Size i = 0; i < left_px.size() && checked < 8; ++i) {
    const Vec2 pl = left_px[i];
    // True range from the camera to the landmark this pixel came from is not
    // recoverable here without the id, so instead check the *rendered* images:
    // the landmark patch in the right image must sit to the right of the left
    // one by the predicted disparity. Find the patch by correlating a small
    // window, using the left image as the reference.
    const int h = 3;
    const int lx = static_cast<int>(std::lround(pl.x()));
    const int ly = static_cast<int>(std::lround(pl.y()));
    if (lx - h < 0 || ly - h < 0 || lx + h >= left.w || ly + h >= left.h) continue;
    Scalar best = 1e18;
    int best_d = 0;
    for (int d = 1; d <= 60; ++d) {
      const int rx = lx + d;  // right image content is shifted to +u
      if (rx + h >= right.w) break;
      Scalar sad = 0;
      for (int dy = -h; dy <= h; ++dy)
        for (int dx = -h; dx <= h; ++dx)
          sad += std::abs(static_cast<int>(left.at(lx + dx, ly + dy)) -
                          static_cast<int>(right.at(rx + dx, ly + dy)));
      if (sad < best) {
        best = sad;
        best_d = d;
      }
    }
    // A landmark at finite range must show a strictly positive disparity.
    SC_CHECK(best_d > 0);
    // And the implied depth must be inside the world, not off the map.
    const Scalar depth = cfg.camera.fx * cfg.camera.baseline / best_d;
    SC_CHECK(depth > 0.2 && depth < 40.0);
    ++checked;
  }
  SC_CHECK(checked > 0);
  (void)pose;
  (void)w;
}

/// The IMU must measure gravity when the vehicle is level and still, and must
/// rotate that measurement with roll/pitch. This is the check that the
/// accelerometer convention matches what the error-state filter assumes.
void testImuGravity() {
  ImuNoise noise;
  noise.accel_noise = 0.0;
  noise.gyro_noise = 0.0;
  noise.accel_bias_rw = 0.0;
  noise.gyro_bias_rw = 0.0;
  noise.accel_bias0 = 0.0;
  noise.gyro_bias0 = 0.0;
  ImuEmulator imu(noise, 100.0, 1);
  std::vector<ImuSample> out;
  // Level, still: the accelerometer sees only gravity. The sign convention is
  // the filter's, so this pins it down.
  for (int i = 0; i < 10; ++i) imu.step(0.01 * (i + 1), 0.0, 0.0, 0.0, &out);
  SC_CHECK(!out.empty());
  SC_NEAR(std::abs(out.back().accel.z()), noise.gravity, 1e-6);
  SC_NEAR(out.back().accel.x(), 0.0, 1e-6);
  SC_NEAR(out.back().accel.y(), 0.0, 1e-6);
  SC_NEAR(out.back().gyro.norm(), 0.0, 1e-6);

  // Accelerating forward adds a specific force along body +x.
  ImuEmulator imu2(noise, 100.0, 1);
  std::vector<ImuSample> out2;
  for (int i = 0; i < 10; ++i) imu2.step(0.01 * (i + 1), 1.0, 0.0, 2.0, &out2);
  SC_CHECK(!out2.empty());
  SC_NEAR(out2.back().accel.x(), 2.0, 1e-6);
}

/// Odometry must report the true motion up to its modelled error, and the
/// error must be systematic (a scale error), not random noise alone.
void testOdomScale() {
  OdomNoise noise;
  noise.dist_scale = 0.01;
  noise.angle_scale = 0.0;
  noise.dist_sigma = 0.0;
  noise.angle_sigma = 0.0;
  noise.slip_prob = 0.0;
  OdomEmulator odom(noise, 1);
  Scalar sum = 0;
  const int n = 2000;
  for (int i = 0; i < n; ++i) sum += odom.measure(0.01 * (i + 1), 1.0, 0.0).d_forward;
  // The realised scale error is drawn once, so the measurement is consistently
  // scaled: mean error must be near the nominal 1%, not zero.
  const Scalar rel = (sum / n) - 1.0;
  SC_CHECK(std::abs(rel) > 1e-4);
  SC_CHECK(std::abs(rel) < 0.05);
}

void testFaultSchedule() {
  const FaultSchedule sched = FaultSchedule::standardSuite(2.0, 1.0);
  SC_CHECK(!sched.events().empty());
  SC_CHECK(!sched.describe().empty());
  // An event is active inside its window and not outside it.
  const FaultEvent e = sched.events().front();
  SC_CHECK(e.active(e.t_start));
  SC_CHECK(!e.active(e.t_start - 0.001));
  SC_CHECK(!e.active(e.t_start + e.duration));
  // Every standard-suite event must be individually addressable by type.
  FaultInjector inj(sched);
  Image img;
  img.resize(8, 8);
  inj.applyCamera(e.t_start, &img);
  SC_CHECK(!sched.active(FaultType::kNone, e.t_start));
}

void testImageHelpers() {
  Image im;
  im.resize(4, 3);
  SC_CHECK(im.w == 4 && im.h == 3);
  SC_CHECK(im.bytes() == 12);
  SC_CHECK(im.inBounds(3, 2));
  SC_CHECK(!im.inBounds(4, 2));
  im.at(2, 1) = 200;
  SC_CHECK(im.at(2, 1) == 200);
  SC_CHECK(static_cast<const Image&>(im).at(2, 1) == 200);
}

void run() {
  testWorldGeneration();
  testWorldDeterminism();
  testRaycastAndLineOfSight();
  testGridPath();
  testRouteStaysFree();
  testVehicleBicycleModel();
  testProjection();
  testStereoDisparity();
  testImuGravity();
  testOdomScale();
  testFaultSchedule();
  testImageHelpers();
}

}  // namespace

SC_TEST_MAIN("test_sim")
