// LUNA -- the end-to-end simulator harness.
//
// Produces a deterministic stream of sensor measurements from a scripted route
// through a generated world. The same seed reproduces the same stream bit for
// bit, which is what makes `luna bench` figures reproducible across
// machines (the portability claim in the contract).
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "luna/sim/faults.hpp"
#include "luna/sim/sensors.hpp"
#include "luna/sim/vehicle.hpp"
#include "luna/sim/world.hpp"

namespace luna {
namespace sim {

/// Everything the SLAM frontend consumes for one timestep.
struct SensorFrame {
  Scalar t = 0;
  Pose2 ground_truth;
  Image image;
  Image image_right;
  std::vector<ImuSample> imu;
  OdomSample odom;
  bool camera_ok = true;
  /// Ground-truth landmark ids visible this frame (evaluator only).
  std::vector<std::uint32_t> visible_landmarks;
  std::vector<Vec2> visible_pixels;
};

struct SimConfig {
  WorldSpec world;
  CameraModel camera;
  CameraNoise camera_noise;
  ImuNoise imu_noise;
  OdomNoise odom_noise;
  VehicleParams vehicle;
  FaultSchedule faults;

  Scalar duration = 60.0;        ///< simulated seconds
  Scalar imu_rate = 200.0;       ///< Hz
  Scalar camera_rate = 10.0;     ///< Hz
  Scalar odom_rate = 20.0;       ///< Hz
  Scalar control_rate = 20.0;    ///< Hz
  /// Target speed for the scripted reference route. Kept moderate on purpose:
  /// the reference run exists to be a clean, drivable baseline, so a collision
  /// during fault injection is attributable to the fault rather than to
  /// tracking error accumulated at speed.
  Scalar target_speed = 0.9;     ///< m/s for the scripted route
  bool stereo = true;
  bool record_images = false;    ///< keep images in the frame log (RAM heavy)
};

struct SimResult {
  WorldStats world_stats;
  Scalar path_length = 0;
  Size frames = 0;
  Size imu_samples = 0;
  Size dropped_frames = 0;
  Size faults_injected = 0;
  Size collisions = 0;
  Scalar wall_time_s = 0;
  std::vector<Vec2> route;       ///< the scripted reference route
  GroundTruthLog truth;          ///< full ground-truth trajectory
};

/// Builds a scripted route: a closed tour of the world's free space that
/// visits every room/hall, computed from a grid path between sampled extremes.
///
/// `clearance_m` is the preferred margin from every obstacle, and
/// `min_clearance_m` is the smallest margin the result may ever contain. Both
/// must exceed the vehicle's half-diagonal footprint, or a scripted run reports
/// collisions that are route defects rather than control failures. When the
/// preferred margin cannot connect the whole tour, the builder relaxes toward
/// `min_clearance_m` and keeps the *widest* route it can connect, never the
/// narrowest one that happens to visit the most anchors.
std::vector<Vec2> buildRoute(const World& world, Rng* rng, Scalar clearance_m = 0.75,
                             Scalar min_clearance_m = 0.45, Scalar min_turn_radius = 0.0);

/// Runs the simulator. When `frame_sink` is non-null it is called for every
/// camera frame, which is how the pipeline consumes the simulator without the
/// simulator knowing anything about SLAM.
class Simulator {
 public:
  explicit Simulator(const SimConfig& config);

  using FrameSink = std::function<void(const SensorFrame&)>;

  SimResult run(const FrameSink& frame_sink);

  const World& world() const { return world_; }
  const SimConfig& config() const { return config_; }

 private:
  SimConfig config_;
  World world_;
  std::unique_ptr<Renderer> renderer_;
  std::unique_ptr<ImuEmulator> imu_;
  std::unique_ptr<OdomEmulator> odom_;
  std::unique_ptr<FaultInjector> faults_;
};

/// Render the world to a PGM so a human can look at the benchmark scene.
bool writePgm(const std::string& path, const GridMap& grid, Scalar scale = 0.25);

}  // namespace sim
}  // namespace luna
