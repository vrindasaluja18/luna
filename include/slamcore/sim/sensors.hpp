// SLAMCore-8G -- sensor emulators.
//
// These are *renderers*, not shortcut models. The camera emulator rasterises
// the occupancy grid into a real grayscale image (walls become shaded vertical
// surfaces by per-column ray casting, landmarks become perspective-projected
// textured patches), and the frontend then runs FAST/BRIEF/KLT on those actual
// pixel buffers. That is what makes the "validated without hardware" claim
// meaningful: the feature pipeline sees images, not pre-digested points.
//
// The IMU and wheel encoders are modelled with the standard error taxonomy
// (white noise + random-walk bias + scale error), because those are precisely
// the terms the error-state filter must estimate.
#pragma once

#include <cstdint>
#include <vector>

#include "slamcore/core/se3.hpp"
#include "slamcore/core/stats.hpp"
#include "slamcore/sim/world.hpp"

namespace sc {
namespace sim {

/// Single-channel 8-bit image with preallocated storage.
struct Image {
  int w = 0, h = 0;
  std::vector<std::uint8_t> px;
  void resize(int width, int height) {
    w = width;
    h = height;
    px.assign(static_cast<Size>(w) * h, 0);
  }
  std::uint8_t& at(int x, int y) { return px[static_cast<Size>(y) * w + x]; }
  std::uint8_t at(int x, int y) const { return px[static_cast<Size>(y) * w + x]; }
  bool inBounds(int x, int y) const { return x >= 0 && y >= 0 && x < w && y < h; }
  Size bytes() const { return px.size(); }
};

/// Pinhole camera intrinsics plus body extrinsics.
struct CameraModel {
  int width = 320;
  int height = 240;
  Scalar fx = 260.0;
  Scalar fy = 260.0;
  Scalar cx = 160.0;
  Scalar cy = 120.0;
  Scalar body_height = 0.6;   ///< optical centre height above floor
  Scalar baseline = 0.12;     ///< stereo baseline (metres, lateral)
  Scalar max_range = 20.0;
  Scalar fov_h() const { return 2.0 * std::atan(0.5 * width / fx); }
};

/// Project a body-frame point into pixel coordinates. Returns false if behind
/// the camera. Camera looks along body +x, with +y to the left and z up.
bool project(const CameraModel& cam, const Vec3& p_body, Vec2* pixel);

/// Band-limited procedural surface texture.
///
/// The renderer needs *2D* texture on walls, not just vertical shading: a
/// feature pipeline cannot track structure that only varies along one axis.
/// This samples a precomputed value-noise field indexed by a wall-local
/// coordinate, so the pattern is attached to the surface and therefore moves
/// with correct parallax when the camera moves.
class Texture {
 public:
  Texture(int size = 64, std::uint64_t seed = 0x7E47u);
  /// Sample at (s, t) in metres. Aperiodic, so the pattern never repeats.
  Scalar sample(Scalar s, Scalar t) const;
  Scalar amplitude() const { return amp_; }

 private:
  int size_;
  Scalar amp_;
  std::uint32_t seed_ = 0;
  std::vector<Scalar> field_;  ///< reserved for a future cached variant
};

/// Per-pixel sensor noise parameters for the renderer.
struct CameraNoise {
  Scalar pixel_sigma = 2.0;      ///< additive Gaussian read noise (grey levels)
  Scalar blur_sigma = 0.8;       ///< Gaussian blur (defocus + motion)
  Scalar exposure = 1.0;         ///< multiplicative gain drift
  Scalar dropout_prob = 0.0;     ///< per-frame probability of a dead frame
  Scalar occlusion_prob = 0.0;   ///< per-landmark probability of being hidden
};

/// Renders the world into images. Buffers are allocated once in the
/// constructor so per-frame rendering performs no heap allocation.
class Renderer {
 public:
  Renderer(const World& world, const CameraModel& cam, const CameraNoise& noise, std::uint64_t seed);

  /// Render the monocular view for a vehicle pose.
  void render(const Pose2& pose, Image* out);

  /// Render the stereo pair. The left camera sits at the vehicle pose and the
  /// right camera is displaced along body -y (to the right) by `baseline`. A
  /// landmark at depth Z therefore projects to a smaller u in the right image,
  /// so with disparity defined as u_left - u_right it equals fx * baseline / Z
  /// and is positive for every landmark in front of the pair.
  void renderStereo(const Pose2& pose, Image* left, Image* right);

  /// Landmark ids whose projection landed inside the image and passed the
  /// occlusion test this frame. Used as SLAM ground truth for the evaluator.
  const std::vector<std::uint32_t>& visibleLandmarks() const { return visible_; }

  /// Camera-frame landmark positions for the visible set, parallel to
  /// visibleLandmarks().
  const std::vector<Vec3>& visiblePoints() const { return visible_pts_; }

  /// Exact per-landmark pixel observations, parallel to visibleLandmarks().
  const std::vector<Vec2>& visiblePixels() const { return visible_px_; }

 private:
  void rasterise(const Pose2& pose, const Vec2& body_offset, Image* out);
  void projectLandmarks(const Pose2& pose, const Vec2& body_offset, bool stamp, Image* out,
                        std::vector<std::uint32_t>* ids, std::vector<Vec3>* pts, std::vector<Vec2>* px);
  void stampLandmark(const CameraModel& cam, const Vec3& p_cam, std::uint32_t seed, Image* out);
  void applyNoise(Image* img);
  void blur(Image* img);

  const World& world_;
  Texture texture_;
  CameraModel cam_;
  CameraNoise noise_;
  Rng rng_;
  Image tmp_;
  std::vector<Scalar> depth_col_;
  std::vector<std::uint32_t> visible_;
  std::vector<Vec3> visible_pts_;
  std::vector<Vec2> visible_px_;
};

/// IMU sample: accelerometer and gyroscope, body frame.
struct ImuSample {
  Scalar t = 0;
  Vec3 accel = Vec3::Zero();
  Vec3 gyro = Vec3::Zero();
};

/// IMU error model. Random-walk biases are the dominant long-horizon error and
/// are exactly what the error-state filter estimates online.
struct ImuNoise {
  Scalar accel_noise = 0.02;        ///< m/s^2/sqrt(Hz)
  Scalar gyro_noise = 0.002;        ///< rad/s/sqrt(Hz)
  Scalar accel_bias_rw = 0.001;     ///< m/s^3/sqrt(Hz)
  Scalar gyro_bias_rw = 0.0001;     ///< rad/s^2/sqrt(Hz)
  Scalar accel_bias0 = 0.05;
  Scalar gyro_bias0 = 0.005;
  Scalar gravity = 9.81;
};

/// Stateful IMU emulator: integrates the true trajectory between samples and
/// injects the modelled errors.
class ImuEmulator {
 public:
  ImuEmulator(const ImuNoise& noise, Scalar rate_hz, std::uint64_t seed);

  /// Integrate from the previous sample to time `t` given the true motion, and
  /// append samples to `out`. `a_true` and `omega_true` are body-frame.
  void step(Scalar t, Scalar v, Scalar omega, Scalar a_body, std::vector<ImuSample>* out);

  const Vec3& accelBias() const { return ba_; }
  const Vec3& gyroBias() const { return bg_; }

 private:
  ImuNoise noise_;
  Scalar dt_ = 0.005;
  Scalar last_t_ = -1.0;
  Vec3 ba_ = Vec3::Zero();
  Vec3 bg_ = Vec3::Zero();
  Rng rng_;
};

/// Wheel-odometry error model: a systematic scale error (the usual culprit for
/// slow drift) plus white noise.
struct OdomNoise {
  Scalar dist_scale = 0.01;      ///< 1% under/over-read of distance
  Scalar angle_scale = 0.01;     ///< 1% yaw scale error
  Scalar dist_sigma = 0.01;      ///< metres per metre travelled
  Scalar angle_sigma = 0.005;    ///< rad per rad turned
  Scalar slip_prob = 0.0;        ///< probability of a wheel-slip outlier
  Scalar slip_magnitude = 5.0;   ///< outlier multiplier
};

/// Relative motion measurement between two consecutive frames.
struct OdomSample {
  Scalar t = 0;
  Scalar d_forward = 0;
  Scalar d_yaw = 0;
  bool valid = true;
};

class OdomEmulator {
 public:
  OdomEmulator(const OdomNoise& noise, std::uint64_t seed);
  OdomSample measure(Scalar t, Scalar d_forward_true, Scalar d_yaw_true);

 private:
  OdomNoise noise_;
  Rng rng_;
  Scalar dist_scale_realised_;
  Scalar angle_scale_realised_;
};

}  // namespace sim
}  // namespace sc
