#include "luna/sim/sensors.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace luna {
namespace sim {

namespace {

/// Hash an integer lattice cell to a value in [-1, 1). Deterministic, no state,
/// and -- crucially -- *aperiodic*: the same pattern never repeats anywhere on
/// the surface.
inline Scalar hashCell(int x, int y, std::uint32_t seed) {
  std::uint32_t h = static_cast<std::uint32_t>(x) * 374761393u +
                    static_cast<std::uint32_t>(y) * 668265263u + seed * 2246822519u;
  h = (h ^ (h >> 13)) * 1274126177u;
  h ^= h >> 16;
  return static_cast<Scalar>(h & 0xFFFFu) / 32768.0 - 1.0;
}

inline Scalar smoothstep(Scalar t) { return t * t * (3.0 - 2.0 * t); }

}  // namespace

Texture::Texture(int size, std::uint64_t seed) : size_(size), amp_(1.0) {
  (void)size_;
  seed_ = static_cast<std::uint32_t>(seed ^ (seed >> 32));
}

Scalar Texture::sample(Scalar s, Scalar t) const {
  // Aperiodic value noise: three octaves of hashed lattice cells, with the
  // period expressed in *metres* so the projected pattern scales correctly
  // with range. Periods of 0.8/0.4/0.2 m project to roughly 20/10/5 px at
  // 10 m and stay above 2 px across the camera's working range, which keeps
  // the finest detail clear of Nyquist.
  //
  // The pattern is deliberately non-repeating. A tiled texture gives KLT a
  // periodic cost surface, and then a coarse pyramid level can converge to a
  // neighbouring period instead of the true minimum -- a silent multi-pixel
  // error. Aperiodic noise has a unique minimum and no such ambiguity.
  // Periods are in metres of surface. The 0.55/0.28/0.14 m band projects to
  // roughly 13/7/3.5 px at 10 m and stays clear of Nyquist out to the camera's
  // 20 m range, so the pyramid can decimate without aliasing.
  const Scalar periods[3] = {0.55, 0.28, 0.14};
  const Scalar weights[3] = {0.45, 0.33, 0.22};
  Scalar acc = 0;
  for (int o = 0; o < 3; ++o) {
    const Scalar inv = 1.0 / periods[o];
    const Scalar fs = s * inv, ft = t * inv;
    const int i = static_cast<int>(std::floor(fs));
    const int j = static_cast<int>(std::floor(ft));
    const Scalar us = smoothstep(fs - i);
    const Scalar ut = smoothstep(ft - j);
    const std::uint32_t sd = seed_ + static_cast<std::uint32_t>(o) * 0x9E3779B9u;
    const Scalar v00 = hashCell(i, j, sd), v10 = hashCell(i + 1, j, sd);
    const Scalar v01 = hashCell(i, j + 1, sd), v11 = hashCell(i + 1, j + 1, sd);
    const Scalar top = v00 * (1 - us) + v10 * us;
    const Scalar bot = v01 * (1 - us) + v11 * us;
    acc += weights[o] * (top * (1 - ut) + bot * ut);
  }
  return acc;
}

bool project(const CameraModel& cam, const Vec3& p_body, Vec2* pixel) {
  // Body frame: x forward, y left, z up. Camera optical frame: z forward,
  // x right, y down. Body +y (left) is therefore optical -x, and body +z (up)
  // is optical -y. Both signs matter: getting one wrong mirrors the image, and
  // a mirrored projection still lands landmarks inside the frame, so the
  // failure is silent until the depth comes out inconsistent.
  const Vec3 p_opt(-p_body.y(), -p_body.z(), p_body.x());
  if (p_opt.z() < 1e-6) return false;
  const Scalar u = cam.cx + cam.fx * p_opt.x() / p_opt.z();
  const Scalar v = cam.cy + cam.fy * p_opt.y() / p_opt.z();
  if (pixel) *pixel = Vec2(u, v);
  return true;
}

// ---------------------------------------------------------------------------
// Renderer
// ---------------------------------------------------------------------------

Renderer::Renderer(const World& world, const CameraModel& cam, const CameraNoise& noise, std::uint64_t seed)
    : world_(world), cam_(cam), noise_(noise), rng_(seed) {
  tmp_.resize(cam.width, cam.height);
  depth_col_.assign(cam.width, cam.max_range);
}

void Renderer::rasterise(const Pose2& pose, const Vec2& body_offset, Image* out) {
  const int W = cam_.width, H = cam_.height;
  out->resize(W, H);

  // Camera eye expressed in the world, accounting for the body-frame offset.
  const Vec2 eye = pose2Apply(pose, body_offset);

  for (int u = 0; u < W; ++u) {
    // Optical x maps to body -y (see project()), so a pixel left of centre
    // looks toward body +y. Optical z (forward) maps to body +x.
    const Scalar x_opt = (u - cam_.cx) / cam_.fx;
    const Vec2 dir_body(1.0, -x_opt);
    const Scalar dirn = dir_body.norm();
    const Vec2 dirw = Vec2(std::cos(pose.th) * dir_body.x() - std::sin(pose.th) * dir_body.y(),
                           std::sin(pose.th) * dir_body.x() + std::cos(pose.th) * dir_body.y()) /
                      dirn;
    const Scalar range = world_.raycast(eye, dirw, cam_.max_range);
    depth_col_[u] = range;

    // Perspective column: a wall at range r spanning floor to ceiling maps to
    // a vertical strip. Wall height 2.0 m. Body z is up and optical y is down,
    // so the floor (z = -body_height) lands below the principal point.
    const Scalar wall_h = 2.0;
    const Scalar inv_r = 1.0 / std::max(range, 1e-3);
    const Scalar y_floor = cam_.cy + cam_.fy * cam_.body_height * inv_r;
    const Scalar y_top = cam_.cy + cam_.fy * (cam_.body_height - wall_h) * inv_r;
    int v0 = std::max(static_cast<int>(std::ceil(y_top)), 0);
    int v1 = std::min(static_cast<int>(std::floor(y_floor)), H - 1);

    // Shading: distance falloff, a vertical gradient, and a band-limited
    // surface texture sampled in wall-local coordinates. The texture is what
    // gives KLT a trackable 2D pattern; the smooth falloffs alone would be
    // nearly edge-like and ambiguous along the wall direction.
    const Scalar base = 160.0 * (1.0 - std::min(range / cam_.max_range, 1.0)) + 45.0;
    // Wall-local horizontal coordinate: arc length along the wall, which for
    // our purposes is the projection of the ray direction onto the wall plane.
    const Scalar wall_s = eye.x() + range * dirw.x();
    const Scalar wall_s2 = eye.y() + range * dirw.y();
    const Scalar tex_u = (std::abs(dirw.x()) > std::abs(dirw.y())) ? wall_s2 : wall_s;
    for (int v = v0; v <= v1; ++v) {
      const Scalar grad = 1.0 - 0.30 * static_cast<Scalar>(v - v0) / std::max(1, v1 - v0 + 1);
      // Vertical surface coordinate: v maps to height on the wall plane.
      const Scalar tex_v = static_cast<Scalar>(v) * 0.01;
      const Scalar tex = texture_.sample(tex_u, tex_v);
      out->at(u, v) = static_cast<std::uint8_t>(clamp(base * grad * (1.0 + 0.45 * tex), 0.0, 255.0));
    }
    for (int v = v1 + 1; v < H; ++v) {
      const Scalar frac = static_cast<Scalar>(v - v1) / std::max(1, H - v1);
      out->at(u, v) = static_cast<std::uint8_t>(clamp(70.0 - 40.0 * frac, 0.0, 255.0));
    }
    for (int v = 0; v < v0; ++v) {
      out->at(u, v) = static_cast<std::uint8_t>(clamp(20.0 + 10.0 * (v % 3), 0.0, 255.0));
    }
  }
}

void Renderer::projectLandmarks(const Pose2& pose, const Vec2& body_offset, bool stamp, Image* out,
                                std::vector<std::uint32_t>* ids, std::vector<Vec3>* pts,
                                std::vector<Vec2>* px) {
  // Transform a landmark into the camera's body frame. The camera's world
  // position is pose (+) body_offset; the camera shares the body orientation.
  const Pose2 cam_pose = pose2Compose(pose, Pose2{body_offset.x(), body_offset.y(), 0.0});
  const Vec2 cam_xy(cam_pose.x, cam_pose.y);
  const Scalar c = std::cos(pose.th), s = std::sin(pose.th);

  // Pass 1: cull, project and collect. Stamping has to happen after this,
  // because the draw order is what resolves overlaps between landmark patches.
  struct Visible {
    const Landmark* lm;
    Vec3 p_body;
    Vec2 p;
  };
  std::vector<Visible> visible;
  for (const Landmark& lm : world_.landmarks()) {
    const Vec2 rel = lm.p - cam_xy;
    const Vec3 p_body(c * rel.x() + s * rel.y(), -s * rel.x() + c * rel.y(),
                      lm.height - cam_.body_height);
    Vec2 p;
    if (!project(cam_, p_body, &p)) continue;
    if (p.x() < 0 || p.y() < 0 || p.x() >= out->w || p.y() >= out->h) continue;
    const int u = static_cast<int>(p.x());
    if (p_body.x() > depth_col_[u] + 0.15) continue;  // occluded by a wall
    if (ids) ids->push_back(lm.id);
    if (pts) pts->push_back(p_body);
    if (px) px->push_back(p);
    if (stamp && lm.textured) visible.push_back({&lm, p_body, p});
  }
  if (!stamp) return;

  // Pass 2: stamp far-to-near so a nearer patch paints over a farther one.
  // Stamping in landmark-id order instead lets an arbitrary far landmark erase
  // the near one it overlaps, and the block matcher then measures the far
  // landmark's disparity at the near one's pixel: a confident, wrong depth.
  std::sort(visible.begin(), visible.end(), [](const Visible& a, const Visible& b) {
    return a.p_body.x() > b.p_body.x();
  });
  for (const Visible& v : visible) stampLandmark(cam_, v.p_body, v.lm->descSeed, out);
}

void Renderer::stampLandmark(const CameraModel& cam, const Vec3& p_cam, std::uint32_t seed, Image* out) {
  Vec2 px;
  if (!project(cam, p_cam, &px)) return;
  // Landmark patch radius shrinks with distance.
  const Scalar depth = std::max(p_cam.x(), 1e-3);
  // The patch must be larger than the stereo matching window at every range,
  // otherwise the window straddles the patch boundary and the block matcher
  // locks onto the *wall* behind the landmark instead. That failure is silent:
  // it reports a confident depth, just the wrong surface's. Clamping the lower
  // bound above the window half-width keeps the window entirely inside the
  // patch, so what is matched is always the landmark.
  const Scalar radius = clamp(cam.fx * 0.35 / depth, 9.0, 16.0);
  const int rad = static_cast<int>(std::ceil(radius));
  // Stamp around the *sub-pixel* centre, not the rounded one. Rounding here
  // would shift the pattern by up to half a pixel, and because the two cameras
  // round independently that offset does not cancel: it becomes a systematic
  // disparity bias, which the block matcher then reads as a depth error.
  const int cu = static_cast<int>(std::floor(px.x()));
  const int cv = static_cast<int>(std::floor(px.y()));
  const Scalar sub_u = px.x() - cu;
  const Scalar sub_v = px.y() - cv;
  // Deterministic per-landmark texture. This is deliberately *aperiodic* noise
  // rather than a checkerboard: a regular checker has a periodic cost surface,
  // so a coarse pyramid level can converge to a neighbouring checker cell and
  // report a displacement off by exactly one period. Aperiodic noise has a
  // single global minimum, which is what lets the pyramid help rather than
  // introduce a silent multi-pixel error.
  Rng tex(seed);
  const Scalar phase_s = tex.uniform(0.0, 100.0);
  const Scalar phase_t = tex.uniform(0.0, 100.0);
  const Scalar contrast = 55.0 + 35.0 * tex.uniform();
  // Coordinate scale for the patch pattern. Chosen so the texture's finest
  // octave (0.2 m) spans roughly 8 px on the image: detail finer than that
  // would sit near Nyquist, and decimating for the pyramid would alias it.
  // Finest octave spans ~3 px, so a 9x9 matching window covers about three
  // periods. That is the regime where the cost curve has a sharp, unique
  // minimum; a smoother patch (period >= 5 px) gives a broad basin and the
  // block matcher can settle a couple of pixels off, which reads as tens of
  // percent of depth error at mid range.
  constexpr Scalar kPatchScale = 0.045;
  for (int dy = -rad; dy <= rad; ++dy) {
    for (int dx = -rad; dx <= rad; ++dx) {
      if (dx * dx + dy * dy > rad * rad) continue;
      const int x = cu + dx, y = cv + dy;
      if (!out->inBounds(x, y)) continue;
      const Scalar du = static_cast<Scalar>(dx) - sub_u;
      const Scalar dv = static_cast<Scalar>(dy) - sub_v;
      const Scalar n = texture_.sample(phase_s + du * kPatchScale, phase_t + dv * kPatchScale);
      out->at(x, y) = static_cast<std::uint8_t>(clamp(150.0 + contrast * n, 0.0, 255.0));
    }
  }
}

void Renderer::blur(Image* img) {
  if (noise_.blur_sigma <= 1e-3) return;
  // Separable 5-tap Gaussian with a sigma-derived kernel.
  const Scalar s = noise_.blur_sigma;
  Scalar k[5];
  Scalar sum = 0;
  for (int i = -2; i <= 2; ++i) {
    k[i + 2] = std::exp(-0.5 * (i * i) / (s * s));
    sum += k[i + 2];
  }
  for (int i = 0; i < 5; ++i) k[i] /= sum;

  const int W = img->w, H = img->h;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      Scalar acc = 0;
      for (int i = -2; i <= 2; ++i) acc += k[i + 2] * img->at(clamp(x + i, 0, W - 1), y);
      tmp_.at(x, y) = static_cast<std::uint8_t>(clamp(acc, 0.0, 255.0));
    }
  }
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      Scalar acc = 0;
      for (int i = -2; i <= 2; ++i) acc += k[i + 2] * tmp_.at(x, clamp(y + i, 0, H - 1));
      img->at(x, y) = static_cast<std::uint8_t>(clamp(acc, 0.0, 255.0));
    }
  }
}

void Renderer::applyNoise(Image* img) {
  const Scalar sigma = noise_.pixel_sigma;
  const Scalar gain = noise_.exposure;
  for (auto& p : img->px) {
    Scalar v = p * gain;
    if (sigma > 0) v += rng_.gaussian(sigma);
    p = static_cast<std::uint8_t>(clamp(v, 0.0, 255.0));
  }
}

void Renderer::render(const Pose2& pose, Image* out) {
  visible_.clear();
  visible_pts_.clear();
  visible_px_.clear();

  rasterise(pose, Vec2::Zero(), out);
  projectLandmarks(pose, Vec2::Zero(), true, out, &visible_, &visible_pts_, &visible_px_);

  blur(out);
  applyNoise(out);
}

void Renderer::renderStereo(const Pose2& pose, Image* left, Image* right) {
  // Left view also defines the reported visible set, in the left camera frame.
  render(pose, left);

  // The right camera is displaced along body -y (to the right, since +y is
  // left). A landmark at depth Z therefore lands at a *smaller* u in the right
  // image by exactly fx * baseline / Z. Defining disparity as u_left - u_right
  // makes it positive and equal to that quantity, which is the convention the
  // block matcher searches over (rx = lx - d).
  rasterise(pose, Vec2(0.0, -cam_.baseline), right);
  projectLandmarks(pose, Vec2(0.0, -cam_.baseline), true, right, nullptr, nullptr, nullptr);
  blur(right);
  applyNoise(right);
}

// ---------------------------------------------------------------------------
// IMU
// ---------------------------------------------------------------------------

ImuEmulator::ImuEmulator(const ImuNoise& noise, Scalar rate_hz, std::uint64_t seed)
    : noise_(noise), dt_(1.0 / rate_hz), rng_(seed) {
  ba_ = Vec3(rng_.gaussian(noise_.accel_bias0), rng_.gaussian(noise_.accel_bias0),
             rng_.gaussian(noise_.accel_bias0));
  bg_ = Vec3(rng_.gaussian(noise_.gyro_bias0), rng_.gaussian(noise_.gyro_bias0),
             rng_.gaussian(noise_.gyro_bias0));
}

void ImuEmulator::step(Scalar t, Scalar v, Scalar omega, Scalar a_body, std::vector<ImuSample>* out) {
  if (last_t_ < 0) last_t_ = t - dt_;
  const Scalar sqrt_dt = std::sqrt(dt_);
  for (Scalar ts = last_t_ + dt_; ts <= t + 1e-12; ts += dt_) {
    // Random-walk bias update.
    ba_ += Vec3(rng_.gaussian(noise_.accel_bias_rw * sqrt_dt), rng_.gaussian(noise_.accel_bias_rw * sqrt_dt),
                rng_.gaussian(noise_.accel_bias_rw * sqrt_dt));
    bg_ += Vec3(rng_.gaussian(noise_.gyro_bias_rw * sqrt_dt), rng_.gaussian(noise_.gyro_bias_rw * sqrt_dt),
                rng_.gaussian(noise_.gyro_bias_rw * sqrt_dt));

    ImuSample s;
    s.t = ts;
    // Specific force in body frame: forward acceleration plus gravity reaction
    // (level vehicle, so gravity acts on +z), then bias and white noise.
    s.accel = Vec3(a_body, 0.0, noise_.gravity) + ba_ +
              Vec3(rng_.gaussian(noise_.accel_noise / sqrt_dt), rng_.gaussian(noise_.accel_noise / sqrt_dt),
                   rng_.gaussian(noise_.accel_noise / sqrt_dt));
    s.gyro = Vec3(0.0, 0.0, omega) + bg_ +
             Vec3(rng_.gaussian(noise_.gyro_noise / sqrt_dt), rng_.gaussian(noise_.gyro_noise / sqrt_dt),
                  rng_.gaussian(noise_.gyro_noise / sqrt_dt));
    out->push_back(s);
  }
  last_t_ = t;
  (void)v;
}

// ---------------------------------------------------------------------------
// Wheel odometry
// ---------------------------------------------------------------------------

OdomEmulator::OdomEmulator(const OdomNoise& noise, std::uint64_t seed) : noise_(noise), rng_(seed) {
  dist_scale_realised_ = 1.0 + rng_.gaussian(noise_.dist_scale);
  angle_scale_realised_ = 1.0 + rng_.gaussian(noise_.angle_scale);
}

OdomSample OdomEmulator::measure(Scalar t, Scalar d_forward_true, Scalar d_yaw_true) {
  OdomSample s;
  s.t = t;
  const Scalar dist = std::abs(d_forward_true);
  const Scalar turn = std::abs(d_yaw_true);
  Scalar slip = 1.0;
  if (noise_.slip_prob > 0 && rng_.uniform() < noise_.slip_prob) {
    slip = noise_.slip_magnitude;
    s.valid = false;
  }
  s.d_forward = d_forward_true * dist_scale_realised_ * slip +
                rng_.gaussian(noise_.dist_sigma * dist) * slip;
  s.d_yaw = d_yaw_true * angle_scale_realised_ * slip + rng_.gaussian(noise_.angle_sigma * turn) * slip;
  return s;
}

}  // namespace sim
}  // namespace luna
