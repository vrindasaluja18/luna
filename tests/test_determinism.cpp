// SLAMCore-8G -- determinism tests.
//
// The contract promises that every benchmark number is reproducible from
// public data on any machine. That claim is only true if the simulator, the
// RNG and the whole sensor pipeline are bit-deterministic for a fixed seed.
// These tests assert exactly that, by running the same configuration twice and
// comparing the produced streams byte for byte, and by running a different seed
// and requiring the streams to differ.
//
// This is also the test that would catch a stray use of std::random_device, a
// time-seeded RNG, or iteration over an unordered container whose order
// depends on pointer values.
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "slamcore/frontend/features.hpp"
#include "slamcore/sim/faults.hpp"
#include "slamcore/sim/simulator.hpp"
#include "test_util.hpp"

using namespace sc;
using namespace sc::sim;
using namespace sc::frontend;

namespace {

/// Everything a run produces that a benchmark figure could depend on.
struct Fingerprint {
  std::uint64_t image_hash = 0;
  std::uint64_t right_hash = 0;
  std::uint64_t truth_hash = 0;
  std::uint64_t imu_hash = 0;
  Size frames = 0;
  Scalar path = 0;
  Scalar final_x = 0, final_y = 0, final_th = 0;
  Size features = 0, tracked = 0, stereo_valid = 0;
};

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  h ^= v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
  return h;
}

Fingerprint run(const SimConfig& cfg, bool with_frontend) {
  Simulator sim(cfg);
  Fingerprint fp;
  BriefPattern pat;
  FastOptions fast;
  fast.threshold = 12;
  BriefOptions brief;
  StereoOptions stereo_opt;
  KltTracker klt{KltOptions()};

  Image prev;
  std::vector<Feature> prev_feats;
  bool have_prev = false;

  const SimResult r = sim.run([&](const SensorFrame& f) {
    ++fp.frames;
    for (std::uint8_t v : f.image.px) fp.image_hash = mix(fp.image_hash, v);
    for (std::uint8_t v : f.image_right.px) fp.right_hash = mix(fp.right_hash, v);
    fp.truth_hash = mix(fp.truth_hash, static_cast<std::uint64_t>(f.ground_truth.x * 1e9));
    fp.truth_hash = mix(fp.truth_hash, static_cast<std::uint64_t>(f.ground_truth.y * 1e9));
    for (const ImuSample& s : f.imu)
      fp.imu_hash = mix(fp.imu_hash, static_cast<std::uint64_t>((s.accel.x() + 10.0) * 1e9));

    if (!with_frontend) {
      prev = f.image;
      return;
    }
    std::vector<Feature> feats;
    detectAndDescribe(f.image, pat, fast, brief, &feats);
    if (!f.image_right.px.empty())
      stereoDepth(f.image, f.image_right, cfg.camera, stereo_opt, &feats);
    fp.features += feats.size();
    for (const Feature& x : feats)
      if (x.depth > 0) ++fp.stereo_valid;
    if (have_prev) {
      std::vector<Feature> out;
      klt.track(prev, f.image, prev_feats, &out);
      fp.tracked += out.size();
    }
    prev = f.image;
    prev_feats = feats;
    have_prev = true;
  });

  fp.path = r.path_length;
  if (!r.truth.frames().empty()) {
    const auto& last = r.truth.frames().back();
    fp.final_x = last.pose.x;
    fp.final_y = last.pose.y;
    fp.final_th = last.pose.th;
  }
  return fp;
}

SimConfig cfg(std::uint64_t seed, WorldType t = WorldType::kRooms) {
  SimConfig c;
  c.world.type = t;
  c.world.seed = seed;
  c.world.size_m = 24.0;
  c.world.rooms = 4;
  c.duration = 6.0;
  c.stereo = true;
  return c;
}

void testSameSeedReproducesExactly() {
  const SimConfig c = cfg(12345);
  const Fingerprint a = run(c, true);
  const Fingerprint b = run(c, true);
  SC_CHECK(a.frames == b.frames);
  SC_CHECK(a.image_hash == b.image_hash);
  SC_CHECK(a.right_hash == b.right_hash);
  SC_CHECK(a.truth_hash == b.truth_hash);
  SC_CHECK(a.imu_hash == b.imu_hash);
  SC_NEAR(a.path, b.path, 0.0);
  SC_NEAR(a.final_x, b.final_x, 0.0);
  SC_NEAR(a.final_y, b.final_y, 0.0);
  SC_NEAR(a.final_th, b.final_th, 0.0);
  // The frontend is deterministic too, which is what makes a reported ATE
  // reproducible rather than merely similar.
  SC_CHECK(a.features == b.features);
  SC_CHECK(a.tracked == b.tracked);
  SC_CHECK(a.stereo_valid == b.stereo_valid);
}

void testDifferentSeedDiffers() {
  const Fingerprint a = run(cfg(1), false);
  const Fingerprint b = run(cfg(2), false);
  SC_CHECK(a.image_hash != b.image_hash);
  SC_CHECK(a.truth_hash != b.truth_hash);
}

/// Same seed and world, different world type: the type must actually change
/// the generated world rather than being ignored.
void testWorldTypeMatters() {
  const Fingerprint a = run(cfg(7, WorldType::kRooms), false);
  const Fingerprint b = run(cfg(7, WorldType::kSaddleHall), false);
  SC_CHECK(a.image_hash != b.image_hash);
  SC_CHECK(a.truth_hash != b.truth_hash);
}

/// Faults must be scheduled deterministically as well, or a fault-injection
/// safety result would not be reproducible.
void testFaultsAreDeterministic() {
  SimConfig c = cfg(99, WorldType::kRooms);
  c.faults = FaultSchedule::standardSuite(2.0, 1.0);
  const Fingerprint a = run(c, false);
  const Fingerprint b = run(c, false);
  SC_CHECK(a.image_hash == b.image_hash);
  SC_CHECK(a.truth_hash == b.truth_hash);
  // And a faulted run must differ from a clean one, i.e. the faults did
  // something rather than being silently inert.
  const Fingerprint clean = run(cfg(99, WorldType::kRooms), false);
  SC_CHECK(a.image_hash != clean.image_hash);
}

/// The frontend over a real simulated run must produce a usable amount of data.
/// This is the end-to-end check that the renderer, the detector, the tracker
/// and the stereo matcher all agree with each other: a renderer that draws
/// nothing trackable, or a matcher that rejects everything, shows up here as a
/// collapsed rate rather than as a wrong number somewhere downstream.
void testFrontendOverSimRun() {
  SimConfig c = cfg(2024, WorldType::kRooms);
  c.duration = 8.0;
  const Fingerprint fp = run(c, true);
  SC_CHECK(fp.frames > 50);
  SC_CHECK(fp.features > 0);
  // Most detected features must survive tracking: the frontend is only useful
  // if it carries a stable set between frames.
  const Scalar track_rate = static_cast<Scalar>(fp.tracked) / static_cast<Scalar>(fp.features);
  SC_CHECK(track_rate > 0.6);
  // And a real fraction of them must acquire a stereo depth. The exact share
  // depends on how much of the view is bare wall, so the bar is set where a
  // broken matcher (near zero) fails but a working one on a sparse scene passes.
  const Scalar depth_rate = static_cast<Scalar>(fp.stereo_valid) / static_cast<Scalar>(fp.features);
  SC_CHECK(depth_rate > 0.1);
}

/// Wall-clock must be far below real time on a single core, which is the
/// headroom the real-time budget needs on the reference hardware.
void testRunsFasterThanRealTime() {
  SimConfig c = cfg(31, WorldType::kRooms);
  c.duration = 10.0;
  Simulator sim(c);
  const SimResult r = sim.run(nullptr);
  SC_CHECK(r.wall_time_s > 0.0);
  // 10 s of simulated time with stereo rendering must fit in well under 10 s.
  SC_CHECK(r.wall_time_s < c.duration);
}

void run() {
  testSameSeedReproducesExactly();
  testDifferentSeedDiffers();
  testWorldTypeMatters();
  testFaultsAreDeterministic();
  testFrontendOverSimRun();
  testRunsFasterThanRealTime();
}

}  // namespace

SC_TEST_MAIN("test_determinism")
