// SLAMCore-8G -- fault injection.
//
// Every safety claim in the project contract is backed by a fault that is
// actually injected, not asserted. Faults act on the sensor streams and on the
// plant, and are scheduled deterministically so the same fault fires at the
// same simulated time on every machine.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "slamcore/core/stats.hpp"
#include "slamcore/sim/sensors.hpp"
#include "slamcore/sim/vehicle.hpp"

namespace sc {
namespace sim {

enum class FaultType {
  kNone,
  kImuFreeze,       ///< IMU output latched at its last value
  kImuSpike,        ///< large single-sample accelerometer outlier
  kImuDropout,      ///< IMU stops producing samples
  kCamBlackout,     ///< camera returns a dead frame
  kCamOcclusion,    ///< landmarks randomly hidden
  kOdomSlip,        ///< wheel slip outliers
  kOdomFreeze,      ///< odometry latched
  kActuatorLag,     ///< steering responds late (plant mismatch)
  kTextureless,     ///< landmarks stop being detected
  kExtremeNoise,    ///< all noise scaled up
};

const char* toString(FaultType t);
bool parseFaultType(const std::string& s, FaultType* out);

/// A scheduled fault window.
struct FaultEvent {
  FaultType type = FaultType::kNone;
  Scalar t_start = 0;
  Scalar duration = 2.0;
  Scalar magnitude = 1.0;

  bool active(Scalar t) const { return t >= t_start && t < t_start + duration; }
};

/// Deterministic fault schedule.
class FaultSchedule {
 public:
  void add(const FaultEvent& e) { events_.push_back(e); }
  void clear() { events_.clear(); }
  const std::vector<FaultEvent>& events() const { return events_; }

  /// True when any event of `type` is active at time `t`.
  bool active(FaultType type, Scalar t) const;
  Scalar magnitude(FaultType type, Scalar t) const;

  /// A standard suite: every fault class, well separated in time.
  static FaultSchedule standardSuite(Scalar start_time, Scalar spacing = 6.0);

  std::string describe() const;

 private:
  std::vector<FaultEvent> events_;
};

/// Applies a fault schedule to the sensor streams. Kept separate from the
/// emulators so the emulators stay pure noise models.
class FaultInjector {
 public:
  explicit FaultInjector(const FaultSchedule& schedule) : schedule_(schedule) {}

  /// Filter IMU samples for time `t`. May drop or corrupt samples.
  void applyImu(Scalar t, std::vector<ImuSample>* samples);

  /// Corrupt a rendered frame. Returns false when the frame should be dropped
  /// entirely (camera blackout), which the frontend must survive.
  bool applyCamera(Scalar t, Image* img);

  /// Corrupt an odometry sample.
  void applyOdom(Scalar t, OdomSample* s);

  /// Whether landmarks are suppressed this frame (textureless fault).
  bool suppressLandmarks(Scalar t) const { return schedule_.active(FaultType::kTextureless, t); }

  /// Extra noise multiplier applied to the frontend's measurement covariances,
  /// so a "noisier" world is reflected in the filter rather than silently
  /// absorbed.
  Scalar noiseScale(Scalar t) const { return 1.0 + 4.0 * schedule_.magnitude(FaultType::kExtremeNoise, t); }

  /// Actuator lag in seconds, applied by the plant.
  Scalar actuatorLag(Scalar t) const { return 0.25 * schedule_.magnitude(FaultType::kActuatorLag, t); }

  Size injectedCount() const { return injected_; }

 private:
  FaultSchedule schedule_;
  Vec3 last_imu_accel_ = Vec3::Zero();
  Vec3 last_imu_gyro_ = Vec3::Zero();
  bool have_last_imu_ = false;
  OdomSample last_odom_;
  Size injected_ = 0;
};

}  // namespace sim
}  // namespace sc
