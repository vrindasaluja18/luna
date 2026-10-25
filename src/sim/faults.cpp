#include "slamcore/sim/faults.hpp"

#include <algorithm>
#include <sstream>

namespace sc {
namespace sim {

const char* toString(FaultType t) {
  switch (t) {
    case FaultType::kNone: return "none";
    case FaultType::kImuFreeze: return "imu_freeze";
    case FaultType::kImuSpike: return "imu_spike";
    case FaultType::kImuDropout: return "imu_dropout";
    case FaultType::kCamBlackout: return "cam_blackout";
    case FaultType::kCamOcclusion: return "cam_occlusion";
    case FaultType::kOdomSlip: return "odom_slip";
    case FaultType::kOdomFreeze: return "odom_freeze";
    case FaultType::kActuatorLag: return "actuator_lag";
    case FaultType::kTextureless: return "textureless";
    case FaultType::kExtremeNoise: return "extreme_noise";
  }
  return "?";
}

bool parseFaultType(const std::string& s, FaultType* out) {
  static const std::pair<const char*, FaultType> kTable[] = {
      {"none", FaultType::kNone},
      {"imu_freeze", FaultType::kImuFreeze},
      {"imu_spike", FaultType::kImuSpike},
      {"imu_dropout", FaultType::kImuDropout},
      {"cam_blackout", FaultType::kCamBlackout},
      {"cam_occlusion", FaultType::kCamOcclusion},
      {"odom_slip", FaultType::kOdomSlip},
      {"odom_freeze", FaultType::kOdomFreeze},
      {"actuator_lag", FaultType::kActuatorLag},
      {"textureless", FaultType::kTextureless},
      {"extreme_noise", FaultType::kExtremeNoise},
  };
  for (const auto& e : kTable) {
    if (s == e.first) {
      *out = e.second;
      return true;
    }
  }
  return false;
}

bool FaultSchedule::active(FaultType type, Scalar t) const {
  for (const auto& e : events_)
    if (e.type == type && e.active(t)) return true;
  return false;
}

Scalar FaultSchedule::magnitude(FaultType type, Scalar t) const {
  Scalar m = 0;
  for (const auto& e : events_)
    if (e.type == type && e.active(t)) m = std::max(m, e.magnitude);
  return m;
}

FaultSchedule FaultSchedule::standardSuite(Scalar start_time, Scalar spacing) {
  FaultSchedule s;
  const FaultType types[] = {
      FaultType::kImuFreeze,  FaultType::kImuSpike,   FaultType::kImuDropout,
      FaultType::kCamBlackout, FaultType::kCamOcclusion, FaultType::kOdomSlip,
      FaultType::kOdomFreeze, FaultType::kActuatorLag, FaultType::kTextureless,
      FaultType::kExtremeNoise,
  };
  Scalar t = start_time;
  for (FaultType ft : types) {
    FaultEvent e;
    e.type = ft;
    e.t_start = t;
    e.duration = 2.0;
    e.magnitude = 1.0;
    s.add(e);
    t += spacing;
  }
  return s;
}

std::string FaultSchedule::describe() const {
  std::ostringstream os;
  for (const auto& e : events_) {
    os << toString(e.type) << " @ " << e.t_start << "s for " << e.duration << "s\n";
  }
  return os.str();
}

void FaultInjector::applyImu(Scalar t, std::vector<ImuSample>* samples) {
  if (samples->empty()) return;

  if (schedule_.active(FaultType::kImuDropout, t)) {
    samples->clear();
    ++injected_;
    return;
  }
  if (schedule_.active(FaultType::kImuFreeze, t)) {
    if (!have_last_imu_) {
      last_imu_accel_ = samples->back().accel;
      last_imu_gyro_ = samples->back().gyro;
      have_last_imu_ = true;
    }
    for (auto& s : *samples) {
      s.accel = last_imu_accel_;
      s.gyro = last_imu_gyro_;
    }
    ++injected_;
    return;
  }
  if (schedule_.active(FaultType::kImuSpike, t)) {
    // One large outlier in the middle of the window.
    const Size mid = samples->size() / 2;
    const Scalar mag = 8.0 * schedule_.magnitude(FaultType::kImuSpike, t);
    (*samples)[mid].accel += Vec3(mag, -mag * 0.5, mag * 0.25);
    ++injected_;
    return;
  }
  if (schedule_.active(FaultType::kExtremeNoise, t)) {
    const Scalar scale = noiseScale(t);
    for (auto& s : *samples) {
      s.accel *= scale;
      s.gyro *= scale;
    }
    ++injected_;
    return;
  }
  have_last_imu_ = true;
  last_imu_accel_ = samples->back().accel;
  last_imu_gyro_ = samples->back().gyro;
}

bool FaultInjector::applyCamera(Scalar t, Image* img) {
  if (schedule_.active(FaultType::kCamBlackout, t)) {
    std::fill(img->px.begin(), img->px.end(), static_cast<std::uint8_t>(0));
    ++injected_;
    return false;
  }
  if (schedule_.active(FaultType::kExtremeNoise, t)) {
    // Salt-and-pepper at high density: destroys feature matching.
    Rng r(static_cast<std::uint64_t>(t * 1000.0) + 7u);
    const Scalar p = 0.05 * noiseScale(t);
    for (auto& px : img->px) {
      if (r.uniform() < p) px = (r.uniform() < 0.5) ? 0 : 255;
    }
    ++injected_;
  }
  return true;
}

void FaultInjector::applyOdom(Scalar t, OdomSample* s) {
  if (schedule_.active(FaultType::kOdomFreeze, t)) {
    if (last_odom_.t == 0 && last_odom_.d_forward == 0 && last_odom_.d_yaw == 0) last_odom_ = *s;
    s->d_forward = last_odom_.d_forward;
    s->d_yaw = last_odom_.d_yaw;
    s->valid = false;
    ++injected_;
    return;
  }
  if (schedule_.active(FaultType::kOdomSlip, t)) {
    s->d_forward *= 3.0 * schedule_.magnitude(FaultType::kOdomSlip, t);
    s->d_yaw += 0.4 * schedule_.magnitude(FaultType::kOdomSlip, t);
    s->valid = false;
    ++injected_;
    return;
  }
  last_odom_ = *s;
}

}  // namespace sim
}  // namespace sc
