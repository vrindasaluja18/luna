// LUNA -- rolling occupancy mapping under a hard RAM budget.
//
// The map is the largest thing the stack owns, and the 8 GB ceiling is a
// *process* budget shared with everything else. That rules out the usual
// "allocate a grid the size of the world" approach: an unbounded world means an
// unbounded allocation, and a 6-hour run has to stay inside the same ceiling as
// a 60-second one.
//
// The design is therefore tiled and paged:
//
//   * The map is a fixed-capacity cache of square tiles. Tiles are allocated
//     on first observation and evicted LRU once the cache is full. The cache
//     size is a constructor parameter, so peak RAM is a property of the
//     configuration rather than of the trajectory.
//   * Each tile is a log-odds occupancy grid. Inverse-sensor-model updates
//     accumulate log-odds and clamp them, which is what makes the map stable
//     under repeated observation instead of flickering.
//   * Evicted tiles are written to a memory-mapped file rather than freed.
//     Writing is what keeps a long run's map complete while its resident set
//     stays bounded; the mmap means a revisit faults the tile back in from the
//     page cache without a read syscall in the hot path.
//
// Everything here is allocation-free after `reserve()`: the tile pool, the
// index, and the per-update scratch are all preallocated.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "luna/core/se2.hpp"
#include "luna/core/types.hpp"
#include "luna/frontend/features.hpp"
#include "luna/sim/sensors.hpp"
#include "luna/sim/world.hpp"

namespace luna {
namespace mapping {

/// Occupancy of a cell in log-odds. 0 = unknown (p = 0.5), positive = occupied.
using LogOdds = float;

struct OccupancyOptions {
  /// Tile side length in cells. 64 cells at 5 cm is a 3.2 m tile: large enough
  /// that the per-tile bookkeeping is negligible, small enough that evicting
  /// one loses little local detail.
  int tile_cells = 64;
  Scalar resolution = 0.05;
  /// Resident tile cap. This, not the world size, is what bounds map RAM:
  /// `max_resident_tiles * tile_cells^2 * 4 bytes`. The default is 4096 tiles
  /// (about 42 MB of log-odds), covering a 200 x 200 m area before eviction.
  Size max_resident_tiles = 4096;
  /// Log-odds increment for a cell the beam passed through, and for the cell
  /// the beam terminated in. The asymmetry is the standard inverse sensor
  /// model: a hit is much stronger evidence than a miss, because a miss can be
  /// caused by a reflective or thin surface that the beam actually hit.
  LogOdds free_delta = -0.40f;
  LogOdds occupied_delta = 0.85f;
  /// Clamp on accumulated log-odds. Without it, a wall observed a thousand
  /// times becomes infinitely certain and can never be revised, so a dynamic
  /// obstacle would be permanently baked into the map.
  LogOdds clamp_min = -4.0f;
  LogOdds clamp_max = 4.0f;
  /// Log-odds above which a cell is reported occupied, below which free.
  LogOdds occupied_threshold = 0.6f;
  LogOdds free_threshold = -0.4f;
  /// Range at which a beam is considered to have hit nothing (max-range
  /// reading). Such a beam marks cells free but the terminal cell unknown.
  Scalar max_range = 12.0;
  /// Cells closer than this to the sensor are not updated: the beam origin is
  /// inside the robot footprint, and marking it free carves a hole around the
  /// vehicle that grows into a trail as it drives.
  Scalar min_range = 0.25;
  /// Path to the mmap backing file. Empty disables paging (evicted tiles are
  /// simply dropped), which is the right behaviour for a bounded world.
  std::string backing_path;
  /// Whether to write evicted tiles to the backing file.
  bool page_to_disk = false;
};

/// A ray to be integrated. Sensor-frame direction plus measured range.
struct RangeReading {
  Vec2 dir = Vec2(1, 0);  ///< unit direction, sensor frame
  Scalar range = 0;       ///< measured distance, metres
};

struct OccupancyStats {
  Size tiles_allocated = 0;
  Size tiles_resident = 0;
  Size tiles_evicted = 0;
  Size tiles_loaded = 0;
  Size cells_updated = 0;
  Size rays_integrated = 0;
  Size bytes_resident = 0;
  Size bytes_peak = 0;
};

/// Tile-keyed occupancy map with LRU eviction and optional mmap paging.
class OccupancyMap {
 public:
  explicit OccupancyMap(const OccupancyOptions& opt = OccupancyOptions());

  /// Integrate one scan. `pose` is the sensor pose in world coordinates; the
  /// readings are in the sensor frame. Beams are traced with a DDA walk.
  void integrateScan(const Pose2& pose, const std::vector<RangeReading>& readings);

  /// Integrate a single beam. Exposed so a caller with its own scan model (a
  /// depth image, a lidar sweep) does not have to materialise a vector.
  void integrateRay(const Pose2& pose, const Vec2& dir, Scalar range);

  /// Occupancy at a world point, as log-odds. Unknown cells return 0.
  LogOdds logOddsAt(const Vec2& p) const;

  /// Classified occupancy at a world point.
  sim::Cell atWorld(const Vec2& p) const;

  /// Extract a dense occupancy grid covering the observed area, for the
  /// topology and planning modules to consume. `res` and `pad_m` control the
  /// raster; the grid is in the same frame as the world.
  sim::GridMap toGrid(Scalar res, Scalar pad_m = 2.0) const;

  /// Extract a dense grid over an explicit window. Cheaper than toGrid when the
  /// caller already knows the region of interest.
  sim::GridMap toGridWindow(const Vec2& lo, const Vec2& hi, Scalar res) const;

  /// Mark a disc as free or occupied directly. Used by the closed-loop harness
  /// to fold the controller's own footprint knowledge back into the map.
  void markDisc(const Vec2& center, Scalar radius, bool occupied);

  /// Drop every resident tile, keeping statistics. Paged tiles stay on disk.
  void clear();

  /// Flush dirty tiles to the backing file and release them.
  void flush();

  const OccupancyStats& stats() const { return stats_; }
  const OccupancyOptions& options() const { return opt_; }

  /// World-space bounds of every tile ever touched. Grows monotonically.
  const Vec2& minSeen() const { return min_seen_; }
  const Vec2& maxSeen() const { return max_seen_; }

  /// Number of cells currently marked occupied in the resident set.
  Size countOccupied() const;
  /// Number of cells currently marked free in the resident set.
  Size countFree() const;

 private:
  struct Tile {
    std::vector<LogOdds> data;
    std::uint64_t last_use = 0;
    bool dirty = false;
    bool paged = false;  ///< has a copy in the backing file
  };

  struct Key {
    int tx = 0, ty = 0;
    bool operator==(const Key& o) const { return tx == o.tx && ty == o.ty; }
  };
  struct KeyHash {
    Size operator()(const Key& k) const {
      return static_cast<Size>(static_cast<std::uint64_t>(static_cast<std::uint32_t>(k.tx)) *
                                   0x9E3779B97F4A7C15ull ^
                               static_cast<std::uint64_t>(static_cast<std::uint32_t>(k.ty)) *
                                   0xC2B2AE3D27D4EB4Full);
    }
  };

  Tile* tile(int tx, int ty, bool create);
  Tile* findTile(int tx, int ty) const;
  void evictIfNeeded();
  void evictOne();
  bool loadPaged(const Key& k, Tile* t) const;
  void writePaged(const Key& k, const Tile& t) const;
  void touch(Tile* t);
  void bumpPeak();

  OccupancyOptions opt_;
  std::unordered_map<Key, Tile, KeyHash> tiles_;
  std::uint64_t clock_ = 0;
  OccupancyStats stats_;
  Vec2 min_seen_ = Vec2::Constant(std::numeric_limits<Scalar>::max());
  Vec2 max_seen_ = Vec2::Constant(-std::numeric_limits<Scalar>::max());
  /// Scratch for the beam walk, reused so integration allocates nothing.
  mutable std::vector<Key> touched_;
};

/// Build the inverse sensor model readings from a simulated stereo depth frame.
///
/// This is the bridge from the frontend to the map: each feature with a valid
/// depth becomes one beam whose range is the feature's depth. Features without
/// depth contribute nothing, which is the honest behaviour -- a bearing-only
/// observation says nothing about where the obstacle is.
std::vector<RangeReading> readingsFromFeatures(const std::vector<frontend::Feature>& feats,
                                               const sim::CameraModel& cam);

}  // namespace mapping
}  // namespace luna
