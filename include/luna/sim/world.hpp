// LUNA -- 2D world model and procedural generator.
//
// The simulator exists because public datasets cannot exercise every failure
// mode the stack must survive: textureless corridors, 1.5x-width bottlenecks,
// dynamic clutter, sensor dropouts. Worlds are generated from an explicit
// seed so every benchmark figure is reproducible.
//
// The generator deliberately produces worlds with *Morse structure*: rooms
// (local maxima of the clearance function) joined by doorways (saddles). That
// is what gives the topology module something real to extract, and it is what
// makes the corridor-graph planner's behaviour checkable against ground truth.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "luna/core/se3.hpp"
#include "luna/core/stats.hpp"
#include "luna/core/types.hpp"

namespace luna {
namespace sim {

/// Occupancy value stored per cell. 0 = free, 100 = occupied, 255 = unknown.
enum class Cell : std::uint8_t { kFree = 0, kOccupied = 100, kUnknown = 255 };

/// Dense 2D occupancy grid. Resolution 5 cm by default: fine enough that a
/// 0.6 m robot has ~12 cells of body width, coarse enough that a 100 m world
/// is 4M cells (~4 MB) instead of something that would blow the RAM budget.
struct GridMap {
  int w = 0, h = 0;
  Scalar res = 0.05;
  Vec2 origin = Vec2::Zero();
  std::vector<std::uint8_t> cells;

  void resize(int width, int height, Scalar resolution, const Vec2& org);
  bool inBounds(int x, int y) const { return x >= 0 && y >= 0 && x < w && y < h; }
  Cell at(int x, int y) const;
  void set(int x, int y, Cell c);
  Cell atWorld(const Vec2& p) const;
  bool isFreeWorld(const Vec2& p) const { return atWorld(p) == Cell::kFree; }
  bool isOccupiedWorld(const Vec2& p) const { return atWorld(p) == Cell::kOccupied; }
  Vec2 cellCenter(int x, int y) const;
  Eigen::Vector2i toCell(const Vec2& p) const;
  Size countFree() const;
  Size bytes() const { return cells.size(); }
  Vec2 minCorner() const { return origin; }
  Vec2 maxCorner() const { return origin + Vec2(w * res, h * res); }

  /// Paint a filled axis-aligned rectangle in world coordinates.
  void fillRect(const Vec2& lo, const Vec2& hi, Cell c);
  /// Paint a disc.
  void fillDisc(const Vec2& center, Scalar radius, Cell c);
  /// Thicken all occupied cells by `radius` (used for clearance/inflation).
  void dilate(Scalar radius);
};

enum class WorldType {
  kRooms,         ///< grid of rectangular rooms joined by doorways
  kCorridorLoop,  ///< closed corridor loop with branches (loop-closure test)
  kSaddleHall,    ///< two halls joined by a narrow gap (Morse saddle test)
  kMaze,          ///< recursive-division maze
  kCluttered,     ///< open hall with random obstacles
  kTextureless,   ///< long straight corridor with no landmarks (stress test)
};

const char* toString(WorldType t);
bool parseWorldType(const std::string& s, WorldType* out);

/// A visual landmark. Descriptors are *generated*, not stored: `descSeed`
/// feeds a deterministic PRNG so the same landmark always yields the same
/// BRIEF-like bit vector. That is what makes loop closure ground-truthable.
struct Landmark {
  Vec2 p = Vec2::Zero();
  std::uint32_t id = 0;
  std::uint32_t descSeed = 0;
  Scalar height = 0.5;  ///< metres above the floor
  bool textured = true;
};

struct WorldSpec {
  WorldType type = WorldType::kRooms;
  Scalar size_m = 40.0;        ///< world is size_m x size_m
  Scalar res = 0.05;
  std::uint64_t seed = 1;
  int rooms = 6;
  /// Doorway width. Must clear the vehicle's diagonal *and* leave room for the
  /// turn the route takes through it. A doorway is where the path changes
  /// direction, so it has to be wide enough for the vehicle's turning circle,
  /// not merely wide enough for the vehicle: at a 0.55 m minimum turn radius a
  /// 90-degree turn sweeps about 1.1 m laterally, so a 1.2 m door cannot be
  /// negotiated at all however well the tracker performs.
  Scalar door_width = 1.8;
  Scalar landmark_density = 2.0;  ///< landmarks per m^2 of free space
  Scalar clutter = 0.0;        ///< fraction of free space filled with obstacles
  bool textureless = false;
};

struct WorldStats {
  Size freeCells = 0;
  Size occupiedCells = 0;
  Scalar freeArea = 0;
  Scalar corridorWidthMin = 0;
  Scalar corridorWidthMedian = 0;
  Size landmarks = 0;
};

/// A generated world: occupancy grid plus landmark field.
class World {
 public:
  static World generate(const WorldSpec& spec);

  const GridMap& grid() const { return grid_; }
  const std::vector<Landmark>& landmarks() const { return landmarks_; }
  const WorldSpec& spec() const { return spec_; }
  const WorldStats& stats() const { return stats_; }

  /// Line-of-sight test by grid ray marching (Amanatides-Woo DDA).
  bool lineOfSight(const Vec2& a, const Vec2& b) const;

  /// Distance to the first occupied cell along `dir`, capped at max_range.
  /// Returns max_range when nothing is hit.
  Scalar raycast(const Vec2& from, const Vec2& dir, Scalar max_range) const;

  /// Clearance field: distance in metres from each cell to the nearest
  /// occupied cell. This is the height function the topology module analyses,
  /// so it is computed once at generation time and shared.
  const std::vector<Scalar>& clearance() const { return clearance_; }
  Scalar clearanceAtWorld(const Vec2& p) const;

  /// Minimum clearance along a path, doubled to give a corridor width. Used to
  /// report the true bottleneck width of a generated route.
  Scalar minCorridorWidth(const std::vector<Vec2>& path) const;

  /// Breadth-first shortest grid path between two world points. Returns an
  /// empty vector when unreachable. Used to build ground-truth routes.
  ///
  /// `clearance_m` keeps the path at least that far from any obstacle, which
  /// is what makes a scripted route drivable by a finite-width vehicle.
  std::vector<Vec2> gridPath(const Vec2& from, const Vec2& to, Scalar clearance_m = 0.0) const;

  /// Round the corners of a dense path into arcs a car can actually follow.
  ///
  /// A 4-connected grid path is a staircase: its turns are 90 degrees, which no
  /// finite-turn-radius vehicle can track, so a scripted run collides at every
  /// corner even though the path itself is collision-free. Smoothing moves each
  /// interior point toward the midpoint of its neighbours, which rounds the
  /// staircase into a polyline with bounded curvature, and rejects any move
  /// that would drop the point below `clearance_m`.
  /// `min_turn_radius` bounds the curvature: smoothing continues until every
  /// corner is at least that gentle, or the iteration budget runs out. Without
  /// it the corners are rounded but arbitrarily tight, and a vehicle whose
  /// minimum turn radius is larger simply leaves the path on the first hard
  /// corner and drives away in a straight line.
  std::vector<Vec2> smoothPath(const std::vector<Vec2>& path, Scalar clearance_m,
                               Scalar min_turn_radius = 0.0, Scalar spacing = 0.15,
                               int iterations = 200) const;

  /// Reduce a dense path to the minimal set of waypoints whose straight-line
  /// connections remain collision-free *with at least `clearance_m` margin*.
  ///
  /// The margin is not cosmetic: plain line-of-sight permits a segment that
  /// grazes a wall corner, which is free of occupied cells but too tight for a
  /// finite-width vehicle to drive. Scripted runs then report collisions that
  /// look like controller failures but are really route defects.
  std::vector<Vec2> shortcut(const std::vector<Vec2>& path, Scalar clearance_m = 0.0) const;

  /// Render the world as an ASCII picture (docs, debugging, test fixtures).
  std::string toAscii(int max_cols = 96) const;

  /// Serialise to a compact PGM-like binary blob for the replay cache.
  std::vector<std::uint8_t> serialise() const;

 private:
  WorldSpec spec_;
  GridMap grid_;
  std::vector<Scalar> clearance_;
  std::vector<Landmark> landmarks_;
  WorldStats stats_;
};

/// Sample `n` points uniformly at random from free space, at least `min_sep`
/// apart. Deterministic for a given RNG state.
///
/// `min_clearance_m` rejects points that are free but too tight to drive
/// through. Sampling anywhere free is not enough: the generator leaves thin
/// slivers between walls, and an anchor that lands in one cannot be connected
/// to any other anchor by a path wide enough for the vehicle, so the whole
/// route collapses to that single point.
std::vector<Vec2> sampleFreePoints(const GridMap& grid, int n, Scalar min_sep, Rng* rng,
                                   Scalar min_clearance_m = 0.0,
                                   const std::vector<Scalar>* clearance = nullptr);

}  // namespace sim
}  // namespace luna
