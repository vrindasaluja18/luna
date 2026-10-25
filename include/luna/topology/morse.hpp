// LUNA -- Morse-theoretic topology extraction.
//
// A metric map ("there is a wall 3 m to the left") is not enough to plan
// through: the planner also needs the *structure* of the free space -- which
// rooms exist, which doorways join them, and which of those connections are
// real rather than noise. This module extracts that structure with Morse
// theory.
//
// The construction, in order:
//
//   1. Height function. h(p) = clearance to the nearest obstacle, smoothed.
//      Rooms are local maxima (far from every wall), doorways are saddles
//      (the ridge you must cross between two rooms), wall tips are minima.
//
//   2. Critical points and their persistence. Cells are flooded in descending
//      h with a union-find over the superlevel sets. Each merge of two
//      components happens at a saddle cell and kills the younger component's
//      maximum; the difference in height is that maximum's *persistence*, the
//      standard measure of how prominent the feature is. This gives the join
//      tree directly, without a separate critical-point detector: a cell that
//      merges two components has two or more lower-link components, which is
//      exactly the discrete saddle criterion, and a cell whose neighbours are
//      all lower is a maximum.
//
//   3. Persistence filtering. A bump 2 cm deep in the clearance field is
//      noise, not a room. Critical points below the persistence threshold are
//      merged away rather than reported, which is what makes the corridor
//      graph stable across seeds. Filtering happens *inside* the flooding, so
//      a filtered maximum simply never appears and its basin is absorbed by
//      its neighbour.
//
//   4. Corridor graph. Nodes are the surviving maxima (rooms) and saddles
//      (doorways); edges are the connections recorded when a component merged
//      with significant persistence. For a simply-connected domain this graph
//      is a tree, and its edges are exactly the Morse-Smale 1-cells that a
//      route must traverse.
//
// The result is exposed as a plain graph so the planner can run A* over it.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "luna/core/se3.hpp"
#include "luna/core/types.hpp"
#include "luna/sim/world.hpp"

namespace luna {
namespace topology {

/// Morse index of a critical point, i.e. the number of negative Hessian
/// eigenvalues. 0 = minimum, 1 = saddle, 2 = maximum.
enum class CriticalType : int { kMinimum = 0, kSaddle = 1, kMaximum = 2 };

const char* toString(CriticalType t);

struct CriticalPoint {
  Vec2 p = Vec2::Zero();
  int cx = 0, cy = 0;     ///< grid cell
  CriticalType type = CriticalType::kSaddle;
  Scalar h = 0;           ///< height (clearance) at the point
  /// Persistence, metres. For a maximum or minimum this is the height
  /// difference to the saddle at which its basin merges with a deeper one; for
  /// a saddle it is the persistence of the pair it forms. A critical point that
  /// survives the filter but has no partner (the global extremum) carries the
  /// height difference to the far end of its component.
  Scalar persistence = 0;
};

struct CorridorEdge {
  int a = 0, b = 0;        ///< node indices into the critical-point list
  int saddle = -1;         ///< node index of the saddle this corridor crosses
  Scalar length = 0;       ///< metres, along the corridor
  /// Clearance at the saddle: half the true width of the doorway. This is the
  /// bottleneck the vehicle must fit through, so the planner can reject a
  /// corridor that is too tight without a second map query.
  Scalar bottleneck = 0;
  Scalar persistence = 0;
};

struct MorseOptions {
  /// Gaussian smoothing applied to the height field before analysis, in cells.
  /// The clearance field is piecewise-smooth but has creases along the medial
  /// axis of every corner, and those creases generate spurious critical points
  /// that persistence alone will not remove. Smoothing first is what makes the
  /// extraction stable.
  Scalar smoothing_sigma = 1.2;
  /// Features with persistence below this (metres) are discarded as noise.
  /// The right value scales with the world, not with the resolution: it is the
  /// smallest room-to-doorway height difference that should count as a room.
  Scalar min_persistence = 0.20;
  /// Cells with clearance below this are treated as obstacles for the
  /// analysis. This is what keeps the corridor graph inside drivable space
  /// instead of threading the 5 cm skin along every wall.
  Scalar min_clearance = 0.22;
  /// A corridor is a doorway only if its bottleneck is this fraction of the
  /// narrower of the two regions it joins. A large open hall has many local
  /// maxima in its clearance field -- the room centre is a broad, nearly flat
  /// plateau -- and their mutual persistence is a genuine metre or so, so a
  /// persistence threshold alone reports them as separate rooms. What actually
  /// distinguishes a doorway is *contrast*: the vehicle must squeeze, so the
  /// gap is much narrower than the space on either side. Regions joined by a
  /// corridor that is nearly as wide as they are get merged into one room,
  /// which is the physically meaningful answer.
  Scalar min_bottleneck_contrast = 0.75;
  /// Include minima (dead ends, wall tips) as graph nodes. Off by default:
  /// they are useful for diagnostics and for a "back out of a dead end"
  /// behaviour, but they are not needed for point-to-point routing.
  bool include_minima = false;
  /// Safety valve. A pathological field can produce a critical point per cell;
  /// beyond this many the analysis stops adding nodes rather than growing
  /// without bound.
  Size max_nodes = 20000;
};

struct MorseStats {
  Size free_cells = 0;
  Size raw_maxima = 0;
  Size raw_saddles = 0;
  Size kept_maxima = 0;
  Size kept_saddles = 0;
  Size kept_minima = 0;
  Size corridors = 0;
  Scalar extract_ms = 0;
  Scalar max_persistence = 0;
};

/// Corridor graph: a plain weighted graph plus spatial lookup.
class CorridorGraph {
 public:
  std::vector<CriticalPoint> nodes;
  std::vector<CorridorEdge> edges;
  /// node -> indices into `edges`.
  std::vector<std::vector<int>> adjacency;

  Size nodeCount() const { return nodes.size(); }
  Size edgeCount() const { return edges.size(); }

  /// Index of the node nearest `p` within `max_dist`, or -1.
  int nearestNode(const Vec2& p, Scalar max_dist) const;

  /// Indices of the nodes of a given type, in order.
  std::vector<int> nodesOfType(CriticalType t) const;

  /// A* over the corridor graph. Returns false when start and goal are in
  /// different components. `route` receives the node indices in order,
  /// including the start and goal nodes.
  bool route(int start_node, int goal_node, std::vector<int>* route) const;

  /// Largest clearance among nodes: the widest room. Used to pick a fallback
  /// target when the goal itself is unreachable.
  int widestNode() const;

  /// Serialise for the benchmark report: one line per node and per edge.
  std::string toDot() const;
};

class MorseComplex {
 public:
  /// Extract from a metric map. The height function is the clearance field of
  /// the grid's free space, computed internally by an exact Euclidean distance
  /// transform, so the caller does not have to supply one.
  static MorseComplex extract(const sim::GridMap& grid, const MorseOptions& opt = MorseOptions());

  /// Extract from an explicit scalar field, one value per cell in row-major
  /// order. `field[i]` must be the height of cell i; cells marked occupied are
  /// forced to -infinity. This is the entry point for analytic tests
  /// (f = x^3 - 3x, monkey saddle) and for a cost map other than clearance.
  static MorseComplex analyseField(const sim::GridMap& grid, const std::vector<Scalar>& field,
                                   const MorseOptions& opt = MorseOptions());

  const CorridorGraph& graph() const { return graph_; }
  /// The smoothed height field actually analysed, row-major over the grid.
  const std::vector<Scalar>& heightField() const { return height_; }
  const sim::GridMap& grid() const { return grid_; }
  const MorseStats& stats() const { return stats_; }
  const MorseOptions& options() const { return opt_; }

  /// Height at a world point, from the analysed field.
  Scalar heightAt(const Vec2& p) const;

  /// ASCII rendering of the corridor graph over the occupancy map, for
  /// debugging and for the docs. '#' wall, ' ' free, 'O' room, 'x' doorway,
  /// '.' minimum, '*' corridor path.
  std::string toAscii(int max_cols = 96) const;

 private:
  sim::GridMap grid_;
  std::vector<Scalar> height_;
  CorridorGraph graph_;
  MorseStats stats_;
  MorseOptions opt_;
};

/// Exact Euclidean distance transform: for every free cell, the distance in
/// metres to the nearest occupied cell. Two-pass Felzenszwalb-Huttenlocher,
/// O(n). Exposed because the planner and the mapping module both need a
/// clearance field and neither should reimplement this.
void clearanceTransform(const sim::GridMap& grid, std::vector<Scalar>* out);

/// Separable Gaussian blur of a scalar field on a grid. Cells outside the grid
/// and occupied cells are excluded from the average (treated as missing), so
/// the blur does not drag the wall value into the free space.
void gaussianBlurField(const sim::GridMap& grid, const std::vector<Scalar>& in, Scalar sigma,
                       std::vector<Scalar>* out);

}  // namespace topology
}  // namespace luna
