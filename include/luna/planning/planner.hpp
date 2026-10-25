// LUNA -- corridor-graph planning and B-spline smoothing.
//
// The planner works on two levels, and the split is the point of the topology
// module:
//
//   * Global: A* over the Morse corridor graph. The graph has a node per room
//     and an edge per doorway, so a 40 m world with nine rooms is a nine-node
//     problem -- solved in microseconds, and, more importantly, *stable*. A
//     grid A* over the same world re-plans a visibly different route every time
//     one occupancy cell flips; the corridor route only changes when a room or
//     a doorway genuinely appears or disappears.
//
//   * Local: the corridor route is expanded into waypoints and smoothed into a
//     B-spline. A corridor centreline is a sequence of room centres joined at
//     doorways, which is a polyline with hard corners; a car cannot follow a
//     corner. The spline interpolates the waypoints with bounded curvature,
//     which is what makes the result drivable.
//
// Replanning is event-driven rather than periodic: it runs when the corridor
// route changes, when the goal changes, or when the vehicle leaves the corridor
// it is following. A planner that re-solves every tick is both wasteful and
// jittery, and its latency claim is meaningless.
#pragma once

#include <chrono>
#include <queue>
#include <vector>

#include "luna/core/se3.hpp"
#include "luna/core/types.hpp"
#include "luna/sim/world.hpp"
#include "luna/topology/morse.hpp"

namespace luna {
namespace planning {

struct PlanOptions {
  /// Preferred clearance from obstacles for the local path, metres. The
  /// planner falls back to a smaller margin when the corridor is tight, down
  /// to `min_clearance`.
  ///
  /// These must exceed the vehicle's *half-diagonal*, not its half-width: the
  /// body is a rectangle and it is the corner that hits first when the vehicle
  /// is not aligned with the wall. For the default 0.34 x 0.62 m body the
  /// half-diagonal is 0.36 m, so a 0.30 m margin produced paths that the
  /// controller tracked perfectly and that still put a corner inside a wall.
  Scalar preferred_clearance = 0.70;
  /// Hard floor on the clearance the plan will accept, metres.
  ///
  /// This is not a preference, it is the number that makes "no collisions" a
  /// property rather than a hope, and it has to be set from the controller's
  /// measured tracking error rather than from the body alone. A plan point with
  /// clearance `c` guarantees the *body centre* is `c - half_diagonal` from the
  /// wall, but the controller does not hold the centre on the plan; it tracks
  /// it to within some RMS error `e`. The body corner therefore comes within
  /// `c - half_diagonal - e` of the wall, and the plan only prevents collisions
  /// if that stays positive:
  ///
  ///     c > half_diagonal + e
  ///
  /// For the default body the half-diagonal is 0.354 m, so a 0.52 m floor
  /// leaves 0.166 m -- less than the loop's own tracking error -- and the
  /// vehicle grazed a wall on 42 consecutive ticks while the plan reported
  /// itself valid. This default applies to analytic maps, where 0.62 m is
  /// available on every route; a reconstructed map cannot supply it (see the
  /// note in the navigation options) and the controller's speed profile, which
  /// carries the same budget, is what keeps the body off the wall there.
  Scalar min_clearance = 0.62;
  /// Waypoint spacing along the expanded corridor route, metres.
  Scalar waypoint_spacing = 0.6;
  /// B-spline samples per waypoint span. Higher is smoother and costs more.
  int spline_samples = 8;
  /// Maximum curvature of the smoothed path, 1/m. The default is the inverse of
  /// the vehicle's minimum turn radius, so a path that satisfies it is
  /// trackable; the smoother relaxes the waypoints toward straight lines until
  /// every corner meets it.
  /// Maximum path curvature, 1/m. The vehicle's minimum turn radius is
  /// L / tan(max_steer) = 0.5 / 0.93 = 0.54 m, i.e. curvature 1.86. Asking the
  /// smoother for 1.8 leaves essentially no steering margin: every corner is
  /// taken at full lock, so there is no authority left to correct tracking
  /// error and the loop drifts wide. A limit of 0.9 (1.1 m radius) keeps half
  /// the steering range in reserve.
  Scalar max_curvature = 0.9;
  /// Steps of clearance ascent used to re-centre the path in corridors.
  int clearance_iterations = 8;
  /// Largest turn a single path point may contain, radians. The clearance
  /// relaxation moves points toward the medial axis one at a time, so it has to
  /// refuse moves that fold the path; this is the bound it checks. It is a
  /// per-point *corner* bound rather than a curvature bound because a discrete
  /// path has no curvature at a point, only an angle between the two segments
  /// that meet there.
  Scalar max_turn_per_point = 0.9;
  /// Iterations of the curvature-relaxation loop.
  int smooth_iterations = 60;
  /// Reject a plan whose route crosses an occupied cell. Off by default: the
  /// topology graph is built from the same map, so a crossing means the map
  /// changed under the plan, which is exactly when replanning is wanted.
  bool verify_collision_free = true;
  /// Radius of the vehicle body, metres: half its diagonal, so the corner is
  /// covered when the body is not aligned with the wall.
  Scalar body_half_diagonal = 0.353;
  /// Slack the controller keeps above `body_half_diagonal` before it refuses to
  /// move, metres. It must match `ControllerOptions`' own stop threshold: the
  /// planner reports `kOk`, the controller sees a margin of zero and safe-stops,
  /// and the run is scored as a success that drove nothing. Keeping the two in
  /// step is what makes `kOk` mean "this can be driven".
  Scalar body_margin = 0.02;
  /// The clearance a plan must hold at *every* point to be accepted, metres.
  ///
  /// This is the physical drivability floor, not the search target. A plan point
  /// with clearance `c` leaves the body centre `c - half_diagonal` from the wall,
  /// and the controller refuses to move while that is at or below `body_margin`,
  /// so below `half_diagonal + body_margin` the vehicle will not drive the plan
  /// at all. Reporting `kOk` there is a false success: the run drives zero
  /// metres and the benchmark counts it as a win.
  ///
  /// `min_clearance` is deliberately *not* folded in here. It is the margin the
  /// grid search routes through, and a reconstructed map often cannot supply it
  /// -- its doorways come out narrower than the walls that produced them. A plan
  /// that dips below `min_clearance` between grid cells but still clears the body
  /// is drivable, and rejecting it would report a route as blocked that the
  /// vehicle can actually follow.
  Scalar acceptanceClearance() const { return body_half_diagonal + body_margin; }
};

enum class PlanStatus {
  kOk,
  kNoCorridorRoute,   ///< start and goal are not connected in the corridor graph
  kNoPath,            ///< the local path could not be built
  kBlocked,           ///< the path crosses an obstacle
  kInvalid,
};

const char* toString(PlanStatus s);

struct PlanStats {
  int corridor_nodes = 0;
  int corridor_edges = 0;
  int waypoints = 0;
  Scalar corridor_route_length = 0;
  Scalar path_length = 0;
  Scalar min_clearance_along = 0;
  Scalar max_curvature_along = 0;
  Scalar plan_ms = 0;
  /// Bottleneck (clearance) of the tightest doorway the route crosses.
  Scalar route_bottleneck = 0;
};

struct Plan {
  PlanStatus status = PlanStatus::kInvalid;
  /// Smoothed, drivable path in world coordinates.
  std::vector<Vec2> path;
  /// The corridor route that produced it: node indices into the corridor graph.
  std::vector<int> corridor_route;
  /// The raw waypoints before smoothing, kept for diagnostics and for the
  /// controller to reason about corners.
  std::vector<Vec2> waypoints;
  /// Distance from each `path` point to the nearest obstacle, metres. Carried
  /// on the plan so the controller can bound its speed by the room it actually
  /// has without needing the occupancy map.
  std::vector<Scalar> clearance;
  PlanStats stats;
  bool ok() const { return status == PlanStatus::kOk; }
};

/// Plans over a corridor graph.
class CorridorPlanner {
 public:
  CorridorPlanner(const topology::CorridorGraph& graph, const sim::GridMap& grid,
                  const PlanOptions& opt = PlanOptions());

  /// Plan from a world point to a world point. Start and goal are snapped to
  /// the nearest corridor node within `snap_radius`; a goal that is not near
  /// any node is connected to the nearest one directly.
  Plan plan(const Vec2& start, const Vec2& goal, Scalar snap_radius = 2.5);

  /// Whether the straight line from `a` to `b` stays at least `clearance`
  /// metres from every obstacle.
  bool segmentClear(const Vec2& a, const Vec2& b, Scalar clearance) const;

  /// Label the corridor nodes by the connectivity the occupancy grid can
  /// actually realise.
  ///
  /// The topology graph and the occupancy grid are built from the same field but
  /// answer different questions. The graph says "these two regions are joined by
  /// a doorway"; the grid says "there is a line of free cells between them". On a
  /// clean map they agree. On a reconstructed one they do not: a doorway that
  /// the topology extractor accepted can be a single row of cells that the
  /// planner's clearance test rejects, and then the graph advertises a route
  /// that A* cannot follow. The planner falls back to a straight segment between
  /// the two nodes, the segment crosses the wall the doorway was supposed to be
  /// in, and the result is reported as blocked.
  ///
  /// Labelling by the grid's own connectivity removes that class of failure at
  /// the source: an edge the grid cannot support does not join two components,
  /// so the caller never asks for a route across it.
  std::vector<int> drivableComponents(Scalar clearance) const;

  /// Plan between two corridor nodes chosen by the caller.
  ///
  /// `plan()` and `planToNode()` snap the start and goal to the *nearest* node,
  /// which is right when the caller has nothing better to say. It is wrong when
  /// the caller does: the nearest node to a point is often an isolated pocket of
  /// high clearance with no corridors at all, and re-snapping a request that was
  /// deliberately aimed at a well-connected node redirects it to a dead end and
  /// reports `no_corridor_route` for a route the caller verified exists. When the
  /// caller has already reasoned about connectivity -- as a navigator that
  /// inspected the graph's components has -- it should be able to say so.
  Plan planViaNodes(const Vec2& start, int start_node, const Vec2& goal, int goal_node);

  /// Plan to a specific corridor node.
  Plan planToNode(const Vec2& start, int goal_node);

  /// True when a point is inside the map's free space with at least
  /// `clearance` metres to the nearest obstacle.
  bool isDrivable(const Vec2& p, Scalar clearance) const;

  /// Distance to the nearest obstacle at a world point, metres. Backed by a
  /// clearance field computed once in the constructor: the planner queries it
  /// thousands of times per plan (every sample of every candidate shortcut,
  /// every path vertex), and a ring search per query made a 253-room maze plan
  /// take over a second.
  Scalar clearanceAt(const Vec2& p) const;

  /// Curvature of a polyline at each interior point, 1/m.
  static std::vector<Scalar> curvature(const std::vector<Vec2>& path);

  const PlanOptions& options() const { return opt_; }

 private:
  /// Shared implementation of the two public plan() entry points.
  /// Direct grid-A* plan, used when the corridor graph cannot supply a route or
  /// its route does not survive collision verification. Returns a plan whose
  /// status is kOk, kBlocked, or kNoCorridorRoute.
  Plan planDirectGrid(const Vec2& start, const Vec2& goal,
                      const std::chrono::steady_clock::time_point& t0);

  Plan planToNodeImpl(const Vec2& start, const Vec2& goal, int start_node, int goal_node,
                      const std::chrono::steady_clock::time_point& t0);

  /// Expand a corridor route into a dense polyline through room centres and
  /// doorways.
  std::vector<Vec2> expandRoute(const std::vector<int>& route) const;
  /// Iteratively relax waypoints until the curvature bound is met.
  void smoothPath(std::vector<Vec2>* path) const;

  /// Re-centre the path by ascending the clearance field. Smoothing pulls
  /// points toward the straight line between waypoints, which in a doorway
  /// means off the centreline; a path that enters a 1.2 m doorway 0.18 m off
  /// centre has 0.08 m of margin against a 0.36 m half-diagonal and is not
  /// drivable no matter how well the controller tracks it. A few steps of
  /// gradient ascent on clearance is the discrete equivalent of snapping to the
  /// medial axis, and it is what turns "geometrically valid" into "drivable".
  void refineClearance(std::vector<Vec2>* path) const;
  /// Whether replacing p[i] with `cand` keeps the turn at i within
  /// `max_turn_per_point`.
  bool drivableTurn(const std::vector<Vec2>& p, Size i, const Vec2& cand) const;

  const topology::CorridorGraph& graph_;
  const sim::GridMap& grid_;
  PlanOptions opt_;
  /// Row-major over the grid: distance to the nearest occupied cell.
  std::vector<Scalar> clearance_;
};

/// Resample a polyline at a fixed arc-length spacing.
std::vector<Vec2> resample(const std::vector<Vec2>& path, Scalar spacing);

/// Arc-length parameterisation: cumulative distance at each vertex.
std::vector<Scalar> arcLength(const std::vector<Vec2>& path);

/// Point at arc length `s` along the path, by linear interpolation.
Vec2 pointAt(const std::vector<Vec2>& path, const std::vector<Scalar>& s_axis, Scalar s);

}  // namespace planning
}  // namespace luna
