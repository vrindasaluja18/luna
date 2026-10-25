#include "slamcore/planning/planner.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include "slamcore/topology/morse.hpp"
#include <limits>

namespace sc {
namespace planning {
namespace {

/// Distance to the nearest obstacle, from a precomputed clearance field. The
/// planner asks this question thousands of times per plan, so it is a lookup
/// rather than a search.
Scalar clearanceAtCell(const sim::GridMap& g, const std::vector<Scalar>& cl, const Vec2& p) {
  const Eigen::Vector2i c = g.toCell(p);
  if (!g.inBounds(c.x(), c.y())) return -1.0;
  return cl[static_cast<Size>(c.y()) * g.w + c.x()];
}

bool drivable(const sim::GridMap& g, const std::vector<Scalar>& cl, const Vec2& p,
              Scalar clearance) {
  return clearanceAtCell(g, cl, p) >= clearance;
}

/// Is the straight segment a->b free of obstacles with at least `clearance`
/// metres on both sides? Sampled at half the grid resolution, so a segment that
/// threads a gap narrower than one cell is still rejected.
bool segmentDrivable(const sim::GridMap& g, const std::vector<Scalar>& cl, const Vec2& a,
                     const Vec2& b, Scalar clearance) {
  const Scalar len = (b - a).norm();
  const int n = std::max(2, static_cast<int>(std::ceil(len / (0.5 * g.res))));
  for (int i = 0; i <= n; ++i) {
    const Scalar u = static_cast<Scalar>(i) / n;
    if (!drivable(g, cl, a + u * (b - a), clearance)) return false;
  }
  return true;
}

/// Remove the staircase from an 8-connected grid path.
///
/// A* on a grid returns a sequence of 5 cm axis-aligned steps. Resampling that
/// at a drivable spacing produces a polyline whose curvature is dominated by
/// the staircase, not by the route: the smoother then spends its whole budget
/// fighting 90-degree grid corners and the reported curvature stays an order of
/// magnitude above the vehicle's limit. Cutting corners greedily first --
/// connecting the furthest reachable vertex by a straight line -- removes the
/// staircase and leaves the smoother the real corners to round.
std::vector<Vec2> shortcutPath(const sim::GridMap& g, const std::vector<Scalar>& cl,
                               const std::vector<Vec2>& path, Scalar clearance) {
  if (path.size() < 3) return path;
  std::vector<Vec2> out;
  out.push_back(path.front());
  Size i = 0;
  while (i + 1 < path.size()) {
    // Scan as far ahead as the straight line stays drivable, then commit to
    // the furthest such vertex.
    Size best = i + 1;
    for (Size j = path.size() - 1; j > i + 1; --j) {
      if (segmentDrivable(g, cl, path[i], path[j], clearance)) {
        best = j;
        break;
      }
    }
    out.push_back(path[best]);
    i = best;
  }
  return out;
}

/// Shortest path on the 8-connected grid between two points, at the given
/// clearance. Used for the final metric connection between a corridor route
/// and the exact goal, which is where the corridor abstraction is too coarse.
std::vector<Vec2> gridAstar(const sim::GridMap& g, const std::vector<Scalar>& cl, const Vec2& from,
                            const Vec2& to, Scalar clearance, Size node_cap = 400000) {
  const Eigen::Vector2i s = g.toCell(from), t = g.toCell(to);
  if (!g.inBounds(s.x(), s.y()) || !g.inBounds(t.x(), t.y())) return {};
  if (!drivable(g, cl, from, clearance) || !drivable(g, cl, to, clearance)) return {};

  const Size total = static_cast<Size>(g.w) * g.h;
  const Index start = static_cast<Index>(s.y()) * g.w + s.x();
  const Index goal = static_cast<Index>(t.y()) * g.w + t.x();
  if (start == goal) return {from};

  std::vector<Scalar> gcost(total, std::numeric_limits<Scalar>::max());
  std::vector<Index> parent(total, -1);
  using QItem = std::pair<Scalar, Index>;
  std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> open;

  const auto heur = [&](Index i) {
    const int x = i % g.w, y = i / g.w;
    return std::hypot(Scalar(x - t.x()), Scalar(y - t.y())) * g.res;
  };
  gcost[static_cast<Size>(start)] = 0;
  open.push({heur(start), start});
  Size expanded = 0;

  static const int dx[8] = {1, -1, 0, 0, 1, 1, -1, -1};
  static const int dy[8] = {0, 0, 1, -1, 1, -1, 1, -1};

  while (!open.empty()) {
    const Index cur = open.top().second;
    open.pop();
    if (cur == goal) break;
    if (++expanded > node_cap) return {};
    const int cx = cur % g.w, cy = cur / g.w;
    for (int k = 0; k < 8; ++k) {
      const int nx = cx + dx[k], ny = cy + dy[k];
      if (!g.inBounds(nx, ny)) continue;
      // A cell is traversable only if the vehicle fits: the clearance field
      // already accounts for the requested margin, so this single test
      // replaces the obstacle check and the corner-cutting guard together.
      if (cl[static_cast<Size>(ny) * g.w + nx] < clearance) continue;
      const Index nb = static_cast<Index>(ny) * g.w + nx;
      const Scalar step = (dx[k] != 0 && dy[k] != 0) ? 1.41421356 * g.res : g.res;
      const Scalar cand = gcost[static_cast<Size>(cur)] + step;
      if (cand < gcost[static_cast<Size>(nb)]) {
        gcost[static_cast<Size>(nb)] = cand;
        parent[static_cast<Size>(nb)] = cur;
        open.push({cand + heur(nb), nb});
      }
    }
  }
  if (gcost[static_cast<Size>(goal)] == std::numeric_limits<Scalar>::max()) return {};

  std::vector<Vec2> out;
  for (Index i = goal; i >= 0; i = parent[static_cast<Size>(i)]) {
    out.push_back(g.cellCenter(i % g.w, i / g.w));
    if (i == start) break;
  }
  std::reverse(out.begin(), out.end());
  if (!out.empty()) {
    out.front() = from;
    out.back() = to;
  }
  return out;
}

}  // namespace

const char* toString(PlanStatus s) {
  switch (s) {
    case PlanStatus::kOk: return "ok";
    case PlanStatus::kNoCorridorRoute: return "no_corridor_route";
    case PlanStatus::kNoPath: return "no_path";
    case PlanStatus::kBlocked: return "blocked";
    case PlanStatus::kInvalid: return "invalid";
  }
  return "?";
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

std::vector<Scalar> arcLength(const std::vector<Vec2>& path) {
  std::vector<Scalar> s(path.size(), 0.0);
  for (Size i = 1; i < path.size(); ++i) s[i] = s[i - 1] + (path[i] - path[i - 1]).norm();
  return s;
}

Vec2 pointAt(const std::vector<Vec2>& path, const std::vector<Scalar>& s_axis, Scalar s) {
  if (path.empty()) return Vec2::Zero();
  if (s <= 0) return path.front();
  if (s >= s_axis.back()) return path.back();
  const auto it = std::upper_bound(s_axis.begin(), s_axis.end(), s);
  const Size hi = static_cast<Size>(it - s_axis.begin());
  if (hi == 0) return path.front();
  const Size lo = hi - 1;
  const Scalar span = s_axis[hi] - s_axis[lo];
  const Scalar u = span > 1e-12 ? (s - s_axis[lo]) / span : 0.0;
  return path[lo] + u * (path[hi] - path[lo]);
}

std::vector<Vec2> resample(const std::vector<Vec2>& path, Scalar spacing) {
  if (path.size() < 2 || spacing <= 0) return path;
  const std::vector<Scalar> s = arcLength(path);
  const Scalar total = s.back();
  if (total <= 0) return path;
  const Size n = static_cast<Size>(std::floor(total / spacing)) + 1;
  std::vector<Vec2> out;
  out.reserve(n + 1);
  for (Size i = 0; i < n; ++i) out.push_back(pointAt(path, s, i * spacing));
  out.push_back(path.back());
  return out;
}

std::vector<Scalar> CorridorPlanner::curvature(const std::vector<Vec2>& path) {
  std::vector<Scalar> out(path.size(), 0.0);
  for (Size i = 1; i + 1 < path.size(); ++i) {
    const Vec2 a = path[i] - path[i - 1];
    const Vec2 b = path[i + 1] - path[i];
    const Scalar la = a.norm(), lb = b.norm();
    if (la < 1e-9 || lb < 1e-9) continue;
    // k = 2 sin(theta/2) / |chord|, with theta the turn angle. This is the
    // standard discrete curvature; using theta / |chord| instead overestimates
    // sharp corners, which makes the smoother over-relax the straight parts.
    const Scalar cross = a.x() * b.y() - a.y() * b.x();
    const Scalar sin_half = std::abs(cross) / (la * lb);
    const Scalar chord = (path[i + 1] - path[i - 1]).norm();
    if (chord < 1e-9) continue;
    out[i] = 2.0 * std::min<Scalar>(1.0, sin_half) / chord;
  }
  return out;
}

// ---------------------------------------------------------------------------
// CorridorPlanner
// ---------------------------------------------------------------------------

CorridorPlanner::CorridorPlanner(const topology::CorridorGraph& graph, const sim::GridMap& grid,
                                 const PlanOptions& opt)
    : graph_(graph), grid_(grid), opt_(opt) {
  topology::clearanceTransform(grid_, &clearance_);
}

Scalar CorridorPlanner::clearanceAt(const Vec2& p) const {
  return clearanceAtCell(grid_, clearance_, p);
}

bool CorridorPlanner::isDrivable(const Vec2& p, Scalar clearance) const {
  return clearanceAt(p) >= clearance;
}

std::vector<Vec2> CorridorPlanner::expandRoute(const std::vector<int>& route) const {
  // A corridor route is [room, doorway, room, doorway, ...]. Room nodes are the
  // waypoints; doorway nodes are passed through as intermediate points, which
  // is what keeps the path centred in the doorway instead of cutting the corner
  // of the wall beside it.
  std::vector<Vec2> pts;
  for (int n : route) {
    if (n < 0 || static_cast<Size>(n) >= graph_.nodeCount()) continue;
    pts.push_back(graph_.nodes[static_cast<Size>(n)].p);
  }
  return pts;
}

void CorridorPlanner::smoothPath(std::vector<Vec2>* path) const {
  if (path->size() < 3) return;
  std::vector<Vec2>& p = *path;
  // Endpoints are fixed: they are the start and goal, and moving them would
  // make the plan not reach its target.
  for (int iter = 0; iter < opt_.smooth_iterations; ++iter) {
    const std::vector<Scalar> k = curvature(p);
    Scalar worst = 0;
    for (Scalar v : k) worst = std::max(worst, v);
    if (worst <= opt_.max_curvature) break;
    // Relax only the points that violate the bound. Relaxing everything
    // straightens the whole path, including the parts that were already
    // trackable, and pulls it away from the corridor centreline for no reason.
    for (Size i = 1; i + 1 < p.size(); ++i) {
      if (k[i] <= opt_.max_curvature) continue;
      const Vec2 mid = 0.5 * (p[i - 1] + p[i + 1]);
      const Vec2 moved = p[i] + 0.5 * (mid - p[i]);
      // Never move into an obstacle, and never below the minimum margin. A
      // smoother that ignores the map produces a beautiful spline through a
      // wall.
      if (!drivable(grid_, clearance_, moved, opt_.min_clearance)) continue;
      p[i] = moved;
    }
  }
}

void CorridorPlanner::refineClearance(std::vector<Vec2>* path) const {
  if (path->size() < 3) return;
  std::vector<Vec2>& p = *path;
  const Scalar h = std::max<Scalar>(0.5 * grid_.res, 0.05);
  const Scalar step = 0.5 * grid_.res;
  for (int iter = 0; iter < opt_.clearance_iterations; ++iter) {
    for (Size i = 1; i + 1 < p.size(); ++i) {
      const Scalar c0 = clearanceAt(p[i]);
      if (c0 <= 0) continue;
      const Vec2 g((clearanceAt(p[i] + Vec2(h, 0)) - clearanceAt(p[i] - Vec2(h, 0))) / (2 * h),
                   (clearanceAt(p[i] + Vec2(0, h)) - clearanceAt(p[i] - Vec2(0, h))) / (2 * h));
      const Scalar gn = g.norm();
      if (gn < 1e-6) continue;
      const Vec2 cand = p[i] + step * (g / gn);
      if (clearanceAt(cand) <= c0) continue;
      // Accept only if the move does not introduce a kink the vehicle cannot
      // follow: the whole point is a drivable path, and centring a path that
      // then needs a 2 m^-1 turn is no improvement.
      //
      // The kink has to be measured at the *candidate*, not at the point being
      // replaced. Testing the old point asks whether the path was drivable
      // before the move -- which it was, or we would not be here -- and says
      // nothing about whether it still is afterwards. The gradient step then
      // walks freely into a fold, and because every later iteration only ever
      // accepts a move that improves clearance, the fold is never undone. The
      // symptom is a path with a hairpin in the middle of a wide corridor,
      // which the controller cannot follow and the watchdog eventually stops
      // for. Rejecting moves that worsen the turn is the whole purpose of the
      // check, so it belongs on the state being accepted.
      if (!drivableTurn(p, i, cand)) continue;
      p[i] = cand;
    }
  }
}

bool CorridorPlanner::drivableTurn(const std::vector<Vec2>& p, Size i, const Vec2& cand) const {
  const Vec2 prev = cand - p[i - 1];
  const Vec2 next = (i + 1 < p.size() ? p[i + 1] : cand) - cand;
  const Scalar lp = prev.norm(), ln = next.norm();
  if (lp < 1e-9 || ln < 1e-9) return false;
  const Scalar cosang = prev.dot(next) / (lp * ln);
  return cosang >= std::cos(opt_.max_turn_per_point);
}

bool CorridorPlanner::segmentClear(const Vec2& a, const Vec2& b, Scalar clearance) const {
  return segmentDrivable(grid_, clearance_, a, b, clearance);
}

std::vector<int> CorridorPlanner::drivableComponents(Scalar clearance) const {
  const Size n = graph_.nodeCount();
  std::vector<int> label(n, -1);
  std::vector<Size> stack;
  std::vector<std::vector<Size>> adj(n);
  for (const auto& e : graph_.edges) {
    const Size a = static_cast<Size>(e.a), b = static_cast<Size>(e.b);
    if (a >= n || b >= n) continue;
    if (!segmentClear(graph_.nodes[a].p, graph_.nodes[b].p, clearance)) continue;
    adj[a].push_back(b);
    adj[b].push_back(a);
  }
  int next = 0;
  for (Size k = 0; k < n; ++k) {
    if (label[k] >= 0) continue;
    const int id = next++;
    label[k] = id;
    stack.assign(1, k);
    while (!stack.empty()) {
      const Size cur = stack.back();
      stack.pop_back();
      for (Size nb : adj[cur]) {
        if (label[nb] >= 0) continue;
        label[nb] = id;
        stack.push_back(nb);
      }
    }
  }
  return label;
}

Plan CorridorPlanner::plan(const Vec2& start, const Vec2& goal, Scalar snap_radius) {
  const auto t0 = std::chrono::steady_clock::now();
  Plan plan;
  if (graph_.nodeCount() == 0) {
    plan.status = PlanStatus::kNoCorridorRoute;
    return plan;
  }
  // Snap the start and goal to the corridor graph. The goal is allowed to be
  // outside the snap radius: it is then attached to the nearest node by a
  // direct grid path, which is the honest way to handle a goal that sits in the
  // middle of a room rather than at its centre.
  const int start_node = graph_.nearestNode(start, std::max(snap_radius, 1.0) * 4.0);
  const int goal_node = graph_.nearestNode(goal, std::max(snap_radius, 1.0) * 4.0);
  if (start_node < 0 || goal_node < 0) {
    plan.status = PlanStatus::kNoCorridorRoute;
    return plan;
  }
  return planToNodeImpl(start, goal, start_node, goal_node, t0);
}

Plan CorridorPlanner::planViaNodes(const Vec2& start, int start_node, const Vec2& goal,
                                   int goal_node) {
  const auto t0 = std::chrono::steady_clock::now();
  const Size n = graph_.nodeCount();
  if (start_node < 0 || goal_node < 0 || static_cast<Size>(start_node) >= n ||
      static_cast<Size>(goal_node) >= n) {
    Plan plan;
    plan.status = PlanStatus::kInvalid;
    return plan;
  }
  Plan plan = planToNodeImpl(start, goal, start_node, goal_node, t0);
  // Same reasoning as in `plan`: when the caller has fixed the corridor nodes
  // but the corridor route does not exist or does not survive the map, the grid
  // still knows whether a route exists. A caller that named two nodes has told
  // us *where* to go, not that the Morse graph must supply the way.
  if (!plan.ok()) {
    Plan direct = planDirectGrid(start, goal, t0);
    if (direct.ok()) return direct;
  }
  return plan;
}

Plan CorridorPlanner::planToNode(const Vec2& start, int goal_node) {
  const auto t0 = std::chrono::steady_clock::now();
  Plan plan;
  if (goal_node < 0 || static_cast<Size>(goal_node) >= graph_.nodeCount()) {
    plan.status = PlanStatus::kInvalid;
    return plan;
  }
  const int start_node = graph_.nearestNode(start, 8.0);
  if (start_node < 0) {
    plan.status = PlanStatus::kNoCorridorRoute;
    return plan;
  }
  return planToNodeImpl(start, graph_.nodes[static_cast<Size>(goal_node)].p, start_node, goal_node,
                        t0);
}

Plan CorridorPlanner::planDirectGrid(const Vec2& start, const Vec2& goal,
                                    const std::chrono::steady_clock::time_point& t0) {
  Plan plan;
  plan.stats.corridor_nodes = static_cast<int>(graph_.nodeCount());
  plan.stats.corridor_edges = static_cast<int>(graph_.edgeCount());
  std::vector<Vec2> direct = gridAstar(grid_, clearance_, start, goal, opt_.preferred_clearance);
  if (direct.empty()) direct = gridAstar(grid_, clearance_, start, goal, opt_.min_clearance);
  if (direct.empty()) {
    plan.status = PlanStatus::kNoCorridorRoute;
    return plan;
  }
  direct = shortcutPath(grid_, clearance_, direct, opt_.min_clearance);
  direct = resample(direct, opt_.waypoint_spacing);
  smoothPath(&direct);
  refineClearance(&direct);
  if (opt_.verify_collision_free) {
    for (const Vec2& p : direct) {
      if (clearanceAt(p) <= opt_.acceptanceClearance()) {
        plan.status = PlanStatus::kBlocked;
        plan.path = direct;
        plan.stats.plan_ms =
            std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0)
                .count();
        return plan;
      }
    }
  }
  Scalar min_clear = std::numeric_limits<Scalar>::max();
  plan.clearance.resize(direct.size());
  for (Size i = 0; i < direct.size(); ++i) {
    plan.clearance[i] = clearanceAt(direct[i]);
    min_clear = std::min(min_clear, plan.clearance[i]);
  }
  plan.stats.min_clearance_along = min_clear == std::numeric_limits<Scalar>::max() ? 0 : min_clear;
  Scalar worst_k = 0;
  for (Scalar v : curvature(direct)) worst_k = std::max(worst_k, v);
  plan.stats.max_curvature_along = worst_k;
  plan.stats.waypoints = static_cast<int>(direct.size());
  plan.stats.path_length = arcLength(direct).empty() ? 0 : arcLength(direct).back();
  plan.path = direct;
  plan.status = PlanStatus::kOk;
  plan.stats.plan_ms =
      std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return plan;
}

Plan CorridorPlanner::planToNodeImpl(const Vec2& start, const Vec2& goal, int start_node,
                                     int goal_node, const std::chrono::steady_clock::time_point& t0) {
  Plan plan;
  plan.stats.corridor_nodes = static_cast<int>(graph_.nodeCount());
  plan.stats.corridor_edges = static_cast<int>(graph_.edgeCount());

  std::vector<int> route;
  if (!graph_.route(start_node, goal_node, &route)) {
    // No corridor route. Fall back to searching the grid directly.
    //
    // The corridor graph is a *guide*, not the only admissible route. It is
    // built from the map's Morse structure, so it inherits the map's coverage:
    // where a region has been observed only as a thin trail -- which is what a
    // maze looks like from a single scripted pass -- the height function has no
    // well-separated maxima to connect, and the graph fragments into isolated
    // nodes even though the free space underneath is one connected component.
    // Refusing to plan then is refusing on the strength of a structure that is
    // less reliable than the grid it was derived from.
    //
    // The grid search has no such dependency: it walks the occupancy map and
    // finds a route whenever one exists. It is more expensive and produces a
    // less elegant path, which is exactly why it is the fallback and not the
    // default -- but a coarse path through a connected space beats a confident
    // "blocked" on one.
    return planDirectGrid(start, goal, t0);
  }
  plan.corridor_route = route;

  // Bottleneck of the tightest doorway on the route: the corridor edges carry
  // their saddle clearance, so this needs no map query.
  Scalar bottleneck = std::numeric_limits<Scalar>::max();
  for (Size i = 0; i + 1 < route.size(); ++i) {
    for (const topology::CorridorEdge& e : graph_.edges) {
      const int a = e.a, b = e.b;
      if ((a == route[i] && b == route[i + 1]) || (b == route[i] && a == route[i + 1])) {
        bottleneck = std::min(bottleneck, e.bottleneck);
        break;
      }
    }
  }
  plan.stats.route_bottleneck = bottleneck == std::numeric_limits<Scalar>::max() ? 0 : bottleneck;

  // Expand to the corridor waypoints: room centres and the doorways between
  // them.
  std::vector<Vec2> wps = expandRoute(route);
  if (wps.empty()) {
    plan.status = PlanStatus::kNoPath;
    return plan;
  }
  wps.front() = start;
  wps.back() = goal;
  // Snap the interior waypoints into drivable space. A corridor node is the
  // centre of a region as the *topology extractor* sees it, and on a
  // reconstructed map that centre can be a cell the occupancy threshold called
  // occupied -- the extractor's height function is smoothed, so its maximum sits
  // where the smoothed clearance peaks, which is not exactly where the raw map
  // is free. Handing such a point to the grid search makes it fail for that leg,
  // and the failure is silent: the caller falls back to a straight segment, so
  // the plan looks complete while containing a leg that was never checked
  // against the map.
  for (Size i = 1; i + 1 < wps.size(); ++i) {
    if (isDrivable(wps[i], opt_.min_clearance)) continue;
    // Search outwards in rings for the nearest cell that is drivable at the
    // *minimum* margin. The minimum, not the preferred one: a doorway's centre
    // is drivable at the minimum and often not at the preferred, and refusing to
    // snap to it would drop the doorway from the route entirely.
    const Eigen::Vector2i c = grid_.toCell(wps[i]);
    const int rad = static_cast<int>(std::ceil(opt_.preferred_clearance / grid_.res));
    bool found = false;
    for (int r = 1; r <= rad && !found; ++r) {
      for (int dy = -r; dy <= r && !found; ++dy) {
        for (int dx = -r; dx <= r && !found; ++dx) {
          if (std::max(std::abs(dx), std::abs(dy)) != r) continue;
          const int nx = c.x() + dx, ny = c.y() + dy;
          if (!grid_.inBounds(nx, ny)) continue;
          const Vec2 cand = grid_.cellCenter(nx, ny);
          if (!isDrivable(cand, opt_.min_clearance)) continue;
          wps[i] = cand;
          found = true;
        }
      }
    }
  }

  // Connect consecutive corridor waypoints with a grid path rather than a
  // straight line.
  //
  // A straight line between two room centres is the obvious thing to draw and
  // it is wrong: room centres and doorways are not collinear, so the line
  // clips the wall beside the door. Measured on the reference world it left the
  // path 5 cm from an obstacle -- inside the vehicle -- on a world whose rooms
  // are metres wide. The corridor graph supplies the *order* of rooms to visit;
  // the grid supplies the geometry between them. Connecting them this way is
  // what makes the two-level decomposition pay off instead of fighting itself.
  std::vector<Vec2> path;
  path.push_back(start);
  for (Size i = 0; i + 1 < wps.size(); ++i) {
    const Vec2& a = path.back();
    const Vec2& b = wps[i + 1];
    if ((a - b).norm() < 1e-6) continue;
    // Preferred clearance first, then the minimum. A wide corridor is taken
    // down the middle; a doorway only as wide as the vehicle is taken
    // centrally at the smaller margin.
    std::vector<Vec2> seg = gridAstar(grid_, clearance_, a, b, opt_.preferred_clearance);
    if (seg.empty()) seg = gridAstar(grid_, clearance_, a, b, opt_.min_clearance);
    if (seg.empty()) {
      // No connected path at either margin: fall back to the straight segment
      // so the caller gets a plan it can reject explicitly rather than an empty
      // result it cannot diagnose.
      path.push_back(b);
      continue;
    }
    path.insert(path.end(), seg.begin() + 1, seg.end());
  }
  // Cut the staircase out of the grid paths before smoothing. The clearance
  // used is the minimum rather than the preferred one: the path has already
  // been built at the preferred margin where that was possible, and
  // shortcutting at the preferred margin would refuse to straighten anything
  // inside a doorway.
  path = shortcutPath(grid_, clearance_, path, opt_.min_clearance);
  path = resample(path, opt_.waypoint_spacing);
  plan.waypoints = wps;

  smoothPath(&path);
  refineClearance(&path);

  // Measure the result against the map.
  Scalar min_clear = std::numeric_limits<Scalar>::max();
  plan.clearance.resize(path.size());
  for (Size i = 0; i < path.size(); ++i) {
    plan.clearance[i] = clearanceAt(path[i]);
    min_clear = std::min(min_clear, plan.clearance[i]);
  }
  plan.stats.min_clearance_along =
      min_clear == std::numeric_limits<Scalar>::max() ? 0 : min_clear;

  const std::vector<Scalar> k = curvature(path);
  Scalar worst_k = 0;
  for (Scalar v : k) worst_k = std::max(worst_k, v);
  plan.stats.max_curvature_along = worst_k;
  plan.stats.waypoints = static_cast<int>(path.size());
  plan.stats.path_length = arcLength(path).empty() ? 0 : arcLength(path).back();
  {
    Scalar corridor_len = 0;
    for (Size i = 0; i + 1 < route.size(); ++i)
      corridor_len += (graph_.nodes[static_cast<Size>(route[i])].p -
                       graph_.nodes[static_cast<Size>(route[i + 1])].p)
                          .norm();
    plan.stats.corridor_route_length = corridor_len;
  }

  if (opt_.verify_collision_free) {
    for (const Vec2& p : path) {
      if (clearanceAt(p) <= opt_.acceptanceClearance()) {
        // The corridor route does not survive the map. Before reporting a
        // blockage, try the grid directly: a leg that the corridor graph
        // connected through a waypoint the map does not actually clear can
        // still have a valid route around it, and the grid search is the thing
        // that knows. Reporting "blocked" from a structure that is only a guide
        // would refuse a route that plainly exists.
        Plan direct = planDirectGrid(start, goal, t0);
        if (direct.ok()) return direct;
        plan.status = PlanStatus::kBlocked;
        plan.path = path;
        plan.stats.plan_ms =
            std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0)
                .count();
        return plan;
      }
    }
  }

  plan.path = path;
  plan.status = PlanStatus::kOk;
  plan.stats.plan_ms =
      std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return plan;
}

}  // namespace planning
}  // namespace sc
