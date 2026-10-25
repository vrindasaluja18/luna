#include "slamcore/sim/world.hpp"

#include <cmath>
#include <limits>

#include <algorithm>
#include <cstring>
#include <deque>
#include <sstream>

namespace sc {
namespace sim {

// ---------------------------------------------------------------------------
// GridMap
// ---------------------------------------------------------------------------

void GridMap::resize(int width, int height, Scalar resolution, const Vec2& org) {
  w = width;
  h = height;
  res = resolution;
  origin = org;
  cells.assign(static_cast<Size>(w) * static_cast<Size>(h), static_cast<std::uint8_t>(Cell::kUnknown));
}

Cell GridMap::at(int x, int y) const {
  if (!inBounds(x, y)) return Cell::kUnknown;
  return static_cast<Cell>(cells[static_cast<Size>(y) * w + x]);
}

void GridMap::set(int x, int y, Cell c) {
  if (!inBounds(x, y)) return;
  cells[static_cast<Size>(y) * w + x] = static_cast<std::uint8_t>(c);
}

Vec2 GridMap::cellCenter(int x, int y) const {
  return origin + Vec2((x + 0.5) * res, (y + 0.5) * res);
}

Eigen::Vector2i GridMap::toCell(const Vec2& p) const {
  return Eigen::Vector2i(static_cast<int>(std::floor((p.x() - origin.x()) / res)),
                         static_cast<int>(std::floor((p.y() - origin.y()) / res)));
}

Cell GridMap::atWorld(const Vec2& p) const {
  const Eigen::Vector2i c = toCell(p);
  return at(c.x(), c.y());
}

Size GridMap::countFree() const {
  Size n = 0;
  for (auto v : cells)
    if (static_cast<Cell>(v) == Cell::kFree) ++n;
  return n;
}

void GridMap::fillRect(const Vec2& lo, const Vec2& hi, Cell c) {
  const Eigen::Vector2i a = toCell(lo);
  const Eigen::Vector2i b = toCell(hi);
  for (int y = a.y(); y <= b.y(); ++y)
    for (int x = a.x(); x <= b.x(); ++x) set(x, y, c);
}

void GridMap::fillDisc(const Vec2& center, Scalar radius, Cell c) {
  const int r = static_cast<int>(std::ceil(radius / res));
  const Eigen::Vector2i cc = toCell(center);
  for (int dy = -r; dy <= r; ++dy)
    for (int dx = -r; dx <= r; ++dx) {
      if (dx * dx + dy * dy > r * r) continue;
      set(cc.x() + dx, cc.y() + dy, c);
    }
}

void GridMap::dilate(Scalar radius) {
  if (radius <= 0) return;
  std::vector<std::uint8_t> out = cells;
  const int r = static_cast<int>(std::ceil(radius / res));
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (static_cast<Cell>(cells[static_cast<Size>(y) * w + x]) != Cell::kOccupied) continue;
      for (int dy = -r; dy <= r; ++dy) {
        for (int dx = -r; dx <= r; ++dx) {
          if (dx * dx + dy * dy > r * r) continue;
          const int nx = x + dx, ny = y + dy;
          if (inBounds(nx, ny)) out[static_cast<Size>(ny) * w + nx] = static_cast<std::uint8_t>(Cell::kOccupied);
        }
      }
    }
  }
  cells.swap(out);
}

// ---------------------------------------------------------------------------
// WorldType helpers
// ---------------------------------------------------------------------------

const char* toString(WorldType t) {
  switch (t) {
    case WorldType::kRooms: return "rooms";
    case WorldType::kCorridorLoop: return "corridor_loop";
    case WorldType::kSaddleHall: return "saddle_hall";
    case WorldType::kMaze: return "maze";
    case WorldType::kCluttered: return "cluttered";
    case WorldType::kTextureless: return "textureless";
  }
  return "?";
}

bool parseWorldType(const std::string& s, WorldType* out) {
  static const std::pair<const char*, WorldType> kTable[] = {
      {"rooms", WorldType::kRooms},
      {"corridor_loop", WorldType::kCorridorLoop},
      {"corridor", WorldType::kCorridorLoop},
      {"saddle_hall", WorldType::kSaddleHall},
      {"saddle", WorldType::kSaddleHall},
      {"maze", WorldType::kMaze},
      {"cluttered", WorldType::kCluttered},
      {"textureless", WorldType::kTextureless},
  };
  for (const auto& e : kTable) {
    if (s == e.first) {
      *out = e.second;
      return true;
    }
  }
  return false;
}

// ---------------------------------------------------------------------------
// generation helpers
// ---------------------------------------------------------------------------

namespace {

/// Punch a doorway through a wall segment between two world points.
///
/// The carved rectangle has to be wider than the wall is thick *and* centred
/// exactly on the wall, or it merely scrapes the surface and leaves the rooms
/// sealed. That failure is invisible in the ASCII map and only shows up as an
/// unreachable world, so the carve width is derived from the caller's wall
/// thickness rather than being a fixed pixel margin.
void carveDoorway(GridMap* g, const Vec2& center, Scalar width, Scalar wall_thickness,
                  bool horizontal_wall) {
  const Scalar half = std::max(width, 2.0 * g->res) * 0.5;
  const Scalar through = wall_thickness * 0.5 + 2.0 * g->res;
  if (horizontal_wall) {
    g->fillRect(Vec2(center.x() - half, center.y() - through),
                Vec2(center.x() + half, center.y() + through), Cell::kFree);
  } else {
    g->fillRect(Vec2(center.x() - through, center.y() - half),
                Vec2(center.x() + through, center.y() + half), Cell::kFree);
  }
}

void borderWalls(GridMap* g, Scalar thickness) {
  const Vec2 lo = g->minCorner();
  const Vec2 hi = g->maxCorner();
  g->fillRect(lo, Vec2(hi.x(), lo.y() + thickness), Cell::kOccupied);
  g->fillRect(Vec2(lo.x(), hi.y() - thickness), hi, Cell::kOccupied);
  g->fillRect(lo, Vec2(lo.x() + thickness, hi.y()), Cell::kOccupied);
  g->fillRect(Vec2(hi.x() - thickness, lo.y()), hi, Cell::kOccupied);
}

void genRooms(GridMap* g, const WorldSpec& spec, Rng* rng) {
  const Scalar L = spec.size_m;
  borderWalls(g, 0.4);
  // Grid of rooms: choose a factorisation of `rooms` into (cols, rows).
  int cols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(spec.rooms))));
  int rows = static_cast<int>(std::ceil(static_cast<double>(spec.rooms) / cols));
  const Scalar cellW = L / cols, cellH = L / rows;
  const Scalar wall = 0.3;

  // Interior walls with doorways. Each shared edge gets a doorway; the
  // resulting free-space graph is a grid, so it contains loops (good for
  // testing loop closure) and bottlenecks (good for Morse saddles).
  for (int cx = 1; cx < cols; ++cx) {
    const Scalar x = cx * cellW;
    g->fillRect(Vec2(x - wall / 2, 0.0), Vec2(x + wall / 2, L), Cell::kOccupied);
  }
  for (int cy = 1; cy < rows; ++cy) {
    const Scalar y = cy * cellH;
    g->fillRect(Vec2(0.0, y - wall / 2), Vec2(L, y + wall / 2), Cell::kOccupied);
  }
  // The doorway jitter is bounded so a door never drifts into a corner, where
  // the carve would eat the perpendicular wall and join three rooms at once.
  const Scalar jitter_y = std::max(0.0, 0.5 * (cellH - spec.door_width) - wall);
  const Scalar jitter_x = std::max(0.0, 0.5 * (cellW - spec.door_width) - wall);
  for (int cx = 1; cx < cols; ++cx) {
    const Scalar x = cx * cellW;
    for (int cy = 0; cy < rows; ++cy) {
      const Scalar y = (cy + 0.5) * cellH + rng->uniform(-jitter_y, jitter_y);
      carveDoorway(g, Vec2(x, y), spec.door_width, wall, false);
    }
  }
  for (int cy = 1; cy < rows; ++cy) {
    const Scalar y = cy * cellH;
    for (int cx = 0; cx < cols; ++cx) {
      const Scalar x = (cx + 0.5) * cellW + rng->uniform(-jitter_x, jitter_x);
      carveDoorway(g, Vec2(x, y), spec.door_width, wall, true);
    }
  }
}

void genCorridorLoop(GridMap* g, const WorldSpec& spec, Rng* rng) {
  const Scalar L = spec.size_m;
  borderWalls(g, 0.4);
  const Scalar inset = 0.25 * L;
  const Scalar w = 2.0;  // corridor half-width
  // Ring corridor: outer square minus inner block.
  const Scalar lo = inset, hi = L - inset;
  // Carve the ring by painting the inner block as wall.
  g->fillRect(Vec2(lo + w, lo + w), Vec2(hi - w, hi - w), Cell::kOccupied);
  // Two cross-branches to create additional cycles.
  const Scalar ymid = 0.5 * L;
  g->fillRect(Vec2(lo + w, ymid - w), Vec2(hi - w, ymid + w), Cell::kFree);
  const Scalar xmid = 0.5 * L;
  g->fillRect(Vec2(xmid - w, lo + w), Vec2(xmid + w, hi - w), Cell::kFree);
  (void)rng;
}

void genSaddleHall(GridMap* g, const WorldSpec& spec, Rng* rng) {
  const Scalar L = spec.size_m;
  borderWalls(g, 0.4);
  // Two large halls separated by a wall with one narrow gap: the gap is a
  // genuine Morse saddle of the clearance function.
  const Scalar wallX = 0.5 * L;
  const Scalar wall = 0.4;
  g->fillRect(Vec2(wallX - wall / 2, 0.0), Vec2(wallX + wall / 2, L), Cell::kOccupied);
  const Scalar gapY = 0.5 * L + rng->uniform(-0.05 * L, 0.05 * L);
  carveDoorway(g, Vec2(wallX, gapY), spec.door_width, wall, false);
  // A couple of pillars inside each hall create secondary critical points.
  for (int i = 0; i < 4; ++i) {
    const Scalar x = rng->uniform(2.0, wallX - 2.0);
    const Scalar y = rng->uniform(2.0, L - 2.0);
    g->fillDisc(Vec2(x, y), rng->uniform(0.6, 1.2), Cell::kOccupied);
    const Scalar x2 = rng->uniform(wallX + 2.0, L - 2.0);
    const Scalar y2 = rng->uniform(2.0, L - 2.0);
    g->fillDisc(Vec2(x2, y2), rng->uniform(0.6, 1.2), Cell::kOccupied);
  }
}

void genMaze(GridMap* g, const WorldSpec& spec, Rng* rng) {
  const Scalar L = spec.size_m;
  const int cells = std::max(4, spec.rooms * 2);
  const Scalar cs = L / cells;
  const Scalar wall = std::max(0.15, cs * 0.12);
  borderWalls(g, 0.4);

  // Start from every cell fully walled off, then knock down walls along a
  // spanning tree. Independent random knockouts would leave the maze in several
  // disconnected pieces, and a route built from sampled free points would then
  // be a single point: the vehicle drives straight into a wall and the run
  // reports collisions that are really a generator defect.
  for (int y = 0; y < cells; ++y) {
    for (int x = 0; x < cells; ++x) {
      const Scalar x0 = x * cs, y0 = y * cs;
      if (x + 1 < cells) g->fillRect(Vec2(x0 + cs - wall, y0), Vec2(x0 + cs, y0 + cs), Cell::kOccupied);
      if (y + 1 < cells) g->fillRect(Vec2(x0, y0 + cs - wall), Vec2(x0 + cs, y0 + cs), Cell::kOccupied);
    }
  }

  // Randomised depth-first search over the cell grid: each visited cell knocks
  // down the wall to one unvisited neighbour. Every cell is reached exactly
  // once, so the free space is a single connected region by construction.
  const int n = cells * cells;
  std::vector<std::uint8_t> visited(static_cast<Size>(n), 0);
  std::vector<std::pair<int, int>> stack;
  const auto cell_index = [cells](int x, int y) { return y * cells + x; };
  stack.push_back({0, 0});
  visited[0] = 1;
  while (!stack.empty()) {
    const auto [x, y] = stack.back();
    // Collect unvisited neighbours in a deterministic but seed-dependent order.
    int cand[4][2];
    int nc = 0;
    const int dx[4] = {1, -1, 0, 0};
    const int dy[4] = {0, 0, 1, -1};
    for (int k = 0; k < 4; ++k) {
      const int nx = x + dx[k], ny = y + dy[k];
      if (nx < 0 || ny < 0 || nx >= cells || ny >= cells) continue;
      if (visited[cell_index(nx, ny)]) continue;
      cand[nc][0] = nx;
      cand[nc][1] = ny;
      ++nc;
    }
    if (nc == 0) {
      stack.pop_back();
      continue;
    }
    const int pick = std::min(nc - 1, static_cast<int>(rng->uniform() * nc));
    const int nx = cand[pick][0], ny = cand[pick][1];
    const Scalar x0 = x * cs, y0 = y * cs;
    const Scalar jitter = cs * 0.15;
    // The wall between two cells is drawn as [edge - wall, edge], so its
    // centre sits half a wall thickness inside the cell, not on the edge.
    // Carving at the edge would only remove the outer half of the wall and
    // leave the cells sealed.
    const Scalar yc = y0 + cs * 0.5 + rng->uniform(-jitter, jitter);
    const Scalar xc = x0 + cs * 0.5 + rng->uniform(-jitter, jitter);
    if (nx > x) {
      carveDoorway(g, Vec2(x0 + cs - wall * 0.5, yc), cs * 0.6, wall, false);
    } else if (nx < x) {
      carveDoorway(g, Vec2(x0 - wall * 0.5, yc), cs * 0.6, wall, false);
    } else if (ny > y) {
      carveDoorway(g, Vec2(xc, y0 + cs - wall * 0.5), cs * 0.6, wall, true);
    } else {
      carveDoorway(g, Vec2(xc, y0 - wall * 0.5), cs * 0.6, wall, true);
    }
    visited[cell_index(nx, ny)] = 1;
    stack.push_back({nx, ny});
  }

  // A spanning tree has no loops, so add a few extra openings. Loops are what
  // make loop closure testable, and a pure tree would make every revisit a new
  // place.
  const int extra = std::max(1, cells);
  for (int i = 0; i < extra; ++i) {
    const int x = std::min(cells - 2, static_cast<int>(rng->uniform() * (cells - 1)));
    const int y = std::min(cells - 2, static_cast<int>(rng->uniform() * (cells - 1)));
    const Scalar x0 = x * cs, y0 = y * cs;
    const Scalar jitter = cs * 0.15;
    if (rng->uniform() < 0.5) {
      carveDoorway(g, Vec2(x0 + cs - wall * 0.5, y0 + cs * 0.5 + rng->uniform(-jitter, jitter)),
                   cs * 0.6, wall, false);
    } else {
      carveDoorway(g, Vec2(x0 + cs * 0.5 + rng->uniform(-jitter, jitter), y0 + cs - wall * 0.5),
                   cs * 0.6, wall, true);
    }
  }
}

void genCluttered(GridMap* g, const WorldSpec& spec, Rng* rng) {
  borderWalls(g, 0.4);
  const Scalar L = spec.size_m;
  const int n = static_cast<int>(spec.clutter * 400.0);
  for (int i = 0; i < n; ++i) {
    const Vec2 c(rng->uniform(1.5, L - 1.5), rng->uniform(1.5, L - 1.5));
    const Scalar r = rng->uniform(0.3, 1.0);
    g->fillDisc(c, r, Cell::kOccupied);
  }
}

void genTextureless(GridMap* g, const WorldSpec& spec, Rng* rng) {
  const Scalar L = spec.size_m;
  borderWalls(g, 0.4);
  // A long S-shaped corridor: long straight sections with no visual texture,
  // which is where feature-based odometry is expected to degrade.
  const Scalar w = 2.0;
  const Scalar y1 = 0.3 * L, y2 = 0.7 * L;
  g->fillRect(Vec2(1.0, y1 - w), Vec2(0.6 * L, y1 + w), Cell::kFree);
  g->fillRect(Vec2(0.6 * L - w, y1), Vec2(0.6 * L + w, y2), Cell::kFree);
  g->fillRect(Vec2(0.4 * L, y2 - w), Vec2(L - 1.0, y2 + w), Cell::kFree);
  (void)rng;
}

/// Interior clearance via a two-pass exact Euclidean distance transform
/// (Felzenszwalb-Huttenlocher). O(n) per row/column pass.
void distanceTransform(const GridMap& g, std::vector<Scalar>* dist) {
  const int w = g.w, h = g.h;
  const Scalar big = 1e12;
  std::vector<Scalar> f(std::max(w, h));
  std::vector<Scalar> d(static_cast<Size>(w) * h, 0.0);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
      d[static_cast<Size>(y) * w + x] = (g.at(x, y) == Cell::kOccupied) ? 0.0 : big;

  const auto dt1d = [&](int n, Scalar* f_in, Scalar* out) {
    std::vector<int> v(n);
    std::vector<Scalar> z(n + 1);
    int k = 0;
    v[0] = 0;
    z[0] = -big;
    z[1] = big;
    for (int q = 1; q < n; ++q) {
      Scalar s = 0;
      while (true) {
        s = ((f_in[q] + q * q) - (f_in[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
        if (s > z[k] || k == 0) break;
        --k;
      }
      ++k;
      v[k] = q;
      z[k] = s;
      z[k + 1] = big;
    }
    k = 0;
    for (int q = 0; q < n; ++q) {
      while (z[k + 1] < q) ++k;
      const Scalar dq = q - v[k];
      out[q] = dq * dq + f_in[v[k]];
    }
  };

  for (int x = 0; x < w; ++x) {
    for (int y = 0; y < h; ++y) f[y] = d[static_cast<Size>(y) * w + x];
    dt1d(h, f.data(), f.data());
    for (int y = 0; y < h; ++y) d[static_cast<Size>(y) * w + x] = f[y];
  }
  for (int y = 0; y < h; ++y) {
    Scalar* row = d.data() + static_cast<Size>(y) * w;
    dt1d(w, row, row);
  }
  dist->resize(d.size());
  for (Size i = 0; i < d.size(); ++i) (*dist)[i] = std::sqrt(d[i]) * g.res;
}

}  // namespace

// ---------------------------------------------------------------------------
// World
// ---------------------------------------------------------------------------

World World::generate(const WorldSpec& spec) {
  World world;
  world.spec_ = spec;
  Rng rng(spec.seed);

  const int n = static_cast<int>(std::round(spec.size_m / spec.res));
  world.grid_.resize(n, n, spec.res, Vec2::Zero());
  world.grid_.fillRect(Vec2::Zero(), world.grid_.maxCorner(), Cell::kFree);

  switch (spec.type) {
    case WorldType::kRooms: genRooms(&world.grid_, spec, &rng); break;
    case WorldType::kCorridorLoop: genCorridorLoop(&world.grid_, spec, &rng); break;
    case WorldType::kSaddleHall: genSaddleHall(&world.grid_, spec, &rng); break;
    case WorldType::kMaze: genMaze(&world.grid_, spec, &rng); break;
    case WorldType::kCluttered: genCluttered(&world.grid_, spec, &rng); break;
    case WorldType::kTextureless: genTextureless(&world.grid_, spec, &rng); break;
  }

  // Statistics, including the clearance distribution that the topology module
  // consumes as its height function.
  std::vector<Scalar> clearance;
  distanceTransform(world.grid_, &clearance);
  world.clearance_ = clearance;
  std::vector<Scalar> free_clear;
  free_clear.reserve(clearance.size());
  for (Size i = 0; i < clearance.size(); ++i) {
    const Cell c = static_cast<Cell>(world.grid_.cells[i]);
    if (c == Cell::kFree) {
      ++world.stats_.freeCells;
      free_clear.push_back(clearance[i]);
    } else if (c == Cell::kOccupied) {
      ++world.stats_.occupiedCells;
    }
  }
  world.stats_.freeArea = world.stats_.freeCells * spec.res * spec.res;
  if (!free_clear.empty()) {
    std::sort(free_clear.begin(), free_clear.end());
    world.stats_.corridorWidthMedian = free_clear[free_clear.size() / 2] * 2.0;
    // Minimum width over the *reachable* free space, ignoring the 1-cell skin
    // against walls (which is always ~1 cell by construction).
    Size idx = static_cast<Size>(free_clear.size() * 0.02);
    world.stats_.corridorWidthMin = free_clear[idx] * 2.0;
  }

  // Landmarks: sample free cells with clearance, which keeps them on visible
  // surfaces rather than inside walls.
  const int nlm = static_cast<int>(spec.landmark_density * world.stats_.freeArea);
  std::vector<std::uint8_t> taken(world.grid_.cells.size(), 0);
  int placed = 0;
  const Scalar min_sep = std::max(2.0 * spec.res, 0.35);
  const int sep_cells = static_cast<int>(std::ceil(min_sep / spec.res));
  int attempts = 0;
  const int max_attempts = nlm * 40 + 1000;
  while (placed < nlm && attempts < max_attempts) {
    ++attempts;
    const int x = static_cast<int>(rng.uniform(0, world.grid_.w));
    const int y = static_cast<int>(rng.uniform(0, world.grid_.h));
    if (world.grid_.at(x, y) != Cell::kFree) continue;
    if (clearance[static_cast<Size>(y) * world.grid_.w + x] < 0.15) continue;
    bool clash = false;
    for (int dy = -sep_cells; dy <= sep_cells && !clash; ++dy)
      for (int dx = -sep_cells; dx <= sep_cells && !clash; ++dx) {
        const int nx = x + dx, ny = y + dy;
        if (world.grid_.inBounds(nx, ny) && taken[static_cast<Size>(ny) * world.grid_.w + nx]) clash = true;
      }
    if (clash) continue;
    taken[static_cast<Size>(y) * world.grid_.w + x] = 1;
    Landmark lm;
    lm.p = world.grid_.cellCenter(x, y);
    lm.id = static_cast<std::uint32_t>(placed);
    lm.descSeed = static_cast<std::uint32_t>(rng.next() & 0xFFFFFFFFu);
    lm.height = rng.uniform(0.2, 1.2);
    lm.textured = !(spec.textureless || spec.type == WorldType::kTextureless);
    world.landmarks_.push_back(lm);
    ++placed;
  }
  world.stats_.landmarks = world.landmarks_.size();
  return world;
}

bool World::lineOfSight(const Vec2& a, const Vec2& b) const {
  const Vec2 d = b - a;
  const Scalar len = d.norm();
  if (len < 1e-9) return grid_.atWorld(a) != Cell::kOccupied;
  const Scalar step = grid_.res * 0.5;
  const int nsteps = static_cast<int>(std::ceil(len / step));
  const Vec2 inc = d / static_cast<Scalar>(nsteps);
  Vec2 p = a;
  for (int i = 0; i <= nsteps; ++i) {
    if (grid_.atWorld(p) == Cell::kOccupied) return false;
    p += inc;
  }
  return true;
}

Scalar World::raycast(const Vec2& from, const Vec2& dir, Scalar max_range) const {
  const Scalar len = dir.norm();
  if (len < 1e-12) return max_range;
  const Vec2 u = dir / len;
  const Scalar step = grid_.res * 0.5;
  for (Scalar t = 0; t <= max_range; t += step) {
    const Vec2 p = from + u * t;
    if (!grid_.inBounds(grid_.toCell(p).x(), grid_.toCell(p).y())) return t;
    if (grid_.atWorld(p) == Cell::kOccupied) return t;
  }
  return max_range;
}

Scalar World::clearanceAtWorld(const Vec2& p) const {
  const Eigen::Vector2i c = grid_.toCell(p);
  if (!grid_.inBounds(c.x(), c.y())) return 0.0;
  if (clearance_.empty()) return 0.0;
  return clearance_[static_cast<Size>(c.y()) * grid_.w + c.x()];
}

Scalar World::minCorridorWidth(const std::vector<Vec2>& path) const {
  Scalar best = std::numeric_limits<Scalar>::max();
  for (const Vec2& p : path) best = std::min(best, clearanceAtWorld(p));
  if (best == std::numeric_limits<Scalar>::max()) return 0.0;
  return 2.0 * best;
}

std::vector<Vec2> World::gridPath(const Vec2& from, const Vec2& to, Scalar clearance_m) const {
  const Eigen::Vector2i s = grid_.toCell(from);
  const Eigen::Vector2i g = grid_.toCell(to);
  if (!grid_.inBounds(s.x(), s.y()) || !grid_.inBounds(g.x(), g.y())) return {};
  if (grid_.at(s.x(), s.y()) == Cell::kOccupied || grid_.at(g.x(), g.y()) == Cell::kOccupied) return {};

  // A cell is traversable when it is free and far enough from any obstacle.
  const auto passable = [&](int x, int y) {
    if (grid_.at(x, y) == Cell::kOccupied) return false;
    if (clearance_m <= 0 || clearance_.empty()) return true;
    return clearance_[static_cast<Size>(y) * grid_.w + x] >= clearance_m;
  };

  const Size total = static_cast<Size>(grid_.w) * grid_.h;
  std::vector<Index> parent(total, -1);
  std::vector<std::uint8_t> seen(total, 0);
  std::deque<Index> q;
  const Index start = s.y() * grid_.w + s.x();
  const Index goal = g.y() * grid_.w + g.x();
  seen[start] = 1;
  q.push_back(start);
  const int dx[4] = {1, -1, 0, 0};
  const int dy[4] = {0, 0, 1, -1};
  bool found = false;
  while (!q.empty()) {
    const Index cur = q.front();
    q.pop_front();
    if (cur == goal) {
      found = true;
      break;
    }
    const int cx = cur % grid_.w, cy = cur / grid_.w;
    for (int k = 0; k < 4; ++k) {
      const int nx = cx + dx[k], ny = cy + dy[k];
      if (!grid_.inBounds(nx, ny)) continue;
      if (!passable(nx, ny)) continue;
      const Index ni = ny * grid_.w + nx;
      if (seen[ni]) continue;
      seen[ni] = 1;
      parent[ni] = cur;
      q.push_back(ni);
    }
  }
  if (!found) return {};
  std::vector<Vec2> path;
  for (Index cur = goal; cur != -1; cur = parent[cur]) {
    path.push_back(grid_.cellCenter(cur % grid_.w, cur / grid_.w));
    if (cur == start) break;
  }
  std::reverse(path.begin(), path.end());
  return path;
}

std::vector<Vec2> World::smoothPath(const std::vector<Vec2>& path, Scalar clearance_m,
                                   Scalar min_turn_radius, Scalar spacing, int iterations) const {
  if (path.size() < 3) return path;

  // Resample to uniform spacing first: smoothing a non-uniform polyline biases
  // the result toward whichever segment happened to be sampled more densely.
  std::vector<Vec2> pts;
  pts.push_back(path.front());
  Scalar carry = 0.0;
  for (Size i = 1; i < path.size(); ++i) {
    const Vec2 a = path[i - 1], b = path[i];
    const Scalar len = (b - a).norm();
    if (len < 1e-9) continue;
    Scalar t = (spacing - carry) / len;
    while (t <= 1.0) {
      pts.push_back(a + (b - a) * t);
      t += spacing / len;
    }
    carry = std::fmod(carry + len, spacing);
  }
  if ((pts.back() - path.back()).norm() > 1e-6) pts.push_back(path.back());
  if (pts.size() < 3) return path;

  // Laplacian smoothing with a clearance guard. The guard is what keeps this
  // safe: without it, smoothing a path around a convex corner cuts the corner
  // and walks the line into the wall.
  const Scalar guard = std::max(clearance_m, 0.0);

  // Circumradius of the turn at each interior vertex; the smallest value is the
  // tightest corner in the path. A collinear triple has infinite radius.
  const auto tightest_radius = [&]() {
    Scalar best = std::numeric_limits<Scalar>::max();
    for (Size i = 1; i + 1 < pts.size(); ++i) {
      const Vec2 a = pts[i] - pts[i - 1], b = pts[i + 1] - pts[i];
      const Scalar la = a.norm(), lb = b.norm();
      if (la < 1e-9 || lb < 1e-9) continue;
      const Scalar cross = std::abs(a.x() * b.y() - a.y() * b.x());
      if (cross < 1e-9) continue;
      const Scalar r = (la * lb * (pts[i + 1] - pts[i - 1]).norm()) / (2.0 * cross);
      best = std::min(best, r);
    }
    return best;
  };

  for (int it = 0; it < iterations; ++it) {
    bool moved = false;
    for (Size i = 1; i + 1 < pts.size(); ++i) {
      const Vec2 mid = 0.5 * (pts[i - 1] + pts[i + 1]);
      const Vec2 cand = pts[i] + 0.5 * (mid - pts[i]);
      if (guard > 0.0 && clearanceAtWorld(cand) < guard) continue;
      if ((cand - pts[i]).norm() > 1e-6) moved = true;
      pts[i] = cand;
    }
    if (!moved) break;
    // Stop as soon as the path is gentle enough for the vehicle to follow.
    if (min_turn_radius > 0.0 && tightest_radius() >= min_turn_radius) break;
  }
  return pts;
}

std::vector<Vec2> World::shortcut(const std::vector<Vec2>& path, Scalar clearance_m) const {
  if (path.size() < 3) return path;

  // A segment is admissible when it is collision-free and every sampled point
  // along it keeps the requested clearance. Sampling at roughly the grid
  // resolution is enough: the clearance field is 1-Lipschitz, so a missed dip
  // between samples is bounded by the sample spacing.
  const Scalar step = std::max(grid_.res, 0.05);
  const auto clear_segment = [&](const Vec2& a, const Vec2& b) {
    if (!lineOfSight(a, b)) return false;
    if (clearance_m <= 0.0 || clearance_.empty()) return true;
    const Scalar len = (b - a).norm();
    const int n = std::max(1, static_cast<int>(std::ceil(len / step)));
    for (int k = 0; k <= n; ++k) {
      const Vec2 p = a + (b - a) * (static_cast<Scalar>(k) / static_cast<Scalar>(n));
      if (clearanceAtWorld(p) < clearance_m) return false;
    }
    return true;
  };

  std::vector<Vec2> out;
  out.push_back(path.front());
  Size i = 0;
  while (i + 1 < path.size()) {
    Size j = path.size() - 1;
    // Greedy farthest admissible point.
    for (; j > i + 1; --j)
      if (clear_segment(path[i], path[j])) break;
    out.push_back(path[j]);
    i = j;
  }
  return out;
}

std::string World::toAscii(int max_cols) const {
  // Block-max downsampling: a thin wall must not disappear because the sample
  // grid happened to miss it.
  const int step = std::max(1, grid_.w / std::max(1, max_cols));
  const int vstep = step * 2;
  std::ostringstream os;
  for (int by = grid_.h - vstep; by >= 0; by -= vstep) {
    for (int bx = 0; bx < grid_.w; bx += step) {
      bool occupied = false;
      for (int y = by; y < by + vstep && !occupied; ++y)
        for (int x = bx; x < bx + step; ++x)
          if (grid_.at(x, y) == Cell::kOccupied) {
            occupied = true;
            break;
          }
      os << (occupied ? '#' : ' ');
    }
    os << '\n';
  }
  return os.str();
}

std::vector<std::uint8_t> World::serialise() const {
  std::vector<std::uint8_t> out;
  out.reserve(grid_.cells.size() + 32);
  const auto push32 = [&](std::uint32_t v) {
    out.push_back(static_cast<std::uint8_t>(v & 0xFF));
    out.push_back(static_cast<std::uint8_t>((v >> 8) & 0xFF));
    out.push_back(static_cast<std::uint8_t>((v >> 16) & 0xFF));
    out.push_back(static_cast<std::uint8_t>((v >> 24) & 0xFF));
  };
  push32(static_cast<std::uint32_t>(grid_.w));
  push32(static_cast<std::uint32_t>(grid_.h));
  push32(static_cast<std::uint32_t>(std::round(grid_.res * 1e6)));
  out.insert(out.end(), grid_.cells.begin(), grid_.cells.end());
  return out;
}

std::vector<Vec2> sampleFreePoints(const GridMap& grid, int n, Scalar min_sep, Rng* rng,
                                   Scalar min_clearance_m, const std::vector<Scalar>* clearance) {
  std::vector<Vec2> out;
  const int sep = static_cast<int>(std::ceil(min_sep / grid.res));
  int attempts = 0;
  while (static_cast<int>(out.size()) < n && attempts < n * 200 + 1000) {
    ++attempts;
    const int x = static_cast<int>(rng->uniform(0, grid.w));
    const int y = static_cast<int>(rng->uniform(0, grid.h));
    if (grid.at(x, y) != Cell::kFree) continue;
    const Vec2 p = grid.cellCenter(x, y);
    if (min_clearance_m > 0.0 && clearance != nullptr &&
        (*clearance)[static_cast<Size>(y) * grid.w + x] < min_clearance_m)
      continue;
    bool clash = false;
    for (const Vec2& q : out) {
      if ((q - p).cwiseAbs().maxCoeff() < min_sep) {
        clash = true;
        break;
      }
    }
    (void)sep;
    if (!clash) out.push_back(p);
  }
  return out;
}

}  // namespace sim
}  // namespace sc
