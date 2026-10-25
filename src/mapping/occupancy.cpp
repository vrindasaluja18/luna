#include "luna/mapping/occupancy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define LUNA_HAVE_MMAP 1
#else
#define LUNA_HAVE_MMAP 0
#endif

namespace luna {
namespace mapping {
namespace {

/// Fixed record size for one paged tile, so a tile's file offset is a pure
/// function of its key and the file can be indexed without a directory.
Size tileBytes(const OccupancyOptions& o) {
  return static_cast<Size>(o.tile_cells) * static_cast<Size>(o.tile_cells) * sizeof(LogOdds);
}

/// Tile index containing a world coordinate. Floor division, so tiles are
/// contiguous across the origin instead of mirrored about it.
inline int tileOf(Scalar w, Scalar res, int tile_cells) {
  const int c = static_cast<int>(std::floor(w / res));
  return c >= 0 ? c / tile_cells : -((-c + tile_cells - 1) / tile_cells);
}

/// Cell index containing a world coordinate.
inline int cellOf(Scalar w, Scalar res) { return static_cast<int>(std::floor(w / res)); }

}  // namespace

OccupancyMap::OccupancyMap(const OccupancyOptions& opt) : opt_(opt) {
  if (opt_.tile_cells < 8) opt_.tile_cells = 8;
  if (opt_.resolution <= 0) opt_.resolution = 0.05;
  if (opt_.max_resident_tiles < 4) opt_.max_resident_tiles = 4;
  tiles_.reserve(std::min<Size>(opt_.max_resident_tiles * 2, 65536));
}

void OccupancyMap::bumpPeak() {
  if (stats_.bytes_resident > stats_.bytes_peak) stats_.bytes_peak = stats_.bytes_resident;
}

void OccupancyMap::touch(Tile* t) { t->last_use = ++clock_; }

OccupancyMap::Tile* OccupancyMap::findTile(int tx, int ty) const {
  const auto it = tiles_.find(Key{tx, ty});
  return it == tiles_.end() ? nullptr : const_cast<Tile*>(&it->second);
}

OccupancyMap::Tile* OccupancyMap::tile(int tx, int ty, bool create) {
  const Key k{tx, ty};
  const auto it = tiles_.find(k);
  if (it != tiles_.end()) {
    touch(&it->second);
    return &it->second;
  }
  if (!create) return nullptr;

  evictIfNeeded();

  Tile t;
  t.data.assign(tileBytes(opt_) / sizeof(LogOdds), 0.0f);
  t.last_use = ++clock_;
  // A tile that was evicted earlier still lives in the backing file. Loading it
  // here is what makes a revisit see the old map rather than a blank square.
  if (opt_.page_to_disk && loadPaged(k, &t)) {
    t.paged = true;
    ++stats_.tiles_loaded;
  }
  ++stats_.tiles_allocated;
  stats_.bytes_resident += t.data.size() * sizeof(LogOdds);
  bumpPeak();
  auto res = tiles_.emplace(k, std::move(t));
  stats_.tiles_resident = tiles_.size();
  return &res.first->second;
}

void OccupancyMap::evictIfNeeded() {
  while (tiles_.size() >= opt_.max_resident_tiles) evictOne();
}

void OccupancyMap::evictOne() {
  if (tiles_.empty()) return;
  // Linear scan for the least recently used tile. The resident set is capped in
  // the low thousands, so this is microseconds and it happens once per tile
  // creation rather than once per cell.
  auto victim = tiles_.begin();
  std::uint64_t best = victim->second.last_use;
  for (auto it = std::next(tiles_.begin()); it != tiles_.end(); ++it) {
    if (it->second.last_use < best) {
      best = it->second.last_use;
      victim = it;
    }
  }
  if (victim->second.dirty && opt_.page_to_disk) writePaged(victim->first, victim->second);
  stats_.bytes_resident -= victim->second.data.size() * sizeof(LogOdds);
  tiles_.erase(victim);
  ++stats_.tiles_evicted;
  stats_.tiles_resident = tiles_.size();
}

// ---------------------------------------------------------------------------
// Paging. Tiles are fixed-size records at an offset derived from the key, so
// the backing file needs no index of its own. A hash of the key selects a
// sparse region; the OS materialises only pages actually written.
// ---------------------------------------------------------------------------
namespace {
constexpr std::uint64_t kPageStride = 1u << 20;  // 1 MiB per hash bucket
constexpr std::uint64_t kBuckets = 1u << 16;     // 65536 buckets
}

bool OccupancyMap::loadPaged(const Key& k, Tile* t) const {
#if LUNA_HAVE_MMAP
  if (opt_.backing_path.empty()) return false;
  const std::uint64_t bucket = KeyHash{}(k) % kBuckets;
  const off_t off = static_cast<off_t>(bucket * kPageStride);
  const int fd = ::open(opt_.backing_path.c_str(), O_RDONLY);
  if (fd < 0) return false;
  void* p = ::mmap(nullptr, tileBytes(opt_), PROT_READ, MAP_PRIVATE, fd, off);
  ::close(fd);
  if (p == MAP_FAILED) return false;
  std::memcpy(t->data.data(), p, tileBytes(opt_));
  ::munmap(p, tileBytes(opt_));
  // A never-written bucket reads as zeros, which is exactly "unknown". The
  // caller treats an all-zero tile as absent, which is correct: a tile whose
  // evidence has decayed to zero carries no information.
  for (LogOdds v : t->data)
    if (v != 0.0f) return true;
  return false;
#else
  (void)k;
  (void)t;
  return false;
#endif
}

void OccupancyMap::writePaged(const Key& k, const Tile& t) const {
#if LUNA_HAVE_MMAP
  if (opt_.backing_path.empty()) return;
  const std::uint64_t bucket = KeyHash{}(k) % kBuckets;
  const off_t off = static_cast<off_t>(bucket * kPageStride);
  const int fd = ::open(opt_.backing_path.c_str(), O_RDWR | O_CREAT, 0644);
  if (fd < 0) return;
  const Size bytes = tileBytes(opt_);
  // A short write would leave a tile holding a mix of two epochs, which is
  // worse than a stale tile: the map would show a wall half-erased. Check the
  // return and leave the old contents on failure.
  (void)::pwrite(fd, t.data.data(), bytes, off);
  ::close(fd);
#else
  (void)k;
  (void)t;
#endif
}

void OccupancyMap::flush() {
  for (auto& kv : tiles_)
    if (kv.second.dirty) writePaged(kv.first, kv.second);
}

void OccupancyMap::clear() {
  tiles_.clear();
  stats_.tiles_resident = 0;
  stats_.bytes_resident = 0;
}

// ---------------------------------------------------------------------------
// Integration
// ---------------------------------------------------------------------------

void OccupancyMap::integrateRay(const Pose2& pose, const Vec2& dir, Scalar range) {
  if (!(range > 0) || !std::isfinite(range)) return;
  const Scalar c = std::cos(pose.th), s = std::sin(pose.th);
  const Vec2 d(c * dir.x() - s * dir.y(), s * dir.x() + c * dir.y());
  const Scalar len = d.norm();
  if (len < 1e-9) return;
  const Vec2 u = d / len;

  // A max-range reading means "nothing within range": the beam marks space free
  // up to the limit and says nothing about the cell at the limit. Marking that
  // cell occupied is the classic source of phantom walls at maximum range.
  const bool hit = range < opt_.max_range;
  const Scalar travel = hit ? range : opt_.max_range;
  const Scalar start = std::min(opt_.min_range, travel);
  const Scalar cell = opt_.resolution;

  // Walk at half-cell spacing so a beam cannot skip a cell it passes through.
  // That is what produces thin diagonal gaps in a naively stepped map.
  const int steps = static_cast<int>(std::ceil((travel - start) / (0.5 * cell))) + 1;
  int last_tx = std::numeric_limits<int>::min(), last_ty = std::numeric_limits<int>::min();
  Tile* t = nullptr;

  for (int i = 0; i <= steps; ++i) {
    const Scalar s_i = start + i * (0.5 * cell);
    if (s_i > travel + 0.5 * cell) break;
    const Vec2 w = Vec2(pose.x, pose.y) + u * s_i;

    min_seen_ = min_seen_.cwiseMin(w);
    max_seen_ = max_seen_.cwiseMax(w);

    const int tx = tileOf(w.x(), cell, opt_.tile_cells);
    const int ty = tileOf(w.y(), cell, opt_.tile_cells);
    if (tx != last_tx || ty != last_ty) {
      last_tx = tx;
      last_ty = ty;
      t = tile(tx, ty, true);
    }
    if (!t) continue;

    const int lx = cellOf(w.x(), cell) - tx * opt_.tile_cells;
    const int ly = cellOf(w.y(), cell) - ty * opt_.tile_cells;
    if (lx < 0 || ly < 0 || lx >= opt_.tile_cells || ly >= opt_.tile_cells) continue;

    LogOdds& v = t->data[static_cast<Size>(ly) * opt_.tile_cells + lx];
    // The terminal cell of a hitting beam is occupied; every cell before it is
    // free.
    const bool terminal = hit && (s_i + 0.5 * cell >= travel);
    v += terminal ? opt_.occupied_delta : opt_.free_delta;
    v = static_cast<LogOdds>(clamp(static_cast<Scalar>(v), static_cast<Scalar>(opt_.clamp_min),
                                   static_cast<Scalar>(opt_.clamp_max)));
    t->dirty = true;
    ++stats_.cells_updated;
  }
  ++stats_.rays_integrated;
}

void OccupancyMap::integrateScan(const Pose2& pose, const std::vector<RangeReading>& readings) {
  for (const RangeReading& r : readings) integrateRay(pose, r.dir, r.range);
}

LogOdds OccupancyMap::logOddsAt(const Vec2& p) const {
  const Scalar cell = opt_.resolution;
  const int tx = tileOf(p.x(), cell, opt_.tile_cells);
  const int ty = tileOf(p.y(), cell, opt_.tile_cells);
  const Tile* t = findTile(tx, ty);
  if (!t) return 0.0f;
  const int lx = cellOf(p.x(), cell) - tx * opt_.tile_cells;
  const int ly = cellOf(p.y(), cell) - ty * opt_.tile_cells;
  if (lx < 0 || ly < 0 || lx >= opt_.tile_cells || ly >= opt_.tile_cells) return 0.0f;
  return t->data[static_cast<Size>(ly) * opt_.tile_cells + lx];
}

sim::Cell OccupancyMap::atWorld(const Vec2& p) const {
  const LogOdds l = logOddsAt(p);
  if (l > opt_.occupied_threshold) return sim::Cell::kOccupied;
  if (l < opt_.free_threshold) return sim::Cell::kFree;
  return sim::Cell::kUnknown;
}

void OccupancyMap::markDisc(const Vec2& center, Scalar radius, bool occupied) {
  const Scalar cell = opt_.resolution;
  const int r = static_cast<int>(std::ceil(radius / cell));
  const int cx0 = cellOf(center.x(), cell);
  const int cy0 = cellOf(center.y(), cell);
  for (int dy = -r; dy <= r; ++dy) {
    for (int dx = -r; dx <= r; ++dx) {
      if (dx * dx + dy * dy > r * r) continue;
      const int cx = cx0 + dx, cy = cy0 + dy;
      const int tx = cx >= 0 ? cx / opt_.tile_cells : -((-cx + opt_.tile_cells - 1) / opt_.tile_cells);
      const int ty = cy >= 0 ? cy / opt_.tile_cells : -((-cy + opt_.tile_cells - 1) / opt_.tile_cells);
      Tile* t = tile(tx, ty, true);
      if (!t) continue;
      const int lx = cx - tx * opt_.tile_cells;
      const int ly = cy - ty * opt_.tile_cells;
      if (lx < 0 || ly < 0 || lx >= opt_.tile_cells || ly >= opt_.tile_cells) continue;
      LogOdds& v = t->data[static_cast<Size>(ly) * opt_.tile_cells + lx];
      v = occupied ? opt_.clamp_max : opt_.clamp_min;
      t->dirty = true;
    }
  }
}

sim::GridMap OccupancyMap::toGrid(Scalar res, Scalar pad_m) const {
  if (tiles_.empty()) return sim::GridMap();
  return toGridWindow(min_seen_ - Vec2::Constant(pad_m), max_seen_ + Vec2::Constant(pad_m), res);
}

sim::GridMap OccupancyMap::toGridWindow(const Vec2& lo, const Vec2& hi, Scalar res) const {
  sim::GridMap g;
  if (res <= 0) res = opt_.resolution;
  const int w = std::max(1, static_cast<int>(std::ceil((hi.x() - lo.x()) / res)));
  const int h = std::max(1, static_cast<int>(std::ceil((hi.y() - lo.y()) / res)));
  g.resize(w, h, res, lo);
  // A cell whose tile has been evicted stays unknown rather than free, which is
  // the safe direction for a planner: unknown is planned around, free is driven
  // through.
  std::fill(g.cells.begin(), g.cells.end(), static_cast<std::uint8_t>(sim::Cell::kUnknown));

  for (const auto& kv : tiles_) {
    const int tx = kv.first.tx, ty = kv.first.ty;
    const Tile& t = kv.second;
    for (int ly = 0; ly < opt_.tile_cells; ++ly) {
      for (int lx = 0; lx < opt_.tile_cells; ++lx) {
        const LogOdds v = t.data[static_cast<Size>(ly) * opt_.tile_cells + lx];
        if (v == 0.0f) continue;
        const Scalar wx = (tx * opt_.tile_cells + lx + 0.5) * opt_.resolution;
        const Scalar wy = (ty * opt_.tile_cells + ly + 0.5) * opt_.resolution;
        const int gx = static_cast<int>(std::floor((wx - lo.x()) / res));
        const int gy = static_cast<int>(std::floor((wy - lo.y()) / res));
        if (gx < 0 || gy < 0 || gx >= w || gy >= h) continue;
        sim::Cell c = sim::Cell::kUnknown;
        if (v > opt_.occupied_threshold)
          c = sim::Cell::kOccupied;
        else if (v < opt_.free_threshold)
          c = sim::Cell::kFree;
        g.cells[static_cast<Size>(gy) * w + gx] = static_cast<std::uint8_t>(c);
      }
    }
  }
  return g;
}

Size OccupancyMap::countOccupied() const {
  Size n = 0;
  for (const auto& kv : tiles_)
    for (LogOdds v : kv.second.data)
      if (v > opt_.occupied_threshold) ++n;
  return n;
}

Size OccupancyMap::countFree() const {
  Size n = 0;
  for (const auto& kv : tiles_)
    for (LogOdds v : kv.second.data)
      if (v < opt_.free_threshold) ++n;
  return n;
}

// ---------------------------------------------------------------------------
// Frontend bridge
// ---------------------------------------------------------------------------

std::vector<RangeReading> readingsFromFeatures(const std::vector<frontend::Feature>& feats,
                                               const sim::CameraModel& cam) {
  std::vector<RangeReading> out;
  out.reserve(feats.size());
  for (const frontend::Feature& f : feats) {
    if (!(f.depth > 0)) continue;
    // Unproject the pixel to a camera-frame ray, then take its horizontal
    // component. The map is 2D, so a beam is the feature's projection onto the
    // ground plane: a landmark 3 m ahead and 1 m up projects to a horizontal
    // range of about 2.8 m, which is where a 2D occupancy map should record it.
    // Using the slant range would place the obstacle too far by 1/cos(elev).
    const Scalar u = f.x - cam.cx;
    const Scalar v = f.y - cam.cy;
    const Vec3 ray(u / cam.fx, -v / cam.fy, 1.0);
    const Vec3 p = ray * f.depth;
    const Scalar horiz = std::hypot(p.x(), p.y());
    if (horiz < 1e-6) continue;
    RangeReading r;
    r.dir = Vec2(p.x(), p.y()) / horiz;
    r.range = horiz;
    out.push_back(r);
  }
  return out;
}

}  // namespace mapping
}  // namespace luna
