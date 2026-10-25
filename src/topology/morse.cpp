#include "slamcore/topology/morse.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <functional>
#include <sstream>

namespace sc {
namespace topology {
namespace {

constexpr Scalar kNegInf = -std::numeric_limits<Scalar>::max();

/// Disjoint-set over cell indices with path compression and union by size.
/// This is the workhorse of the flooding: it is what turns "which components
/// have merged" into an O(alpha) query, so the whole extraction is O(n log n)
/// dominated by the initial sort rather than by the union-find.
class DisjointSet {
 public:
  explicit DisjointSet(Size n) : parent_(n), size_(n, 1) {
    std::iota(parent_.begin(), parent_.end(), 0);
  }
  Index find(Index x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }
  /// Returns the surviving root. `merged`, when non-null, reports whether a
  /// real merge happened (two distinct components) rather than a redundant
  /// union. It is optional so the common "just attach this cell" case does not
  /// have to pass a throwaway.
  Index unite(Index a, Index b, bool* merged) {
    Index ra = find(a), rb = find(b);
    if (ra == rb) {
      if (merged) *merged = false;
      return ra;
    }
    if (size_[ra] < size_[rb]) std::swap(ra, rb);
    parent_[rb] = ra;
    size_[ra] += size_[rb];
    if (merged) *merged = true;
    return ra;
  }
  bool same(Index a, Index b) { return find(a) == find(b); }

 private:
  std::vector<Index> parent_;
  std::vector<Index> size_;
};

/// 4-connected neighbour offsets. 8-connectivity is deliberately not used:
/// with diagonal links the superlevel sets of a clearance field connect
/// through single-cell pinches at wall corners, which merges two rooms that a
/// real doorway does not join and destroys exactly the saddle the analysis is
/// looking for.
constexpr int kDx[4] = {1, -1, 0, 0};
constexpr int kDy[4] = {0, 0, 1, -1};

}  // namespace

const char* toString(CriticalType t) {
  switch (t) {
    case CriticalType::kMinimum: return "minimum";
    case CriticalType::kSaddle: return "saddle";
    case CriticalType::kMaximum: return "maximum";
  }
  return "?";
}

// ---------------------------------------------------------------------------
// Field utilities
// ---------------------------------------------------------------------------

void clearanceTransform(const sim::GridMap& grid, std::vector<Scalar>* out) {
  const int w = grid.w, h = grid.h;
  if (w <= 0 || h <= 0) {
    out->clear();
    return;
  }
  const Scalar big = 1e12;
  out->assign(static_cast<Size>(w) * h, big);
  // Seed: occupied cells are distance 0; everything else starts at big.
  for (Size i = 0; i < out->size(); ++i)
    if (static_cast<sim::Cell>(grid.cells[i]) == sim::Cell::kOccupied) (*out)[i] = 0.0;

  const auto dt1d = [big](int n, Scalar* f, Scalar* res) {
    std::vector<int> v(static_cast<Size>(n));
    std::vector<Scalar> z(static_cast<Size>(n) + 1);
    int k = 0;
    v[0] = 0;
    z[0] = -big;
    z[1] = big;
    for (int q = 1; q < n; ++q) {
      Scalar s = 0;
      while (true) {
        s = ((f[q] + Scalar(q) * q) - (f[v[k]] + Scalar(v[k]) * v[k])) /
            (2.0 * q - 2.0 * v[k]);
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
      res[q] = dq * dq + f[v[k]];
    }
  };

  std::vector<Scalar> f(static_cast<Size>(std::max(w, h)));
  for (int x = 0; x < w; ++x) {
    for (int y = 0; y < h; ++y) f[y] = (*out)[static_cast<Size>(y) * w + x];
    dt1d(h, f.data(), f.data());
    for (int y = 0; y < h; ++y) (*out)[static_cast<Size>(y) * w + x] = f[y];
  }
  for (int y = 0; y < h; ++y) {
    Scalar* row = out->data() + static_cast<Size>(y) * w;
    dt1d(w, row, row);
  }
  // The transform above yields squared distance in cell units; convert to
  // metres so every threshold in this module can be stated in physical units.
  const Scalar r = grid.res;
  for (Scalar& v : *out) v = std::sqrt(v) * r;
}

void gaussianBlurField(const sim::GridMap& grid, const std::vector<Scalar>& in, Scalar sigma,
                       std::vector<Scalar>* out) {
  const int w = grid.w, h = grid.h;
  out->assign(in.size(), 0.0);
  if (sigma <= 0) {
    *out = in;
    return;
  }
  const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));
  std::vector<Scalar> kernel(static_cast<Size>(2 * radius + 1));
  Scalar ksum = 0;
  for (int i = -radius; i <= radius; ++i) {
    const Scalar v = std::exp(-0.5 * (i / sigma) * (i / sigma));
    kernel[static_cast<Size>(i + radius)] = v;
    ksum += v;
  }
  for (Scalar& v : kernel) v /= ksum;

  // Valid (non-missing) mask: occupied cells and out-of-grid cells are missing
  // and are excluded from both numerator and denominator. Including them would
  // pull the height down near every wall, which flattens the very ridges the
  // saddle detection depends on.
  std::vector<Scalar> valid(in.size(), 0.0);
  for (Size i = 0; i < in.size(); ++i)
    valid[i] = (static_cast<sim::Cell>(grid.cells[i]) == sim::Cell::kOccupied) ? 0.0 : 1.0;

  std::vector<Scalar> tmp(in.size(), 0.0), num(in.size(), 0.0);
  std::vector<Scalar> wsum(in.size(), 0.0);
  // Horizontal pass.
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      Scalar acc = 0, wacc = 0;
      for (int k = -radius; k <= radius; ++k) {
        const int xx = x + k;
        if (xx < 0 || xx >= w) continue;
        const Size idx = static_cast<Size>(y) * w + xx;
        const Scalar kv = kernel[static_cast<Size>(k + radius)];
        acc += kv * in[idx] * valid[idx];
        wacc += kv * valid[idx];
      }
      const Size idx = static_cast<Size>(y) * w + x;
      num[idx] = acc;
      wsum[idx] = wacc;
    }
  }
  for (Size i = 0; i < in.size(); ++i)
    tmp[i] = wsum[i] > 1e-12 ? num[i] / wsum[i] : in[i];
  // Vertical pass.
  std::fill(num.begin(), num.end(), 0.0);
  std::fill(wsum.begin(), wsum.end(), 0.0);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      Scalar acc = 0, wacc = 0;
      for (int k = -radius; k <= radius; ++k) {
        const int yy = y + k;
        if (yy < 0 || yy >= h) continue;
        const Size idx = static_cast<Size>(yy) * w + x;
        const Scalar kv = kernel[static_cast<Size>(k + radius)];
        acc += kv * tmp[idx] * valid[idx];
        wacc += kv * valid[idx];
      }
      const Size idx = static_cast<Size>(y) * w + x;
      num[idx] = acc;
      wsum[idx] = wacc;
    }
  }
  for (Size i = 0; i < in.size(); ++i)
    (*out)[i] = wsum[i] > 1e-12 ? num[i] / wsum[i] : in[i];
}

// ---------------------------------------------------------------------------
// Morse analysis
// ---------------------------------------------------------------------------

MorseComplex MorseComplex::analyseField(const sim::GridMap& grid, const std::vector<Scalar>& field,
                                        const MorseOptions& opt) {
  MorseComplex mc;
  mc.grid_ = grid;
  mc.opt_ = opt;
  const int w = grid.w, h = grid.h;
  const Size n = static_cast<Size>(w) * h;
  if (n == 0 || field.size() != n) return mc;

  // Missing cells (occupied, or below the drivable clearance) are removed from
  // the analysis by being pushed to -infinity. They then never join a
  // superlevel set and never merge anything, which is the discrete equivalent
  // of deleting them from the domain.
  std::vector<Scalar> hf(n, kNegInf);
  for (Size i = 0; i < n; ++i) {
    if (static_cast<sim::Cell>(grid.cells[i]) == sim::Cell::kOccupied) continue;
    if (field[i] < opt.min_clearance) continue;
    hf[i] = field[i];
    ++mc.stats_.free_cells;
  }
  mc.height_ = hf;

  // ---------------------------------------------------------------------
  // Join-tree flooding.
  //
  // Cells are visited in descending height, and each is unioned with its
  // already-visited neighbours. The subtlety is equal heights: the clearance
  // field is full of flat plateaus (the centreline of a straight corridor is
  // exactly flat), and two *disconnected* plateaus can share the same height by
  // symmetry. Grouping cells by height alone therefore merges rooms that are
  // nowhere near each other, which destroys their local maxima -- measured as
  // 3 rooms extracted from a 9-room world. The plateau is therefore identified
  // by connectivity, not by value:
  //
  //   * a cell with no strictly-higher neighbour and no already-visited
  //     equal-height neighbour *starts* a plateau, and is its candidate
  //     maximum. Exactly one cell per plateau qualifies, because every other
  //     cell of the plateau has the first one as an equal-height neighbour by
  //     the time it is visited;
  //   * equal-height neighbours are unioned first, so the plateau is one
  //     component before it meets anything higher;
  //   * the first strictly-higher component it touches simply absorbs it
  //     (that is an ordinary descent); each *further* higher component it
  //     touches is a genuine saddle, because two previously separate basins
  //     meet there. Counting the first attachment would make every cell in the
  //     map a saddle.
  // ---------------------------------------------------------------------
  std::vector<Index> order;
  order.reserve(mc.stats_.free_cells);
  for (Size i = 0; i < n; ++i)
    if (hf[i] != kNegInf) order.push_back(static_cast<Index>(i));
  std::sort(order.begin(), order.end(), [&](Index a, Index b) { return hf[a] > hf[b]; });

  DisjointSet ds(n);
  std::vector<std::uint8_t> processed(n, 0);
  std::vector<Index> comp_max(n, -1);   ///< per root: the plateau's candidate maximum
  std::vector<std::uint8_t> is_max(n, 0);
  /// For every cell, the maximum of the basin it drained into. Recorded as the
  /// flood proceeds, because at the end everything in a connected region shares
  /// one component and the per-basin assignment is no longer recoverable from
  /// the union-find alone.
  std::vector<Index> cell_peak(n, -1);

  /// One merge of two basins, in the order it happened.
  struct Merge {
    Index winner = -1;  ///< cell of the surviving maximum
    Index loser = -1;   ///< cell of the maximum whose basin was absorbed
    Index saddle = -1;  ///< cell where the merge happened
    Scalar persistence = 0;
  };
  std::vector<Merge> merges;

  const auto neighbours = [&](Index c, Index out[4]) {
    const int cx = static_cast<int>(c % w), cy = static_cast<int>(c / w);
    int k = 0;
    for (int d = 0; d < 4; ++d) {
      const int nx = cx + kDx[d], ny = cy + kDy[d];
      if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
      const Index nb = static_cast<Index>(ny) * w + nx;
      if (hf[nb] == kNegInf) continue;
      out[k++] = nb;
    }
    return k;
  };

  Index nbr[4];
  for (Size oi = 0; oi < order.size(); ++oi) {
    const Index c = order[oi];
    const Scalar hc = hf[c];
    const int cnt = neighbours(c, nbr);
    processed[c] = 1;

    // Distinct roots of the processed neighbours. Descending order means every
    // processed neighbour is at height >= h(c), and equal-height neighbours are
    // included deliberately: they are how a flat plateau becomes one component
    // instead of a chain of cells. Because only *adjacent* cells are ever
    // united, two disconnected plateaus that happen to share a height stay
    // separate -- grouping by value instead of connectivity is what made a
    // nine-room world extract as a single room.
    Index roots[4];
    int nroots = 0;
    for (int t = 0; t < cnt; ++t) {
      if (!processed[nbr[t]]) continue;
      const Index r = ds.find(nbr[t]);
      bool seen = false;
      for (int q = 0; q < nroots; ++q)
        if (roots[q] == r) {
          seen = true;
          break;
        }
      if (!seen && nroots < 4) roots[nroots++] = r;
    }

    if (nroots == 0) {
      // No neighbour at or above this height: the highest cell of a new
      // plateau, and therefore a local maximum. Exactly one cell per plateau
      // qualifies, because every other cell of it has the first one as an
      // already-processed equal-height neighbour.
      const Index r = ds.find(c);
      comp_max[r] = c;
      is_max[c] = 1;
      cell_peak[c] = c;
      ++mc.stats_.raw_maxima;
      continue;
    }

    // Absorb into the neighbouring components. The first attachment is an
    // ordinary descent -- this cell is simply lower than what is above it. Each
    // *further* distinct component is a saddle: two basins that were separate
    // above this height become one here, and the younger one's maximum dies.
    //
    // Roots are re-resolved through `find` on every use. They were captured
    // before any of these unions, so after the first union the later entries
    // are stale; uniting against a stale root is a no-op that silently drops
    // the merge, and reading a stale `comp_max` slot reads a value that has
    // since moved. Together those made a 236-maximum world record five merges
    // and extract as one room.
    Index proot = ds.find(c);
    Index cur_max = comp_max[proot];
    int attached = 0;
    for (int t = 0; t < nroots; ++t) {
      const Index r = ds.find(roots[t]);
      if (r == proot) continue;
      const Index other_max = comp_max[r];
      // Order the two sides by peak height. A component whose peak slot is
      // still empty is treated as lower than any real peak, which is what it
      // is: `comp_max` is filled in as the flood descends, so the cell that
      // creates a component is its first peak and every later cell inherits
      // one. Comparing the slots directly without this guard keeps a null
      // slot on the "winner" side whenever the two heights tie, and a null
      // winner makes the merge meaningless.
      const Scalar h_cur = cur_max >= 0 ? hf[cur_max] : -1.0;
      const Scalar h_oth = other_max >= 0 ? hf[other_max] : -1.0;
      const bool cur_higher = h_cur >= h_oth;
      const Index hi = cur_higher ? cur_max : other_max;
      const Index lo = cur_higher ? other_max : cur_max;
      const Scalar h_hi = std::max(h_cur, h_oth);
      const Scalar h_lo = std::min(h_cur, h_oth);
      proot = ds.unite(proot, r, nullptr);
      if (attached++ == 0) {
        // Ordinary descent: this cell is simply lower than what is above it.
        cur_max = hi;
        continue;
      }
      Merge m;
      m.winner = hi;
      m.loser = lo;
      m.saddle = c;
      m.persistence = (lo >= 0) ? std::max<Scalar>(0.0, h_hi - hc) : 0.0;
      // Equal-height peaks are one feature, not two: a flat-topped room under
      // a symmetric smoothing kernel produces mirrored maxima of exactly equal
      // height. The merge is still recorded so the room decomposition can
      // union the two basins -- dropping it leaves them as separate rooms,
      // which is how a nine-room world extracted as two hundred -- but it is
      // never a doorway and is not counted as a saddle, because nothing
      // meaningful meets there.
      if (lo >= 0 && std::abs(h_hi - h_lo) > 1e-9) {
        ++mc.stats_.raw_saddles;
        if (m.persistence > mc.stats_.max_persistence) mc.stats_.max_persistence = m.persistence;
      }
      merges.push_back(m);
      cur_max = hi;
    }
    comp_max[proot] = cur_max;
    cell_peak[c] = cur_max;
  }

  // ---------------------------------------------------------------------
  // Persistence bookkeeping.
  //
  // A maximum's persistence is the height difference to the saddle at which
  // its basin was absorbed. It is recorded for reporting and for the corridor
  // metadata, but -- importantly -- it is NOT used to decide which maxima are
  // rooms. See the room decomposition below for why.
  // ---------------------------------------------------------------------
  std::vector<Scalar> max_persist(n, 0.0);
  for (const Merge& m : merges)
    if (m.loser >= 0) max_persist[m.loser] = m.persistence;

  Scalar field_min = std::numeric_limits<Scalar>::max();
  Scalar field_max = kNegInf;
  for (Size i = 0; i < n; ++i) {
    if (hf[i] == kNegInf) continue;
    field_min = std::min(field_min, hf[i]);
    field_max = std::max(field_max, hf[i]);
  }
  if (field_min > field_max) field_min = field_max;

  std::vector<Index> parent_max(n, -1);
  for (const Merge& m : merges)
    if (m.loser >= 0) parent_max[m.loser] = m.winner;

  // ---------------------------------------------------------------------
  // Room decomposition by bottleneck contrast.
  //
  // Persistence cannot separate rooms here, and it is worth being explicit
  // about why, because the obvious design (threshold the persistence) is
  // wrong on exactly the worlds this module exists for. The generator produces
  // rooms of similar size, so their clearance maxima sit at nearly the *same*
  // height; the basins therefore merge at a height barely below the peaks, and
  // every room-to-room merge has near-zero persistence. A persistence
  // threshold deletes all of them and the world extracts as a single room.
  //
  // What actually distinguishes a doorway is contrast: the vehicle has to
  // squeeze through, so the gap is far narrower than the space on either side.
  // A merge is a doorway when its saddle is below `min_bottleneck_contrast`
  // times the narrower of the two peaks. Every other merge joins two peaks of
  // the same open space, so those peaks are unioned into one room.
  //
  // This also absorbs the two other failure modes seen in practice: the
  // fragmented ripple of a large flat hall (peaks merge at almost their own
  // height -> interior -> one room), and mirrored equal-height peak pairs from
  // a symmetric smoothing kernel (same thing).
  // ---------------------------------------------------------------------
  std::vector<Index> peak_comp(n, -1);
  Index npeaks = 0;
  for (Size i = 0; i < n; ++i)
    if (is_max[i]) peak_comp[i] = npeaks++;

  std::vector<Index> room_parent(static_cast<Size>(npeaks));
  for (Index i = 0; i < npeaks; ++i) room_parent[static_cast<Size>(i)] = i;
  std::function<Index(Index)> room_find = [&](Index x) {
    while (room_parent[static_cast<Size>(x)] != x) {
      room_parent[static_cast<Size>(x)] =
          room_parent[static_cast<Size>(room_parent[static_cast<Size>(x)])];
      x = room_parent[static_cast<Size>(x)];
    }
    return x;
  };

  /// A doorway: the merge's saddle is much narrower than the regions it joins.
  std::vector<std::uint8_t> is_doorway(merges.size(), 0);
  for (Size k = 0; k < merges.size(); ++k) {
    const Merge& m = merges[k];
    if (m.winner < 0 || m.loser < 0) continue;
    const Scalar narrower = std::min(hf[m.winner], hf[m.loser]);
    is_doorway[k] =
        (narrower > 1e-9 && hf[m.saddle] < opt.min_bottleneck_contrast * narrower) ? 1 : 0;
  }
  // Union every non-doorway merge. Done as a second pass so the classification
  // is not order-dependent: a peak can be merged by several edges and the
  // interior ones must all be applied before the rooms are read off.
  for (Size k = 0; k < merges.size(); ++k) {
    if (is_doorway[k]) continue;
    const Index a = room_find(peak_comp[merges[k].winner]);
    const Index b = room_find(peak_comp[merges[k].loser]);
    if (a != b) room_parent[static_cast<Size>(a)] = b;
  }

  // Node per room: the most open cell in it, which is where a route should aim.
  //
  // It is NOT the representative plateau cell. A flat plateau's representative
  // is whichever cell the flood happened to visit first, and that cell can sit
  // anywhere on the plateau -- including its edge, against a wall. Aiming there
  // put the planned path within 5 cm of an obstacle on worlds whose rooms are
  // metres wide. Every cell is therefore assigned to its room through the basin
  // it drained into, and the room's node is its maximum-clearance cell.
  std::vector<Index> room_best(static_cast<Size>(npeaks), -1);
  for (Size i = 0; i < n; ++i) {
    if (hf[i] == kNegInf || cell_peak[i] < 0) continue;
    const Index r = room_find(peak_comp[cell_peak[i]]);
    if (room_best[static_cast<Size>(r)] < 0 || hf[i] > hf[room_best[static_cast<Size>(r)]])
      room_best[static_cast<Size>(r)] = static_cast<Index>(i);
  }
  std::vector<int> node_of_room(static_cast<Size>(npeaks), -1);
  for (Index r = 0; r < npeaks; ++r) {
    const Index best = room_best[static_cast<Size>(r)];
    if (best < 0) continue;
    if (mc.graph_.nodes.size() >= opt.max_nodes) break;
    CriticalPoint cp;
    cp.cx = static_cast<int>(best % w);
    cp.cy = static_cast<int>(best / w);
    cp.p = grid.cellCenter(cp.cx, cp.cy);
    cp.h = hf[best];
    cp.type = CriticalType::kMaximum;
    cp.persistence = max_persist[best];
    node_of_room[static_cast<Size>(r)] = static_cast<int>(mc.graph_.nodes.size());
    mc.graph_.nodes.push_back(cp);
    ++mc.stats_.kept_maxima;
  }

  // Corridors: the doorway merges, remapped onto rooms. Two rooms joined by
  // two doorways get one edge carrying the wider door, because the route
  // through them should use the roomier one.
  std::vector<CorridorEdge> edges;
  std::vector<Index> edge_saddle_cell;
  for (Size k = 0; k < merges.size(); ++k) {
    if (!is_doorway[k]) continue;
    const Index ra = room_find(peak_comp[merges[k].winner]);
    const Index rb = room_find(peak_comp[merges[k].loser]);
    if (ra == rb) continue;  // the door leads back into the same room
    const int na = node_of_room[static_cast<Size>(ra)];
    const int nb2 = node_of_room[static_cast<Size>(rb)];
    if (na < 0 || nb2 < 0) continue;
    CorridorEdge e;
    e.a = std::min(na, nb2);
    e.b = std::max(na, nb2);
    e.bottleneck = hf[merges[k].saddle];
    e.persistence = merges[k].persistence;
    e.length = (mc.graph_.nodes[static_cast<Size>(e.a)].p -
                mc.graph_.nodes[static_cast<Size>(e.b)].p)
                   .norm();
    bool found = false;
    for (Size q = 0; q < edges.size(); ++q) {
      if (edges[q].a == e.a && edges[q].b == e.b) {
        if (e.bottleneck > edges[q].bottleneck) {
          edges[q] = e;
          edge_saddle_cell[q] = merges[k].saddle;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      edges.push_back(e);
      edge_saddle_cell.push_back(merges[k].saddle);
    }
  }

  // Turn each corridor into a path *through* its doorway: two edges, room ->
  // saddle -> room, rather than one edge straight from room to room.
  //
  // This is what makes the corridor graph usable for planning. With a direct
  // room-to-room edge the route is a straight line between two room centres,
  // which cuts the corner of every doorway it crosses -- the path then clips
  // the wall beside the door and the planner (correctly) rejects its own plan
  // as blocked. Routing through the saddle node puts the doorway itself in the
  // route, so the expanded path is centred in the opening.
  std::vector<CorridorEdge> through;
  for (Size q = 0; q < edges.size(); ++q) {
    const Index cell = edge_saddle_cell[q];
    if (cell < 0) continue;
    if (mc.graph_.nodes.size() >= opt.max_nodes) break;
    CriticalPoint cp;
    cp.cx = static_cast<int>(cell % w);
    cp.cy = static_cast<int>(cell / w);
    cp.p = grid.cellCenter(cp.cx, cp.cy);
    cp.h = hf[cell];
    cp.type = CriticalType::kSaddle;
    cp.persistence = edges[q].persistence;
    const int sid = static_cast<int>(mc.graph_.nodes.size());
    mc.graph_.nodes.push_back(cp);
    ++mc.stats_.kept_saddles;

    CorridorEdge e1 = edges[q];
    e1.b = sid;
    e1.saddle = sid;
    e1.length = (mc.graph_.nodes[static_cast<Size>(e1.a)].p - cp.p).norm();
    CorridorEdge e2 = edges[q];
    e2.a = sid;
    e2.b = edges[q].b;
    e2.saddle = sid;
    e2.length = (cp.p - mc.graph_.nodes[static_cast<Size>(edges[q].b)].p).norm();
    through.push_back(e1);
    through.push_back(e2);
  }
  for (const CorridorEdge& e : through) mc.graph_.edges.push_back(e);
  mc.stats_.corridors = mc.graph_.edges.size();

  // Minima, when asked for: a cell lower than every valid neighbour. They are
  // wall corners and dead-end tips -- diagnostic rather than navigational, so
  // they are added as isolated nodes with no corridors.
  if (opt.include_minima) {
    for (Size i = 0; i < n; ++i) {
      if (hf[i] == kNegInf || is_max[i]) continue;
      Index nb[4];
      const int cnt = neighbours(static_cast<Index>(i), nb);
      if (cnt == 0) continue;
      bool lowest = true;
      for (int t = 0; t < cnt; ++t)
        if (hf[nb[t]] < hf[i]) {
          lowest = false;
          break;
        }
      if (!lowest) continue;
      if (mc.graph_.nodes.size() >= opt.max_nodes) break;
      CriticalPoint cp;
      cp.cx = static_cast<int>(i % w);
      cp.cy = static_cast<int>(i / w);
      cp.p = grid.cellCenter(cp.cx, cp.cy);
      cp.h = hf[i];
      cp.type = CriticalType::kMinimum;
      cp.persistence = field_max - hf[i];
      mc.graph_.nodes.push_back(cp);
      ++mc.stats_.kept_minima;
    }
  }

  mc.graph_.adjacency.assign(mc.graph_.nodes.size(), {});
  for (Size i = 0; i < mc.graph_.edges.size(); ++i) {
    const CorridorEdge& e = mc.graph_.edges[i];
    mc.graph_.adjacency[static_cast<Size>(e.a)].push_back(static_cast<int>(i));
    mc.graph_.adjacency[static_cast<Size>(e.b)].push_back(static_cast<int>(i));
  }
  return mc;
}

MorseComplex MorseComplex::extract(const sim::GridMap& grid, const MorseOptions& opt) {
  std::vector<Scalar> clearance;
  clearanceTransform(grid, &clearance);
  std::vector<Scalar> smoothed;
  gaussianBlurField(grid, clearance, opt.smoothing_sigma, &smoothed);
  MorseComplex mc = analyseField(grid, smoothed, opt);
  mc.opt_ = opt;
  return mc;
}

Scalar MorseComplex::heightAt(const Vec2& p) const {
  if (height_.empty()) return 0;
  const Eigen::Vector2i c = grid_.toCell(p);
  if (!grid_.inBounds(c.x(), c.y())) return 0;
  const Scalar v = height_[static_cast<Size>(c.y()) * grid_.w + c.x()];
  return v == kNegInf ? 0 : v;
}

// ---------------------------------------------------------------------------
// CorridorGraph
// ---------------------------------------------------------------------------

int CorridorGraph::nearestNode(const Vec2& p, Scalar max_dist) const {
  int best = -1;
  Scalar best_d = max_dist * max_dist;
  for (Size i = 0; i < nodes.size(); ++i) {
    const Scalar d = (nodes[i].p - p).squaredNorm();
    if (d < best_d) {
      best_d = d;
      best = static_cast<int>(i);
    }
  }
  return best;
}

std::vector<int> CorridorGraph::nodesOfType(CriticalType t) const {
  std::vector<int> out;
  for (Size i = 0; i < nodes.size(); ++i)
    if (nodes[i].type == t) out.push_back(static_cast<int>(i));
  return out;
}

int CorridorGraph::widestNode() const {
  int best = -1;
  Scalar best_h = -1;
  for (Size i = 0; i < nodes.size(); ++i) {
    if (nodes[i].type != CriticalType::kMaximum) continue;
    if (nodes[i].h > best_h) {
      best_h = nodes[i].h;
      best = static_cast<int>(i);
    }
  }
  return best;
}

bool CorridorGraph::route(int start_node, int goal_node, std::vector<int>* route) const {
  route->clear();
  if (start_node < 0 || goal_node < 0) return false;
  if (static_cast<Size>(start_node) >= nodes.size()) return false;
  if (static_cast<Size>(goal_node) >= nodes.size()) return false;
  if (start_node == goal_node) {
    route->push_back(start_node);
    return true;
  }

  // A* with the straight-line distance as the heuristic. The graph is a tree
  // for a simply-connected domain, so the heuristic is admissible and the
  // search expands essentially the single connecting path; on a domain with
  // holes (an island obstacle) it is still admissible because the straight
  // line lower-bounds any corridor route.
  const Size n = nodes.size();
  std::vector<Scalar> g(n, std::numeric_limits<Scalar>::max());
  std::vector<int> parent(n, -1);
  using QItem = std::pair<Scalar, int>;
  std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> open;

  const Vec2 goal_p = nodes[static_cast<Size>(goal_node)].p;
  const auto heur = [&](int i) { return (nodes[static_cast<Size>(i)].p - goal_p).norm(); };

  g[static_cast<Size>(start_node)] = 0;
  open.push({heur(start_node), start_node});
  std::vector<std::uint8_t> closed(n, 0);

  while (!open.empty()) {
    const int cur = open.top().second;
    open.pop();
    if (closed[static_cast<Size>(cur)]) continue;
    closed[static_cast<Size>(cur)] = 1;
    if (cur == goal_node) {
      for (int i = cur; i >= 0; i = parent[static_cast<Size>(i)]) route->push_back(i);
      std::reverse(route->begin(), route->end());
      return true;
    }
    for (int ei : adjacency[static_cast<Size>(cur)]) {
      const CorridorEdge& e = edges[static_cast<Size>(ei)];
      const int other = (e.a == cur) ? e.b : e.a;
      // Cost is the corridor length divided by its bottleneck clearance: a
      // tight doorway is more expensive than a wide one of the same length,
      // which is what makes the planner prefer the roomy route when two exist.
      const Scalar cost = e.length / std::max<Scalar>(0.05, e.bottleneck);
      const Scalar cand = g[static_cast<Size>(cur)] + cost;
      if (cand < g[static_cast<Size>(other)]) {
        g[static_cast<Size>(other)] = cand;
        parent[static_cast<Size>(other)] = cur;
        open.push({cand + heur(other), other});
      }
    }
  }
  return false;
}

std::string CorridorGraph::toDot() const {
  std::ostringstream os;
  os << "graph corridors {\n";
  for (Size i = 0; i < nodes.size(); ++i) {
    os << "  n" << i << " [label=\"" << toString(nodes[i].type) << " p="
       << nodes[i].persistence << "\"];\n";
  }
  for (Size i = 0; i < edges.size(); ++i) {
    const CorridorEdge& e = edges[i];
    os << "  n" << e.a << " -- n" << e.b << " [label=\"w=" << e.bottleneck << "\"];\n";
  }
  os << "}\n";
  return os.str();
}

std::string MorseComplex::toAscii(int max_cols) const {
  if (grid_.w <= 0 || grid_.h <= 0) return "";
  const int step = std::max(1, (grid_.w + max_cols - 1) / max_cols);
  std::vector<std::string> rows;
  for (int y = 0; y < grid_.h; y += step * 2) {
    std::string line;
    for (int x = 0; x < grid_.w; x += step) {
      char ch = '.';
      const sim::Cell c = grid_.at(x, y);
      if (c == sim::Cell::kOccupied) ch = '#';
      else if (c == sim::Cell::kUnknown) ch = ' ';
      line.push_back(ch);
    }
    rows.push_back(line);
  }
  // Overlay the critical points, largest first so a room label is not hidden by
  // a later doorway mark.
  for (const CriticalPoint& cp : graph_.nodes) {
    const int ry = cp.cy / (step * 2);
    const int rx = cp.cx / step;
    if (ry < 0 || ry >= static_cast<int>(rows.size())) continue;
    if (rx < 0 || rx >= static_cast<int>(rows[static_cast<Size>(ry)].size())) continue;
    rows[static_cast<Size>(ry)][static_cast<Size>(rx)] =
        cp.type == CriticalType::kMaximum ? 'O'
                                          : (cp.type == CriticalType::kSaddle ? 'x' : '.');
  }
  std::string out;
  for (const std::string& r : rows) {
    out += r;
    out += '\n';
  }
  return out;
}

}  // namespace topology
}  // namespace sc
