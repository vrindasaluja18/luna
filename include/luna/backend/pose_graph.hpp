// LUNA -- pose-graph optimizer over SE(2).
//
// The backend's job is to turn a chain of noisy relative-pose measurements
// (odometry) plus a handful of long-range ones (loop closures) into a globally
// consistent trajectory. It is a least-squares problem on a graph:
//
//   nodes       = poses x_i, each carrying 3 unknowns
//   edges       = constraints (i, j, z_ij, Omega_ij) meaning
//                     x_j = x_i * Exp(eps),  eps ~ N(0, Omega^{-1})
//
// We minimise the Mahalanobis residual over all edges,
//
//   F(x) = sum_ij  r_ij(x)^T Omega_ij r_ij(x),
//   r_ij(x) = Log( Z_ij^{-1} * (x_i^{-1} x_j) )        (with Z_ij = Exp(z_ij))
//
// by Gauss-Newton with Levenberg-Marquardt damping and a Huber robust kernel.
// The normal equations H dx = -b are assembled block-sparsely and solved with
// the sparse Cholesky in core/sparse.hpp.
//
// Why SE(2) and not SE(3): the simulator's ground truth is planar and every
// measurement it emits is planar. Six extra unobservable degrees of freedom
// would make the reported ATE meaningless.
#pragma once

#include <string>
#include <vector>

#include "luna/core/se2.hpp"
#include "luna/core/se3.hpp"
#include "luna/core/sparse.hpp"

namespace luna {
namespace backend {

/// One relative-pose measurement: pose `j` as seen from pose `i`.
///
/// `z` is the measured relative transform in tangent form (the [rho; phi]
/// returned by se2Log). `info` is the 3x3 information matrix, i.e. the inverse
/// of the measurement covariance. It must be symmetric positive definite; a
/// singular information matrix would make the corresponding edge meaningless,
/// so it is validated on insertion.
struct PoseEdge {
  Index i = 0;
  Index j = 0;
  Vec3 z = Vec3::Zero();
  Mat3 info = Mat3::Identity();
  bool robust = true;      ///< apply the Huber kernel to this edge
  Scalar chi2 = 0;         ///< residual at the last evaluation (diagnostics)
};

struct GraphOptions {
  int max_iterations = 30;
  Scalar lambda0 = 1e-4;          ///< initial LM damping
  Scalar lambda_up = 10.0;        ///< damping multiplier on a rejected step
  Scalar lambda_down = 0.1;       ///< damping multiplier on an accepted step
  Scalar huber_delta = 3.0;       ///< Huber threshold in sqrt(chi2) units
  Scalar convergence_grad = 1e-8;  ///< stop when |b|_inf falls below this
  Scalar convergence_step = 1e-8;  ///< stop when |dx|_inf falls below this
  bool fix_first = true;          ///< gauge freedom: pin node 0
  /// Cap on iterations that actually improve the cost. A pose graph that needs
  /// more than a few dozen accepted steps is diverging, not converging.
  int max_accepted = 20;
};

struct GraphResult {
  bool ok = false;
  int iterations = 0;
  Scalar chi2_initial = 0;
  Scalar chi2_final = 0;
  Scalar grad_norm = 0;
  Scalar solve_ms = 0;          ///< total time in assemble+solve
  Scalar last_solve_ms = 0;     ///< last iteration only
  int nonzeros = 0;
  /// Per-edge chi2 at the final estimate. An edge far out in the tail is the
  /// signature of a false loop closure, which is what the caller inspects to
  /// decide whether to roll the graph back.
  std::vector<Scalar> edge_chi2;
  /// Rejected loop closures, in descending chi2 order.
  std::vector<Index> outliers;
};

/// Pose graph. Node poses are owned here; edges are appended as measurements
/// arrive. Node 0 is created by the constructor and is the gauge anchor.
class PoseGraph {
 public:
  PoseGraph() = default;

  /// Append a node with the given initial guess and return its index.
  Index addNode(const Pose2& p);

  /// Append an edge. Rejects (and reports) an edge whose information matrix is
  /// not positive definite, or that references an unknown node, because a bad
  /// edge silently corrupts every pose in the graph rather than just its own.
  bool addEdge(const PoseEdge& e);

  Size nodes() const { return static_cast<Size>(poses_.size()); }
  Size edgeCount() const { return static_cast<Size>(edges_.size()); }
  const Pose2& pose(Index i) const { return poses_[i]; }
  std::vector<Pose2>& poses() { return poses_; }
  const std::vector<PoseEdge>& edges() const { return edges_; }

  /// Optimise in place. Returns false when the problem is structurally broken
  /// (fewer than two nodes, or a failed factorisation).
  GraphResult optimize(const GraphOptions& opt = GraphOptions());

  /// Chi2 of every edge at the current estimate, without optimising.
  std::vector<Scalar> residuals() const;

  /// Largest chi2 among edges and the index of that edge. This is the
  /// loop-closure gate: a single edge far outside the rest of the distribution
  /// means a geometrically plausible but wrong match.
  Scalar maxResidual(Index* edge_index) const;

  /// Drop the edge at `index`, used to recover from a rejected closure. The
  /// edge list is rebuilt, so indices shift; the caller must re-read edges.
  void removeEdge(Index index);

  /// Re-seed node `i` from its current neighbours, used to recover a pose that
  /// a bad closure dragged far off. Returns false when the node has no usable
  /// neighbours.
  bool reinitialiseNode(Index i);

  /// Marginal covariance of node `i` at the current estimate, by the Laplace
  /// approximation: Sigma_i = (H^{-1})_{ii}, where H is the information matrix
  /// of the linearised problem.
  ///
  /// The inverse is never formed. Every other variable is eliminated with a
  /// Schur complement first, which leaves a 3x3 block to invert; forming the
  /// full H^{-1} of a graph with a thousand nodes is 9 million entries of which
  /// three are wanted. The result is the *marginal* covariance and not the
  /// diagonal of H^{-1} read off directly, and the difference matters: the
  /// diagonal of a sparse information matrix is local, while a trajectory's
  /// uncertainty accumulates along the chain and a marginal carries that
  /// accumulation. The controller's chance constraint depends on the marginal.
  ///
  /// Returns false when the system is singular or the node is out of range.
  bool poseCovariance(Index i, Mat3* cov) const;

  /// Covariance of the *relative* pose of node `j` as seen from node `i`, i.e.
  /// of the transform x_i^{-1} x_j.
  ///
  /// This is a different quantity from `poseCovariance`, and for a controller it
  /// is the right one. A marginal covariance accumulates drift over the whole
  /// trajectory: the longer the run, the larger it grows, without bound. A
  /// chance constraint built on it would make the vehicle slower and slower and
  /// eventually stop, purely because it had been driving for a long time --
  /// which is not what "is this gap wide enough for me, now" is asking.
  ///
  /// The relative covariance measures how well the vehicle knows its pose
  /// *relative to where it recently was*, and that stays bounded because
  /// consecutive observations are correlated. It is what the vehicle's own
  /// recent odometry constrains, and it is what a local planning horizon can
  /// act on.
  ///
  /// Computed to first order as
  ///   Sigma_rel = Sigma_jj + Ad Sigma_ii Ad^T - Sigma_ji Ad^T - Ad Sigma_ij,
  ///   Ad = Ad((x_i^{-1} x_j)^{-1}),
  /// which is the BCH expansion of the relative-pose error truncated after the
  /// linear term. It is exact when the two nodes are joined by a single
  /// constraint, and a good approximation while the heading change across the
  /// window is small; the neglected terms are second order in that heading
  /// change times the per-step noise, so the estimate errs on the conservative
  /// side for a window of a few keyframes and degrades for a very long one.
  /// Keep the window short -- a few metres -- which is also what a local
  /// planning horizon needs.
  bool relativeCovariance(Index i, Index j, Mat3* cov) const;

  /// Total information-matrix dimension (3 per node).
  Size variables() const { return 3 * poses_.size(); }

  /// Dense information matrix at the current estimate. Test/diagnostic use
  /// only: it is O(n^2) and must never be called on a production graph.
  MatX infoMatrixForTest() const;

 private:
  /// Block (a, b) of the information matrix's inverse, 3x3, in the nodes'
  /// tangent spaces. The building block of every covariance query.
  bool covarianceBlock(Index a, Index b, Mat3* blk) const;

  /// Evaluate the total robust cost.
  Scalar evaluate(bool robust, const GraphOptions& opt) const;
  /// Assemble H and the right-hand side b into `sys`. `sys` owns both, so the
  /// solve that follows reads the same vector that was assembled -- passing a
  /// separate b is how a graph silently converges to nothing.
  void assemble(SparseSystem* sys, bool robust, const GraphOptions& opt) const;

  std::vector<Pose2> poses_;
  std::vector<PoseEdge> edges_;
};

}  // namespace backend
}  // namespace luna
