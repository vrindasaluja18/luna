#include "slamcore/backend/pose_graph.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include <Eigen/Cholesky>

namespace sc {
namespace backend {

namespace {

/// Huber weight for a residual of squared norm `chi2`.
///
/// The kernel is applied in sqrt(chi2) units, so `delta` is a Mahalanobis
/// distance rather than a squared one. That is what makes the threshold
/// comparable to a sigma level: a false loop closure sits many sigmas out and
/// is down-weighted to almost nothing, while an ordinary edge is untouched.
Scalar huberWeight(Scalar chi2, Scalar delta) {
  const Scalar r = std::sqrt(std::max(chi2, 0.0));
  if (r <= delta || r <= 0.0) return 1.0;
  return delta / r;
}

/// Is `info` symmetric positive definite?
///
/// A non-PD information matrix is not a measurement: it either encodes a
/// negative variance (a sign error) or is rank deficient (a constraint that
/// says nothing about one direction). Either one poisons the whole graph, so
/// such an edge is rejected at insertion instead of being silently used.
bool isPositiveDefinite(const Mat3& info) {
  if (!info.allFinite()) return false;
  Eigen::LLT<Mat3> llt(info);
  return llt.info() == Eigen::Success;
}

/// Residual and Jacobians for a single edge at the current node estimates.
///
/// Derivation of the Jacobians, because getting these wrong is silent: with
/// right perturbations x_i -> x_i Exp(dx_i) and x_j -> x_j Exp(dx_j), let
/// e = T_i^{-1} T_j and r = Log(Z^{-1} e) with Z = Exp(z). Then
///
///     dr / ddx_j =  J_r^{-1}(r)
///     dr / ddx_i = -J_r^{-1}(r) Ad_{Z^{-1}}
///
/// The second follows from conjugating the perturbation through the residual:
/// with r = Log(Z^{-1} T_i^{-1} T_j) and T_i -> T_i Exp(dx_i), the increment
/// enters as Exp(-Ad_{Z^{-1}} dx_i) inside the log. The Adjoint argument is
/// Z^{-1} = Exp(-z), NOT the measured transform e -- using e there is a real
/// error that still converges, because the linearisation is only used to pick
/// a step, but it biases every turning segment. Both Jacobians are checked
/// against central finite differences of the residual in the tests.
struct EdgeEval {
  Vec3 r = Vec3::Zero();
  Mat3 Ji = Mat3::Zero();
  Mat3 Jj = Mat3::Zero();
  Scalar chi2 = 0;
  bool ok = false;
};

EdgeEval evaluateEdge(const PoseEdge& e, const std::vector<Pose2>& poses, Scalar huber_delta,
                      bool robust) {
  EdgeEval out;
  const Mat3 Ti = pose2ToSE2(poses[e.i]);
  const Mat3 Tj = pose2ToSE2(poses[e.j]);
  const Mat3 Z = se2Exp(e.z);

  const Mat3 einv = se2Inverse(Ti) * Tj;
  const Vec3 r = se2Log(se2Inverse(Z) * einv);
  if (!r.allFinite()) return out;

  // Node j's perturbation lands on the right of the argument inside the log,
  // so it takes the right Jacobian inverse. Node i's lands on the left, so it
  // takes the *left* one. Using the right inverse for both is wrong by a term
  // of order phi*r and biases every turning segment.
  const Mat3 JjInv = se2RightJacobianMatInverse(r);
  const Mat3 JiInv = se2LeftJacobianMatInverse(r);

  out.r = r;
  out.Jj = JjInv;
  out.Ji = -JiInv * se2Adjoint(se2Inverse(Z));

  out.chi2 = r.dot(e.info * r);
  if (robust) {
    const Scalar w = huberWeight(out.chi2, huber_delta);
    out.chi2 *= w;
  }
  out.ok = true;
  return out;
}

}  // namespace

Index PoseGraph::addNode(const Pose2& p) {
  poses_.push_back(p);
  return static_cast<Index>(poses_.size()) - 1;
}

bool PoseGraph::addEdge(const PoseEdge& e) {
  if (e.i < 0 || e.j < 0) return false;
  if (e.i >= static_cast<Index>(poses_.size()) || e.j >= static_cast<Index>(poses_.size())) {
    return false;
  }
  if (e.i == e.j) return false;
  if (!isPositiveDefinite(e.info)) return false;
  edges_.push_back(e);
  return true;
}

void PoseGraph::removeEdge(Index index) {
  if (index < 0 || index >= static_cast<Index>(edges_.size())) return;
  edges_.erase(edges_.begin() + index);
}

Scalar PoseGraph::evaluate(bool robust, const GraphOptions& opt) const {
  Scalar cost = 0;
  for (const PoseEdge& e : edges_) {
    const EdgeEval ev = evaluateEdge(e, poses_, opt.huber_delta, robust && e.robust);
    if (ev.ok) cost += ev.chi2;
  }
  return cost;
}

void PoseGraph::assemble(SparseSystem* sys, bool robust, const GraphOptions& opt) const {
  for (const PoseEdge& e : edges_) {
    const EdgeEval ev = evaluateEdge(e, poses_, opt.huber_delta, robust && e.robust);
    if (!ev.ok) continue;

    Scalar w = 1.0;
    if (robust && e.robust) {
      const Scalar raw = ev.r.dot(e.info * ev.r);
      w = huberWeight(raw, opt.huber_delta);
    }
    const Mat3 Omega = w * e.info;
    const Mat3 JiO = ev.Ji.transpose() * Omega;
    const Mat3 JjO = ev.Jj.transpose() * Omega;

    sys->addH(3 * e.i, 3 * e.i, JiO * ev.Ji);
    sys->addH(3 * e.i, 3 * e.j, JiO * ev.Jj);
    sys->addH(3 * e.j, 3 * e.j, JjO * ev.Jj);
    // The gradient of the cost is 2 J^T Omega r; the normal equations solve
    // for the step, so the right-hand side is the negated gradient.
    sys->addB(3 * e.i, -JiO * ev.r);
    sys->addB(3 * e.j, -JjO * ev.r);
  }
}

GraphResult PoseGraph::optimize(const GraphOptions& opt) {
  GraphResult res;
  const Size n = poses_.size();
  // A local copy, because the robust threshold is annealed below. See the note
  // where huber_delta is updated.
  GraphOptions go = opt;
  if (n < 2 || edges_.empty()) {
    res.ok = false;
    return res;
  }

  const Size nv = n * 3;
  const Size nnz_hint = edges_.size() * 40 + nv;
  SparseSystem sys(nv, nnz_hint);

  res.chi2_initial = evaluate(true, go);
  Scalar cost = res.chi2_initial;
  res.chi2_final = cost;

  Scalar lambda = opt.lambda0;
  int accepted = 0;
  VecX dx;

  const auto t_start = std::chrono::steady_clock::now();

  for (int iter = 0; iter < opt.max_iterations; ++iter) {
    res.iterations = iter + 1;

    // Graduated robustness. A genuine loop closure starts life with a large
    // residual -- that residual *is* the accumulated drift it exists to remove.
    // Applying the Huber kernel from the first iteration therefore
    // down-weights the very measurement that should dominate, and the graph
    // settles in a local minimum only a fraction of the way to the fix. The
    // threshold is annealed from effectively no suppression down to the
    // configured value, so the closure acts like a plain least-squares term
    // while the gross correction is made and only then guards against outliers.
    go.huber_delta =
        std::max(opt.huber_delta, 1e3 * std::pow(0.1, static_cast<double>(iter)));
    sys.reset(nv, nnz_hint);
    assemble(&sys, true, go);

    // The assembled right-hand side is the negated gradient, so its infinity
    // norm is the stopping test: zero means a stationary point.
    res.grad_norm = sys.rhs().cwiseAbs().maxCoeff();
    if (res.grad_norm < opt.convergence_grad) {
      res.ok = true;
      break;
    }

    // Gauge: pin node 0. Implemented as a stiff prior on its own variables
    // rather than by deleting rows, which keeps H symmetric positive definite
    // and leaves the sparsity pattern uniform across iterations.
    if (opt.fix_first) {
      sys.addH(0, 0, 1e8 * Mat3::Identity(), false);
    }

    sys.addDamping(lambda);

    const auto t_iter = std::chrono::steady_clock::now();
    if (!sys.factorise(true)) {
      // H is not positive definite even with damping. Usually transient, so
      // raise the damping and retry rather than aborting.
      lambda *= opt.lambda_up;
      if (lambda > 1e12) break;
      continue;
    }
    if (!sys.solve(&dx)) break;
    res.last_solve_ms =
        std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t_iter)
            .count();
    res.solve_ms += res.last_solve_ms;
    res.nonzeros = static_cast<int>(sys.nonzeros());

    if (!dx.allFinite()) {
      lambda *= opt.lambda_up;
      if (lambda > 1e12) break;
      continue;
    }
    if (dx.cwiseAbs().maxCoeff() < opt.convergence_step && lambda <= opt.lambda0) {
      res.ok = true;
      break;
    }

    // Tentative step, using the right-perturbation update the Jacobians were
    // derived for. Node 0 is skipped when the gauge is pinned.
    const std::vector<Pose2> saved = poses_;
    const Index first = opt.fix_first ? 1 : 0;
    for (Index i = first; i < static_cast<Index>(n); ++i) {
      poses_[i] = pose2Update(poses_[i], dx.segment<3>(3 * i));
    }

    const Scalar trial = evaluate(true, go);
    if (trial < cost) {
      cost = trial;
      ++accepted;
      lambda = std::max(lambda * opt.lambda_down, 1e-10);
      if (accepted >= opt.max_accepted) {
        res.ok = true;
        break;
      }
    } else {
      poses_ = saved;
      lambda *= opt.lambda_up;
      if (lambda > 1e12) {
        res.ok = true;
        break;
      }
    }
  }

  res.solve_ms =
      std::chrono::duration<Scalar, std::milli>(std::chrono::steady_clock::now() - t_start).count();
  res.chi2_final = cost;
  if (!res.ok) res.ok = std::isfinite(cost);
  res.edge_chi2 = residuals();

  std::vector<Index> order(res.edge_chi2.size());
  for (Size i = 0; i < order.size(); ++i) order[i] = static_cast<Index>(i);
  std::sort(order.begin(), order.end(),
            [&](Index a, Index b2) { return res.edge_chi2[a] > res.edge_chi2[b2]; });
  const Size n_out = std::min<Size>(order.size(), 16);
  res.outliers.assign(order.begin(), order.begin() + n_out);
  return res;
}


namespace {
/// Adjoint of an SE(2) transform, acting on the tangent vector [x, y, theta].
///
///   Ad_T = [ R(theta)   J t ]      J = [ 0  1 ]
///          [    0       1 ]              [-1  0 ]
///
/// Used to express a covariance block that lives in one node's tangent space in
/// another's, which is what makes a relative covariance meaningful rather than a
/// sum of unrelated blocks. The rotation part is a real rotation and not the
/// identity: dropping it -- which is tempting, because for a vehicle whose
/// heading is near zero the two agree -- makes the formula wrong by the heading
/// change over the window, and a heading change is exactly what a vehicle does
/// while driving.
Mat3 se2Adjoint(const Pose2& p) {
  Mat3 ad = Mat3::Zero();
  const Scalar c = std::cos(p.th), s = std::sin(p.th);
  ad(0, 0) = c;
  ad(0, 1) = -s;
  ad(1, 0) = s;
  ad(1, 1) = c;
  ad(0, 2) = p.y;
  ad(1, 2) = -p.x;
  ad(2, 2) = 1.0;
  return ad;
}
}  // namespace

bool PoseGraph::poseCovariance(Index i, Mat3* cov) const {
  const Size n = poses_.size();
  if (i < 0 || static_cast<Size>(i) >= n || cov == nullptr) return false;
  Mat3 blk;
  if (!covarianceBlock(i, i, &blk)) return false;
  *cov = blk;
  return true;
}

MatX PoseGraph::infoMatrixForTest() const {
  const Size n = poses_.size();
  SparseSystem sys(3 * n, edges_.size() * 40 + 3 * n);
  const GraphOptions go;
  assemble(&sys, false, go);
  sys.factorise(false);
  return sys.denseH();
}

bool PoseGraph::covarianceBlock(Index a, Index b, Mat3* blk) const {
  const Size n = poses_.size();
  if (a < 0 || b < 0 || static_cast<Size>(a) >= n || static_cast<Size>(b) >= n) return false;

  // Assemble the information matrix at the current estimate. The robust kernel
  // is deliberately *off*: the Laplace approximation is a statement about a
  // Gaussian posterior, and a Huber weight is a function of the residual, so
  // including it would make the "covariance" depend on the very errors it is
  // meant to describe. The kernel's job is to keep an outlier from moving the
  // estimate, not to redefine the uncertainty of the estimate it produced.
  SparseSystem sys(3 * n, edges_.size() * 40 + 3 * n);
  const GraphOptions go;
  assemble(&sys, /*robust=*/false, go);
  // Gauge. The information matrix of a pose graph is singular along the global
  // translation: shifting every pose costs nothing. Pin node 0 with a strong
  // prior, which is what the optimiser does and what makes the trajectory
  // meaningful only *relative to* node 0. Node 0's own covariance comes out
  // near zero, which is the correct statement of that convention rather than an
  // artefact.
  // `addH` places a whole *block* starting at (i, j), so the pin must be a 1x1
  // block. Handing it a 3x3 identity here adds the same block three times at
  // three overlapping origins: the corner entries land on top of each other and
  // the diagonal ends up 1x, 2x, 3x the intended value, which quietly
  // mis-weights the gauge and every covariance derived from it.
  MatX pin(1, 1);
  pin(0, 0) = 1e6;
  for (Size k = 0; k < 3; ++k) {
    sys.addH(static_cast<Index>(k), static_cast<Index>(k), pin, /*symmetric=*/false);
  }
  if (!sys.factorise(true)) return false;

  // Extract a 3x3 block of H^{-1} without forming H^{-1}. Each column of H^{-1}
  // is one sparse solve against a unit vector, reusing the same factorisation,
  // so a block costs three triangular solves.
  for (Size c = 0; c < 3; ++c) {
    VecX e = VecX::Zero(3 * n);
    e(static_cast<Index>(3 * static_cast<Size>(b) + c)) = 1.0;
    if (!sparseCholeskySolve(sys.matrix(), e, &e)) return false;
    for (Size r = 0; r < 3; ++r) {
      (*blk)(static_cast<Index>(r), static_cast<Index>(c)) =
          e(static_cast<Index>(3 * static_cast<Size>(a) + r));
    }
  }
  return true;
}

bool PoseGraph::relativeCovariance(Index i, Index j, Mat3* cov) const {
  Mat3 sii, sjj, sij, sji;
  if (!covarianceBlock(i, i, &sii)) return false;
  if (!covarianceBlock(j, j, &sjj)) return false;
  if (!covarianceBlock(i, j, &sij)) return false;
  if (!covarianceBlock(j, i, &sji)) return false;

  const Pose2 rel = pose2Compose(pose2Inverse(poses_[i]), poses_[j]);
  const Mat3 ad = se2Adjoint(pose2Inverse(rel));
  // Right-perturbation error propagation: xi_ij = xi_j - Ad xi_i, so
  //   Sigma_rel = Sigma_jj + Ad Sigma_ii Ad^T - Sigma_ji Ad^T - Ad Sigma_ij.
  Mat3 r = sjj + ad * sii * ad.transpose() - sji * ad.transpose() - ad * sij;
  r = 0.5 * (r + r.transpose());
  // A covariance must have a non-negative diagonal. The solves are accurate up
  // to the factorisation's conditioning, and a negative variance would make
  // every downstream chance constraint nonsense, so it is rejected here rather
  // than propagated.
  for (Size k = 0; k < 3; ++k) {
    if (!std::isfinite(r(static_cast<Index>(k), static_cast<Index>(k))) ||
        r(static_cast<Index>(k), static_cast<Index>(k)) < 0) {
      return false;
    }
  }
  *cov = r;
  return true;
}

std::vector<Scalar> PoseGraph::residuals() const {
  std::vector<Scalar> out;
  out.reserve(edges_.size());
  const GraphOptions opt;
  for (const PoseEdge& e : edges_) {
    const EdgeEval ev = evaluateEdge(e, poses_, opt.huber_delta, false);
    out.push_back(ev.ok ? ev.r.dot(e.info * ev.r) : kInf);
  }
  return out;
}

Scalar PoseGraph::maxResidual(Index* edge_index) const {
  const std::vector<Scalar> r = residuals();
  Scalar best = -1;
  Index best_i = -1;
  for (Size i = 0; i < r.size(); ++i) {
    if (r[i] > best) {
      best = r[i];
      best_i = static_cast<Index>(i);
    }
  }
  if (edge_index) *edge_index = best_i;
  return best_i < 0 ? 0.0 : best;
}

bool PoseGraph::reinitialiseNode(Index i) {
  if (i < 0 || i >= static_cast<Index>(poses_.size())) return false;
  Scalar sx = 0, sy = 0, cx = 0, sn = 0;
  int count = 0;
  for (const PoseEdge& e : edges_) {
    // z encodes x_j = x_i Exp(z). Solve it for whichever endpoint is `i`.
    if (e.j == i) {
      const Pose2 cand = matToPose2(pose2ToSE2(poses_[e.i]) * se2Exp(e.z));
      sx += cand.x; sy += cand.y;
      cx += std::cos(cand.th); sn += std::sin(cand.th);
      ++count;
    } else if (e.i == i) {
      const Pose2 cand = matToPose2(pose2ToSE2(poses_[e.j]) * se2Inverse(se2Exp(e.z)));
      sx += cand.x; sy += cand.y;
      cx += std::cos(cand.th); sn += std::sin(cand.th);
      ++count;
    }
  }
  if (count == 0) return false;
  const Scalar inv = 1.0 / static_cast<Scalar>(count);
  // Average yaw through its sine and cosine rather than the raw angle, so a
  // cluster straddling +/-pi averages the short way round instead of to 0.
  poses_[i] = Pose2{sx * inv, sy * inv, std::atan2(sn, cx)};
  return true;
}

}  // namespace backend
}  // namespace sc
