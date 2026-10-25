#include "slamcore/core/sparse.hpp"

#include <algorithm>
#include <chrono>

#include <Eigen/Cholesky>
#include <Eigen/LU>

namespace sc {

SparseSystem::SparseSystem(Size n_vars, Size nnz_hint) { reset(n_vars, nnz_hint); }

void SparseSystem::reset(Size n_vars, Size nnz_hint) {
  n_ = n_vars;
  triplets_.clear();
  triplets_.reserve(nnz_hint ? nnz_hint : 64);
  H_.resize(n_vars, n_vars);
  b_.setZero(n_vars);
  dx_.setZero(n_vars);
  factorised_ = false;
}

void SparseSystem::reserve(Size nnz_hint) {
  if (nnz_hint > triplets_.capacity()) triplets_.reserve(nnz_hint);
}

void SparseSystem::addH(Index i, Index j, const MatX& blk, bool symmetric) {
  const Index r0 = i, c0 = j;
  for (Index r = 0; r < blk.rows(); ++r) {
    for (Index c = 0; c < blk.cols(); ++c) {
      const Scalar v = blk(r, c);
      if (v == 0.0) continue;
      triplets_.emplace_back(r0 + r, c0 + c, v);
      if (symmetric && i != j) triplets_.emplace_back(c0 + c, r0 + r, v);
    }
  }
  factorised_ = false;
}

void SparseSystem::addB(Index i, const VecX& v) { b_.segment(i, v.size()) += v; }

void SparseSystem::addDamping(Scalar lambda) {
  for (Index i = 0; i < static_cast<Index>(n_); ++i) triplets_.emplace_back(i, i, lambda);
  factorised_ = false;
}

void SparseSystem::buildMatrix() {
  H_.resize(n_, n_);
  H_.setFromTriplets(triplets_.begin(), triplets_.end());
  H_.makeCompressed();
}

bool SparseSystem::factorise(bool use_ldlt) {
  const auto t0 = std::chrono::steady_clock::now();
  buildMatrix();
  use_ldlt_ = use_ldlt;
  bool ok = false;
  if (use_ldlt) {
    ldlt_.compute(H_);
    ok = ldlt_.info() == Eigen::Success && (ldlt_.vectorD().array().abs() > 0).all();
  } else {
    llt_.compute(H_);
    ok = llt_.info() == Eigen::Success;
  }
  const auto t1 = std::chrono::steady_clock::now();
  factor_ms_ = std::chrono::duration<Scalar, std::milli>(t1 - t0).count();
  factorised_ = ok;
  return ok;
}

bool SparseSystem::solve(VecX* dx) {
  if (!factorised_) return false;
  const auto t0 = std::chrono::steady_clock::now();
  if (use_ldlt_) {
    dx_ = ldlt_.solve(b_);
  } else {
    dx_ = llt_.solve(b_);
  }
  const auto t1 = std::chrono::steady_clock::now();
  solve_ms_ = std::chrono::duration<Scalar, std::milli>(t1 - t0).count();
  if (!dx_.allFinite()) return false;
  if (dx) *dx = dx_;
  return true;
}

MatX SparseSystem::denseH() const { return MatX(H_); }

MatX SparseSystem::diagonalBlock(Index i, Size block_size) const {
  return MatX(H_.block(i, i, block_size, block_size));
}

SchurResult schurComplement(const SpMat& H, const VecX& b, const std::vector<Index>& marginal,
                            const std::vector<Index>& keep) {
  const auto t0 = std::chrono::steady_clock::now();
  SchurResult out;
  const Size nm = marginal.size();
  const Size nk = keep.size();
  if (nm == 0 || nk == 0) return out;

  // Gather rows into a dense matrix: the marginalised block is bounded by the
  // lag window, so this stays small by construction.
  const auto gather = [&](const std::vector<Index>& rows, const std::vector<Index>& cols) {
    MatX out(static_cast<Index>(rows.size()), static_cast<Index>(cols.size()));
    for (Size r = 0; r < rows.size(); ++r)
      for (Size c = 0; c < cols.size(); ++c) out(static_cast<Index>(r), static_cast<Index>(c)) = H.coeff(rows[r], cols[c]);
    return out;
  };
  const auto gatherVec = [&](const std::vector<Index>& idx) {
    VecX out(static_cast<Index>(idx.size()));
    for (Size i = 0; i < idx.size(); ++i) out(static_cast<Index>(i)) = b(idx[i]);
    return out;
  };
  const MatX Hmm = gather(marginal, marginal);
  const MatX Hmr = gather(marginal, keep);
  const MatX Hrm = gather(keep, marginal);
  const MatX Hrr = gather(keep, keep);
  const VecX bm = gatherVec(marginal);
  const VecX br = gatherVec(keep);

  MatX Hmm_inv;
  if (!invertSymmetric(Hmm, &Hmm_inv)) return out;

  out.Hmm_inv_Hmr = Hmm_inv * Hmr;
  out.Hmm_inv_bm = Hmm_inv * bm;
  out.S = Hrr - Hrm * out.Hmm_inv_Hmr;
  out.b = br - Hrm * out.Hmm_inv_bm;
  out.S = 0.5 * (out.S + out.S.transpose()).eval();
  out.ok = true;
  const auto t1 = std::chrono::steady_clock::now();
  out.eliminate_ms = std::chrono::duration<Scalar, std::milli>(t1 - t0).count();
  return out;
}

bool solveDense(const MatX& A, const VecX& b, VecX* x) {
  Eigen::LDLT<MatX> ldlt(A);
  if (ldlt.info() != Eigen::Success) {
    const VecX sol = A.fullPivLu().solve(b);
    if (!sol.allFinite()) return false;
    if (x) *x = sol;
    return true;
  }
  const VecX sol = ldlt.solve(b);
  if (!sol.allFinite()) return false;
  if (x) *x = sol;
  return true;
}

bool invertSymmetric(const MatX& H, MatX* inv) {
  Eigen::LLT<MatX> llt(H);
  if (llt.info() != Eigen::Success) return false;
  const MatX I = MatX::Identity(H.rows(), H.cols());
  const MatX sol = llt.solve(I);
  if (!sol.allFinite()) return false;
  *inv = 0.5 * (sol + sol.transpose()).eval();
  return true;
}

bool sparseCholeskySolve(const SpMat& H, const VecX& b, VecX* x) {
  Eigen::SimplicialLDLT<SpMat> ldlt(H);
  if (ldlt.info() != Eigen::Success) return false;
  const VecX sol = ldlt.solve(b);
  if (!sol.allFinite()) return false;
  if (x) *x = sol;
  return true;
}

}  // namespace sc
