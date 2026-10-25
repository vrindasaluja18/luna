// SLAMCore-8G -- sparse-only linear algebra.
//
// The 8 GB / 4 core budget is what makes this module exist: every normal
// equation in the backend is assembled as a block-sparse system and solved by
// a sparse Cholesky factorisation. Dense N x N solves are banned by the code
// review checklist; the only dense paths allowed here are the small block
// operations (6x6, 3x3) that feed the triplets.
#pragma once

#include <vector>

#include <Eigen/Sparse>

#include "slamcore/core/se3.hpp"
#include "slamcore/core/types.hpp"

namespace sc {

using SpMat = Eigen::SparseMatrix<Scalar>;
using Triplet = Eigen::Triplet<Scalar>;

/// Block-sparse normal-equation accumulator.
///
/// Usage:
///     SparseSystem sys(n_vars);
///     sys.addH(i, j, JtWJ);   // symmetric: only i <= j must be added
///     sys.addB(i, JtWr);
///     sys.solve(&dx);
///
/// Storage is a triplet list, sorted once per solve. Reserving the triplet
/// capacity up front keeps assembly allocation-free in the optimize loop.
class SparseSystem {
 public:
  explicit SparseSystem(Size n_vars = 0, Size nnz_hint = 0);

  void reset(Size n_vars, Size nnz_hint = 0);
  void reserve(Size nnz_hint);

  /// Add a block at (row_block, col_block). If `symmetric` is true the
  /// transpose is added automatically (pass symmetric=false when both halves
  /// are supplied explicitly, as the Schur module does).
  void addH(Index i, Index j, const MatX& blk, bool symmetric = true);
  void addB(Index i, const VecX& v);
  /// Levenberg-Marquardt style damping on the diagonal.
  void addDamping(Scalar lambda);

  /// Assemble and factorise. Returns false if the factorisation failed, in
  /// which case the caller must fall back to the degraded path.
  bool factorise(bool use_ldlt = true);
  bool solve(VecX* dx);

  /// The assembled matrix (valid after factorise). Exposed for Schur
  /// marginalisation and for the Laplace covariance H^{-1}.
  const SpMat& matrix() const { return H_; }
  SpMat& matrix() { return H_; }
  const VecX& rhs() const { return b_; }

  Size variables() const { return n_; }
  Size nonzeros() const { return H_.nonZeros(); }
  Scalar factoriseMs() const { return factor_ms_; }
  Scalar solveMs() const { return solve_ms_; }

  /// Dense view of H, for tests and for the tiny systems in the topology
  /// module. Never call this on a production-sized graph.
  MatX denseH() const;

  /// Extract the diagonal block for variable `i` as a dense matrix.
  MatX diagonalBlock(Index i, Size block_size) const;

 private:
  void buildMatrix();

  Size n_ = 0;
  std::vector<Triplet> triplets_;
  SpMat H_;
  VecX b_;
  Eigen::SimplicialLDLT<SpMat> ldlt_;
  Eigen::SimplicialLLT<SpMat> llt_;
  VecX dx_;
  Scalar factor_ms_ = 0, solve_ms_ = 0;
  bool factorised_ = false;
  bool use_ldlt_ = true;
};

/// Schur complement of a block system.
///
/// Given
///     [ Hmm  Hmr ] [ dxm ]   [ bm ]
///     [ Hrm  Hrr ] [ dxr ] = [ br ]
/// eliminate the `marginal` variables:
///     S = Hrr - Hrm Hmm^{-1} Hmr,      br' = br - Hrm Hmm^{-1} bm
/// The reduced matrix S is returned as a dense matrix because the marginalised
/// block is the fixed-lag window itself (bounded by the lag length), which is
/// exactly the regime where dense S is faster than sparse S.
struct SchurResult {
  MatX S;          ///< reduced (dense) information matrix
  VecX b;          ///< reduced rhs
  MatX Hmm_inv_Hmr;  ///< cached for back-substitution
  VecX Hmm_inv_bm;
  bool ok = false;
  Scalar eliminate_ms = 0;
};

SchurResult schurComplement(const SpMat& H, const VecX& b, const std::vector<Index>& marginal,
                            const std::vector<Index>& keep);

/// Solve a small dense symmetric system with the same robustness contract as
/// SparseSystem::factorise (LDLT with a pivoted fallback).
bool solveDense(const MatX& A, const VecX& b, VecX* x);

/// Cholesky-based symmetric inverse for covariance extraction (Laplace
/// approximation H^{-1}). Returns false when H is not positive definite.
bool invertSymmetric(const MatX& H, MatX* inv);

/// Sparse Cholesky solve of a raw matrix, used by the fixed-lag smoother when
/// it owns the assembly itself.
bool sparseCholeskySolve(const SpMat& H, const VecX& b, VecX* x);

}  // namespace sc
