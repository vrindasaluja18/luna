// LUNA -- convex QP solver (OSQP-style ADMM), pure CPU, no dependencies.
//
// Solves
//     min 0.5 x'Px + q'x   s.t.  l <= Ax <= u
// using the scaled ADMM splitting of OSQP (Stellato et al., 2020):
//
//     x^{k+1} = (P + sigma I + rho A'A)^-1 (sigma x^k - q + rho A'(z^k - y^k))
//     z^{k+1} = proj_[l,u](A x^{k+1} + y^k)
//     y^{k+1} = y^k + A x^{k+1} - z^{k+1}
//
// The KKT matrix depends only on (P, rho), so it is factorised once per
// problem and reused across every MPC tick: the per-iteration cost is then a
// triangular solve, which is what keeps a 25-step horizon inside the 20 Hz
// budget on a 4-core CPU. Warm starting across ticks keeps iterations in the
// low single digits in steady state.
//
// rho adaptation (the same rule OSQP uses) is enabled by default because a
// fixed rho is genuinely slow on badly scaled problems; the factorisation is
// rebuilt only when rho actually changes.
#pragma once

#include <vector>

#include "luna/core/se3.hpp"
#include "luna/core/types.hpp"

namespace luna {

enum class QpStatus { kSolved, kMaxIter, kPrimalInfeasible, kDualInfeasible, kInvalid };

const char* toString(QpStatus s);

struct QpOptions {
  Size max_iter = 4000;
  Scalar eps_abs = 1e-6;
  Scalar eps_rel = 1e-6;
  /// Proximal (Tikhonov) regularisation added to the KKT diagonal. It must be
  /// large enough that the KKT matrix stays well conditioned when P is only
  /// positive *semi*-definite, which is the normal case for an MPC whose state
  /// cost cannot see some input directions within the horizon. The default is
  /// fine for a strictly convex problem; a condensed MPC should raise it.
  Scalar sigma = 1e-6;
  Scalar rho = 0.1;                 ///< initial penalty
  bool adaptive_rho = true;         ///< OSQP-style residual balancing
  Scalar rho_min = 1e-6;
  Scalar rho_max = 1e6;
  Size rho_interval = 25;           ///< iterations between rho updates
  bool polish = true;               ///< final active-set equality polish
  Size polish_iter = 3;
};

struct QpResult {
  QpStatus status = QpStatus::kInvalid;
  Size iterations = 0;
  Scalar primal_residual = 0;
  Scalar dual_residual = 0;
  Scalar objective = 0;
  Scalar rho_final = 0;
  Size factorisations = 0;
  Scalar solve_ms = 0;
};

/// Dense convex QP. Problems seen here (MPC with n < 200) are far below the
/// size where sparse storage pays for its indirection on cache-resident data,
/// so dense row-major is the faster and simpler choice.
class QpSolver {
 public:
  QpSolver() = default;

  /// Pre-allocate for a problem of the given shape so that solve() is
  /// allocation-free after the first call with that shape.
  void reserve(Size n, Size m);

  /// Set the problem. P must be symmetric PSD; only the upper triangle is read.
  /// A may have zero rows when m == 0.
  void setProblem(const MatX& P, const VecX& q, const MatX& A, const VecX& l, const VecX& u);

  QpResult solve(VecX* x, const QpOptions& opt = QpOptions());
  QpResult solve(VecX* x, const VecX& warm_start, const QpOptions& opt = QpOptions());

  const VecX& primalSolution() const { return x_; }
  const VecX& dualSolution() const { return y_; }

 private:
  void setupFactorisation(Scalar rho);
  void admmIteration();
  void evaluateResiduals();
  bool converged(Scalar eps_pri, Scalar eps_dual) const;
  void polish(Size iterations);
  Scalar objective() const;

  Size n_ = 0, m_ = 0;
  MatX P_, A_, At_, AtA_, KKT_;
  VecX q_, l_, u_, x_, z_, z_prev_, y_, tmp_m_, d_;
  Eigen::LLT<MatX> llt_;
  bool factorised_ = false;
  Scalar cached_rho_ = -1.0;
  Scalar sigma_ = 1e-6;
  Scalar rho_ = 0.1;
  Scalar primal_res_ = 0, dual_res_ = 0;
  Scalar eps_pri_ = 0, eps_dual_ = 0;
  Size factorisations_ = 0;
};

}  // namespace luna
