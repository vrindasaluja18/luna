#include "luna/core/qp.hpp"

#include <chrono>

namespace luna {

const char* toString(QpStatus s) {
  switch (s) {
    case QpStatus::kSolved: return "solved";
    case QpStatus::kMaxIter: return "max_iter";
    case QpStatus::kPrimalInfeasible: return "primal_infeasible";
    case QpStatus::kDualInfeasible: return "dual_infeasible";
    case QpStatus::kInvalid: return "invalid";
  }
  return "?";
}

void QpSolver::reserve(Size n, Size m) {
  n_ = n;
  m_ = m;
  P_.setZero(n, n);
  A_.setZero(m, n);
  At_.setZero(n, m);
  AtA_.setZero(n, n);
  KKT_.setZero(n, n);
  q_.setZero(n);
  l_.setZero(m);
  u_.setZero(m);
  x_.setZero(n);
  z_.setZero(m);
  z_prev_.setZero(m);
  y_.setZero(m);
  tmp_m_.setZero(m);
  d_.setZero(n);
  factorised_ = false;
  cached_rho_ = -1.0;
}

void QpSolver::setProblem(const MatX& P, const VecX& q, const MatX& A, const VecX& l, const VecX& u) {
  const Size n = static_cast<Size>(q.size());
  const Size m = static_cast<Size>(A.rows());
  if (n != n_ || m != m_) reserve(n, m);
  P_.setZero();
  P_.template triangularView<Eigen::Upper>() = P.template triangularView<Eigen::Upper>();
  P_ = P_.template selfadjointView<Eigen::Upper>();
  q_ = q;
  if (m_ > 0) {
    A_ = A;
    l_ = l;
    u_ = u;
    At_ = A_.transpose();
    AtA_.noalias() = At_ * A_;
  }
  factorised_ = false;
}

void QpSolver::setupFactorisation(Scalar rho) {
  KKT_.noalias() = P_;
  KKT_.diagonal().array() += sigma_;
  if (m_ > 0) KKT_.noalias() += rho * AtA_;
  KKT_ = 0.5 * (KKT_ + KKT_.transpose()).eval();
  llt_.compute(KKT_);
  cached_rho_ = rho;
  factorised_ = true;
  ++factorisations_;
}

QpResult QpSolver::solve(VecX* x, const QpOptions& opt) {
  return solve(x, VecX::Zero(n_), opt);
}

QpResult QpSolver::solve(VecX* x, const VecX& warm_start, const QpOptions& opt) {
  const auto t0 = std::chrono::steady_clock::now();
  QpResult res;
  if (n_ == 0) {
    res.status = QpStatus::kInvalid;
    return res;
  }
  x_ = warm_start;
  if (static_cast<Size>(x_.size()) != n_) x_.setZero();
  // The proximal term is taken from the options on every solve: it is part of
  // the problem's conditioning, not of the warm-start state, and a caller that
  // changes it must get a fresh factorisation.
  if (std::abs(opt.sigma - sigma_) > 1e-15) {
    sigma_ = opt.sigma;
    factorised_ = false;
  }
  rho_ = opt.rho;
  factorisations_ = 0;
  if (m_ > 0) {
    z_ = (A_ * x_).cwiseMax(l_).cwiseMin(u_);
    z_prev_ = z_;
    y_.setZero();
  }
  if (!factorised_ || cached_rho_ != rho_) setupFactorisation(rho_);

  const Scalar sqrt_n = std::sqrt(static_cast<Scalar>(n_));
  Size it = 0;
  for (; it < opt.max_iter; ++it) {
    admmIteration();
    evaluateResiduals();

    // OSQP termination: scaled residuals against the norm of the iterate.
    eps_pri_ = sqrt_n * opt.eps_abs + opt.eps_rel * std::max((A_ * x_).norm(), z_.norm());
    eps_dual_ = sqrt_n * opt.eps_abs + opt.eps_rel * (m_ > 0 ? (At_ * y_).norm() : 0.0);
    if (converged(eps_pri_, eps_dual_)) {
      ++it;
      break;
    }

    if (opt.adaptive_rho && m_ > 0 && it > 0 && (it % opt.rho_interval) == 0) {
      // Balance primal and dual residuals by rescaling rho; only rebuild the
      // factorisation when the change is material.
      const Scalar ratio = primal_res_ / std::max(dual_res_, 1e-12);
      if (ratio > 10.0 || ratio < 0.1) {
        const Scalar rho_new = clamp(rho_ * std::sqrt(ratio), opt.rho_min, opt.rho_max);
        if (std::abs(rho_new - rho_) > 1e-9 * std::max(1.0, rho_)) {
          const Scalar scale = rho_ / rho_new;
          y_ *= scale;
          rho_ = rho_new;
          setupFactorisation(rho_);
        }
      }
    }
  }
  if (opt.polish && m_ > 0) polish(opt.polish_iter);

  res.status = (it < opt.max_iter) ? QpStatus::kSolved : QpStatus::kMaxIter;
  res.iterations = it;
  res.primal_residual = primal_res_;
  res.dual_residual = dual_res_;
  res.objective = objective();
  res.rho_final = rho_;
  res.factorisations = factorisations_;
  const auto t1 = std::chrono::steady_clock::now();
  res.solve_ms = std::chrono::duration<Scalar, std::milli>(t1 - t0).count();
  if (x) *x = x_;
  return res;
}

void QpSolver::admmIteration() {
  // x-update: (P + sigma I + rho A'A) x = sigma*x - q + rho A'(z - y)
  d_.noalias() = sigma_ * x_ - q_;
  if (m_ > 0) d_.noalias() += rho_ * (At_ * (z_ - y_));
  x_ = llt_.solve(d_);
  if (m_ == 0) return;

  // z-update: Euclidean projection of (A x + y) onto the box [l, u].
  tmp_m_.noalias() = A_ * x_;
  tmp_m_ += y_;
  z_ = tmp_m_.cwiseMax(l_).cwiseMin(u_);

  // y-update: y += A x - z. Note this uses the *pre-update* y inside tmp_m_,
  // so the increment is (A x + y_old - z), i.e. y_new = y_old + Ax - z.
  y_ = tmp_m_ - z_;
}

void QpSolver::evaluateResiduals() {
  if (m_ == 0) {
    primal_res_ = 0;
    dual_res_ = 0;
    return;
  }
  primal_res_ = (A_ * x_ - z_).norm();
  dual_res_ = (rho_ * (At_ * (z_ - z_prev_))).norm();
  z_prev_ = z_;
}

bool QpSolver::converged(Scalar eps_pri, Scalar eps_dual) const {
  return primal_res_ <= eps_pri && dual_res_ <= eps_dual;
}

void QpSolver::polish(Size iterations) {
  // Equality-constrained polish on the active set: freeze the tight rows as
  // equalities and solve the KKT system exactly. This removes the ADMM tail
  // error so the control loop sees a crisp active set.
  for (Size k = 0; k < iterations; ++k) {
    std::vector<Index> act;
    act.reserve(m_);
    for (Size i = 0; i < m_; ++i) {
      const Scalar ax = A_.row(i).dot(x_);
      if (ax <= l_(i) + 1e-8 || ax >= u_(i) - 1e-8) act.push_back(static_cast<Index>(i));
    }
    if (act.empty()) break;
    const Size na = act.size();
    MatX K = MatX::Zero(n_ + na, n_ + na);
    VecX rhs = VecX::Zero(n_ + na);
    K.topLeftCorner(n_, n_) = P_;
    K.topLeftCorner(n_, n_).diagonal().array() += 1e-9;
    rhs.head(n_) = -q_;
    for (Size a = 0; a < na; ++a) {
      const Index i = act[a];
      K.block(0, n_ + a, n_, 1) = A_.row(i).transpose();
      K.block(n_ + a, 0, 1, n_) = A_.row(i);
      const Scalar ax = A_.row(i).dot(x_);
      rhs(n_ + a) = (ax <= l_(i) + 1e-8) ? l_(i) : u_(i);
    }
    const VecX sol = K.fullPivLu().solve(rhs);
    if (!sol.allFinite()) break;
    const Scalar step = (sol.head(n_) - x_).norm();
    x_ = sol.head(n_);
    if (step < 1e-11) break;
  }
}

Scalar QpSolver::objective() const {
  return 0.5 * x_.dot(P_ * x_) + q_.dot(x_);
}

}  // namespace luna
