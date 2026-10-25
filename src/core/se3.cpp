#include "slamcore/core/se3.hpp"

namespace sc {
namespace {

/// Coefficients (A, B, C) with
///   A = sin(t)/t, B = (1-cos t)/t^2, C = (t-sin t)/t^3
/// evaluated with Taylor series below 1e-4 to avoid catastrophic cancellation.
void rodriguesCoeffs(Scalar t, Scalar* A, Scalar* B, Scalar* C) {
  if (t < 1e-4) {
    const Scalar t2 = t * t;
    *A = 1.0 - t2 / 6.0 + t2 * t2 / 120.0;
    *B = 0.5 - t2 / 24.0 + t2 * t2 / 720.0;
    *C = 1.0 / 6.0 - t2 / 120.0 + t2 * t2 / 5040.0;
  } else {
    const Scalar s = std::sin(t), c = std::cos(t);
    *A = s / t;
    *B = (1.0 - c) / (t * t);
    *C = (t - s) / (t * t * t);
  }
}

}  // namespace

Mat3 skew(const Vec3& v) {
  Mat3 m;
  m << 0.0, -v.z(), v.y(),
       v.z(), 0.0, -v.x(),
       -v.y(), v.x(), 0.0;
  return m;
}

Mat3 so3Exp(const Vec3& phi) {
  const Scalar t = phi.norm();
  if (t < 1e-12) {
    // Exp is I + [phi]_x to first order; second order term vanishes.
    return Mat3::Identity() + skew(phi);
  }
  Scalar A, B, C;
  rodriguesCoeffs(t, &A, &B, &C);
  const Mat3 K = skew(phi);
  return Mat3::Identity() + A * K + B * (K * K);
}

Vec3 so3Log(const Mat3& R) {
  const Scalar c = clamp((R.trace() - 1.0) * 0.5, -1.0, 1.0);
  const Scalar theta = std::acos(c);
  if (theta < 1e-8) {
    // Small-angle: use the antisymmetric part directly.
    return Vec3(0.5 * (R(2, 1) - R(1, 2)), 0.5 * (R(0, 2) - R(2, 0)), 0.5 * (R(1, 0) - R(0, 1)));
  }
  if (theta > kPi - 1e-6) {
    // Near pi: extract the axis from the largest diagonal element of R + I.
    Vec3 axis;
    const Mat3 A = (R + Mat3::Identity()) * 0.5;
    const int i = (A(0, 0) > A(1, 1)) ? ((A(0, 0) > A(2, 2)) ? 0 : 2) : ((A(1, 1) > A(2, 2)) ? 1 : 2);
    axis[i] = std::sqrt(std::max(A(i, i), 0.0));
    if (axis[i] > 1e-8) {
      const int j = (i + 1) % 3, k = (i + 2) % 3;
      axis[j] = A(j, i) / axis[i];
      axis[k] = A(k, i) / axis[i];
    }
    const Scalar n = axis.norm();
    if (n > 1e-12) axis /= n;
    return axis * theta;
  }
  const Scalar k = theta / (2.0 * std::sin(theta));
  return Vec3(k * (R(2, 1) - R(1, 2)), k * (R(0, 2) - R(2, 0)), k * (R(1, 0) - R(0, 1)));
}

Mat3 so3RightJacobian(const Vec3& phi) {
  const Scalar t = phi.norm();
  if (t < 1e-8) return Mat3::Identity() - 0.5 * skew(phi);
  Scalar A, B, C;
  rodriguesCoeffs(t, &A, &B, &C);
  const Mat3 K = skew(phi);
  // J_r = I - B [phi]_x + C [phi]_x^2
  return Mat3::Identity() - B * K + C * (K * K);
}

Mat3 so3RightJacobianInverse(const Vec3& phi) {
  const Scalar t = phi.norm();
  const Mat3 K = skew(phi);
  if (t < 1e-8) return Mat3::Identity() + 0.5 * K;
  const Scalar half = 0.5 * t;
  // J_r^{-1} = I + 0.5 [phi]_x + (1/t^2)(1 - (t/2) cot(t/2)) [phi]_x^2
  const Scalar cot = std::cos(half) / std::sin(half);
  const Scalar coef = (1.0 - half * cot) / (t * t);
  return Mat3::Identity() + 0.5 * K + coef * (K * K);
}

Mat4 se3Exp(const Vec6& xi) {
  const Vec3 rho = xi.head<3>();
  const Vec3 phi = xi.tail<3>();
  const Scalar t = phi.norm();
  Scalar A, B, C;
  rodriguesCoeffs(t, &A, &B, &C);
  const Mat3 K = skew(phi);
  const Mat3 R = Mat3::Identity() + A * K + B * (K * K);
  // V = I + B [phi]_x + C [phi]_x^2   (left Jacobian of SO(3))
  const Mat3 V = Mat3::Identity() + B * K + C * (K * K);
  return makeSE3(R, V * rho);
}

Vec6 se3Log(const Mat4& T) {
  const Mat3 R = rotationOf(T);
  const Vec3 t = translationOf(T);
  const Vec3 phi = so3Log(R);
  const Scalar tn = phi.norm();
  Scalar A, B, C;
  rodriguesCoeffs(tn, &A, &B, &C);
  const Mat3 K = skew(phi);
  const Mat3 V = Mat3::Identity() + B * K + C * (K * K);
  Vec6 xi;
  xi.head<3>() = V.inverse() * t;
  xi.tail<3>() = phi;
  return xi;
}

Mat6 se3Adjoint(const Mat4& T) {
  const Mat3 R = rotationOf(T);
  const Vec3 t = translationOf(T);
  Mat6 Ad = Mat6::Zero();
  Ad.template block<3, 3>(0, 0) = R;
  Ad.template block<3, 3>(0, 3) = skew(t) * R;
  Ad.template block<3, 3>(3, 3) = R;
  return Ad;
}

Mat6 se3RightJacobian(const Vec6& xi) {
  const Vec3 phi = xi.tail<3>();
  const Vec3 rho = xi.head<3>();
  const Scalar t = phi.norm();
  const Mat3 K = skew(phi);
  Scalar A, B, C;
  rodriguesCoeffs(t, &A, &B, &C);
  // Q = 0.5 [rho]_x + C/2 ([phi]_x[rho]_x + [rho]_x[phi]_x + [phi]_x[rho]_x[phi]_x)
  //     - (B - 3C)/6 ([phi]_x^2[rho]_x + [rho]_x[phi]_x^2 - 3[phi]_x[rho]_x[phi]_x)
  //     - (1/2 - A/4 - 3C/2) ... higher-order translation terms are dropped
  //     below 1e-6 rad, where the block reduces to 0.5 [rho]_x.
  Mat6 J = Mat6::Identity();
  J.template block<3, 3>(3, 3) = so3RightJacobian(phi);
  Mat3 Q = Mat3::Zero();
  if (t > 1e-8) {
    const Mat3 K2 = K * K;
    const Mat3 rK = skew(rho);
    const Mat3 K3 = K2 * K;
    const Mat3 K4 = K2 * K2;
    const Scalar t2 = t * t;
    const Scalar t3 = t2 * t;
    const Scalar t4 = t2 * t2;
    const Scalar t5 = t4 * t;
    const Scalar c1 = (t3 - 6.0 * t + 6.0 * std::sin(t)) / (6.0 * t5);
    const Scalar c2 = (t2 - 2.0 + 2.0 * std::cos(t)) / (2.0 * t4);
    const Scalar c3 = (t3 - 6.0 * t + 6.0 * std::sin(t)) / (6.0 * t5);
    const Scalar c4 = (2.0 * t2 - 3.0 - 3.0 * t * std::sin(t) + 6.0 * std::cos(t)) / (6.0 * t4);
    const Scalar c5 = (t2 + 3.0 - 3.0 * t * std::sin(t) - 6.0 * std::cos(t)) / (6.0 * t5);
    Q = c1 * rK + c2 * (K * rK + rK * K + K * rK * K) + c3 * (K2 * rK + rK * K2) +
        c4 * (K3 * rK + rK * K3) + c5 * (K4 * rK + rK * K4);
  } else {
    Q = 0.5 * skew(rho);
  }
  J.template block<3, 3>(0, 0) = so3RightJacobian(phi);
  J.template block<3, 3>(0, 3) = Q;
  return J;
}

Mat6 se3RightJacobianInverse(const Vec6& xi) {
  const Vec6 neg = -xi;
  const Mat6 J = se3RightJacobian(neg);
  const Vec3 rho = xi.head<3>();
  // Use the analytic block structure: J^{-1} = [[Jr^{-1}, -Jr^{-1} Q Jr^{-1}],
  //                                           [0, Jr^{-1}]]
  Mat6 Jinv = Mat6::Zero();
  const Mat3 Jr = so3RightJacobian(xi.tail<3>());
  const Mat3 Jri = so3RightJacobianInverse(xi.tail<3>());
  Jinv.template block<3, 3>(0, 0) = Jri;
  Jinv.template block<3, 3>(3, 3) = Jri;
  Jinv.template block<3, 3>(0, 3) = -Jri * J.template block<3, 3>(0, 3) * Jri;
  (void)rho;
  (void)Jr;
  return Jinv;
}

Vec3 so3Residual(const Mat3& R_meas, const Mat3& R_pred) {
  return so3Log(R_meas * R_pred.transpose());
}

}  // namespace sc
