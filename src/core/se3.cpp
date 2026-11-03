#include "luna/core/se3.hpp"

namespace luna {
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

namespace {

/// The little adjoint ad_xi of SE(3), in the module's tangent ordering
/// [rho; phi]:
///
///     ad = [ [phi^,  rho^] ]
///          [ [  0,   phi^] ]
///
/// It is the matrix of the map y -> xi x y on the Lie algebra. Everything the
/// right Jacobian is built from is a polynomial in this single matrix.
Mat6 se3LittleAdjoint(const Vec6& xi) {
  const Mat3 P = skew(xi.tail<3>());
  const Mat3 R = skew(xi.head<3>());
  Mat6 ad = Mat6::Zero();
  ad.template block<3, 3>(0, 0) = P;
  ad.template block<3, 3>(0, 3) = R;
  ad.template block<3, 3>(3, 3) = P;
  return ad;
}

}  // namespace

Mat6 se3RightJacobian(const Vec6& xi) {
  // The right Jacobian is the matrix series
  //
  //     J_r(xi) = sum_{n>=0} (-1)^n / (n+1)! ad_xi^n
  //
  // which is the SE(3) generalisation of J_r = I - B[phi]_x + C[phi]_x^2 and
  // reduces to it exactly when rho = 0 (then ad is block diagonal and the
  // translation block decouples). The series is the definition the rest of the
  // module relies on: it is what makes
  //
  //     Exp(xi + d) = Exp(xi) Exp(J_r(xi) d) + O(|d|^2)
  //
  // hold. The previous closed form here computed the (0,3) block from a set of
  // coefficients that did not reproduce the series -- the block was wrong by
  // O(1) at |phi| ~ 0.5 -- so the function failed its own defining identity
  // while every test that only touched SO(3) still passed. Building it from the
  // series is both exact and checkable against that identity.
  //
  // The series is summed to a fixed 20 terms. Because the rotational part
  // satisfies ad^3 = -|phi|^2 ad, each term is bounded by |phi|^n / (n+1)!, so
  // the truncation is far below double precision for every |phi| <= pi and
  // every finite rho. The bound is fixed rather than adaptive so the result is
  // bit-reproducible, which the seeded replay depends on.
  const Mat6 ad = se3LittleAdjoint(xi);
  Mat6 J = Mat6::Identity();
  Mat6 term = Mat6::Identity();   // ad^n
  Scalar fact = 1.0;              // (n+1)!
  for (int n = 1; n <= 20; ++n) {
    term = term * ad;
    fact *= static_cast<Scalar>(n + 1);
    J += ((n % 2) ? -1.0 : 1.0) / fact * term;
  }
  return J;
}

Mat6 se3RightJacobianInverse(const Vec6& xi) {
  // The inverse is taken of the Jacobian that is actually computed, so the two
  // cannot drift apart. A hand-derived block inverse (what this used to be) is
  // only correct if it agrees with the forward Jacobian in every block, and the
  // forward block was wrong -- the inverse was then wrong in a *different* way
  // and the product was off by O(10). A 6x6 inverse is a few hundred flops and
  // this is called once per edge, not per pixel.
  return se3RightJacobian(xi).inverse();
}

Vec3 so3Residual(const Mat3& R_meas, const Mat3& R_pred) {
  return so3Log(R_meas * R_pred.transpose());
}

}  // namespace luna
