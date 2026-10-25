#include "slamcore/core/se2.hpp"

namespace sc {

namespace {

/// Series expansion of the SE(2) left Jacobian, valid for small |phi|.
///
///     V(phi) = sum_{k>=0} (-1)^k phi^k / (k+2)! * [ [0, -1], [1, 0] ]^k
/// Expanding the first few terms avoids the 0/0 cancellation in the closed
/// form, which otherwise loses most of its digits below |phi| ~ 1e-4.
Mat2 leftJacobianSeries(Scalar phi) {
  // [c, -s; s, c] coefficients from the even/odd split of the series.
  //   sin(phi)/phi  = 1 - phi^2/6 + phi^4/120 - ...
  //   (1-cos(phi))/phi = phi/2 - phi^3/24 + phi^5/720 - ...
  const Scalar p2 = phi * phi;
  const Scalar a = 1.0 - p2 / 6.0 + p2 * p2 / 120.0;          // sin(phi)/phi
  const Scalar b = phi * (0.5 - p2 / 24.0 + p2 * p2 / 720.0);  // (1-cos phi)/phi
  Mat2 V;
  V << a, -b, b, a;
  return V;
}

}  // namespace

Mat2 se2LeftJacobian(Scalar phi) {
  if (std::abs(phi) < 1e-4) return leftJacobianSeries(phi);
  const Scalar s = std::sin(phi), c = std::cos(phi);
  Mat2 V;
  V << s / phi, -(1.0 - c) / phi, (1.0 - c) / phi, s / phi;
  return V;
}

Mat2 se2LeftJacobianInverse(Scalar phi) {
  if (std::abs(phi) < 1e-4) {
    // Inverse of the series, expanded to the same order:
    //     V^{-1} = [ 1 - phi^2/12,  phi/2 ]
    //              [ -phi/2,       1 - phi^2/12 ]
    // Note the diagonal coefficient is *negative*; the closed form below
    // divides by 2 - 2cos(phi), and expanding that gives 1 - phi^2/12. Using
    // +phi^2/12 here instead is a sign error that is invisible at |phi| ~ 1e-2
    // (where it is ~1e-5) but grows to ~1e-3 by |phi| ~ 1 and silently biases
    // every small rotation in the pose graph.
    const Scalar p2 = phi * phi;
    const Scalar d = 1.0 - p2 / 12.0;
    Mat2 Vinv;
    Vinv << d, phi / 2.0, -phi / 2.0, d;
    return Vinv;
  }
  const Scalar s = std::sin(phi), c = std::cos(phi);
  // V^{-1} = 1/(2(1-c)) * [ phi sin f,  phi (1-c) ]
  //                       [ -phi(1-c), phi sin f  ]
  // The denominator is written as 2 - 2c rather than 2(1-c) to keep the
  // relative error small when phi is near pi, where 1-c is close to 2 and
  // the naive form subtracts two nearly equal numbers.
  const Scalar denom = 2.0 - 2.0 * c;
  const Scalar k = phi * s / denom;
  const Scalar m = phi * (1.0 - c) / denom;
  Mat2 Vinv;
  Vinv << k, m, -m, k;
  return Vinv;
}

Mat3 se2Exp(const Vec3& xi) {
  const Scalar phi = xi(2);
  const Mat2 V = se2LeftJacobian(phi);
  const Scalar c = std::cos(phi), s = std::sin(phi);
  Mat3 T = Mat3::Identity();
  T(0, 0) = c; T(0, 1) = -s;
  T(1, 0) = s; T(1, 1) = c;
  const Vec2 t = V * Vec2(xi(0), xi(1));
  T(0, 2) = t.x();
  T(1, 2) = t.y();
  return T;
}

Vec3 se2Log(const Mat3& T) {
  const Scalar phi = std::atan2(T(1, 0), T(0, 0));
  const Vec2 t(T(0, 2), T(1, 2));
  const Vec2 rho = se2LeftJacobianInverse(phi) * t;
  return Vec3(rho.x(), rho.y(), phi);
}

Mat3 se2AdMatrix(const Vec3& xi) {
  Mat3 ad = Mat3::Zero();
  ad(0, 1) = -xi(2);
  ad(0, 2) = xi(1);
  ad(1, 0) = xi(2);
  ad(1, 2) = -xi(0);
  return ad;
}

namespace {

/// Closed form, valid for moderate |phi|, from ad^3 = -phi^2 ad:
///     J_l = I + a1 ad + a2 ad^2,  J_r = I - a1 ad + a2 ad^2
///     a1 = (1-cos phi)/phi^2,     a2 = (phi - sin phi)/phi^3
void jacobianCoefficients(Scalar phi, Scalar* a1, Scalar* a2) {
  const Scalar p2 = phi * phi;
  *a1 = (1.0 - std::cos(phi)) / p2;
  *a2 = (phi - std::sin(phi)) / (p2 * phi);
}

/// Series coefficients for small |phi|.
///     a1 = 1/2 - phi^2/24 + phi^4/720
///     a2 = 1/6 - phi^2/120 + phi^4/5040
/// The closed form subtracts nearly equal numbers here and loses most of its
/// digits well before the series does.
void jacobianCoefficientsSeries(Scalar phi, Scalar* a1, Scalar* a2) {
  const Scalar p2 = phi * phi;
  const Scalar p4 = p2 * p2;
  *a1 = 0.5 - p2 / 24.0 + p4 / 720.0;
  *a2 = 1.0 / 6.0 - p2 / 120.0 + p4 / 5040.0;
}

}  // namespace

Mat3 se2LeftJacobianMat(const Vec3& xi) {
  const Scalar phi = xi(2);
  Scalar a1, a2;
  if (std::abs(phi) < 1e-4) {
    jacobianCoefficientsSeries(phi, &a1, &a2);
  } else {
    jacobianCoefficients(phi, &a1, &a2);
  }
  const Mat3 ad = se2AdMatrix(xi);
  return Mat3::Identity() + a1 * ad + a2 * ad * ad;
}

Mat3 se2RightJacobianMat(const Vec3& xi) {
  const Scalar phi = xi(2);
  Scalar a1, a2;
  if (std::abs(phi) < 1e-4) {
    jacobianCoefficientsSeries(phi, &a1, &a2);
  } else {
    jacobianCoefficients(phi, &a1, &a2);
  }
  const Mat3 ad = se2AdMatrix(xi);
  return Mat3::Identity() - a1 * ad + a2 * ad * ad;
}

Mat3 se2RightJacobianMatInverse(const Vec3& xi) {
  // Computed as an explicit 3x3 inverse rather than by a truncated series.
  // The series for J_r^{-1} is sum ad^n/(n+1)! = I + ad/2 + ad^2/12 - ad^4/720
  // + ... and because ad^3 = -phi^2 ad it never terminates: every truncation
  // leaves an error of order phi^2. That error is ~1e-2 at phi ~ 1, which is
  // large enough to make the pose graph converge to a visibly biased answer
  // while every finite-difference check at a *small* angle still passes.
  // A 3x3 inverse is a handful of flops per edge and is exact.
  return se2RightJacobianMat(xi).inverse();
}

Mat3 se2LeftJacobianMatInverse(const Vec3& xi) { return se2LeftJacobianMat(xi).inverse(); }

Mat3 se2Adjoint(const Mat3& T) {
  const Scalar phi = std::atan2(T(1, 0), T(0, 0));
  const Scalar c = std::cos(phi), s = std::sin(phi);
  Mat3 Ad = Mat3::Zero();
  Ad(0, 0) = c; Ad(0, 1) = -s;
  Ad(1, 0) = s; Ad(1, 1) = c;
  // (I - R(phi)) t / phi == (t_y, -t_x) for every phi != 0.
  Ad(0, 2) = T(1, 2);
  Ad(1, 2) = -T(0, 2);
  Ad(2, 2) = 1.0;
  return Ad;
}

}  // namespace sc
