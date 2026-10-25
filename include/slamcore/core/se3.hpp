// SLAMCore-8G -- SO(3)/SE(3) Lie group operations, analytic Jacobians.
//
// Conventions (used consistently across every module):
//   * Right-perturbation update:  T (+) d = T * Exp(d)
//   * Rotation matrix -> vector via Log, branch-cut at theta = pi handled by
//     the standard largest-diagonal fallback.
//   * SE(3) tangent ordering is [rho; phi] (translation first), matching
//     exp/log in the closed form below.
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "slamcore/core/types.hpp"

namespace sc {

using Vec2 = Eigen::Matrix<Scalar, 2, 1>;
using Vec3 = Eigen::Matrix<Scalar, 3, 1>;
using Vec4 = Eigen::Matrix<Scalar, 4, 1>;
using Vec6 = Eigen::Matrix<Scalar, 6, 1>;
using Vec7 = Eigen::Matrix<Scalar, 7, 1>;
using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using Mat2 = Eigen::Matrix<Scalar, 2, 2>;
using Mat3 = Eigen::Matrix<Scalar, 3, 3>;
using Mat4 = Eigen::Matrix<Scalar, 4, 4>;
using Mat6 = Eigen::Matrix<Scalar, 6, 6>;
using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using Quat = Eigen::Quaternion<Scalar>;

/// Skew-symmetric matrix [v]_x such that [v]_x w = v x w.
Mat3 skew(const Vec3& v);

/// SO(3) exponential map (Rodrigues). Stable near 0.
Mat3 so3Exp(const Vec3& phi);
/// SO(3) logarithm; returns rotation vector with |phi| <= pi.
Vec3 so3Log(const Mat3& R);
/// Right Jacobian of SO(3): J_r(phi) = I - (1-cos t)/t^2 [phi]_x + (t-sin t)/t^3 [phi]_x^2
Mat3 so3RightJacobian(const Vec3& phi);
/// Inverse right Jacobian, used when converting covariance parameterisations.
Mat3 so3RightJacobianInverse(const Vec3& phi);

/// SE(3) exponential map of [rho; phi].
Mat4 se3Exp(const Vec6& xi);
/// SE(3) logarithm of a rigid transform, inverse of se3Exp.
Vec6 se3Log(const Mat4& T);
/// Adjoint of SE(3): Ad_T = [[R, [t]_x R], [0, R]].
Mat6 se3Adjoint(const Mat4& T);
/// Right Jacobian of SE(3) (used for covariance propagation).
Mat6 se3RightJacobian(const Vec6& xi);
/// Inverse of se3RightJacobian.
Mat6 se3RightJacobianInverse(const Vec6& xi);

/// Compose a rigid transform from rotation and translation.
inline Mat4 makeSE3(const Mat3& R, const Vec3& t) {
  Mat4 T = Mat4::Identity();
  T.template block<3, 3>(0, 0) = R;
  T.template block<3, 1>(0, 3) = t;
  return T;
}

inline Mat3 rotationOf(const Mat4& T) { return T.template block<3, 3>(0, 0); }
inline Vec3 translationOf(const Mat4& T) { return T.template block<3, 1>(0, 3); }

/// Inverse of a rigid transform exploiting its structure.
inline Mat4 se3Inverse(const Mat4& T) {
  const Mat3 R = rotationOf(T).transpose();
  Mat4 out = Mat4::Identity();
  out.template block<3, 3>(0, 0) = R;
  out.template block<3, 1>(0, 3) = -R * translationOf(T);
  return out;
}

/// 2D pose (x, y, theta) helpers -- used by the simulator, planner and
/// controller, which live in the plane.
struct Pose2 {
  Scalar x = 0, y = 0, th = 0;
};

inline Mat3 pose2ToSE2(const Pose2& p) {
  Mat3 T = Mat3::Identity();
  const Scalar c = std::cos(p.th), s = std::sin(p.th);
  T(0, 0) = c; T(0, 1) = -s; T(0, 2) = p.x;
  T(1, 0) = s; T(1, 1) = c;  T(1, 2) = p.y;
  return T;
}

inline Vec2 pose2Apply(const Pose2& p, const Vec2& q) {
  const Scalar c = std::cos(p.th), s = std::sin(p.th);
  return Vec2(c * q.x() - s * q.y() + p.x, s * q.x() + c * q.y() + p.y);
}

inline Pose2 pose2Compose(const Pose2& a, const Pose2& b) {
  const Scalar c = std::cos(a.th), s = std::sin(a.th);
  return Pose2{a.x + c * b.x - s * b.y, a.y + s * b.x + c * b.y, wrapPi(a.th + b.th)};
}

inline Pose2 pose2Inverse(const Pose2& p) {
  const Scalar c = std::cos(p.th), s = std::sin(p.th);
  return Pose2{-(c * p.x + s * p.y), -(-s * p.x + c * p.y), wrapPi(-p.th)};
}

/// Normalise an angle-difference weighted residual vector in place.
Vec3 so3Residual(const Mat3& R_meas, const Mat3& R_pred);

}  // namespace sc
