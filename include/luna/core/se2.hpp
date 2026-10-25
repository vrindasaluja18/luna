// LUNA -- SE(2) Lie group operations for the planar backend.
//
// Why SE(2) and not SE(3): the simulator's ground truth is a planar pose
// (x, y, theta) and every sensor it emulates is a 2D measurement. A 3D
// estimator over 2D data would add six degrees of freedom that the data cannot
// observe, so ATE measured against it would be meaningless. The 3D group lives
// in se3.hpp for the modules that need it; this one is what the pose graph and
// the trajectory evaluator use.
//
// Conventions match se3.hpp exactly:
//   * Right-perturbation update:  T (+) d = T * Exp(d)
//   * Tangent ordering [rho; phi] = [dx; dy; dtheta] (translation first)
//   * Homogeneous 3x3 matrices, so composition is plain matrix product
#pragma once

#include "luna/core/se3.hpp"

namespace luna {

/// SE(2) exponential map of [rho; phi] = [dx; dy; dtheta].
///
/// Closed form via the left Jacobian V(phi):
///     V(phi) = (1/phi) [ sin f, -(1 - cos f) ]
///                      [ 1 - cos f, sin f    ]
///     Exp(xi) = [ R(phi)  V(phi) rho ; 0 0 1 ]
/// Near phi = 0 the closed form loses precision to cancellation, so the series
/// expansion of V is used instead.
Mat3 se2Exp(const Vec3& xi);

/// SE(2) logarithm, inverse of se2Exp. Returns phi in (-pi, pi].
Vec3 se2Log(const Mat3& T);

/// Adjoint of SE(2): Ad_T xi is the tangent that satisfies
///     Exp(Ad_T xi) = T Exp(xi) T^{-1}
/// Derivation: conjugation preserves the rotation angle, and the translation
/// picks up (I - R(phi)) t = phi (t_y, -t_x) to first order, giving
///     Ad_T = [ R  (t_y, -t_x) ; 0 0 1 ].
Mat3 se2Adjoint(const Mat3& T);

/// Left Jacobian V(phi) and its inverse, exposed because se2Exp/se2Log need
/// them directly. These are the 2x2 translation blocks only.
Mat2 se2LeftJacobian(Scalar phi);
Mat2 se2LeftJacobianInverse(Scalar phi);

/// The 3x3 adjoint matrix ad_xi, i.e. the generator such that
///     Ad_{Exp(xi)} = exp(ad_xi).
/// For SE(2) it is
///     ad = [ 0  -phi   y ]
///          [ phi  0   -x ]
///          [ 0    0    0 ]
Mat3 se2AdMatrix(const Vec3& xi);

/// 3x3 left and right Jacobians in the [rho; phi] tangent ordering.
///
/// These are NOT block diagonal. SE(2) is non-abelian, so the rotation tangent
/// couples into the translation tangent, and the coupling column is exactly
/// what makes the pose-graph Jacobians correct. Padding the 2x2 V(phi) with a 1
/// -- the obvious-looking shortcut -- produces a Jacobian that is wrong by a
/// term of order phi*rho, which on a turn of 0.8 rad and an arm of 0.6 m is
/// already ~0.5. That error is invisible on a straight-line test and biases
/// every turning segment.
///
/// Closed form, from ad^3 = -phi^2 ad:
///     J_l = I + (1-cos phi)/phi^2 * ad + (phi - sin phi)/phi^3 * ad^2
///     J_r = I - (1-cos phi)/phi^2 * ad + (phi - sin phi)/phi^3 * ad^2
/// Both are evaluated by their series below |phi| ~ 1e-4, where the closed
/// form cancels catastrophically.
Mat3 se2LeftJacobianMat(const Vec3& xi);
Mat3 se2RightJacobianMat(const Vec3& xi);
/// Inverses, used to convert a residual into a correction in the pose-graph
/// linearisation. Computed by explicit 3x3 inversion: the series form does not
/// terminate for SE(2), so any truncated series carries an O(phi^2) error.
Mat3 se2LeftJacobianMatInverse(const Vec3& xi);
Mat3 se2RightJacobianMatInverse(const Vec3& xi);

/// Build an SE(2) matrix from a planar pose.
inline Mat3 pose2ToMat(const Pose2& p) { return pose2ToSE2(p); }

/// Extract a planar pose from an SE(2) matrix.
inline Pose2 matToPose2(const Mat3& T) {
  return Pose2{T(0, 2), T(1, 2), std::atan2(T(1, 0), T(0, 0))};
}

/// Inverse of an SE(2) matrix, exploiting its structure.
inline Mat3 se2Inverse(const Mat3& T) {
  const Scalar c = T(0, 0), s = T(1, 0);
  Mat3 out = Mat3::Identity();
  out(0, 0) = c; out(0, 1) = s;
  out(1, 0) = -s; out(1, 1) = c;
  out(0, 2) = -(c * T(0, 2) + s * T(1, 2));
  out(1, 2) = -(-s * T(0, 2) + c * T(1, 2));
  return out;
}

/// Compose a planar pose with a tangent increment: p (+) d = p * Exp(d).
inline Pose2 pose2Update(const Pose2& p, const Vec3& d) {
  return matToPose2(pose2ToSE2(p) * se2Exp(d));
}

/// Relative pose of `b` expressed in the frame of `a`: a^{-1} b.
inline Pose2 pose2Relative(const Pose2& a, const Pose2& b) {
  return matToPose2(se2Inverse(pose2ToSE2(a)) * pose2ToSE2(b));
}

}  // namespace luna
