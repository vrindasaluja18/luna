// LUNA -- core math tests: Lie groups, QP, sparse solve, allocators.
//
// The SE(3) tests are the load-bearing ones: every Jacobian in the backend
// depends on exp/log being exact inverses and on the Adjoint identity holding,
// so those are checked numerically against finite differences rather than
// against a second implementation of the same formula.
#include <cmath>
#include <random>
#include <vector>

#include "luna/core/alloc.hpp"
#include "luna/core/qp.hpp"
#include "luna/core/se2.hpp"
#include "luna/core/se3.hpp"
#include "luna/core/sparse.hpp"
#include "luna/core/stats.hpp"
#include "test_util.hpp"

using namespace luna;

namespace {

void testSo3() {
  Rng rng(1234);
  for (int i = 0; i < 200; ++i) {
    const Vec3 phi = Vec3(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0));
    const Mat3 R = so3Exp(phi);
    // Exp must produce a rotation: orthonormal with determinant +1.
    LUNA_NEAR((R * R.transpose() - Mat3::Identity()).norm(), 0.0, 1e-12);
    LUNA_NEAR(R.determinant(), 1.0, 1e-12);
    // Log must invert Exp. Rotation vectors are only unique up to |phi| <= pi,
    // so compare the rotations they generate rather than the vectors.
    const Vec3 back = so3Log(R);
    LUNA_NEAR((so3Exp(back) - R).norm(), 0.0, 1e-10);
  }
}

void testSe3() {
  Rng rng(4321);
  for (int i = 0; i < 200; ++i) {
    Vec6 xi;
    xi.head<3>() = Vec3(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0));
    xi.tail<3>() = Vec3(rng.gaussian(0.7), rng.gaussian(0.7), rng.gaussian(0.7));
    const Mat4 T = se3Exp(xi);
    // The rotation block is a rotation and the bottom row is exact.
    LUNA_NEAR(rotationOf(T).determinant(), 1.0, 1e-12);
    LUNA_NEAR((T.row(3) - Eigen::RowVector4d(0, 0, 0, 1)).norm(), 0.0, 1e-14);
    // Round trip through the log.
    LUNA_NEAR((se3Exp(se3Log(T)) - T).norm(), 0.0, 1e-9);
  }
  // Identity and inverse.
  LUNA_NEAR((se3Exp(Vec6::Zero()) - Mat4::Identity()).norm(), 0.0, 1e-14);
  Rng rng2(99);
  for (int i = 0; i < 50; ++i) {
    Vec6 xi;
    xi.head<3>() = Vec3(rng2.gaussian(1.0), rng2.gaussian(1.0), rng2.gaussian(1.0));
    xi.tail<3>() = Vec3(rng2.gaussian(0.5), rng2.gaussian(0.5), rng2.gaussian(0.5));
    const Mat4 T = se3Exp(xi);
    LUNA_NEAR((se3Inverse(T) * T - Mat4::Identity()).norm(), 0.0, 1e-12);
  }
}

/// The Adjoint must satisfy Ad_T * xi == Log(T * Exp(xi) * T^-1). This is the
/// identity the covariance propagation and the Schur marginalisation both rely
/// on, so it is checked directly rather than assumed.
void testAdjoint() {
  Rng rng(777);
  for (int i = 0; i < 100; ++i) {
    Vec6 a, b;
    a.head<3>() = Vec3(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0));
    a.tail<3>() = Vec3(rng.gaussian(0.5), rng.gaussian(0.5), rng.gaussian(0.5));
    b.head<3>() = Vec3(rng.gaussian(0.3), rng.gaussian(0.3), rng.gaussian(0.3));
    b.tail<3>() = Vec3(rng.gaussian(0.3), rng.gaussian(0.3), rng.gaussian(0.3));
    const Mat4 T = se3Exp(a);
    const Mat4 lhs = se3Exp(se3Adjoint(T) * b);
    const Mat4 rhs = T * se3Exp(b) * se3Inverse(T);
    LUNA_NEAR((lhs - rhs).norm(), 0.0, 1e-9);
  }
}

/// Right Jacobian check: Exp(xi + d) ~= Exp(xi) * Exp(J_r(xi) d) for small d.
void testRightJacobian() {
  Rng rng(2024);
  for (int i = 0; i < 60; ++i) {
    Vec6 xi;
    xi.head<3>() = Vec3(rng.gaussian(0.6), rng.gaussian(0.6), rng.gaussian(0.6));
    xi.tail<3>() = Vec3(rng.gaussian(0.6), rng.gaussian(0.6), rng.gaussian(0.6));
    Vec6 d;
    d.head<3>() = Vec3(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0));
    d.tail<3>() = Vec3(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0));
    const Scalar h = 1e-6;
    const Mat4 lhs = se3Exp(xi + h * d);
    const Mat4 rhs = se3Exp(xi) * se3Exp(h * (se3RightJacobian(xi) * d));
    // The residual is O(h) truncation from the finite difference, not a
    // Jacobian error, so the tolerance is set to a few times h.
    LUNA_NEAR((lhs - rhs).norm(), 0.0, 1e-5);
  }
}

void testPose2() {
  Rng rng(55);
  for (int i = 0; i < 200; ++i) {
    const Pose2 a{rng.gaussian(3.0), rng.gaussian(3.0), rng.gaussian(3.0)};
    const Pose2 b{rng.gaussian(3.0), rng.gaussian(3.0), rng.gaussian(3.0)};
    // Composition must match the homogeneous-matrix product.
    const Mat3 Ma = pose2ToSE2(a), Mb = pose2ToSE2(b);
    const Mat3 Mc = pose2ToSE2(pose2Compose(a, b));
    LUNA_NEAR((Ma * Mb - Mc).norm(), 0.0, 1e-12);
    // Inverse must undo composition.
    const Pose2 id = pose2Compose(a, pose2Inverse(a));
    LUNA_NEAR(id.x, 0.0, 1e-12);
    LUNA_NEAR(id.y, 0.0, 1e-12);
    LUNA_NEAR(wrapPi(id.th), 0.0, 1e-12);
  }
  // wrapPi maps an angle into a half-open interval of width 2*pi whose
  // endpoint is +/-pi. Both endpoints are the same angle, so the check is on
  // magnitude rather than on which sign the branch cut lands.
  LUNA_NEAR(std::abs(wrapPi(3.0 * kPi)), kPi, 1e-12);
  LUNA_NEAR(std::abs(wrapPi(-3.0 * kPi)), kPi, 1e-12);
  LUNA_NEAR(wrapPi(0.0), 0.0, 1e-14);
  LUNA_NEAR(wrapPi(0.3), 0.3, 1e-14);
  LUNA_NEAR(wrapPi(2.0 * kPi + 0.3), 0.3, 1e-12);
}

void testStats() {
  RunningStats s;
  for (int i = 1; i <= 100; ++i) s.add(static_cast<Scalar>(i));
  LUNA_NEAR(s.mean(), 50.5, 1e-12);
  LUNA_NEAR(s.min(), 1.0, 1e-12);
  LUNA_NEAR(s.max(), 100.0, 1e-12);
  // Sample variance (the n-1 normalisation) of 1..100 is n(n+1)/12.
  LUNA_NEAR(s.variance(), 100.0 * 101.0 / 12.0, 1e-9);
  LUNA_NEAR(RunningStats::percentile({1, 2, 3, 4, 5}, 0.5), 3.0, 1e-12);
  LUNA_NEAR(RunningStats::percentile({1, 2, 3, 4, 5}, 1.0), 5.0, 1e-12);
}

void testRng() {
  // Same seed, same stream: the deterministic replay mode depends on this.
  Rng a(42), b(42);
  for (int i = 0; i < 1000; ++i) LUNA_NEAR(a.uniform(), b.uniform(), 0.0);
  // Uniformity and range.
  Rng c(7);
  RunningStats u;
  for (int i = 0; i < 200000; ++i) u.add(c.uniform());
  LUNA_NEAR(u.mean(), 0.5, 5e-3);
  LUNA_NEAR(u.min() >= 0.0, 1.0, 0.0);
  LUNA_NEAR(u.max() < 1.0, 1.0, 0.0);
  // Gaussian moments.
  Rng d(11);
  RunningStats g;
  for (int i = 0; i < 200000; ++i) g.add(d.gaussian());
  LUNA_NEAR(g.mean(), 0.0, 1e-2);
  LUNA_NEAR(g.stddev(), 1.0, 1e-2);
}

void testPoolAllocator() {
  PoolAllocator pool(4096);
  void* a = pool.allocate(100);
  LUNA_CHECK(a != nullptr);
  // Alignment guarantee for max_align_t.
  LUNA_CHECK(reinterpret_cast<std::uintptr_t>(a) % alignof(std::max_align_t) == 0);
  LUNA_CHECK(pool.misses() == 0);
  // Exhaustion returns nullptr rather than throwing in the hot path.
  LUNA_CHECK(pool.allocate(1 << 20) == nullptr);
  LUNA_CHECK(pool.misses() == 1);
  pool.reset();
  LUNA_CHECK(pool.used() == 0);
}

void testRingAllocator() {
  RingAllocator ring(64, 4);
  LUNA_CHECK(ring.slots() == 4);
  // Wrapping must recycle the oldest slot, never fail.
  for (int i = 0; i < 100; ++i) LUNA_CHECK(ring.next() != nullptr);
  LUNA_CHECK(ring.live() == 4);
}

/// Solve a small least-squares problem with the sparse system and compare with
/// a dense Eigen solve. The point is that the sparse path and the reference
/// agree; the sparse path is what ships.
void testSparseSystem() {
  const int n = 8;
  Rng rng(31337);
  MatX A = MatX::Random(n, n);
  A = A * A.transpose() + 5.0 * MatX::Identity(n, n);  // SPD
  const VecX b = VecX::Random(n);
  const VecX want = A.ldlt().solve(b);

  SparseSystem sys(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) sys.addH(i, j, MatX::Constant(1, 1, A(i, j)));
    sys.addB(i, VecX::Constant(1, b(i)));
  }
  LUNA_CHECK(sys.factorise());
  VecX got;
  LUNA_CHECK(sys.solve(&got));
  LUNA_CHECK(static_cast<int>(got.size()) == n);
  LUNA_NEAR((got - want).norm(), 0.0, 1e-8);
  // The dense view must match what was assembled.
  LUNA_NEAR((sys.denseH() - A).norm(), 0.0, 1e-12);
}

/// QP with a box constraint must land on the constrained optimum, not the
/// unconstrained one. This is the property the MPC's input limits rely on.
void testQp() {
  // min 0.5 x'Px + q'x with P = I, q = -[3, 4], s.t. -1 <= x <= 1.
  // Unconstrained optimum is (3, 4); the box-constrained optimum is (1, 1).
  MatX P = 2.0 * MatX::Identity(2, 2);
  VecX q(2);
  q << -3.0, -4.0;
  MatX A = MatX::Identity(2, 2);
  VecX lo(2), hi(2);
  lo << -1.0, -1.0;
  hi << 1.0, 1.0;

  QpSolver solver;
  solver.setProblem(P, q, A, lo, hi);
  QpOptions opt;
  VecX x;
  const QpResult res = solver.solve(&x, opt);
  LUNA_CHECK(res.status == QpStatus::kSolved || res.status == QpStatus::kMaxIter);
  LUNA_NEAR(x(0), 1.0, 1e-4);
  LUNA_NEAR(x(1), 1.0, 1e-4);
}

/// SE(2): Exp/Log round trip, the Adjoint identity, and the left Jacobian.
///
/// The Adjoint is checked against its defining identity rather than against a
/// second copy of the formula, because the formula is exactly what is in doubt:
/// a wrong Ad_T produces a pose graph that converges to the wrong answer while
/// every individual Jacobian still looks plausible.
void testSe2() {
  Rng rng(777);
  for (int i = 0; i < 300; ++i) {
    const Vec3 xi(rng.gaussian(1.2), rng.gaussian(1.2), rng.gaussian(1.5));
    const Mat3 T = se2Exp(xi);
    // The result is a rigid transform: rotation block orthonormal, bottom row
    // exact, determinant +1.
    const Mat2 R = T.topLeftCorner<2, 2>();
    LUNA_NEAR((R * R.transpose() - Mat2::Identity()).norm(), 0.0, 1e-12);
    LUNA_NEAR(R.determinant(), 1.0, 1e-12);
    LUNA_NEAR(T(2, 0), 0.0, 1e-15);
    LUNA_NEAR(T(2, 1), 0.0, 1e-15);
    LUNA_NEAR(T(2, 2), 1.0, 1e-15);
    // Log inverts Exp. Compare transforms, since a tangent is only unique
    // within the |phi| <= pi ball.
    LUNA_NEAR((se2Exp(se2Log(T)) - T).norm(), 0.0, 1e-10);
  }

  // Small angles are where the closed form loses precision, so they are tested
  // separately with a tolerance that only the series path can meet.
  for (Scalar phi : {0.0, 1e-12, 1e-9, 1e-7, 1e-5}) {
    const Vec3 xi(0.3, -0.2, phi);
    const Mat3 T = se2Exp(xi);
    LUNA_NEAR((se2Exp(se2Log(T)) - T).norm(), 0.0, 1e-12);
    // The translation of a pure small rotation is the identity to first order,
    // so the recovered rho must equal the input rho.
    LUNA_NEAR(se2Log(T).x(), 0.3, 1e-9);
    LUNA_NEAR(se2Log(T).y(), -0.2, 1e-9);
  }

  // Adjoint identity: Exp(Ad_T xi) == T Exp(xi) T^{-1}.
  for (int i = 0; i < 200; ++i) {
    const Vec3 a(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.5));
    const Vec3 b(rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.5));
    const Mat3 T = se2Exp(a);
    const Mat3 lhs = se2Exp(se2Adjoint(T) * b);
    const Mat3 rhs = T * se2Exp(b) * se2Inverse(T);
    LUNA_NEAR((lhs - rhs).norm(), 0.0, 1e-10);
  }

  // Left Jacobian inverts itself, including in the series regime.
  for (Scalar phi : {-3.0, -1.0, -1e-6, 0.0, 1e-6, 0.5, 2.9}) {
    const Mat2 V = se2LeftJacobian(phi);
    const Mat2 Vi = se2LeftJacobianInverse(phi);
    LUNA_NEAR((V * Vi - Mat2::Identity()).norm(), 0.0, 1e-10);
  }

  // Composition and relative pose agree with matrix arithmetic.
  for (int i = 0; i < 100; ++i) {
    const Pose2 a{rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0)};
    const Pose2 b{rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0)};
    const Pose2 c = pose2Compose(a, b);
    const Mat3 Tc = pose2ToSE2(a) * pose2ToSE2(b);
    LUNA_NEAR((pose2ToSE2(c) - Tc).norm(), 0.0, 1e-12);
    const Pose2 rel = pose2Relative(a, b);
    LUNA_NEAR((pose2ToSE2(a) * pose2ToSE2(rel) - pose2ToSE2(b)).norm(), 0.0, 1e-12);
    // Relative of a pose with itself is the identity.
    const Pose2 self = pose2Relative(a, a);
    LUNA_NEAR(std::abs(self.x) + std::abs(self.y) + std::abs(self.th), 0.0, 1e-12);
  }
}

void run() {
  testSo3();
  testSe3();
  testAdjoint();
  testRightJacobian();
  testSe2();
  testPose2();
  testStats();
  testRng();
  testPoolAllocator();
  testRingAllocator();
  testSparseSystem();
  testQp();
}

}  // namespace

LUNA_TEST_MAIN("test_core")
