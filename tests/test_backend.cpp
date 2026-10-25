// LUNA -- backend tests: pose-graph Jacobians, convergence, robustness,
// loop-closure verification, and trajectory metrics.
//
// The Jacobian test is the important one. A wrong analytic Jacobian still
// converges on a small graph with good initialisation, so the pose graph can
// look correct while biasing every result. Here the analytic Jacobians are
// compared against central finite differences of the *actual residual
// function*, which is the only check that cannot be satisfied by a matching
// pair of errors.
#include <cmath>
#include <vector>

#include "luna/backend/loop_close.hpp"
#include "luna/backend/pose_graph.hpp"
#include "luna/bench/trajectory.hpp"
#include "luna/core/se2.hpp"
#include "luna/core/stats.hpp"
#include "test_util.hpp"

using namespace luna;

namespace {

/// Recompute one edge's residual at perturbed node poses.
Vec3 edgeResidual(const backend::PoseEdge& e, const std::vector<Pose2>& poses) {
  const Mat3 Ti = pose2ToSE2(poses[e.i]);
  const Mat3 Tj = pose2ToSE2(poses[e.j]);
  return se2Log(se2Inverse(se2Exp(e.z)) * (se2Inverse(Ti) * Tj));
}

/// The analytic Jacobians, recomputed here from the same formula the library
/// uses, so the test can difference the residual independently.
void analyticJacobians(const backend::PoseEdge& e, const std::vector<Pose2>& poses, Mat3* Ji,
                       Mat3* Jj) {
  const Mat3 Ti = pose2ToSE2(poses[e.i]);
  const Mat3 Tj = pose2ToSE2(poses[e.j]);
  const Mat3 einv = se2Inverse(Ti) * Tj;
  const Vec3 r = se2Log(se2Inverse(se2Exp(e.z)) * einv);
  (void)einv;
  // Node j takes the right Jacobian inverse, node i the left one; see the
  // derivation in pose_graph.cpp.
  *Jj = se2RightJacobianMatInverse(r);
  *Ji = -se2LeftJacobianMatInverse(r) * se2Adjoint(se2Inverse(se2Exp(e.z)));
}

/// Finite-difference check of both edge Jacobians.
///
/// The perturbation convention must match the library's right-perturbation
/// update: x_i -> x_i Exp(dx_i). Getting that convention wrong here would make
/// the test pass against a consistently wrong implementation, so it is stated
/// explicitly and checked at two step sizes.
void testEdgeJacobians() {
  Rng rng(2024);
  // Central differences have truncation error O(eps^2) and rounding error
  // O(machine_eps/eps). 1e-5 sits near the minimum of the sum; 1e-6 is
  // already inside the rounding floor and produces noise an order of magnitude
  // larger, which would be misread as a wrong Jacobian.
  const Scalar eps = 1e-5;
  for (int trial = 0; trial < 40; ++trial) {
    std::vector<Pose2> poses(2);
    poses[0] = Pose2{rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0)};
    poses[1] = Pose2{rng.gaussian(1.0), rng.gaussian(1.0), rng.gaussian(1.0)};

    backend::PoseEdge e;
    e.i = 0;
    e.j = 1;
    e.z = Vec3(rng.gaussian(0.4), rng.gaussian(0.4), rng.gaussian(0.5));
    e.info = Mat3::Identity();

    Mat3 Ji, Jj;
    analyticJacobians(e, poses, &Ji, &Jj);

    for (int k = 0; k < 3; ++k) {
      // Central difference in node i.
      std::vector<Pose2> p = poses;
      Vec3 d = Vec3::Zero();
      d(k) = eps;
      p[0] = pose2Update(poses[0], d);
      const Vec3 rp = edgeResidual(e, p);
      d(k) = -eps;
      p[0] = pose2Update(poses[0], d);
      const Vec3 rm = edgeResidual(e, p);
      const Vec3 fd = (rp - rm) / (2.0 * eps);
      LUNA_NEAR((fd - Ji.col(k)).norm(), 0.0, 1e-6);

      // Central difference in node j.
      p = poses;
      d(k) = eps;
      p[1] = pose2Update(poses[1], d);
      const Vec3 rp2 = edgeResidual(e, p);
      d(k) = -eps;
      p[1] = pose2Update(poses[1], d);
      const Vec3 rm2 = edgeResidual(e, p);
      const Vec3 fd2 = (rp2 - rm2) / (2.0 * eps);
      LUNA_NEAR((fd2 - Jj.col(k)).norm(), 0.0, 1e-6);
    }
  }
}

/// A consistent loop must converge to near-zero cost, and the recovered poses
/// must match the poses the measurements were generated from.
void testPoseGraphRecovery() {
  Rng rng(99);
  // Generate a ground-truth trajectory, then synthesise exact relative
  // measurements from it. With exact measurements the optimum is exactly the
  // generating trajectory (up to the pinned gauge), so this test has a known
  // right answer rather than merely "the cost went down".
  std::vector<Pose2> truth;
  Pose2 p{0, 0, 0};
  for (int i = 0; i < 60; ++i) {
    truth.push_back(p);
    const Pose2 step{rng.uniform(0.2, 0.6), rng.uniform(-0.1, 0.1), rng.uniform(-0.3, 0.3)};
    p = pose2Compose(p, step);
  }

  // A perfect closure from the last node back to the first.
  backend::PoseEdge loop;
  loop.i = static_cast<Index>(truth.size() - 1);
  loop.j = 0;
  loop.z = se2Log(se2Inverse(pose2ToSE2(truth.back())) * pose2ToSE2(truth[0]));
  loop.info = Mat3::Identity() * 100.0;

  // Start from a deliberately wrong initial guess: dead reckoning with drift.
  backend::PoseGraph noisy;
  noisy.addNode(truth[0]);
  for (Size i = 1; i < truth.size(); ++i) {
    const Pose2 n{rng.gaussian(0.15), rng.gaussian(0.15), rng.gaussian(0.05)};
    noisy.addNode(pose2Compose(truth[i], n));
  }
  for (Size i = 1; i < truth.size(); ++i) {
    backend::PoseEdge e;
    e.i = static_cast<Index>(i - 1);
    e.j = static_cast<Index>(i);
    e.z = se2Log(se2Inverse(pose2ToSE2(truth[i - 1])) * pose2ToSE2(truth[i]));
    e.info = Mat3::Identity() * 100.0;
    noisy.addEdge(e);
  }
  noisy.addEdge(loop);

  const backend::GraphResult res = noisy.optimize();
  LUNA_CHECK(res.ok);
  LUNA_CHECK(res.chi2_final < res.chi2_initial);

  // Every node must land back on the generating trajectory.
  Scalar worst = 0;
  for (Size i = 0; i < truth.size(); ++i) {
    const Pose2& got = noisy.pose(static_cast<Index>(i));
    worst = std::max(worst, std::hypot(got.x - truth[i].x, got.y - truth[i].y));
  }
  LUNA_NEAR(worst, 0.0, 1e-4);
}

/// A single grossly wrong loop closure must be down-weighted by the robust
/// kernel. The test compares against the same graph with robustness disabled,
/// because the meaningful claim is "the kernel suppresses the outlier", not
/// "the outlier has no effect at all": one false edge against 39 good ones
/// still tugs the chain a little, and pretending otherwise would only be true
/// for a kernel so aggressive it also suppressed real closures.
void testRobustness() {
  const int n_nodes = 40;

  auto build = [&](bool robust_bad) {
    backend::PoseGraph graph;
    std::vector<Pose2> truth;
    Pose2 p{0, 0, 0};
    for (int i = 0; i < n_nodes; ++i) {
      truth.push_back(p);
      p = pose2Compose(p, Pose2{0.5, 0, 0});
    }
    graph.addNode(truth[0]);
    for (Size i = 1; i < truth.size(); ++i) {
      graph.addNode(pose2Compose(truth[i], Pose2{0.05, 0.05, 0}));
    }
    for (Size i = 1; i < truth.size(); ++i) {
      backend::PoseEdge e;
      e.i = static_cast<Index>(i - 1);
      e.j = static_cast<Index>(i);
      e.z = Vec3(0.5, 0, 0);
      e.info = Mat3::Identity() * 400.0;
      graph.addEdge(e);
    }
    // A false closure claiming the far end is back at the start.
    backend::PoseEdge bad;
    bad.i = static_cast<Index>(truth.size() - 1);
    bad.j = 0;
    bad.z = Vec3(0, 0, 0);
    bad.info = Mat3::Identity() * 400.0;
    bad.robust = robust_bad;
    graph.addEdge(bad);
    return graph;
  };

  backend::PoseGraph robust = build(true);
  const backend::GraphResult res = robust.optimize();
  LUNA_CHECK(res.ok);
  // The bad edge must be flagged as the worst residual.
  LUNA_CHECK(!res.outliers.empty());
  LUNA_CHECK(res.edge_chi2[res.outliers[0]] > 100.0);

  backend::PoseGraph plain = build(false);
  const backend::GraphResult res_plain = plain.optimize();
  LUNA_CHECK(res_plain.ok);

  // Measure how far each trajectory sags away from the straight line it should
  // be. The robust run must sag substantially less.
  auto sag = [&](backend::PoseGraph& g) {
    Scalar worst = 0;
    for (Index i = 1; i < static_cast<Index>(n_nodes); ++i) {
      worst = std::max(worst, std::abs(g.pose(i).y));
    }
    return worst;
  };
  const Scalar sag_robust = sag(robust);
  const Scalar sag_plain = sag(plain);
  LUNA_CHECK(sag_robust < 0.5 * sag_plain);
  // And the robust trajectory must stay close to the line.
  LUNA_CHECK(sag_robust < 0.4);
}

/// 3D-3D RANSAC must recover a known transform and reject a scrambled one.
void testRansacVerification() {
  Rng rng(7);
  const Mat3 R = so3Exp(Vec3(0.05, -0.1, 0.6));
  const Vec3 t(1.5, -2.0, 0.3);
  std::vector<Vec3> a, b;
  for (int i = 0; i < 40; ++i) {
    const Vec3 p(rng.uniform(-3, 3), rng.uniform(-3, 3), rng.uniform(-1, 1));
    a.push_back(p);
    b.push_back(R * p + t);
  }
  backend::LoopCloseOptions opt;
  Mat4 T;
  std::vector<int> inliers;
  Size attempts = 0;
  LUNA_CHECK(backend::verifyRansac(a, b, opt, &T, &inliers, &attempts));
  LUNA_CHECK(inliers.size() >= 38);
  LUNA_NEAR((T.topLeftCorner<3, 3>() - R).norm(), 0.0, 1e-6);
  LUNA_NEAR((T.topRightCorner<3, 1>() - t).norm(), 0.0, 1e-6);

  // With a third of the correspondences randomised, RANSAC must still find
  // the true transform and report it through the inlier set.
  std::vector<Vec3> b2 = b;
  for (int i = 0; i < 13; ++i) {
    b2[i * 3] = Vec3(rng.uniform(-8, 8), rng.uniform(-8, 8), rng.uniform(-8, 8));
  }
  std::vector<int> inliers2;
  LUNA_CHECK(backend::verifyRansac(a, b2, opt, &T, &inliers2, &attempts));
  LUNA_CHECK(inliers2.size() >= 24);
  LUNA_NEAR((T.topRightCorner<3, 1>() - t).norm(), 0.0, 1e-3);
}

/// ATE against a trajectory that is exactly the ground truth must be zero, and
/// against a known offset must equal that offset.
void testAte() {
  bench::Trajectory truth;
  for (int i = 0; i <= 100; ++i) {
    truth.push_back(bench::TrajectoryPoint{i * 0.1, Pose2{i * 0.1, 0, 0}});
  }

  // Exact copy: zero error.
  bench::AteResult exact = bench::absoluteTrajectoryError(truth, truth);
  LUNA_CHECK(exact.ok);
  LUNA_NEAR(exact.rmse, 0.0, 1e-12);
  LUNA_NEAR(exact.percent_of_distance, 0.0, 1e-12);

  // A constant 0.3 m offset in y: RMSE must be exactly 0.3.
  bench::Trajectory off = truth;
  for (auto& p : off) p.pose.y += 0.3;
  bench::AteResult shifted = bench::absoluteTrajectoryError(off, truth);
  LUNA_CHECK(shifted.ok);
  LUNA_NEAR(shifted.rmse, 0.3, 1e-9);
  LUNA_NEAR(shifted.max, 0.3, 1e-9);
  // Path length is 10 m, so 0.3 m is 3 percent of distance.
  LUNA_NEAR(shifted.path_length, 10.0, 1e-9);
  LUNA_NEAR(shifted.percent_of_distance, 3.0, 1e-6);

  // An estimate outside the ground-truth time span must be counted as
  // unmatched, never extrapolated.
  bench::Trajectory late = truth;
  late.push_back(bench::TrajectoryPoint{100.0, Pose2{0, 0, 0}});
  bench::AteResult l = bench::absoluteTrajectoryError(late, truth);
  LUNA_CHECK(l.unmatched == 1);
}

/// RPE must detect a pure scale error that ATE at the endpoints would miss.
void testRpe() {
  bench::Trajectory truth;
  for (int i = 0; i <= 100; ++i) truth.push_back(bench::TrajectoryPoint{i * 0.1, Pose2{i * 0.1, 0, 0}});
  // Estimate travels 1.1x as far: a 10 percent scale error.
  bench::Trajectory est;
  for (int i = 0; i <= 100; ++i) est.push_back(bench::TrajectoryPoint{i * 0.1, Pose2{i * 0.11, 0, 0}});

  bench::EvalOptions opt;
  opt.rpe_interval = 1.0;
  bench::RpeResult r = bench::relativePoseError(est, truth, opt);
  LUNA_CHECK(r.ok);
  // Over a 1 s interval the truth moves 1 m and the estimate 1.1 m, so the
  // relative error is 0.1 m.
  LUNA_NEAR(r.mean_trans, 0.1, 1e-6);
  LUNA_CHECK(r.pairs > 80);
}

/// The descriptor database must stay bounded.
void testLoopDbBound() {
  backend::LoopCloseOptions opt;
  opt.max_keyframes = 10;
  backend::LoopCloser lc(opt);
  for (int i = 0; i < 100; ++i) {
    backend::Keyframe kf;
    kf.id = static_cast<Index>(i);
    kf.pose = Pose2{i * 1.0, 0, 0};
    backend::MapPoint mp;
    mp.p = Vec3(1, 0, 0);
    kf.points.push_back(mp);
    lc.addKeyframe(kf);
  }
  LUNA_CHECK(lc.size() == 10);
  LUNA_CHECK(lc.stats().db_evictions == 90);
}

/// A closure between two identical keyframes must be found, and one between
/// unrelated keyframes must not.
void testLoopClosureMatch() {
  backend::LoopCloseOptions opt;
  opt.min_index_gap = 5;
  opt.min_inliers = 8;
  backend::LoopCloser lc(opt);

  auto makeKf = [](Index id, Scalar x, Scalar y, Scalar th, Rng* rng) {
    backend::Keyframe kf;
    kf.id = id;
    kf.pose = Pose2{x, y, th};
    for (int i = 0; i < 30; ++i) {
      backend::MapPoint mp;
      mp.p = Vec3(rng->uniform(-2, 2), rng->uniform(-2, 2), rng->uniform(0.5, 4.0));
      for (int w = 0; w < 4; ++w) mp.desc.w[w] = rng->next();
      kf.points.push_back(mp);
    }
    return kf;
  };

  Rng rng(3);
  // Two keyframes from the same place with the same view, separated in index.
  const backend::Keyframe kf0 = makeKf(0, 0, 0, 0, &rng);
  backend::Keyframe kf1 = makeKf(1, 0.2, 0.0, 0.05, &rng);
  kf1.points = kf0.points;  // same structure seen from a nearby pose
  lc.addKeyframe(kf0);
  lc.addKeyframe(kf1);

  backend::Keyframe kf2 = makeKf(20, 0, 0, 0, &rng);
  kf2.points = kf0.points;
  backend::LoopClosure lc_out;
  LUNA_CHECK(lc.query(kf2, &lc_out));
  LUNA_CHECK(lc_out.inliers >= 8);

  // An unrelated keyframe with random descriptors must not close a loop.
  backend::Keyframe kf3 = makeKf(30, 40, 40, 1.0, &rng);
  backend::LoopClosure bad;
  LUNA_CHECK(!lc.query(kf3, &bad));
}

/// A descriptor match must be rejected by the ratio test when the two best
/// candidates are equally close, which is the repetitive-corridor failure.
void testDescriptorRatio() {
  backend::LoopCloseOptions opt;
  opt.ratio = 0.8;
  backend::LoopCloser lc(opt);

  std::vector<backend::MapPoint> a(1), b(2);
  a[0].desc.w[0] = 0xFFFFFFFFFFFFFFFFull;
  // Two database descriptors equidistant from the query: ambiguous, so the
  // ratio test must reject it rather than pick arbitrarily.
  b[0].desc.w[0] = 0xFFFFFFFFFFFFFFFEull;
  b[1].desc.w[0] = 0xFFFFFFFFFFFFFFFDull;
  std::vector<std::pair<int, int>> matches;
  LUNA_CHECK(lc.matchDescriptors(a, b, &matches) == 0);

  // With one clearly better candidate the match must be accepted.
  b[1].desc.w[0] = 0x0000000000000000ull;
  LUNA_CHECK(lc.matchDescriptors(a, b, &matches) == 1);
}

/// The accepted closure's transform must be the one *measured from the matched
/// points*, expressed in SE(2) tangent coordinates.
///
/// Two bugs lived here and both were invisible without this test:
///
///   1. `z` was set from the two keyframes' current *estimated* poses. Those
///      estimates already contain the drift the closure exists to correct, so
///      the edge residual came out at zero and the loop added no information at
///      all. The pipeline reported loop closures while the graph did nothing.
///   2. Even once `z` came from RANSAC, it was stored as a raw rigid transform
///      while a graph edge expects an SE(2) tangent vector. The two differ by
///      the left Jacobian, so the error grew with the relative rotation --
///      exactly where closures matter most.
///
/// The test builds two keyframes whose points are related by a known transform
/// with a large rotation, gives the keyframes deliberately *wrong* poses to
/// stand in for drift, and checks that se2Exp(z) maps the query points onto the
/// stored ones.
void testLoopClosureTransformConvention() {
  backend::LoopCloseOptions opt;
  opt.min_index_gap = 5;
  opt.min_inliers = 6;
  backend::LoopCloser lc(opt);

  // The true relationship between the two views, with a rotation large enough
  // that a raw-transform/tangent mix-up is unmistakable.
  const Pose2 T_rel{3.0, -1.5, 2.6};

  Rng rng(11);
  backend::Keyframe stored;
  stored.id = 0;
  stored.pose = Pose2{0, 0, 0};
  for (int i = 0; i < 40; ++i) {
    backend::MapPoint mp;
    mp.p = Vec3(rng.uniform(1.0, 8.0), rng.uniform(-3.0, 3.0), rng.uniform(-1.0, 1.0));
    for (int w = 0; w < 4; ++w) mp.desc.w[w] = rng.next();
    stored.points.push_back(mp);
  }

  // The query sees the same structure from a different place. A graph edge
  // stores z = Log(x_from^-1 x_to), and the same world point observed from the
  // two cameras gives p_to = (x_from^-1 x_to)^-1 p_from, so the query's points
  // are the stored ones mapped by the inverse of the transform under test.
  backend::Keyframe query;
  query.id = 40;
  // The estimate is far from the truth, standing in for accumulated drift.
  query.pose = Pose2{12.0, 7.0, -1.0};
  const Mat3 T_rel_mat = pose2ToSE2(T_rel);
  const Mat3 T_inv_mat = se2Inverse(T_rel_mat);
  const Eigen::Matrix<Scalar, 2, 2> R2 = T_inv_mat.topLeftCorner<2, 2>();
  const Vec2 t2 = T_inv_mat.topRightCorner<2, 1>();
  for (const backend::MapPoint& mp : stored.points) {
    backend::MapPoint q;
    const Vec2 xy = R2 * mp.p.head<2>() + t2;
    q.p = Vec3(xy.x(), xy.y(), mp.p.z());
    q.desc = mp.desc;
    query.points.push_back(q);
  }

  lc.addKeyframe(stored);
  backend::LoopClosure out;
  LUNA_CHECK(lc.query(query, &out));

  // z must reproduce the relative pose the two views actually have.
  const Mat3 measured = se2Exp(out.z);
  const Vec3 err = se2Log(se2Inverse(T_rel_mat) * measured);
  LUNA_CHECK(std::hypot(err.x(), err.y()) < 0.02);
  LUNA_CHECK(std::abs(err.z()) < 0.02);

  // Guard against the specific regression: using the estimated poses would give
  // a transform nowhere near T_rel.
  const Pose2 est_rel = pose2Relative(stored.pose, query.pose);
  LUNA_CHECK(std::hypot(est_rel.x - T_rel.x, est_rel.y - T_rel.y) > 5.0);
}

void run() {
  testEdgeJacobians();
  testPoseGraphRecovery();
  testRobustness();
  testRansacVerification();
  testAte();
  testRpe();
  testLoopDbBound();
  testLoopClosureMatch();
  testDescriptorRatio();
  testLoopClosureTransformConvention();
}

}  // namespace

LUNA_TEST_MAIN("test_backend")
