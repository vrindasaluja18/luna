# AGENTS.md

Repository notes for future sessions. Keep this factual and current.

## What this repo is

SLAMCore-8G: a CPU-only, no-CUDA C++17 SLAM and navigation stack targeting
4 CPU cores and 8 GB RAM, with no ROS. Eigen is the only external dependency.

## Build, test, bench

```sh
make            # configure + build into build/
make test       # CTest, 7 suites: core, sim, frontend, backend, estimator,
                #   navigation, determinism
make bench      # regenerate docs/benchmark_report.md from real runs
```

The CMake build also works directly: `cmake -S . -B build && cmake --build build -j4`.

Build warnings are treated as real. Keep the build warning-free; there were
`-Wunused-but-set-variable` and `-Wunused-function` warnings at one point that
were genuine leftovers.

## Module status (do not overstate)

Implemented and tested: `core`, `sim`, `frontend`, `backend`, `mapping`,
`topology`, `planning`, `control`, `bench`, `app`.

The pipeline is end to end: the `navigate` CLI command runs SLAM -> map ->
Morse topology -> corridor planning -> MPC control against the simulator's
world, and reports ATE, plan status, tracking RMS, and world collisions. The
closed loop is what `tests/test_navigation.cpp` covers.

Honest limits, measured, not aspirational:

- Navigation succeeds on most worlds and seeds, but not all. On a fragmented or
  thinly-observed map the planner can report `blocked` or `no_corridor_route`
  and drive nothing. That is a real coverage outcome, not a bug: see the maze
  note below. Do not "fix" it by loosening clearance until the plan appears.
- `textureless` worlds produce one Morse node (a featureless image has no
  corners), so there is no topology to plan on. Expected.
- The scripted-route tracker still collides on some seeds (see the Known
  remaining issue section); the bench table reports this rather than hiding it.

`app::runEstimator` (in `src/app/estimator.cpp`, exposed as the `estimate` CLI
command) wires odometry -> pose graph -> loop closure end to end and reports
ATE/RPE against the simulator's ground truth, with dead reckoning printed
alongside for comparison. `tests/test_estimator.cpp` covers it.

Current measured accuracy (simulator, `--type rooms`, 300 s, ~110-180 m):

| estimate | ATE | notes |
|---|---|---|
| dead reckoning | 0.1-0.6 % | odometry only |
| pose graph | 0.1-0.6 % | improves on DR when closures fire |

Loop closures only fire in worlds whose route genuinely returns to a previous
pose. Several worlds do not, at the durations used in CI: check the revisit
distance before treating "0 closures" as a matcher bug.

## Simulator invariants (these were all real bugs)

1. **Route clearance must exceed the vehicle's half-diagonal footprint**
   (`0.5 * hypot(length, track_width)`). Below that, every wall corner is a
   collision. The margin is derived from the vehicle and `VehicleParams`, never
   hard-coded.

2. **Never cap the planning margin by the world's minimum corridor width.**
   That minimum is set by thin slivers between walls, which are free but
   undrivable. Capping to it lets the planner route straight through a slit.
   The preferred margin should be *wide*; `gridPath`'s clearance test does the
   filtering.

3. **Anchor sampling must require drivable clearance.** Sampling anywhere free
   (`sampleFreePoints` with `min_clearance_m = 0`) puts anchors in slivers that
   cannot be connected to anything, collapsing the route to one point. Pass the
   world's clearance field.

4. **Route points must not repeat.** A duplicate creates a zero-length segment;
   the initial heading is taken from the first segment, so a duplicate start
   makes the heading `atan2(0, 0) = 0` and the vehicle leaves in an arbitrary
   direction it can never correct from. `buildRoute` dedupes; keep it.

5. **The tracker must see corners before entering them.** Speed is scheduled
   from curvature previewed over the braking distance plus several lookaheads.
   Sampling curvature only over the lookahead reacts when the corner is already
   underneath the vehicle.

6. **The route cursor search window must be bounded** (currently +1.5 m). A
   route that doubles back can otherwise project the vehicle onto a *later*
   pass that is nearer, and the controller steers toward a point behind it.

7. **Doorways must fit the turning circle, not just the body.** A 90-degree
   turn at the minimum radius sweeps about 1.1 m laterally, so `door_width`
   needs to be well above the body diagonal.

8. **The route follower must hold station at the end of the route.** The route
   is a closed tour, so finishing it means arriving back at the start. Without
   an explicit stop, the arc-length projection clamps to the final waypoint, the
   cross-track error keeps growing as the vehicle drives past it, and the
   steering saturates into a circle that carries the vehicle off the map. This
   destroys the revisit that loop closure is supposed to detect, so it looks
   like a matcher failure.

The tracker is Stanley, not pure pursuit. Pure pursuit converges only
asymptotically and cuts the inside of corners, which in a sub-metre corridor is
the wall. It was replaced; do not reintroduce it.

## Estimator invariants (these were all real bugs)

1. **The odometry increment must be the relative transform between two
   dead-reckoned poses.** Taking the previous pose from the graph and the
   current one from dead reckoning mixes two frames and injects the optimiser's
   own correction as a bogus measurement. The graph then degrades the odometry
   it was given, which is the opposite of its job.

2. **The optimised trajectory must carry keyframe timestamps, one per node and
   in node order.** ATE matches the estimate to ground truth *by time*. Stamping
   nodes with keyframe indices compares the tenth keyframe against the tenth
   second and reports a large error for a perfect estimate.

3. **The simulator must accumulate odometry ticks across the frame interval.**
   Odometry runs at 20 Hz and camera frames at 10 Hz. Handing each frame only
   the last tick's delta loses roughly half the distance and presents as a
   dead-reckoning scale error of almost exactly 2x.

4. **Blend the yaw sources once per window, over the same interval.** Odometry
   yaw is accurate but drifts with wheel scale error; the gyro is drift-free
   over short windows but biased. Both must describe the *same* interval.
   Blending a windowed aggregate against a single sub-tick gyro reading halves
   the rotation.

5. **Noise must scale with the distance travelled.** A constant sigma makes a
   long interval look as trustworthy as a short one, so the graph over-trusts
   the segments where odometry is weakest.

6. **`SparseSystem::solve` must use the accumulated right-hand side.** Solving
   `H dx = 0` converges to the initial guess and hides every other bug in the
   backend. `PoseGraph::evaluate` assembles the residual through `addB`; keep
   that path.

7. **The FAST threshold must suit the image statistics.** The simulator renders
   band-limited procedural texture, which has far less high-frequency energy
   than a real camera. The default threshold tuned for real imagery finds only a
   handful of corners per frame, and a handful of descriptors cannot support
   loop-closure retrieval. `EstimatorOptions::fast_threshold` (default 10) is
   the knob; the failure otherwise looks like a broken matcher.

8. **TUM output must be written at full precision.** The default 6 significant
   digits quantises a 20 m coordinate to ~1e-4 m, the same order as the RPE the
   evaluator resolves, so the file cannot round-trip.

9. **A loop-closure measurement must be the RANSAC-fitted transform, stored as
   an SE(2) tangent.** Copying the current estimates' relative pose makes the
   edge residual zero, so the graph looks like it is working while doing nothing
   and every accepted closure scores as false against ground truth. Storing the
   raw fitted translation instead of `se2Log` of it leaves an error proportional
   to the relative rotation, largest exactly where the closure matters.

10. **A loop closure's information must come from its residual RMS, not the
    RANSAC inlier threshold.** The threshold describes how tightly the points
    were required to agree, not how accurate the fitted transform is. Deriving
    sigma from it (`inlier_dist / sqrt(n)`) makes the closure look ~40x more
    precise than the reconstruction it came from; the graph then over-trusts it
    and drags a trajectory whose odometry was already good *off* course --
    measured at nearly 2x the odometry error. Dividing the RMS by `sqrt(n)` is
    the same mistake for the same reason: the residual is dominated by the
    stereo-lite depth bias, which is systematic across the whole keyframe and
    does not average away. Use the RMS as-is, floored.

11. **The closure acceptance gate must be loose enough for the drift it is
    checking against.** The gate compares a candidate closure with the graph's
    prediction, so a *correct* closure at the end of a long loop disagrees by
    roughly the accumulated drift. A tight bound rejects precisely the closures
    that are most needed and does so silently (they count as rejected, not
    false). Measured on the reference run: true closures disagree by under
    0.2 m, aliased ones by over 0.8 m, so the default sits at 0.6 m.

12. **Robust kernels and loop closures are a poor match.** A Huber kernel
    assumes outliers have large residuals. A loop closure inverts that: its
    residual *is* the drift it exists to remove. The kernel is still applied
    (it bounds the damage from a wrong match) but it is annealed inside
    `PoseGraph::optimize` so it does not suppress the honest correction. A fixed
    small delta measurably suppressed a 4 m correction to a 2.4 % contribution.

13. **Dead reckoning and the pose graph must be scored on the same samples.**
    The graph has one node per keyframe, so the optimised trajectory is
    keyframe-rate; recording the odometry baseline at full camera rate and
    comparing the two compares different sample sets, and the odometry number is
    the one that shifts. `res.odom_trajectory` is therefore appended at keyframe
    instants, next to `kf_times`, not per frame.

## Loop closure: where the accuracy actually comes from

All eight closures on the reference run are correct (100 % precision) and the
graph is still slightly *worse* than dead reckoning there (0.224 m vs 0.165 m
ATE). This is not a bug in the optimiser -- the final solve converges to
`grad = 1.9e-7` -- it is the measurement being the weakest link. The odometry on
that world is accurate to 0.094 % of distance, so the total drift is ~0.2 m,
about the error of one stereo-lite loop measurement.

The consequence for tuning: raising `LoopCloseOptions::min_sigma` above the
measurement's real error makes the graph monotonically approach dead reckoning,
which looks like "fixing" the regression but is really discarding the closures.
Do not tune that floor to make the reference number look better. The floor is set
to 0.10 m, the measured closure error; the drift sweep in the README shows the
closures helping once the odometry is worse than they are. The real fix is
upstream -- a visual-inertial frontend whose features constrain the graph -- not
a reweighting here.

## Navigation invariants (mapping -> topology -> planning -> control)

1. **The free-space stamp radius must exceed the planner's minimum clearance.**
   `stamp_free_radius_m` was 0.45 against `NavigationOptions::plan.min_clearance`
   of 0.46, so the *centreline of the trail the vehicle actually drove* sat below
   the planner's threshold. On a room that is invisible -- stereo depth fills the
   space anyway -- but in a maze, where almost all observed free space *is* the
   trail, the map fragmented into 20 drivable islands and the planner reported
   `blocked` on a corridor the vehicle had just driven. The default is now
   derived (`plan.min_clearance + 0.04`); keep the two tied together.

2. **The corridor graph is a guide, not the only admissible route.** A corridor
   edge is an adjacency in the Morse complex, whose height function is smoothed;
   it can therefore join two regions whose connecting cells are not drivable in
   the raw grid. Planning must fall back to grid A* when the corridor route does
   not exist or does not survive collision verification (`planDirectGrid`), or a
   map the grid can traverse is reported impassable.

3. **A plan that goes nowhere is not a plan.** The goal search over candidate
   nodes accepts the start node itself as a trivially-valid route, so a map with
   no real route reported `ok` with a zero-length path and the controller -- 
   correctly -- refused to move. Reject candidates shorter than a metre.

4. **`graph_components` is a count of components, not the chosen component's
   node count.** The two were conflated; on a fragmented map (many 1-2 node
   components) it overstated connectivity by an order of magnitude.

5. **Report the honest status when there is no reachable corridor node.** Leaving
   the `PlanStatus` default (`kInvalid`) reads as a malformed request; the real
   outcome is `kNoCorridorRoute`.

## Controller sees the SLAM estimate, not the truth

`runNavigation` propagates the controller's pose estimate by the plant's
odometry increments, starting from the plan's own front (which is SLAM's pose
in the world frame). It does **not** copy `vehicle.state().pose` into the
estimate each tick. Two earlier variants were both wrong, and both are easy to
reintroduce:

* Feeding the controller the plant state is perfect state feedback. The
  controller tracks a plan it can see it is exactly on, so localization error
  never reaches it and every run reports ~0 tracking error.
* Starting the plant at the estimate teleports the body onto its own error.
  A 0.9 m estimate error then puts the body 0.9 m from where it really is --
  sometimes inside a wall -- so collisions become a setup artifact.

With the plant at truth and the estimate dead-reckoned from SLAM's start, the
initial separation between them *is* the run's localization error and reaches
the vehicle as a plan the controller believes it follows while the body sits
elsewhere. Measured separation on `maze` seed 6 is 0.131 m; the estimate then
drifts at the odometry rate. Tracking RMS of ~0.19 m is this offset, not a
controller defect.

### Residual `navigate` collisions are localization error, not a planner defect

Sweep of 6 worlds x 8 seeds (48 runs): 22 plan+drive+finish clean at 0.17-0.21 m
RMS, 19 are refused, 7 accept a plan the body then grazes a wall on. The 7 are
not the old "planner says `kOk` for an undrivable plan" bug: the plan's own
clearance is above the body radius, and the body leaves the corridor because the
*world-frame* estimate error (0.13-0.2 m) exceeds the plan's margin (0.46 -
0.353 = 0.107 m). The controller cannot detect this from inside the loop because
its chance constraint uses the graph's relative covariance (0.018 m), which is
bounded by construction.

These runs are reported honestly -- non-zero `COLLISIONS`, non-zero exit code --
and that is deliberate. Raising the acceptance floor to `min_clearance` was
measured: it turns several clean runs into refusals and still leaves collisions
on other seeds, so the floor stays at the physical body radius.

## Known remaining issue

Reference sweep (6 worlds x 16 seeds = 96 runs): 77 collision-free. The
remainder are hairpin reversals at tour anchors, where a car physically cannot
complete the turn. 2-opt was tried to shorten the tour and made it *worse*, so
it was reverted; do not re-add it without measuring.

The route follower also still loses the route at some sharp corners and then
drives straight off it (observed excursions of tens of metres on several seeds).
The vehicle does not crash, but the resulting run is not a useful benchmark: the
sensor stream stops containing new structure. Worlds and seeds that stay on
route (`rooms` seed 1 and 6, `corridor_loop` seed 6 and 8) are the ones used for
estimator validation. Do not benchmark on a run without checking the excursion
first.

### `navigate` on the maze (honest, not fixed)

Most maze seeds report `blocked` and drive nothing, and this is the correct
answer for what the map contains. A 60 s scripted pass through a 40 m maze
observes roughly 7 % of it -- the reconstructed map is a thin trail one cell
wide, and its drivable area fragments into ~20 islands. There is no route across
the maze in that map because the map has not seen the maze, only the path
through it. The topology extractor then has no well-separated maxima to connect,
which is why the corridor graph and the grid agree that no route exists.

Seed 6 is the exception, and it is the useful one: the map *does* support a
9.6 m route, so the planner accepts it and the loop drives 9.3 m -- but the body
grazes a wall on 198 ticks. This is the localization-error-vs-margin case
described in the README, not a maze-specific defect. Do not read a `blocked`
result as a planner defect, and do not close it by lowering `min_clearance` or
growing `free_completion_m` until a path appears -- that manufactures a route
through unobserved space. The fix, if one is wanted, is exploration: a longer or
multi-pass route that observes more of the maze before planning.

## Determinism

A seeded RNG drives world generation, sensor noise, and fault injection. The
same seed gives byte-identical TUM output. Verified by comparing two runs of
`sim --out`. Keep new randomness routed through the seeded `Rng`.

## Style

- Comments explain *why*, not *what*. Several comments in `simulator.cpp` and
  `world.cpp` record why a plausible-looking alternative is wrong; keep that
  kind of note when changing those paths.
- No allocation in real-time loops after init; use the `core/alloc` pools.
- Sparse linear algebra only; no dense ops in the estimator paths.
