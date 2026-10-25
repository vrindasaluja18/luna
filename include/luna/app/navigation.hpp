// LUNA -- the end-to-end navigation run: SLAM -> map -> plan -> control.
//
// This is the file that makes the project's claim checkable. Every other module
// can be tested in isolation and still fail to compose; the interesting failure
// modes only appear when they are chained, and they are mostly *frame* errors
// rather than algorithmic ones:
//
//   * The SLAM map is not the simulator's world. It is a probabilistic
//     reconstruction with unknown cells where nothing was observed. Planning on
//     it and colliding against the world measures the whole chain, including
//     the map, which is the only honest way to score it.
//   * The plan is expressed in the SLAM frame. The controller drives the
//     vehicle in the world frame. The transform between them is the estimator's
//     own starting offset, and getting it wrong produces a plan that is
//     perfectly drivable and in the wrong place.
//
// The run therefore does exactly this: estimate, build the map from the
// estimate, extract topology from the map, plan on the map, drive the real
// plant, and score collisions against the real world. Nothing in the plan or
// the control path sees ground truth.
#pragma once

#include <string>
#include <vector>

#include "luna/app/estimator.hpp"
#include "luna/control/mpc.hpp"
#include "luna/mapping/occupancy.hpp"
#include "luna/planning/planner.hpp"
#include "luna/sim/world.hpp"
#include "luna/topology/morse.hpp"

namespace luna {
namespace app {

struct NavigationOptions {
  EstimatorOptions estimator;
  planning::PlanOptions plan = [] {
    planning::PlanOptions p;
    // Looser than the analytic-world defaults. A SLAM map's doorways are
    // reconstructed, so they come out narrower than the walls that produced
    // them: the observed gap between two wall faces is smaller than the real
    // one by roughly the map's own resolution. Demanding the ground-truth
    // margin from a reconstructed map rejects every doorway and reports the
    // route as blocked, which reads as a planner failure and is a map
    // resolution effect.
    //
    // "Looser" still has a floor, and it is set by the controller rather than
    // by the body: the plan constrains the body *centre*, and the controller
    // tracks it with a finite error, so the margin has to cover the
    // half-diagonal plus that error. Where the reconstructed map cannot supply
    // that, the planner must not silently plan an optimistic path -- the
    // controller's speed profile carries the same budget and slows the vehicle
    // there instead. Raising this floor to 0.62 does not help: it blocks every
    // doorway in the reconstructed rooms and saddle-hall maps and the route
    // comes back infeasible with zero ticks driven.
    p.min_clearance = 0.46;
    p.preferred_clearance = 0.62;
    return p;
  }();
  control::ControllerOptions control;
  topology::MorseOptions morse;
  mapping::OccupancyOptions map;
  /// Cells below this occupancy probability are treated as obstacles. The map
  /// is probabilistic, so "occupied" is a threshold and the planner's behaviour
  /// depends on it: too low and sensor noise becomes walls, too high and real
  /// walls disappear.
  Scalar occupancy_threshold = 0.60;
  /// Treat unobserved cells as obstacles when planning. See the note in
  /// `buildMap`: this is what keeps the plan inside the region the sensor has
  /// actually seen.
  bool unknown_is_obstacle = true;
  /// Inflate obstacles before planning, metres.
  ///
  /// Kept at zero by default, and the reason is a double-counting trap. The
  /// planner already refuses corridors whose clearance is below its own
  /// `min_clearance`, which is set from the vehicle's footprint; the topology
  /// extractor already has its own `min_clearance` for deciding what counts as
  /// drivable. Inflating the grid as well shrinks every doorway a second time,
  /// and the failure is not graceful: corridors that the topology graph still
  /// believes in become narrower than the map's free space, so a graph edge
  /// connects two regions with no traversable cells between them, and the
  /// planner reports `no_corridor_route` for a route its own graph says exists.
  ///
  /// A small value is still useful when the map is noisier than the planner's
  /// margins assume. It is a knob, not a default.
  Scalar map_inflation_m = 0.0;
  /// Occupied connected components smaller than this many cells are discarded
  /// as sensor noise. A wall is a long run of cells; a stereo depth outlier is
  /// an isolated one. See `removeObstacleSpeckles`.
  Size min_obstacle_component_cells = 12;
  /// Radius around each driven pose that is stamped free, metres. The sensor
  /// model drops returns closer than its own minimum range, so without this the
  /// vehicle's own position is unknown space. See `stampTrajectoryFree`.
  ///
  /// Derived from the plan's minimum clearance rather than fixed, because the
  /// two are not independent: the stamp is what makes the corridor the vehicle
  /// drove drivable, and a stamp narrower than the clearance the planner
  /// demands leaves the *centreline of its own trail* below the threshold. The
  /// value was 0.45 against a `min_clearance` of 0.46 -- a 1 cm shortfall that
  /// is invisible on a wide room, where depth returns fill the space anyway,
  /// and decisive in a maze, where the observed free space is almost entirely
  /// the trail. There the map fragmented into 20 drivable islands, the planner
  /// reported `blocked`, and no route existed even though the vehicle had just
  /// driven one. The margin above the floor accounts for the pose estimate's
  /// own error, so a stamped disc is drivable at the *estimated* pose.
  Scalar stamp_free_radius_m = plan.min_clearance + 0.04;
  /// Radius over which observed free space is grown into unobserved cells.
  /// Stereo depth marks free space only along the beams that returned a depth,
  /// leaving the gaps between beams unknown; without this the map's free space
  /// is a fan of thin streaks and no plan can follow it. See the note in
  /// `buildMap`.
  Scalar free_completion_m = 0.30;
  /// Drive the closed loop after planning.
  bool drive = true;
  /// Stop after this many simulated seconds of driving.
  Scalar drive_duration = 120.0;
};

struct NavigationStats {
  // SLAM
  Size keyframes = 0;
  Scalar slam_ate = 0;
  Scalar slam_ate_percent = 0;
  Size loop_closures = 0;
  Scalar loop_precision = 0;
  /// Position and heading standard deviation the pose graph reports for its
  /// final pose, and whether it was able to compute them at all.
  Scalar slam_sigma_xy = 0;
  Scalar slam_sigma_th = 0;
  bool slam_covariance_ok = false;
  /// Position sigma of the final pose's *marginal* covariance, reported for
  /// comparison with the local figure the controller used. The two differ by an
  /// order of magnitude on a long run and the gap is the point.
  Scalar slam_marginal_sigma_xy = 0;
  // Map
  Size map_occupied_cells = 0;
  Scalar map_known_fraction = 0;
  Size map_tiles_resident = 0;
  Size map_bytes_resident = 0;
  // Topology
  int morse_nodes = 0;
  int morse_edges = 0;
  int maxima = 0;
  int saddles = 0;
  /// Morse nodes that lie in the vehicle's own free-space component. Nodes
  /// outside it are unreachable given what the map has observed, and a plan to
  /// one of them is a plan through unknown space.
  int reachable_nodes = 0;
  /// Number of connected components in the corridor graph.
  int graph_components = 0;
  Scalar topology_ms = 0;
  // Plan
  bool plan_ok = false;
  /// Why the planner returned what it did, for diagnosis.
  planning::PlanStatus plan_status = planning::PlanStatus::kInvalid;
  Scalar plan_length = 0;
  Scalar plan_clearance = 0;
  Scalar plan_ms = 0;
  // Control
  /// Final control status, for diagnosis. A safe stop has several causes and
  /// they call for different fixes; reporting only "safe stop" hides which one.
  control::ControlStatus control_status = control::ControlStatus::kNoPlan;
  Scalar control_speed_cap = 0;
  Scalar control_det_sigma = 0;
  Scalar control_inflation = 0;
  Scalar plan_route_bottleneck = 0;
  Scalar control_sigma_xy = 0;
  Scalar control_sigma_th = 0;
  control::StopReason control_stop_reason = control::StopReason::kNone;
  Scalar drive_distance = 0;
  Scalar tracking_rms = 0;
  Scalar tracking_max = 0;
  Size collisions = 0;
  bool safe_stopped = false;
  Scalar control_ms_p99 = 0;
  Scalar control_ms_max = 0;
  Size control_ticks = 0;
  Scalar goal_error = 0;
};

struct NavigationResult {
  NavigationStats stats;
  /// The occupancy grid the planner actually consumed: SLAM map, thresholded,
  /// unknown space resolved, and inflated. Kept so a failed plan can be
  /// diagnosed from the map that caused it rather than from the trajectory it
  /// produced.
  sim::GridMap planning_grid;
  /// Planned path in world coordinates, for plotting.
  std::vector<Vec2> planned_path;
  /// The driven trajectory in world coordinates.
  std::vector<Vec2> driven_path;
  bool ok = false;
};

/// Runs the whole stack over one simulated session.
NavigationResult runNavigation(const sim::SimConfig& cfg,
                               const NavigationOptions& opt = NavigationOptions());

}  // namespace app
}  // namespace luna
