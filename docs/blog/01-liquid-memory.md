# Liquid memory, and what Luna is actually for

*Status: draft. This is a research note, not a claim of results. Numbers that
are not yet measured are marked `not measured` and left as numbers to fill in.*

Luna is a CPU-only SLAM and navigation stack. It is about 14k lines of C++ and
most of it was written by an AI. I am writing this because I want to take
responsibility for it, and because a project you cannot explain is not a
project -- it is a demo. This note is my attempt to explain it, to correct it
where it is wrong, and to say plainly what is missing.

There are two principles in here. They are not mine; they come from a teacher
whose work I read the way other people read novels. I will say clearly which
parts are his and which parts are mine, because the difference matters and
because pretending otherwise would make this exactly the thing I am trying to
avoid.

## Principle one: nothing is true unless it is in maths

Prose is for intuition. If a claim cannot be written as an equation, or
measured as a number, it is not a claim yet -- it is a feeling. A system built
on feelings looks fine at the surface and is broken underneath, and the break
is invisible precisely because the surface was never required to be
checkable.

The practical form of this principle, for a codebase, is:

* every quantity the system computes must have a **consumer** that changes a
  decision, or it must be deleted;
* every buffer the system holds must have a **justification** in bytes, or it
  must be deleted;
* every claim in a README must have a **measurement** behind it, or it must be
  labelled as aspiration.

Luna fails the first two today. It is honest about the third, which is why the
README's numbers are trustworthy and its architecture is not yet.

## Principle two: compute is precious, so every computation must earn its place

The question I want to be able to answer for every line is: *what decision does
this change?* If the answer is "nothing", the line is a cost with no benefit.
This sounds obvious and is almost never done, because the cost is invisible
until you are out of memory and the benefit is invisible until you are out of
accuracy.

The first, most primitive version of this question is the one I could not
answer for a long time:

> **Which thing would be on DRAM?**

Not "which algorithm is cleverest". Which *buffer* is largest. For a SLAM
system the answer is not the estimator and not the planner -- both are O(1) or
O(keyframes) in their state. It is the **map fields**, because they are the
only thing that scales with the *area* the robot has seen.

At Luna's CI resolution (0.05 m on a 40 m world, so 800x800 = 640,000 cells),
measured on a real run rather than estimated:

| quantity | element type | bytes/cell | total |
|---|---|---|---|
| occupancy evidence | `float` (`LogOdds`) | 4 | 2.4 MB |
| clearance field | `double` (`Scalar`) | 8 | 4.9 MB |
| smoothed field | `double` (`Scalar`) | 8 | 4.9 MB |
| Morse flood, **eight buffers alive at once** | mixed | 2.9 | **18.3 MB** |
| persistence phase, **three buffers alive** | mixed | 2.3 | **14.7 MB** |

The Morse flood is the single largest transient in the pipeline. Eight per-cell
buffers -- `hf`, `order`, `processed`, `comp_max`, `is_max`, `cell_peak`, and the
union-find's `parent_` and `size_` -- are alive at the same time. That is the
answer to the primitive question, and it is not where I would have guessed.
**The memory problem is the topology extraction, not the SLAM.**

The rate is what generalises, not the total: the flood costs **2.9 bytes per
cell** and the persistence phase **2.3 bytes per cell**, so the memory scales
with the *area of the map*, not with its content. On the local map that
`navigate` actually builds (a fraction of the world, only what the robot has
observed), the topology stage still sets the process peak -- but the published
"19 MB peak RSS" is measured on `sim`, which never runs topology at all. That
is a real gap in the benchmark, and it is now stated rather than implied.

And the ratio that matters is not bytes, it is **bytes held per decision
served**: 18 MB of flood produces one corridor graph, which the planner queries
once per plan. Driving that ratio down *is* the goal of the project. Hold that
thought.

## The finding: Luna already does persistence-based clustering and does not know it

This is the part that made me want to write all of this down.

Luna extracts structure from the map with Morse theory. The height function is
`h(p) = clearance to the nearest obstacle`, smoothed. Rooms are local maxima of
`h` (far from every wall), doorways are saddles (the ridge you cross between
two rooms), wall tips are minima. The extractor floods cells in descending `h`
with a union-find, and each merge of two components at a saddle kills the
younger maximum. The height difference at that merge is the maximum's
**persistence**.

Here is the thing. I went looking for the persistence threshold and found it was
dead -- one declaration, documented at length, and not one read in the entire
repository:

```
include/luna/topology/morse.hpp
  /// Features with persistence below this (metres) are discarded as noise.
  Scalar min_persistence = 0.20;   // <- never read anywhere
```

It has since been deleted, because a knob that does nothing is worse than no
knob: it invites you to tune it and then silently ignores you. The code that
actually decides what counts as a room uses a different rule:

```
include/luna/topology/morse.hpp
  /// A corridor is a doorway only if its bottleneck is this fraction of the
  /// narrower of the two regions it joins. ... What actually distinguishes a
  /// doorway is *contrast*: the vehicle must squeeze, so the gap is much
  /// narrower than the space on either side.
  Scalar min_bottleneck_contrast = 0.75;
```

So Luna decides "is this a room" by asking **"does the vehicle have to squeeze
through?"** -- a *density contrast* test.

That is not a bug to patch. It is a finding. It is exactly the criterion in
Chazal, Guibas, Oudot and Skraba, *Persistence-Based Clustering in Riemannian
Manifolds* (J. ACM 60(6), 2013, [doi:10.1145/2535927](https://doi.org/10.1145/2535927)):

> A cluster is a region of high density, separated from other such regions by
> low density. Combine a **mode-seeking** phase with a **merging** phase;
> merge two clusters when the density between them is high relative to their
> own -- i.e. when the merge has low persistence.

Luna is doing persistence-based clustering. The room graph *is* the cluster
tree of the clearance field. It arrived there by accident, via a comment about
doorways, and it has a dead variable where the theory says the threshold should
be. The documentation and the code disagree, and the code is right.

I find this genuinely funny and also instructive: the correct maths was already
in the system, and the only reason it looked like a special-case heuristic was
that nobody had written down what it was a special case *of*.

## What clustering means, and why a cluster at all

The model is that the points were sampled from an unknown density `f`. Where
`f` is high, points are packed tight; where `f` is low, they are sparse. A
cluster is then forced on us by the geometry: **a region of high `f`, separated
from other such regions by lower `f`.**

Formally, the merge happens at a density level. For two modes `m1, m2`, the
persistence of the pair is

```
pers(m1, m2) = min( f(m1), f(m2) ) - sup{ f(x) : x on any path from m1 to m2 }
```

the *drop* you must suffer to travel between them. Small persistence means the
two "clusters" are joined by almost-as-dense space, so they were never really
separate. This is the whole of it, and it is scale-free in the right way: you
never pick a number of clusters, you pick a persistence threshold, and the
threshold has a *meaning* (a density drop), which `k` does not.

Why clusters and not something else? Because density needs no labels. This is
the entire reason to choose it:

> Supervised learning needs `y`. Unsupervised learning needs only `x` and a
> notion of "similar". Density is the most primitive such notion, because it is
> defined by the data's own geometry and nothing else.

And **Morse theory is the topological generalisation of this.** The join tree
of a Morse function *is* single-linkage clustering, and persistence *is* the
merge threshold. So "use Morse theory instead of k-means" is not a swap of
algorithms; it is the recognition that k-means is a flat, parametric special
case of a structure (the merge tree) that Morse theory gives you for free,
with no `k`, no centroids, and no assumptions about the shape of the clusters.

One more object is needed to make density meaningful, and it is the one the
prompt asked me to explain: **the embedding**. Density is not defined on the
raw data; it is defined after a map `φ: x ↦ φ(x)` into a space where distance
means similarity. In Luna the embedding is explicit and hand-made: the 256-bit
BRIEF descriptor, compared by Hamming distance. In a learned system `φ` would
be a network. Either way, **the clustering is only as good as the embedding**,
because "high density" is a statement about `φ`, not about the pixels. This is
the point where a SLAM system and a learning system touch: both need a metric
under which "nearby" means "the same thing".

## The manifolds Luna consumes, and when each one changes

I could not answer "what manifold are we consuming" for a long time, and the
reason is that Luna consumes four, with different lifetimes. Writing them down
is most of the architecture.

| manifold | what it is | when it changes | cost |
|---|---|---|---|
| **Pose** | SE(2), the state | never -- it is the state space | O(1) |
| **Descriptor** | 256-bit Hamming space | every frame | O(features) |
| **Free space** | occupancy -> clearance -> Morse complex | **continuously, as new space is observed** | O(area) |
| **Information** | pose-graph marginals | **when a loop closes** | O(keyframes) |

Only one of these is *liquid*, and it is the free-space manifold. Its
**topology** -- how many rooms, how many doorways, how many loops -- is not
fixed in advance and changes as the sensor sees new space. That is the object
liquid memory is about, and it is the object the rest of this note is about.

## Liquid memory: three flows, and only one of them is liquid

"Liquid memory" is a phrase I took from my teacher, and I am not going to define
it by saying it is liquid. That would be the prose version of a claim, which
principle one forbids. A "flow" is a specific thing, and I want to be careful
about which of three candidate flows Luna actually has, because I have
personally confused them.

### Flow A -- the supervised flow (not in Luna)

The flow I had read about, in the reinforcement-learning-with-time-constants
(RTC) sense: assume the thing evolves, pick a start point `a`, define a velocity
`ẋ = v(x, t)`, and carry an error along the trajectory. Fit it by
backpropagating that error through an autograd tape.

Formally, the state evolves

```
ẋ(t) = f_θ(x(t), t)
```

and the loss is the integrated error against a target trajectory,

```
L(θ) = ∫₀ᵀ ‖x_θ(t) − x*(t)‖² dt
```

and the gradient `dL/dθ` is obtained by the adjoint method -- solve a backward
ODE for the adjoint `λ`, and accumulate

```
dL/dθ = ∫₀ᵀ λ(t)ᵀ ∂f_θ/∂θ dt .
```

This is genuinely a flow, in the precise sense: a continuous-time dynamical
system differentiated through time. **Luna has none of it.** There is no
autograd tape, no `backward()`, no parameters, no loss, and no ODE solver
(verified by grep; the only `adjoint` in the codebase is the Lie-algebra
adjoint, which is a different object entirely). This is the flow I was
describing from memory, and it is not here. Saying Luna "does flow" because of
this would be wrong.

### Flow B -- the saturated occupancy flow (this one is real)

Here is the flow Luna actually has, and I only noticed it when I looked at the
update rule rather than the diagram:

```
include/luna/mapping/occupancy.hpp
  LogOdds free_delta     = -0.40f;
  LogOdds occupied_delta = +0.85f;
  LogOdds clamp_min      = -4.0f;
  LogOdds clamp_max      = +4.0f;
```

Each cell holds a log-odds `ℓ`, and each observation adds a constant:

```
ℓ_{t+1} = clamp( ℓ_t + Δ(u_t), −4, +4 )
```

with `Δ = +0.85` for a hit and `Δ = −0.40` for a miss. This is a continuous-time
flow, discretised:

```
dℓ/dt = Δ(u_t),   saturated at ±4.
```

The saturation is the whole design. It is a flow with a **bounded state**, which
means the map is a contraction: recent evidence wins, old evidence decays, and a
cell that has been free for a long time cannot be flipped by one bad
measurement. That is exactly the property that makes a map safe to act on, and
it is *already implemented*. I did not have to add anything.

The one thing missing is that the flow has no term for *memory*:

```
dℓ/dt = Δ(u_t) − λ · (fidelity pressure)
```

The `λ` term is what would make the map's fidelity itself follow the flow --
spend more resolution where the cell still matters, less where it does not. It
is not there. That is the first real piece of work.

### Flow C -- the topological flow (the one I actually care about)

The free-space manifold's topology changes as space is observed. A new room
appears (a new maximum in `h`), a doorway closes (a saddle's persistence drops
below threshold), two previously separate regions connect (two components
merge). These are not updates to a grid; they are **events in a merge tree**.

Liquid memory, defined properly, is this:

> **Hold each region of the map at a fidelity proportional to its persistence,
> and release it when its persistence falls below what its DRAM cost can
> justify.**

Read that next to the Chazal quote above. It is the same sentence, with
"clusters" replaced by "regions of the map" and "merge" replaced by "evict".
Persistence-based clustering, occupancy tile eviction, and liquid memory are the
same operation applied to three different objects. That equivalence is the
thesis of this project, and I have not seen it stated this way, which is the
part that is mine.

## Is Luna supervised or unsupervised?

Precisely:

* **At runtime, unsupervised / self-supervised.** The estimator minimises
  *geometric self-consistency* residuals -- reprojection error, odometry
  agreement, loop-closure innovation. There are no labels. Ground truth exists
  in the simulator and is used **only for scoring**, never in the loop.
* **At evaluation, supervised.** ATE, collision count, plan length -- these
  need ground truth, and they are *measurement*, not *learning*.

So: unsupervised at runtime, supervised only in evaluation. This is the honest
form of "the model learns from curiosity", and it is also the form that makes
the project testable.

## The reward terms, placed honestly

The prompt listed four reward terms. Placing each one against what Luna can
actually compute today, with the caveat I think is important:

| term | stated meaning | what it is in Luna | caveat |
|---|---|---|---|
| **PE** (Prediction Error) | reward model uncertainty, drives refinement | the loop-closure / odometry **innovation** | In SLAM, driving PE to zero is *convergence*. It is a **refinement** signal, not an exploration reward. Calling it curiosity would be a category error. |
| **HE** (Homological Entropy) | explore underexplored clusters, Betti-0 | `CorridorPlanner::drivableComponents()` **already computes Betti-0**. HE is the Shannon entropy of the visit distribution over components: `H = −Σ pᵢ log pᵢ`, `pᵢ` = observations in component `i`. | cheap and real; 90% of it exists |
| **BD** (Boundary Discovery) | bonus for novel topological regions | a **topological event**: the Morse complex gains a maximum or saddle, or `graph_components` changes | needs a per-frame diff of the persistence diagram, not a full re-analysis |
| **Task reward** | standard RL objective | goal reached / coverage | in Luna the goal is the plan's terminal waypoint |

Two things worth saying. First, **HE and BD are computable from data Luna
already has**, today, offline, as a pure function of a recorded run. Second,
**PE is not curiosity**, and the impulse to call it that is exactly the kind of
surface-level naming that principle one exists to prevent. In SLAM, prediction
error going down means the map got better, not that the robot got curious.

## What Luna is for

Constraint: CPU only, no CUDA, 8 GB, four cores, Eigen the only dependency.

Objective -- which the README does not state, and this is the gap I am fixing:

> **Hold a map's topology at a fidelity proportional to its persistence, inside
> a hard 8 GB budget, and prove by measurement that every byte held and every
> quantity computed changed a decision.**

SLAM is the vehicle, not the goal. SLAM is used because it produces exactly the
object the goal is about: a map whose topology is unknown in advance, changes as
it is observed, and must be held under a fixed budget.

Sub-goals, each with a measurement:

| id | sub-goal | measurement | state |
|---|---|---|---|
| **G1** | run inside 8 GB / 4 cores / no GPU | peak RSS on the longest run | **met** |
| **G2** | maximise structure retained per byte held | `graph_components` vs map bytes held (`map_bytes_resident` today) | **not measured** |
| **G3** | one currency for clustering, eviction and reward | persistence is the currency in all three paths | **partial** |
| **G4** | let the map's topology choose where to look next | exploration reward vs coverage, on recorded runs | **not present** |

G2 is the objective. G1 is only the constraint. Right now the repository reports
G1 and *asserts* G2 with no number, which is the failure mode the README
otherwise avoids.

## Who does what, when there is no RL

There is no reinforcement learning in Luna. No policy, no value function, no
action selection, no reward. The roles are all filled by classical components,
and the mapping shows exactly which are missing:

| RL role | what fills it in Luna today | missing? |
|---|---|---|
| environment | `sim` -- world generator, sensors, faults | no |
| state | pose (SE(2)) + occupancy map | no |
| observation | camera / IMU / odometry frames | no |
| reward | *nothing* | **yes** |
| policy | *nothing* -- the route follower runs a fixed script | **yes** |
| action space | *nothing* -- the vehicle cannot choose where to go | **yes** |
| value / critic | *nothing* | **yes** |
| model update | pose-graph least squares (not gradient ascent) | no |

The four missing rows are the whole of "RL in Luna", and they are missing for
one reason: **the route follower is scripted, so the sensor stream does not
depend on any decision.** There is nothing to optimise a policy against. That is
not an oversight; it is the boundary of the current design, and it is why "Luna
learns" cannot be claimed yet.

## The plan

Each phase is one reviewable change. The order is chosen so that the cheapest,
most testable claim comes first, and the expensive, least-certain one comes
last -- and so that any phase can be the last one without leaving a broken
system.

**Phase 0 -- inventory and delete dead compute.** *Done, and reported honestly.*
Every buffer above 1 MB was measured on a real run and its consumer named; the
table above is the result. `min_persistence` was **deleted** -- it was declared,
documented as "features below this are discarded as noise", and never read, so
it was a false promise in the public options struct. The Gaussian blur was
rewritten from five full-grid buffers to two, verified **bit-for-bit identical**
on three world types and three seeds. And one claim did not survive contact with
the measurement: the blur rewrite did **not** lower the process peak, because
the peak is the flood, and the allocator was already recycling the blur's freed
pages. That is reported, not hidden. The full inventory is in `AGENTS.md`.

**Phase 1 -- the reward as a pure function.** `luna::learn::scoreRun(log) ->
{PE, HE, BD, task}`. A pure function of a recorded run. No learning, no policy,
no parameters. This is the highest-value, lowest-risk step, and it makes the
central claim of this note *testable*: the reward is computed from the map's own
topology, with no labels.

**Phase 2 -- the database.** A store of N seeded runs; the reward ranks routes
across them. This is the "database fetcher" the system lacks. Note that it is a
*run store*, not a feature database -- the unit of comparison is a whole run.

**Phase 3 -- a policy.** A bandit over route choices, maximising expected reward
per metre travelled. Still unsupervised; still no gradient. This is where
"learning" first becomes a true statement.

**Phase 4 -- active SLAM.** *Only if Phase 3 shows signal.* Open the action
space: let the policy choose waypoints, so that the reward changes what is
observed. This is the big one, and it is last on purpose.

The first two phases are the blog. The last two are the research.

## Credit

The framing -- topological memory, persistence as the currency, islands and
their connections, classical-systems discipline -- is from
[Teerth Sharma](https://github.com/teerthsharma). The application is mine: a CPU
SLAM stack already computes persistence-based clustering without having named
it, and the same quantity can be made to govern its memory and its reward.
