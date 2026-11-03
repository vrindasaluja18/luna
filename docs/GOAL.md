# What Luna is for

This file exists because the repository has a *constraint* and no *objective*.
The README states the constraint precisely -- CPU only, no CUDA, 8 GB of RAM,
four cores, Eigen the only dependency -- and then lists what is implemented.
Nowhere does it say what the thing is *for*. A system with a constraint and no
objective can only be scored against the constraint, which is why every
conversation about Luna so far has been a conversation about memory and never a
conversation about results. This file fixes that.

## The goal

> **Hold a map's topology at a fidelity proportional to its persistence, inside
> a hard 8 GB budget, and prove by measurement -- not by claim -- that every
> byte held and every quantity computed changed a decision.**

Three clauses, and each one is a separate obligation.

1. **Topology at fidelity proportional to persistence.** The map is not a grid
   of numbers to be stored. It is a *shape*: rooms, doorways, loops, dead ends.
   Persistence is the number that says how real each feature of that shape is.
   The goal is to spend memory in proportion to that number and no more.

2. **Inside a hard 8 GB budget.** Not "fits comfortably". Hard. The budget is
   the thing that forces the fidelity decision to be made rather than avoided.
   A system with unlimited memory never has to decide what a room is worth.

3. **Prove by measurement that every byte and every computation changed a
   decision.** This is the clause that separates a project from a demonstration.
   A buffer that nothing reads is a cost with no benefit. A quantity computed
   and discarded is compute spent for nothing. Both are measurable and both are
   currently unmeasured.

SLAM is the *vehicle* for that goal, not the goal. SLAM is used because it
produces exactly the object the goal is about: a map whose topology is unknown
in advance, changes as it is observed, and must be held under a fixed budget.

## Why this goal and not another

The alternative goals are all either already solved or unmeasurable.

| candidate goal | why it is not the goal |
|---|---|
| "accurate SLAM" | solved to the limit of the sensors; Luna's ATE is bounded by its odometry, not its algorithm |
| "real-time on a robot" | a constraint, already met; it does not tell you what to spend memory on |
| "a product" | needs users, has none, and would make the honest-numbers discipline a liability |
| "a deep-learning SLAM" | would add a dependency and a failure mode; the interesting part is the memory, not the network |
| "an RL agent" | there is no action space yet, so there is nothing to learn a policy *over* |

The chosen goal is the only one that is (a) not already solved, (b) measurable
with data the repository already produces, and (c) honest about what is hard.

## Sub-goals, each with a measurement

| id | sub-goal | how it is measured | state |
|---|---|---|---|
| **G1** | run inside 8 GB / 4 cores / no GPU | peak RSS on the longest benchmark run | **met** |
| **G2** | maximise structure retained per byte held | `graph_components` and corridor count vs map bytes held (today only `map_bytes_resident` is reported; a peak field is a Phase-0 addition) | **not measured** |
| **G3** | one currency for clustering, eviction and reward | persistence is the currency in all three paths | **partial** -- clustering yes, eviction no, reward no |
| **G4** | let the map's own topology choose where to look next | exploration reward vs coverage, on recorded runs | **not present** |

G2 is the objective. G1 is only the constraint. Today the repository reports G1
and *asserts* G2 without a number, which is the same failure mode the README is
otherwise careful to avoid.

## Non-goals

* Not a product, not a service, not a benchmark entrant.
* Not a deep-learning system. There is no autograd tape, no parameters, no
  loss, and none is planned (see `docs/blog/01-liquid-memory.md`, section
  "Three flows").
* Not real-time guarantees beyond the simulator's own clock.
* Not a general SLAM library. It is a testbed for one idea: topological memory
  management under a hard budget.

## What "who does what" means when there is no learning

The honest statement: **Luna contains no reinforcement learning.** There is no
policy, no value function, no action selection, no reward. The roles one would
name in an RL system are all filled by *classical* components, and the mapping
is worth writing down because it shows exactly which roles are missing.

| RL role | what fills it in Luna today | missing? |
|---|---|---|
| environment | `sim` -- world generator, sensors, fault injection | no |
| state | pose (SE(2)) + occupancy map | no |
| observation | camera / IMU / odometry frames | no |
| reward | *nothing* -- there is no scalar objective in the loop | **yes** |
| policy | *nothing* -- the route follower executes a fixed scripted path | **yes** |
| action space | *nothing* -- the vehicle cannot choose where to go | **yes** |
| value / critic | *nothing* | **yes** |
| model update | pose-graph optimisation (least squares, not gradient ascent) | no |

The four missing rows are the whole of "RL in Luna", and they are missing for
one reason: **the route follower is scripted, so the sensor stream does not
depend on any decision.** There is nothing to optimise a policy against. This
is not an oversight to patch; it is the boundary of the current design.

Making the reward real is possible today (it is a pure function of a recorded
run). Making the *policy* real requires opening the action space, which is the
last phase, not the first. See the plan at the end of the blog post.

## Credit

The framing -- topological memory, persistence as the currency, islands and
their connections, the classical-systems discipline -- is not original to this
project. It comes from a teacher whose public work is cited in
`docs/blog/01-liquid-memory.md`. The contribution here is the *application*:
showing that a CPU SLAM stack already computes persistence-based clustering and
can be made to spend its memory on the same quantity.
