# [Routing four-vectors](@id kinematic_tasks)

```@meta
CurrentModule = CascadeDecays
EditURL = "../src/kinematic-task.md"
```

[`KinematicTask`](@ref) is a reusable specification for turning one event,
represented by final-state four-vectors, into topology-local
[`CascadeKinematics`](@ref). It stores the topologies, the initial-frame
convention, the generated angle-measurement programs, and optional helicity-axis
alignment paths.

## Setup

This page uses a three-body example with physical masses only to get realistic
four-vectors.

```@example kin_task
using FourVectors
using InstructionalDecayTrees
using ThreeBodyDecays:
    ThreeBodyMasses,
    aligned_four_vectors,
    cosθij,
    x2σs

using CascadeDecays

fourvector(p) = FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))

ms = ThreeBodyMasses(0.93827208816, 0.493677, 0.13957039; m0 = 2.28646)
σs = x2σs([0.42, 0.31], ms; k = 3);
```

## From topologies to a task

A task can be generated directly from topologies. The first topology is used as
the reference unless `reference_topology` is supplied explicitly.

```@example kin_task
topologies = (
    DecayTopology(((1, 2), 3)),
    DecayTopology(((3, 1), 2)),
    DecayTopology(((2, 3), 1)),
)

task = KinematicTask(
    topologies;
    reference_topology = topologies[1],
    wigner_finals = (1, 3),
    initial_frame = CurrentFrame(),
)

bracket_notation.(task.topologies)
```

The `wigner_finals` option is only needed for final-state particles with spin.
For those particles, relative Wigner alignment angles are part of the evaluated
kinematics and are stored in the [`KinematicPoint`](@ref).

There are three kinematic objects in this workflow:

- [`KinematicTask`](@ref) is the reusable specification generated from the
  topologies. It owns the topology list, initial-frame convention, generated
  vertex programs, and requested Wigner-alignment paths.
- [`CascadeKinematics`](@ref) is the evaluated kinematics for one topology and
  one event. It stores the line masses squared and the local helicity angles
  used by amplitude evaluation.
- [`KinematicPoint`](@ref) is the evaluated event for the whole task. It stores
  one `CascadeKinematics` object per topology, plus the requested relative
  Wigner rotations.

The generated programs are ordinary `InstructionalDecayTrees.jl` instruction
tuples. They are stored per topology and per vertex.

```@example kin_task
task.programs[1].vertex_programs
```

The public path helper exposes the frame path for one final-state particle.
Internally that particle is a topology line, but the API takes the particle
index. Comparing these paths is exactly how the relative Wigner alignment angles
are obtained.

```@example kin_task
particle_index = 1
(
    reference = helicity_frame_path(topologies[1], particle_index; initial_frame = CurrentFrame()),
    target = helicity_frame_path(topologies[2], particle_index; initial_frame = CurrentFrame()),
)
```

## Current-frame four-vectors

Use [`CurrentFrame`](@ref) when the event is already expressed in the desired
mother-system axes. In this example the aligned event is rotated as
`Rz |> Ry |> Rz`, and the root vertex angles recover the last two rotations.

```@example kin_task
aligned = Tuple(fourvector(p) for p in aligned_four_vectors(σs, ms; k = 3))
current_event = Tuple(p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4) for p in aligned)

point = KinematicPoint(task, current_event)
x12_3 = kinematics_at(point, topologies[1])

(
    root = vertex_angles(x12_3, topologies[1], ((1, 2), 3)),
    isobar = vertex_angles(x12_3, topologies[1], (1, 2)),
    σ12 = line_invariant(x12_3, topologies[1], (1, 2)),
)
```

The `point` object is the evaluated event in the task convention. It stores one
[`CascadeKinematics`](@ref) object per topology and one set of requested
relative Wigner alignment angles per topology. [`kinematics_at`](@ref) is only a
retrieval convenience: it selects the `CascadeKinematics` entry associated with
the requested topology.

```@example kin_task
x31_2 = kinematics_at(point, topologies[2])
(
    topology = bracket_notation(topologies[2]),
    root_angles = vertex_angles(x31_2, topologies[2], ((3, 1), 2)),
    isobar_masses2 = vertex_masses2(x31_2, topologies[2], (3, 1)),
)
```

## Helicity-root four-vectors

Use [`HelicityRootFrame`](@ref) for fully general four-vectors. The event is
first rotated in the mother rest frame and then boosted and reoriented as
`Rz |> Ry |> Rz |> Bz |> Ry |> Rz`.

```@example kin_task
helicity_event = Tuple(
    p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4) |> Bz(1.2) |> Ry(0.1) |> Rz(0.2)
    for p in aligned
)

helicity_task = KinematicTask(topologies; initial_frame = HelicityRootFrame())
helicity_point = KinematicPoint(helicity_task, helicity_event)
helicity_x = kinematics_at(helicity_point, topologies[1])

vertex_angles(helicity_x, topologies[1], ((1, 2), 3))
```

If the summed mother momentum is already numerically at rest, the task falls
back to the current axes so that aligned rest-frame events keep their visible
orientation.

## Relative paths and requested alignments

Relative Wigner angles can be inspected directly for a chosen external line.

```@example kin_task
relative_wigner_angles(
    topologies[1],
    topologies[2],
    1,
    current_event;
    initial_frame = CurrentFrame(),
)
```

Here `wigner_finals = (1, 3)` requested spin-frame alignments for final
particles 1 and 3. The result is indexed by final-particle order; unrequested
final particles carry the identity rotation.

```@example kin_task
alignment_angles_at(point, topologies[2])
```
