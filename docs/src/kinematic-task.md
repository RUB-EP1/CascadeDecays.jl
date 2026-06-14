# Kinematic Tasks And Points

```@meta
CurrentModule = CascadeDecays
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
    ThreeBodySpins,
    aligned_four_vectors,
    cosθij,
    x2σs

using CascadeDecays

fourvector(p) = FourVector(Float64(p[1]), Float64(p[2]), Float64(p[3]); E = Float64(p[4]))

ms = ThreeBodyMasses(0.93827208816, 0.493677, 0.13957039; m0 = 2.28646)
σs = x2σs([0.42, 0.31], ms; k = 3)
system = CascadeSystem(SystemSpins(0, 0, 0; two_h0 = 0), SystemMasses(ms))
nothing
```

## From Topologies To A Task

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

bracket.(task.topologies)
```

The generated programs are ordinary `InstructionalDecayTrees.jl` instruction
tuples. They are stored per topology and per vertex.

```@example kin_task
task.programs[1].vertex_programs
```

The public path helper exposes the frame path for one external line. This is
often the simplest way to check why two topology conventions need a relative
rotation.

```@example kin_task
line = final_line_inds(topologies[2])[1]
(
    reference = helicity_frame_path(topologies[1], line; initial_frame = CurrentFrame()),
    target = helicity_frame_path(topologies[2], line; initial_frame = CurrentFrame()),
)
```

## Current-Frame Four-Vectors

Use [`CurrentFrame`](@ref) when the event is already expressed in the desired
mother-system axes. In this example the aligned event is rotated as
`Rz |> Ry |> Rz`, and the root vertex angles recover the last two rotations.

```@example kin_task
aligned = Tuple(fourvector(p) for p in aligned_four_vectors(σs, ms; k = 3))
current_event = Tuple(p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4) for p in aligned)

point = evaluate(task, current_event, system)
x12_3 = kinematics_at(point, topologies[1])

(
    root = vertex_angles(topologies[1], x12_3, ((1, 2), 3)),
    isobar = vertex_angles(topologies[1], x12_3, (1, 2)),
    σ12 = line_invariant(topologies[1], x12_3, (1, 2)),
)
```

The point stores one [`CascadeKinematics`](@ref) object per topology. Vertex
addresses are usually the most readable way to explore the same event in a
different topology.

```@example kin_task
x31_2 = kinematics_at(point, topologies[2])
(
    topology = bracket(topologies[2]),
    root_lines = vertex_line_inds(topologies[2], vertex_ind_for(topologies[2], ((3, 1), 2))),
    root_angles = vertex_angles(topologies[2], x31_2, ((3, 1), 2)),
    isobar_masses2 = vertex_masses2(topologies[2], x31_2, (3, 1)),
)
```

## Helicity-Root Four-Vectors

Use [`HelicityRootFrame`](@ref) for fully general four-vectors. The event is
first rotated in the mother rest frame and then boosted and reoriented as
`Rz |> Ry |> Rz |> Bz |> Ry |> Rz`.

```@example kin_task
helicity_event = Tuple(
    p |> Rz(0.5) |> Ry(0.3) |> Rz(0.4) |> Bz(1.2) |> Ry(0.1) |> Rz(0.2)
    for p in aligned
)

helicity_task = KinematicTask(topologies; initial_frame = HelicityRootFrame())
helicity_point = evaluate(helicity_task, helicity_event, system)
helicity_x = kinematics_at(helicity_point, topologies[1])

vertex_angles(topologies[1], helicity_x, ((1, 2), 3))
```

If the summed mother momentum is already numerically at rest, the task falls
back to the current axes so that aligned rest-frame events keep their visible
orientation.

## Relative Paths And Requested Alignments

Relative Wigner angles can be inspected directly for a chosen external line.

```@example kin_task
relative_wigner_angles(
    topologies[1],
    topologies[2],
    final_line_inds(topologies[2])[1],
    current_event;
    initial_frame = CurrentFrame(),
)
```

When `wigner_finals` is set, [`evaluate`](@ref) stores alignment angles for only
those requested final particles. Other external axes, including the root, remain
the identity rotation.

```@example kin_task
alignment_angles_at(point, topologies[2])
```

