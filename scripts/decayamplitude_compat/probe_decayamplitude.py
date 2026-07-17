#!/usr/bin/env python3
"""Probe decayamplitude v0.1.9 against the checked one-event reference."""

from __future__ import annotations

import math
import pathlib
import tomllib

import numpy as np
from decayamplitude.chain import AlignedChain, DecayChain
from decayamplitude.kinematics_helpers import mass_from_node
from decayamplitude.resonance import Resonance
from decayamplitude.rotation import QN, wigner_capital_d
from decayangle.decay_topology import Node, Topology


ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "docs" / "fixtures" / "decayamplitude_compat"


def load_toml(name: str) -> dict:
    with (FIXTURES / name).open("rb") as stream:
        return tomllib.load(stream)


def constant_lineshape(_two_l: int, _two_s: int, *_args: float) -> float:
    return 1.0


event = load_toml("event.toml")
reference = load_toml("reference.toml")
momenta = {
    row["particle"]: np.array([row["px"], row["py"], row["pz"], row["E"]])
    for row in event["momenta"]
}

topology = Topology(0, decay_topology=((1, 2), 3))
root_resonance = Resonance(
    Node(0),
    quantum_numbers=QN(2, 1),
    lineshape=constant_lineshape,
    argnames=[],
    preserve_partity=False,
    name="M",
)
isobar_resonance = Resonance(
    Node((1, 2)),
    quantum_numbers=QN(1, 1),
    lineshape=constant_lineshape,
    argnames=[],
    preserve_partity=False,
    name="R",
)
final_state_qn = {1: QN(0, 1), 2: QN(1, 1), 3: QN(1, 1)}
chain = DecayChain(
    topology,
    {0: root_resonance, (1, 2): isobar_resonance},
    momenta,
    final_state_qn,
)

alternate_topology = Topology(0, decay_topology=((1, 3), 2))
alternate_root = Resonance(
    Node(0),
    quantum_numbers=QN(2, 1),
    lineshape=constant_lineshape,
    argnames=[],
    preserve_partity=False,
    name="M_alternate",
)
alternate_isobar = Resonance(
    Node((1, 3)),
    quantum_numbers=QN(1, 1),
    lineshape=constant_lineshape,
    argnames=[],
    preserve_partity=False,
    name="R13",
)
aligned_chain = AlignedChain(
    alternate_topology,
    {0: alternate_root, (1, 3): alternate_isobar},
    momenta,
    final_state_qn,
    topology,
)

# Constructing the recursive root attaches daughter quantum numbers to Resonance.
_ = chain.root
arguments = {
    root_resonance.id: {"couplings": {(0, 2): 1.0}},
    isobar_resonance.id: {"couplings": {(0, 1): 1.0}},
}
aligned_arguments = {
    alternate_root.id: {"couplings": {(0, 2): 1.0}},
    alternate_isobar.id: {"couplings": {(0, 1): 1.0}},
}

angles = chain.helicity_angles
root_angle = angles[((1, 2), 3)]
isobar_angle = angles[(1, 2)]


def mass_squared(indices: tuple[int, ...]) -> float:
    vector = sum((momenta[i] for i in indices), start=np.zeros(4))
    return float(vector[3] ** 2 - np.dot(vector[:3], vector[:3]))


line_masses2 = [
    mass_squared((1,)),
    mass_squared((2,)),
    mass_squared((3,)),
    mass_squared((1, 2)),
    mass_squared((1, 2, 3)),
]
kinematic_errors = [
    *(abs(a - b) for a, b in zip(line_masses2, reference["line_masses2"])),
    abs(float(root_angle.phi_rf) - reference["root_phi"]),
    abs(math.cos(float(root_angle.theta_rf)) - reference["root_cos_theta"]),
    abs(float(isobar_angle.phi_rf) - reference["isobar_phi"]),
    abs(math.cos(float(isobar_angle.theta_rf)) - reference["isobar_cos_theta"]),
]

amplitude_errors: list[float] = []
for row in reference["amplitudes"]:
    lambdas = {1: row["two_h1"], 2: row["two_h2"], 3: row["two_h3"]}
    value = complex(chain.chain_function(row["two_h0"], lambdas, arguments))
    amplitude_errors.append(abs(value - complex(row["re"], row["im"])))

alignment_errors: list[float] = []
for row in reference["alignments"]:
    value = aligned_chain.wigner_rotation[row["particle"]]
    alignment_errors.extend(
        [
            abs(float(value.phi_rf) - row["alpha"]),
            abs(math.cos(float(value.theta_rf)) - row["cos_beta"]),
            abs(float(value.psi_rf) - row["gamma"]),
        ]
    )

aligned_amplitude_errors: list[float] = []
aligned_matrices = {
    two_h0: aligned_chain.aligned_matrix(two_h0, aligned_arguments)
    for two_h0 in (-2, 0, 2)
}
for row in reference["aligned_amplitudes"]:
    key = (row["two_h1"], row["two_h2"], row["two_h3"])
    value = complex(aligned_matrices[row["two_h0"]][key])
    aligned_amplitude_errors.append(abs(value - complex(row["re"], row["im"])))

term_errors: list[float] = []
for row in reference["terms"]:
    two_lambda_r = row["two_lambda_internal"]
    root_coupling = complex(
        root_resonance.amplitude(
            0,
            two_lambda_r,
            1,
            arguments,
            mass_from_node(Node((1, 2)), momenta),
            mass_from_node(Node(3), momenta),
        )
    )
    isobar_coupling = complex(
        isobar_resonance.amplitude(
            two_lambda_r,
            0,
            -1,
            arguments,
            mass_from_node(Node(1), momenta),
            mass_from_node(Node(2), momenta),
        )
    )
    root_rotation = complex(
        np.conj(
            wigner_capital_d(
                root_angle.phi_rf,
                root_angle.theta_rf,
                0,
                2,
                0,
                two_lambda_r - 1,
            )
        )
    )
    isobar_rotation = complex(
        np.conj(
            wigner_capital_d(
                isobar_angle.phi_rf,
                isobar_angle.theta_rf,
                0,
                1,
                two_lambda_r,
                1,
            )
        )
    )
    values = {
        "root_coupling": root_coupling,
        "isobar_coupling": isobar_coupling,
        "root_amplitude": root_coupling * root_rotation,
        "isobar_amplitude": isobar_coupling * isobar_rotation,
    }
    values["product"] = values["root_amplitude"] * values["isobar_amplitude"]
    for name, value in values.items():
        expected = complex(row[f"{name}_re"], row[f"{name}_im"])
        term_errors.append(abs(value - expected))

max_kinematic_error = max(kinematic_errors)
max_term_error = max(term_errors)
max_amplitude_error = max(amplitude_errors)
max_alignment_error = max(alignment_errors)
max_aligned_amplitude_error = max(aligned_amplitude_errors)

assert max_kinematic_error < 1.0e-14
assert max_term_error < 1.0e-14
assert max_amplitude_error < 1.0e-14
assert max_alignment_error < 1.0e-13
assert max_aligned_amplitude_error < 1.0e-14

print("decayamplitude one-event probe")
print(f"  topology: {topology.tuple}")
print(f"  max kinematic error: {max_kinematic_error}")
print(f"  max computation-graph term error: {max_term_error}")
print(f"  max external-amplitude error: {max_amplitude_error}")
print(f"  max alignment-angle error: {max_alignment_error}")
print(f"  max aligned-amplitude error: {max_aligned_amplitude_error}")
