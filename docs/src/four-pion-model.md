```@meta
CurrentModule = CascadeDecays
EditURL = "../four-pion-model.qmd"
```

# [Four-pion model-building catalogue](@id four_pion_model)


This note builds a readable catalogue for

``` math
\chi_{cJ} \to \pi^+_1 \pi^-_2 \pi^+_3 \pi^-_4,
\qquad J = 0,2.
```

The first half is purely kinematical: it identifies which topology
copies are needed before any resonance model is attached. The second
half adds dynamics: toy Breit-Wigner lineshapes, LS-coupling expansion,
a `KinematicTask`, RAMBO phase-space events, and small timing checks.

The two parent hypotheses are

``` math
\chi_{c0}(0^{++}), \qquad \chi_{c2}(2^{++}).
```

The topology catalogue is common to both. The LS lists differ because
the production vertex knows the parent spin and parity.

Regeneration note: this tutorial is intentionally checked in as
generated artifacts and is not rendered by `docs/make.jl` in CI. After
editing this notebook, run

``` sh
julia --project=docs docs/render-four-pion-model.jl
```

and commit the updated `docs/four-pion-model.md` and
`docs/src/four-pion-model.md` files.

## Setup

<details class="code-fold">
<summary>Code</summary>

``` julia
using CascadeDecays
using DataFrames
using HadronicLineshapes
using Random
using RemboOnDiet
using ThreeBodyDecays: @jp_str

const M_PI = 0.13957
const M_CHIC0 = 3.4147
const M_CHIC2 = 3.5562

const M_F0 = 0.98
const Γ_F0 = 0.07
const M_RHO = 0.775
const Γ_RHO = 0.149
const M_F2 = 1.275
const Γ_F2 = 0.185

const M_A1 = 1.23
const Γ_A1 = 0.42
const M_A2 = 1.32
const Γ_A2 = 0.11
const M_PI0R = 1.30
const Γ_PI0R = 0.25
const M_PI1 = 1.60
const Γ_PI1 = 0.30
const M_PI2 = 1.67
const Γ_PI2 = 0.26

const CHARGE = Dict(1 => +1, 2 => -1, 3 => +1, 4 => -1)

particle_label(i) = CHARGE[i] == +1 ? "pi$(i)+" : "pi$(i)-"
tree_label(tree) = replace(bracket_notation(DecayTopology(tree)), "," => ", ")
show_table(df) = (show(IOContext(stdout, :limit => false), MIME("text/html"), df); println())
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
parents = DataFrame(
    parent = ["chi_c0", "chi_c2"],
    jp_label = ["0+", "2+"],
    jp = [jp"0+", jp"2+"],
    mass = [M_CHIC0, M_CHIC2],
)

transform!(
    parents,
    [:jp_label, :mass] => ByRow() do jp_label, mass
        CascadeSystem(
            SystemSpinParities("0-", "0-", "0-", "0-"; jp0=jp_label),
            SystemMasses(M_PI, M_PI, M_PI, M_PI; m0=mass),
        )
    end => :system,
)

two_body = DataFrame(
    r2 = ["f0", "rho", "f2"],
    jp = [jp"0+", jp"1-", jp"2+"],
    mass = [M_F0, M_RHO, M_F2],
    width = [Γ_F0, Γ_RHO, Γ_F2],
)

three_body = DataFrame(
    r3 = ["a1", "a2", "pi0", "pi1", "pi2"],
    jp = [jp"1+", jp"2+", jp"0-", jp"1-", jp"2-"],
    mass = [M_A1, M_A2, M_PI0R, M_PI1, M_PI2],
    width = [Γ_A1, Γ_A2, Γ_PI0R, Γ_PI1, Γ_PI2],
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(select(parents, :parent, :jp, :mass))
```

</details>

<div style="float: left;">

2×3 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

parent
</th>

<th style="text-align: left;">

jp
</th>

<th style="text-align: left;">

mass
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="ThreeBodyDecays.SpinParity" style="text-align: left;">

SpinPari…
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

SpinParity(0, '+')
</td>

<td style="text-align: right;">

3.4147
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

SpinParity(4, '+')
</td>

<td style="text-align: right;">

3.5562
</td>

</tr>

</tbody>

</table>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(select(two_body, :r2 => :name, :jp, :mass, :width))
show_table(select(three_body, :r3 => :name, :jp, :mass, :width))
```

</details>

<div style="float: left;">

3×4 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

name
</th>

<th style="text-align: left;">

jp
</th>

<th style="text-align: left;">

mass
</th>

<th style="text-align: left;">

width
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="ThreeBodyDecays.SpinParity" style="text-align: left;">

SpinPari…
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

SpinParity(0, '+')
</td>

<td style="text-align: right;">

0.98
</td>

<td style="text-align: right;">

0.07
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

SpinParity(2, '-')
</td>

<td style="text-align: right;">

0.775
</td>

<td style="text-align: right;">

0.149
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

SpinParity(4, '+')
</td>

<td style="text-align: right;">

1.275
</td>

<td style="text-align: right;">

0.185
</td>

</tr>

</tbody>

</table>

<div style="float: left;">

5×4 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

name
</th>

<th style="text-align: left;">

jp
</th>

<th style="text-align: left;">

mass
</th>

<th style="text-align: left;">

width
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="ThreeBodyDecays.SpinParity" style="text-align: left;">

SpinPari…
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

SpinParity(2, '+')
</td>

<td style="text-align: right;">

1.23
</td>

<td style="text-align: right;">

0.42
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

SpinParity(4, '+')
</td>

<td style="text-align: right;">

1.32
</td>

<td style="text-align: right;">

0.11
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

SpinParity(0, '-')
</td>

<td style="text-align: right;">

1.3
</td>

<td style="text-align: right;">

0.25
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

SpinParity(2, '-')
</td>

<td style="text-align: right;">

1.6
</td>

<td style="text-align: right;">

0.3
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

SpinParity(4, '-')
</td>

<td style="text-align: right;">

1.67
</td>

<td style="text-align: right;">

0.26
</td>

</tr>

</tbody>

</table>

## Kinematic Topology Catalogue

The cascade class starts from one positive seed,

``` math
(((1,2),3),4)
\equiv ((\pi^+_1\pi^-_2)\pi^+_3)\pi^-_4 .
```

Same-charge permutations generate the positive copies. Charge
conjugation gives one negative seed, and applying the same
charge-preserving permutations to that seed generates the negative
copies. The order of children is preserved in this class, so the
conjugate of `(1,2)` is `(2,1)`.

<details class="code-fold">
<summary>Code</summary>

``` julia
identity_map = Dict(i => i for i in 1:4)
plus_swap = Dict(1 => 3, 3 => 1, 2 => 2, 4 => 4)
minus_swap = Dict(1 => 1, 3 => 3, 2 => 4, 4 => 2)
charge_conjugation = Dict(1 => 2, 2 => 1, 3 => 4, 4 => 3)

bose_maps = DataFrame(
    operation = ["identity", "swap pi+", "swap pi-", "swap pi+ and pi-"],
    pmap = [
        identity_map,
        plus_swap,
        minus_swap,
        Dict(i => plus_swap[minus_swap[i]] for i in 1:4),
    ],
)

map_tree(i::Integer, pmap) = pmap[Int(i)]
map_tree(tree::Tuple, pmap) = map(x -> map_tree(x, pmap), tree)

cascade_seed = (((1, 2), 3), 4)
conjugate_seed = map_tree(cascade_seed, charge_conjugation)

positive_steps = transform(
    bose_maps,
    :operation => ByRow(op -> "positive: " * op) => :stage,
    :pmap => ByRow(pmap -> map_tree(cascade_seed, pmap)) => :tree,
)

conjugate_steps = transform(
    bose_maps,
    :operation => ByRow(op -> "conjugate: " * op) => :stage,
    :pmap => ByRow(pmap -> map_tree(conjugate_seed, pmap)) => :tree,
)

cascade_steps = vcat(positive_steps, conjugate_steps)
cascade_topologies = unique(select(cascade_steps, :stage, :tree))
transform!(
    cascade_topologies,
    :tree => ByRow() do tree
        pair = tree[1][1]
        bachelor = tree[1][2]
        spectator = tree[2]
        triple = Tuple(vcat(collect(pair), bachelor))
        (;
            topology = DecayTopology(tree),
            R3_charge = sum(CHARGE[i] for i in triple) == +1 ? "+" : "-",
            R3_pions = join(particle_label.(triple), " "),
            R2_pair = join(particle_label.(pair), " "),
            bachelor = particle_label(bachelor),
            spectator = particle_label(spectator),
            bracket = "`" * tree_label(tree) * "`",
        )
    end => AsTable,
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(
    select(
        cascade_topologies,
        :stage,
        :R3_charge,
        :R3_pions,
        :R2_pair,
        :bachelor,
        :spectator,
        :bracket,
    ),
)
```

</details>

<div style="float: left;">

8×7 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

stage
</th>

<th style="text-align: left;">

R3_charge
</th>

<th style="text-align: left;">

R3_pions
</th>

<th style="text-align: left;">

R2_pair
</th>

<th style="text-align: left;">

bachelor
</th>

<th style="text-align: left;">

spectator
</th>

<th style="text-align: left;">

bracket
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

positive: identity
</td>

<td style="text-align: left;">

\+
</td>

<td style="text-align: left;">

pi1+ pi2- pi3+
</td>

<td style="text-align: left;">

pi1+ pi2-
</td>

<td style="text-align: left;">

pi3+
</td>

<td style="text-align: left;">

pi4-
</td>

<td style="text-align: left;">

`(((1, 2), 3), 4)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

positive: swap pi+
</td>

<td style="text-align: left;">

\+
</td>

<td style="text-align: left;">

pi3+ pi2- pi1+
</td>

<td style="text-align: left;">

pi3+ pi2-
</td>

<td style="text-align: left;">

pi1+
</td>

<td style="text-align: left;">

pi4-
</td>

<td style="text-align: left;">

`(((3, 2), 1), 4)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

positive: swap pi-
</td>

<td style="text-align: left;">

\+
</td>

<td style="text-align: left;">

pi1+ pi4- pi3+
</td>

<td style="text-align: left;">

pi1+ pi4-
</td>

<td style="text-align: left;">

pi3+
</td>

<td style="text-align: left;">

pi2-
</td>

<td style="text-align: left;">

`(((1, 4), 3), 2)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

positive: swap pi+ and pi-
</td>

<td style="text-align: left;">

\+
</td>

<td style="text-align: left;">

pi3+ pi4- pi1+
</td>

<td style="text-align: left;">

pi3+ pi4-
</td>

<td style="text-align: left;">

pi1+
</td>

<td style="text-align: left;">

pi2-
</td>

<td style="text-align: left;">

`(((3, 4), 1), 2)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

conjugate: identity
</td>

<td style="text-align: left;">

\-
</td>

<td style="text-align: left;">

pi2- pi1+ pi4-
</td>

<td style="text-align: left;">

pi2- pi1+
</td>

<td style="text-align: left;">

pi4-
</td>

<td style="text-align: left;">

pi3+
</td>

<td style="text-align: left;">

`(((2, 1), 4), 3)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

6
</td>

<td style="text-align: left;">

conjugate: swap pi+
</td>

<td style="text-align: left;">

\-
</td>

<td style="text-align: left;">

pi2- pi3+ pi4-
</td>

<td style="text-align: left;">

pi2- pi3+
</td>

<td style="text-align: left;">

pi4-
</td>

<td style="text-align: left;">

pi1+
</td>

<td style="text-align: left;">

`(((2, 3), 4), 1)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

7
</td>

<td style="text-align: left;">

conjugate: swap pi-
</td>

<td style="text-align: left;">

\-
</td>

<td style="text-align: left;">

pi4- pi1+ pi2-
</td>

<td style="text-align: left;">

pi4- pi1+
</td>

<td style="text-align: left;">

pi2-
</td>

<td style="text-align: left;">

pi3+
</td>

<td style="text-align: left;">

`(((4, 1), 2), 3)`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

8
</td>

<td style="text-align: left;">

conjugate: swap pi+ and pi-
</td>

<td style="text-align: left;">

\-
</td>

<td style="text-align: left;">

pi4- pi3+ pi2-
</td>

<td style="text-align: left;">

pi4- pi3+
</td>

<td style="text-align: left;">

pi2-
</td>

<td style="text-align: left;">

pi1+
</td>

<td style="text-align: left;">

`(((4, 3), 2), 1)`
</td>

</tr>

</tbody>

</table>

The second topology class has two two-body isobars,

``` math
\chi_{cJ} \to R_2(\pi^+\pi^-)\,R_2(\pi^+\pi^-).
```

Here reversed child order inside both two-body pairs is a bookkeeping
duplicate. This is the useful cross-check:

``` math
((1,2),(3,4)) \sim ((2,1),(4,3)),
\qquad
A_{(12)(34)} + A_{(21)(43)} = 2 A_{(12)(34)} .
```

<details class="code-fold">
<summary>Code</summary>

``` julia
pair_pair_seed = ((1, 2), (3, 4))

pair_pair_checks = transform(
    bose_maps,
    :operation => ByRow(op -> "Bose: " * op) => :stage,
    :pmap => ByRow() do pmap
        map_tree(pair_pair_seed, pmap)
    end => :raw_tree,
)

push!(
    pair_pair_checks,
    (
        operation = "charge conjugation",
        pmap = charge_conjugation,
        stage = "conjugate seed",
        raw_tree = map_tree(pair_pair_seed, charge_conjugation),
    ),
)

transform!(
    pair_pair_checks,
    :raw_tree => ByRow() do tree
        pairs = map(tree) do pair
            a, b = pair
            CHARGE[a] == +1 ? (a, b) : (b, a)
        end
        Tuple(sort(collect(pairs); by=string))
    end => :tree,
    :raw_tree => ByRow(tree -> "`" * string(tree) * "`") => :raw,
)
transform!(
    pair_pair_checks,
    :tree => ByRow(tree -> "`" * tree_label(tree) * "`") => :canonical,
)

pair_pair_topologies = unique(select(pair_pair_checks, :tree))
transform!(
    pair_pair_topologies,
    :tree => ByRow() do tree
        (;
            topology = DecayTopology(tree),
            R2a_pair = join(particle_label.(tree[1]), " "),
            R2b_pair = join(particle_label.(tree[2]), " "),
            bracket = "`" * tree_label(tree) * "`",
        )
    end => AsTable,
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(select(pair_pair_topologies, :R2a_pair, :R2b_pair, :bracket))
show_table(select(pair_pair_checks, :stage, :raw, :canonical))
```

</details>

<div style="float: left;">

2×3 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

R2a_pair
</th>

<th style="text-align: left;">

R2b_pair
</th>

<th style="text-align: left;">

bracket
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

pi1+ pi2-
</td>

<td style="text-align: left;">

pi3+ pi4-
</td>

<td style="text-align: left;">

`((1, 2), (3, 4))`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

pi1+ pi4-
</td>

<td style="text-align: left;">

pi3+ pi2-
</td>

<td style="text-align: left;">

`((1, 4), (3, 2))`
</td>

</tr>

</tbody>

</table>

<div style="float: left;">

5×3 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

stage
</th>

<th style="text-align: left;">

raw
</th>

<th style="text-align: left;">

canonical
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

Bose: identity
</td>

<td style="text-align: left;">

`((1, 2), (3, 4))`
</td>

<td style="text-align: left;">

`((1, 2), (3, 4))`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

Bose: swap pi+
</td>

<td style="text-align: left;">

`((3, 2), (1, 4))`
</td>

<td style="text-align: left;">

`((1, 4), (3, 2))`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

Bose: swap pi-
</td>

<td style="text-align: left;">

`((1, 4), (3, 2))`
</td>

<td style="text-align: left;">

`((1, 4), (3, 2))`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

Bose: swap pi+ and pi-
</td>

<td style="text-align: left;">

`((3, 4), (1, 2))`
</td>

<td style="text-align: left;">

`((1, 2), (3, 4))`
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

conjugate seed
</td>

<td style="text-align: left;">

`((2, 1), (4, 3))`
</td>

<td style="text-align: left;">

`((1, 2), (3, 4))`
</td>

</tr>

</tbody>

</table>

## LS Catalogue

The following tables attach resonance quantum numbers and ask
`CascadeDecays.jl` for the allowed LS couplings. The construction is
still kinematic plus spin-parity bookkeeping; no lineshape parameters
enter yet.

<details class="code-fold">
<summary>Code</summary>

``` julia
ls_label((two_l, two_s)) =
    "(L=$(iseven(two_l) ? string(two_l ÷ 2) : string(two_l, "/2")), " *
    "S=$(iseven(two_s) ? string(two_s ÷ 2) : string(two_s, "/2")))"
ls_list_label(couplings) = isempty(couplings) ? "forbidden" : join(ls_label.(couplings), ", ")

cascade_ls = crossjoin(
    select(parents, :parent, :system),
    select(three_body, :r3, :jp => :r3_jp),
    select(two_body, :r2, :jp => :r2_jp),
)

transform!(
    cascade_ls,
    [:system, :r3_jp, :r2_jp] => ByRow() do system, r3_jp, r2_jp
        row = first(cascade_topologies)
        propagators = (
            (row.tree[1][1], row.tree[1][2]) => Propagator(r3_jp, ConstantLineshape(:R3)),
            row.tree[1][1] => Propagator(r2_jp, ConstantLineshape(:R2)),
        )
        couplings = possible_vertex_couplings(row.topology, system, propagators)
        (;
            parent_ls = couplings[1].second,
            r3_ls = couplings[2].second,
            r2_ls = couplings[3].second,
            n_ls = prod(length(c.second) for c in couplings),
        )
    end => AsTable,
)

transform!(
    cascade_ls,
    :parent_ls => ByRow(ls_list_label) => Symbol("parent -> R3 pi"),
    :r3_ls => ByRow(ls_list_label) => Symbol("R3 -> R2 pi"),
    :r2_ls => ByRow(ls_list_label) => Symbol("R2 -> pi pi"),
    :n_ls => ByRow(n -> n * nrow(cascade_topologies)) => :explicit_chains,
)

pair_pair_ls = crossjoin(
    select(parents, :parent, :system),
    select(two_body, :r2 => :r2a, :jp => :r2a_jp),
    select(two_body, :r2 => :r2b, :jp => :r2b_jp),
)

transform!(
    pair_pair_ls,
    [:system, :r2a_jp, :r2b_jp] => ByRow() do system, r2a_jp, r2b_jp
        row = first(pair_pair_topologies)
        propagators = (
            row.tree[1] => Propagator(r2a_jp, ConstantLineshape(:R2a)),
            row.tree[2] => Propagator(r2b_jp, ConstantLineshape(:R2b)),
        )
        couplings = possible_vertex_couplings(row.topology, system, propagators)
        (;
            parent_ls = couplings[1].second,
            r2a_ls = couplings[2].second,
            r2b_ls = couplings[3].second,
            n_ls = prod(length(c.second) for c in couplings),
        )
    end => AsTable,
)

transform!(
    pair_pair_ls,
    :parent_ls => ByRow(ls_list_label) => Symbol("parent -> R2a R2b"),
    :r2a_ls => ByRow(ls_list_label) => Symbol("R2a -> pi pi"),
    :r2b_ls => ByRow(ls_list_label) => Symbol("R2b -> pi pi"),
    :n_ls => ByRow(n -> n * nrow(pair_pair_topologies)) => :explicit_chains,
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(
    select(
        cascade_ls,
        :parent,
        :r3,
        :r2,
        Symbol("parent -> R3 pi"),
        Symbol("R3 -> R2 pi"),
        Symbol("R2 -> pi pi"),
        :n_ls,
        :explicit_chains,
    ),
)
```

</details>

<div style="float: left;">

30×8 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

parent
</th>

<th style="text-align: left;">

r3
</th>

<th style="text-align: left;">

r2
</th>

<th style="text-align: left;">

parent -\> R3 pi
</th>

<th style="text-align: left;">

R3 -\> R2 pi
</th>

<th style="text-align: left;">

R2 -\> pi pi
</th>

<th style="text-align: left;">

n_ls
</th>

<th style="text-align: left;">

explicit_chains
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=0, S=1), (L=2, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

16
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

16
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=2, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

6
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

7
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

8
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

9
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

10
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

11
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

12
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

13
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

14
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

16
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

15
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

3
</td>

<td style="text-align: right;">

24
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

16
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

16
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

17
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=0, S=1), (L=2, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

4
</td>

<td style="text-align: right;">

32
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

18
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

4
</td>

<td style="text-align: right;">

32
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

19
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

20
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=2, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

16
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

21
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

a2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=1, S=2), (L=3, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

4
</td>

<td style="text-align: right;">

32
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

22
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

23
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

24
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

25
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=2, S=1)
</td>

<td style="text-align: left;">

forbidden
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

0
</td>

<td style="text-align: right;">

0
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

26
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=2, S=1)
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

27
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi1
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=2, S=1)
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

8
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

28
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

3
</td>

<td style="text-align: right;">

24
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

29
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

6
</td>

<td style="text-align: right;">

48
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

30
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

pi2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

9
</td>

<td style="text-align: right;">

72
</td>

</tr>

</tbody>

</table>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(
    select(
        pair_pair_ls,
        :parent,
        :r2a,
        :r2b,
        Symbol("parent -> R2a R2b"),
        Symbol("R2a -> pi pi"),
        Symbol("R2b -> pi pi"),
        :n_ls,
        :explicit_chains,
    ),
)
```

</details>

<div style="float: left;">

18×8 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

parent
</th>

<th style="text-align: left;">

r2a
</th>

<th style="text-align: left;">

r2b
</th>

<th style="text-align: left;">

parent -\> R2a R2b
</th>

<th style="text-align: left;">

R2a -\> pi pi
</th>

<th style="text-align: left;">

R2b -\> pi pi
</th>

<th style="text-align: left;">

n_ls
</th>

<th style="text-align: left;">

explicit_chains
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=1, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=0, S=0), (L=2, S=2)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

6
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=3)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

7
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=2, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

8
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=3)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

9
</td>

<td style="text-align: left;">

chi_c0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=0, S=0), (L=2, S=2), (L=4, S=4)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

3
</td>

<td style="text-align: right;">

6
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

10
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

1
</td>

<td style="text-align: right;">

2
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

11
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

12
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

3
</td>

<td style="text-align: right;">

6
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

13
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=1, S=1), (L=3, S=1)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

2
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

14
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=0), (L=2, S=1), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

5
</td>

<td style="text-align: right;">

10
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

15
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=1, S=1), (L=1, S=2), (L=1, S=3), (L=3, S=1), (L=3, S=2), (L=3, S=3),
(L=5, S=3)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

7
</td>

<td style="text-align: right;">

14
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

16
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

f0
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=2), (L=4, S=2)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=0, S=0)
</td>

<td style="text-align: right;">

3
</td>

<td style="text-align: right;">

6
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

17
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: left;">

(L=1, S=1), (L=1, S=2), (L=1, S=3), (L=3, S=1), (L=3, S=2), (L=3, S=3),
(L=5, S=3)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=1, S=0)
</td>

<td style="text-align: right;">

7
</td>

<td style="text-align: right;">

14
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

18
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

f2
</td>

<td style="text-align: left;">

(L=0, S=2), (L=2, S=0), (L=2, S=1), (L=2, S=2), (L=2, S=3), (L=2, S=4),
(L=4, S=2), (L=4, S=3), (L=4, S=4), (L=6, S=4)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: left;">

(L=2, S=0)
</td>

<td style="text-align: right;">

10
</td>

<td style="text-align: right;">

20
</td>

</tr>

</tbody>

</table>

## Dynamics

The dynamic model below is deliberately small. It uses the `chi_c2`
parent and the `a1 rho` cascade family, expanded over all eight
charge/Bose topology copies and all allowed LS products. The lineshapes
are toy Breit-Wigners; the numbers only make the benchmark finite and
reproducible.

<details class="code-fold">
<summary>Code</summary>

``` julia
transform!(
    two_body,
    [:mass, :width] => ByRow((m, Γ) -> BreitWigner(m, Γ)) => :lineshape,
)
transform!(
    three_body,
    [:mass, :width] => ByRow((m, Γ) -> BreitWigner(m, Γ)) => :lineshape,
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
dynamic_rows = crossjoin(
    select(subset(parents, :parent => ByRow(==("chi_c2"))), :parent, :system),
    select(cascade_topologies, :stage, :tree, :topology, :bracket),
    select(subset(three_body, :r3 => ByRow(==("a1"))), :r3, :jp => :r3_jp, :lineshape => :r3_lineshape),
    select(subset(two_body, :r2 => ByRow(==("rho"))), :r2, :jp => :r2_jp, :lineshape => :r2_lineshape),
)

transform!(
    dynamic_rows,
    [:tree, :r3_jp, :r3_lineshape, :r2_jp, :r2_lineshape] => ByRow() do tree, r3_jp, r3_ls, r2_jp, r2_ls
        (
            (tree[1][1], tree[1][2]) => Propagator(r3_jp, r3_ls),
            tree[1][1] => Propagator(r2_jp, r2_ls),
        )
    end => :propagators,
)

transform!(
    dynamic_rows,
    [:topology, :system, :propagators] => ByRow() do topology, system, propagators
        all_ls_decay_chains(topology, system, propagators)
    end => :chains,
)
transform!(
    dynamic_rows,
    :chains => ByRow(length) => :n_chains,
)

toy_chains = Tuple(reduce(vcat, dynamic_rows.chains))
toy_system = only(unique(dynamic_rows.system))
toy_topologies = Tuple(unique(chain.topology for chain in toy_chains))
toy_model = CascadeDecay(
    toy_chains,
    toy_system,
    first(toy_topologies);
    couplings=ntuple(_ -> 1.0 + 0.0im, length(toy_chains)),
    names=Tuple("a1rho_$(i)" for i in eachindex(toy_chains)),
)

task = KinematicTask(toy_topologies; reference_topology=first(toy_topologies))
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
show_table(select(dynamic_rows, :parent, :stage, :bracket, :r3, :r2, :n_chains))
```

</details>

<div style="float: left;">

8×6 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

parent
</th>

<th style="text-align: left;">

stage
</th>

<th style="text-align: left;">

bracket
</th>

<th style="text-align: left;">

r3
</th>

<th style="text-align: left;">

r2
</th>

<th style="text-align: left;">

n_chains
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

positive: identity
</td>

<td style="text-align: left;">

`(((1, 2), 3), 4)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

positive: swap pi+
</td>

<td style="text-align: left;">

`(((3, 2), 1), 4)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

positive: swap pi-
</td>

<td style="text-align: left;">

`(((1, 4), 3), 2)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

4
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

positive: swap pi+ and pi-
</td>

<td style="text-align: left;">

`(((3, 4), 1), 2)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

5
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

conjugate: identity
</td>

<td style="text-align: left;">

`(((2, 1), 4), 3)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

6
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

conjugate: swap pi+
</td>

<td style="text-align: left;">

`(((2, 3), 4), 1)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

7
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

conjugate: swap pi-
</td>

<td style="text-align: left;">

`(((4, 1), 2), 3)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

8
</td>

<td style="text-align: left;">

chi_c2
</td>

<td style="text-align: left;">

conjugate: swap pi+ and pi-
</td>

<td style="text-align: left;">

`(((4, 3), 2), 1)`
</td>

<td style="text-align: left;">

a1
</td>

<td style="text-align: left;">

rho
</td>

<td style="text-align: right;">

4
</td>

</tr>

</tbody>

</table>

## Phase Space and Timing

The timing cells use only 100 flat phase-space events. The goal is a
local sanity check: how much time is spent turning four-vectors into
topology-local kinematics, and how much time is spent evaluating the toy
coherent amplitude.

<details class="code-fold">
<summary>Code</summary>

``` julia
const N_BENCH_EVENTS = 100

rng = MersenneTwister(20260707)
generator = PhaseSpaceGenerator(fill(M_PI, 4), M_CHIC2)

events = DataFrame(event = [rand(rng, generator) for _ in 1:N_BENCH_EVENTS])
transform!(
    events,
    :event => ByRow(event -> Tuple(event.momenta)) => :momenta,
    :event => ByRow(event -> event.weight) => :weight,
)
nothing
```

</details>

<details class="code-fold">
<summary>Code</summary>

``` julia
KinematicPoint(task, first(events.momenta))
first_point = KinematicPoint(task, first(events.momenta))
amplitude(toy_model, first_point)

kinematics_seconds = @elapsed begin
    points = [KinematicPoint(task, momenta) for momenta in events.momenta]
end

amplitude_seconds = @elapsed begin
    amplitudes = [amplitude(toy_model, point) for point in points]
end

intensity_seconds = @elapsed begin
    intensities = [sum(abs2, amp) for amp in amplitudes]
end

benchmark = DataFrame(
    step = ["KinematicPoint construction", "coherent amplitude", "sum(abs2, amplitude)"],
    events = N_BENCH_EVENTS,
    total_ms = 1000 .* [kinematics_seconds, amplitude_seconds, intensity_seconds],
)

transform!(
    benchmark,
    [:total_ms, :events] => ByRow((total_ms, n) -> total_ms / n) => :ms_per_event,
)

show_table(benchmark)
```

</details>

<div style="float: left;">

3×4 DataFrame

</div>

<div style="clear: both;">

</div>

<table class="data-frame" style="margin-bottom: 6px;">

<thead>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

Row
</th>

<th style="text-align: left;">

step
</th>

<th style="text-align: left;">

events
</th>

<th style="text-align: left;">

total_ms
</th>

<th style="text-align: left;">

ms_per_event
</th>

</tr>

<tr class="columnLabelRow">

<th class="stubheadLabel" style="font-weight: bold; text-align: right;">

</th>

<th title="String" style="text-align: left;">

String
</th>

<th title="Int64" style="text-align: left;">

Int64
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

<th title="Float64" style="text-align: left;">

Float64
</th>

</tr>

</thead>

<tbody>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

1
</td>

<td style="text-align: left;">

KinematicPoint construction
</td>

<td style="text-align: right;">

100
</td>

<td style="text-align: right;">

1226.32
</td>

<td style="text-align: right;">

12.2632
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

2
</td>

<td style="text-align: left;">

coherent amplitude
</td>

<td style="text-align: right;">

100
</td>

<td style="text-align: right;">

1737.83
</td>

<td style="text-align: right;">

17.3783
</td>

</tr>

<tr class="dataRow">

<td class="rowLabel" style="font-weight: bold; text-align: right;">

3
</td>

<td style="text-align: left;">

sum(abs2, amplitude)
</td>

<td style="text-align: right;">

100
</td>

<td style="text-align: right;">

27.3835
</td>

<td style="text-align: right;">

0.273835
</td>

</tr>

</tbody>

</table>

The benchmark model is intentionally not the full physics model. It
exercises the same kinematic routing and amplitude machinery on a
compact, auditable subset. The full catalogue above tells which
additional resonance and LS families should be added when moving from
this toy benchmark to the analysis model.
