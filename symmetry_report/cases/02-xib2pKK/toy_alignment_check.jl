# Toy numerical check for Xib(0) -> K-(1) p(2) K-(3):
# reconcile the fixed-frame Bose-symmetry rule with the (-1)^(...) phase of the
# invariant-level (aligned) symmetrized amplitude.
#
# Model: manifestly Bose-symmetric covariant amplitude (weak decay, no P or C imposed)
#   M_{λ0,λ2}(p1,p2,p3) = T1 + T2,
#   T1 = f(σ1,σ3) ubar(p2,λ2) Γ p̸1 u(P,λ0),   T2 = T1 with (p1<->p3), i.e. f(σ3,σ1) ... p̸3
# with f a generic asymmetric complex function and Γ = 1 + 0.3 γ5.
# Symmetry under the kaon slot swap (p1<->p3) holds by construction (QFT guarantees
# it for any physical model); what we LEARN from the test is the phase pattern that
# relates aligned amplitudes at swapped invariants.
#
# Conventions (registry):
#   metric (+,-,-,-), Dirac basis, helicity spinors via R(φ,θ,0),
#   standard orientation: all momenta in xz-plane, proton (2) along +z, kaon-1 with px>0.
#   Parent at rest, λ0 quantized along the space-fixed z-axis.

using LinearAlgebra
using Printf

# ---------- Dirac algebra ----------
const I2 = Matrix{ComplexF64}(I, 2, 2)
const sx = ComplexF64[0 1; 1 0]
const sy = ComplexF64[0 -im; im 0]
const sz = ComplexF64[1 0; 0 -1]

blk(A, B, C, D) = [A B; C D]
const g0 = blk(I2, 0I2, 0I2, -I2)
const g1 = blk(0I2, sx, -sx, 0I2)
const g2 = blk(0I2, sy, -sy, 0I2)
const g3 = blk(0I2, sz, -sz, 0I2)
const g5 = blk(0I2, I2, I2, 0I2)

slash(p) = p[1] * g0 - p[2] * g1 - p[3] * g2 - p[4] * g3   # p = (E,px,py,pz)

# ---------- helicity spinors, R(φ,θ,0) convention ----------
Rz(a) = ComplexF64[exp(-im * a / 2) 0; 0 exp(im * a / 2)]
Ry(b) = ComplexF64[cos(b / 2) -sin(b / 2); sin(b / 2) cos(b / 2)]

function uspinor(p, m, λ)  # λ = ±1/2
    P = sqrt(p[2]^2 + p[3]^2 + p[4]^2)
    θ = P < 1e-12 ? 0.0 : acos(clamp(p[4] / P, -1, 1))
    φ = atan(p[3], p[2])
    χ = (Rz(φ) * Ry(θ)) * (λ > 0 ? ComplexF64[1, 0] : ComplexF64[0, 1])
    E = p[1]
    vcat(sqrt(E + m) .* χ, (2λ) * sqrt(max(E - m, 0.0)) .* χ)
end

# ---------- kinematics ----------
const M0, mK, mp = 5.797, 0.4937, 0.938
const m1, m2, m3 = mK, mp, mK

inv2(pa, pb) = (pa[1] + pb[1])^2 - (pa[2] + pb[2])^2 - (pa[3] + pb[3])^2 - (pa[4] + pb[4])^2

# standard-orientation configuration x*(σ1,σ3):
# σ1 = (p2+p3)^2, σ3 = (p1+p2)^2; plane = xz, p2 along +z, p1 with px>0
function config(σ1, σ3)
    E1 = (M0^2 + m1^2 - σ1) / (2M0)
    E3 = (M0^2 + m3^2 - σ3) / (2M0)
    E2 = M0 - E1 - E3
    @assert E1 > m1 && E2 > m2 && E3 > m3 "outside Dalitz region"
    p1m, p2m = sqrt(E1^2 - m1^2), sqrt(E2^2 - m2^2)
    d12 = E1 * E2 - (σ3 - m1^2 - m2^2) / 2        # = p⃗1·p⃗2
    cθ = d12 / (p1m * p2m)
    @assert abs(cθ) < 1 "outside Dalitz region (angle)"
    sθ = sqrt(1 - cθ^2)
    p2 = [E2, 0.0, 0.0, p2m]
    p1 = [E1, p1m * sθ, 0.0, p1m * cθ]
    p3 = [E3, -p1[2], 0.0, -p1[4] - p2m]
    @assert abs(inv2(p3, [0, 0, 0, 0]) + 0 - (p3[1]^2 - p3[2]^2 - p3[3]^2 - p3[4]^2)) < 1e-9
    @assert abs((p3[1]^2 - p3[2]^2 - p3[3]^2 - p3[4]^2) - m3^2) < 1e-9 "p3 off-shell"
    (p1, p2, p3)
end

rotmat(α, β, γ) = [cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1] *
                  [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)] *
                  [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1]
rotx(a) = [1.0 0.0 0.0; 0.0 cos(a) -sin(a); 0.0 sin(a) cos(a)]
rot(R, p) = vcat(p[1], R * p[2:4])

# ---------- the model ----------
f(σa, σb) = (1.0 + 0.7im * σb) / (σa - (2.5 + 0.4im))   # asymmetric in its two arguments
const Γex = Matrix{ComplexF64}(I, 4, 4) + 0.3 * g5

function term(pK, p2, λ0, λ2, fval)   # ubar(p2,λ2) Γ p̸K u(P,λ0) * fval
    P0 = [M0, 0.0, 0.0, 0.0]
    ubar = adjoint(uspinor(p2, m2, λ2)) * g0
    u0 = uspinor(P0, M0, λ0)
    fval * (ubar*Γex*slash(pK)*u0)[1]
end

T1(p1, p2, p3, λ0, λ2) = term(p1, p2, λ0, λ2, f(inv2(p2, p3), inv2(p1, p2)))
T2(p1, p2, p3, λ0, λ2) = term(p3, p2, λ0, λ2, f(inv2(p1, p2), inv2(p2, p3)))
Mfull(p1, p2, p3, λ0, λ2) = T1(p1, p2, p3, λ0, λ2) + T2(p1, p2, p3, λ0, λ2)

# ---------- tests ----------
const hels = [(1 / 2, 1 / 2), (1 / 2, -1 / 2), (-1 / 2, 1 / 2), (-1 / 2, -1 / 2)]
σ1v, σ3v = 8.0, 14.0

println("=== Test A: fixed-frame Bose symmetry  M(p1,p2,p3) = M(p3,p2,p1) ===")
println("(generic orientation: standard config rotated by random Euler angles)")
R = rotmat(0.7, 1.1, -2.3)
x = config(σ1v, σ3v)
q1, q2, q3 = rot(R, x[1]), rot(R, x[2]), rot(R, x[3])
for (λ0, λ2) in hels
    a = Mfull(q1, q2, q3, λ0, λ2)
    b = Mfull(q3, q2, q1, λ0, λ2)   # kaon slots swapped
    @printf("  λ0=%+.1f λ2=%+.1f   |M - M_swap| = %.2e   (M = %+.4f%+.4fim)\n",
        λ0, λ2, abs(a - b), real(a), imag(a))
end

println("\n=== Test B: aligned amplitudes at swapped invariants ===")
println("ratio  T2(x*(σ1,σ3)) / T1(x*(σ3,σ1))  — the cross-chain alignment phase")
xs = config(σ1v, σ3v)   # standard orientation, invariants (σ1,σ3)
ys = config(σ3v, σ1v)   # standard orientation, invariants swapped
for (λ0, λ2) in hels
    r = T2(xs[1], xs[2], xs[3], λ0, λ2) / T1(ys[1], ys[2], ys[3], λ0, λ2)
    @printf("  λ0=%+.1f λ2=%+.1f   ratio = %+.6f %+.6f im\n", λ0, λ2, real(r), imag(r))
end
println("  compare: (-1)^(λ0-λ2) pattern = (+1, -1, -1, +1);  (-1)^(λ0+λ2) pattern = (-1, +1, +1, -1)")

println("\n=== Test C: the full aligned amplitude A(σ1,σ3) := M(x*(σ1,σ3)) at swapped invariants ===")
println("ratio  A(σ3,σ1) / A(σ1,σ3)  — NOT +1, despite exact Bose symmetry (Test A)")
for (λ0, λ2) in hels
    r = Mfull(ys[1], ys[2], ys[3], λ0, λ2) / Mfull(xs[1], xs[2], xs[3], λ0, λ2)
    @printf("  λ0=%+.1f λ2=%+.1f   ratio = %+.6f %+.6f im\n", λ0, λ2, real(r), imag(r))
end

println("\n=== Test D: geometric origin — (x3,x2,x1) = Rz(π) · (y1,y2,y3) ===")
Rpi = rotmat(Float64(π), 0.0, 0.0)
for (a, b) in zip((xs[3], xs[2], xs[1]), (rot(Rpi, ys[1]), rot(Rpi, ys[2]), rot(Rpi, ys[3])))
    @printf("  max component diff: %.2e\n", maximum(abs.(a .- b)))
end

println("\n=== Test E: rotation covariance — where the λ2 phase lives ===")
# Claim: helicity does NOT mix under rotations; the proton contributes only the
# little-group phase φ_W(R, p̂2):
#   M_{λ0λ2}(R·x) = e^{-i φ_W λ2} Σ_{ν0} conj(D^{1/2}_{λ0 ν0}(R)) M_{ν0 λ2}(x)
# For R = Rz(π) and p̂2 = ẑ (spectator-2 standard orientation): φ_W = π,
# i.e. the λ2 part of the test-B/C phase is exactly this little-group phase.
SU2(α, β, γ) = Rz(α) * Ry(β) * Rz(γ)    # = D^{1/2}_{m'm}(α,β,γ), rows/cols ordered (+1/2, -1/2)

function wigner_phase(Ur, R3, p)   # φ_W for a helicity state at momentum p, R(φ,θ,0) convention
    v = p[2:4]; P = norm(v)
    θ, φ = acos(clamp(v[3] / P, -1, 1)), atan(v[2], v[1])
    q = R3 * v
    θq, φq = acos(clamp(q[3] / P, -1, 1)), atan(q[2], q[1])
    Uh(t, f) = Rz(f) * Ry(t)
    W = adjoint(Uh(θq, φq)) * Ur * Uh(θ, φ)
    @assert abs(W[1, 2]) + abs(W[2, 1]) < 1e-10 "little-group element not diagonal"
    2 * angle(W[2, 2])               # W = diag(e^{-iφ_W/2}, e^{+iφ_W/2})
end

let (α, β, γ) = (0.7, 1.1, -2.3)
    R3, Ur = rotmat(α, β, γ), SU2(α, β, γ)
    x = config(σ1v, σ3v)
    q1, q2, q3 = rot(R3, x[1]), rot(R3, x[2]), rot(R3, x[3])
    φw = wigner_phase(Ur, R3, x[2])
    idx(λ) = λ > 0 ? 1 : 2
    for (λ0, λ2) in hels
        lhs = Mfull(q1, q2, q3, λ0, λ2)
        rec(s) = exp(s * im * φw * λ2) *
                 sum(conj(Ur[idx(λ0), idx(ν0)]) * Mfull(x..., ν0, λ2) for ν0 in (1 / 2, -1 / 2))
        @printf("  λ0=%+.1f λ2=%+.1f   |M(Rx) - reco(e^{-iφwλ2})| = %.2e   (with e^{+iφwλ2}: %.2e)\n",
            λ0, λ2, abs(lhs - rec(-1)), abs(lhs - rec(+1)))
    end
    # specialization R = Rz(π), proton along z:
    φpi = wigner_phase(SU2(Float64(π), 0.0, 0.0), rotmat(Float64(π), 0.0, 0.0), x[2] .* 0 .+ [x[2][1], 0.0, 0.0, x[2][4]])
    @printf("  φ_W(Rz(π), p̂2=ẑ) = %.6f  (= π: the λ2 part of the test-B phase)\n", φpi)
end

println("\n=== Test F: closure of the orientation decomposition under p1<->p3 ===")
# All ingredients extracted from the momenta alone. Claims:
#   (ii)  R(x̃) = R(x)·Rz(π)        (body frame flips about the proton axis)
#   (iii) φ_W(x̃) = φ_W(x) + π
#   reconstruction:  M(x) = e^{-iφ_W λ2} Σ_ν conj(D_{λ0ν}) A_{νλ2}(σ1,σ3)
#   closure:         same formula at x̃ (swapped invariants, shifted R, φ_W)
#                    reproduces the SAME number: the two (-1)^(ν-λ2) cancel.
cross3(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]

function orientation(p1, p2, p3)  # R3 such that (p1,p2,p3) = R3 · standard config
    zax = normalize(p2[2:4])
    yax = normalize(cross3(p2[2:4], p1[2:4]))
    R3 = hcat(cross3(yax, zax), yax, zax)
    β = acos(clamp(R3[3, 3], -1, 1))
    α = atan(R3[2, 3], R3[1, 3])
    γ = atan(R3[3, 2], -R3[3, 1])
    R3, SU2(α, β, γ)
end

let (α0, β0, γ0) = (1.9, 0.6, 0.8)   # generic event, orientation to be re-extracted
    R0 = rotmat(α0, β0, γ0)
    xstd = config(σ1v, σ3v)
    x = map(p -> rot(R0, p), xstd)            # generic event x
    xt = (x[3], x[2], x[1])                    # slot-swapped event x̃
    R3x, Ux = orientation(x...)
    R3t, Ut = orientation(xt...)
    @printf("  (ii)  |R(x̃) - R(x)·Rz(π)|            = %.2e\n",
        maximum(abs.(R3t - R3x * rotmat(Float64(π), 0.0, 0.0))))
    p2std = [x[2][1], 0.0, 0.0, norm(x[2][2:4])]   # proton of the standard config (along z)
    φx = wigner_phase(Ux, R3x, p2std)
    φt = wigner_phase(Ut, R3t, p2std)
    @printf("  (iii) φ_W(x̃) - φ_W(x) mod 2π         = %.6f  (= π)\n",
        mod(φt - φx, 2π))
    # aligned amplitude from the model itself
    A(σa, σb, ν, λ2) = Mfull(config(σa, σb)..., ν, λ2)
    idx(λ) = λ > 0 ? 1 : 2
    reco(U, φ, σa, σb, λ0, λ2) = exp(-im * φ * λ2) *
        sum(conj(U[idx(λ0), idx(ν)]) * A(σa, σb, ν, λ2) for ν in (1 / 2, -1 / 2))
    for (λ0, λ2) in hels
        m = Mfull(x..., λ0, λ2)
        r1 = reco(Ux, φx, σ1v, σ3v, λ0, λ2)   # reconstruction at x
        r2 = reco(Ut, φt, σ3v, σ1v, λ0, λ2)   # reconstruction at x̃: everything shifted
        @printf("  λ0=%+.1f λ2=%+.1f   |M(x)-reco(x)| = %.2e   |M(x)-reco(x̃)| = %.2e\n",
            λ0, λ2, abs(m - r1), abs(m - r2))
    end
    println("  closure: reco(x̃) uses (σ3,σ1), γ+π, φ_W+π — and lands on the same M.")
end

println("\n=== Test G: convention map — proton along -z gives the published + exponent ===")
# The LHCb/DPD formula uses (-1)^(M_Xib + λp). Relative to the toy standard
# configuration above, rotating the aligned configuration by Rx(π) puts the
# proton along -z and flips the λp contribution to the cross-chain phase.
let Rminus = rotx(Float64(π))
    xs = map(p -> rot(Rminus, p), config(σ1v, σ3v))
    ys = map(p -> rot(Rminus, p), config(σ3v, σ1v))
    for (λ0, λ2) in hels
        r = T2(xs[1], xs[2], xs[3], λ0, λ2) / T1(ys[1], ys[2], ys[3], λ0, λ2)
        @printf("  λ0=%+.1f λ2=%+.1f   ratio = %+.6f %+.6f im\n",
            λ0, λ2, real(r), imag(r))
    end
    println("  compare: (-1)^(λ0+λ2) pattern = (-1, +1, +1, -1)")
end
