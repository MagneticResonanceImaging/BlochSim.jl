struct BlochMcConnellWorkspace{T<:Real,N}
    A::BlochMcConnellDynamicsMatrix{T,N}
    expmworkspace::MatrixExponentialWorkspace{T,N}

    BlochMcConnellWorkspace(T::Type{<:Real}, N) =
        new{T,N}(BlochMcConnellDynamicsMatrix{T}(N),
                 MatrixExponentialWorkspace{T}(N))
end

BlochMcConnellWorkspace(::Type{SpinMC{T,N}}) where {T,N} = BlochMcConnellWorkspace(T, N)
BlochMcConnellWorkspace(spin::SpinMC) = BlochMcConnellWorkspace(typeof(spin))

"""
    freeprecess(t, M0, T1, T2, Δf)

Simulate free-precession, i.e., relaxation and off-resonance precession. Returns
`(A, B)` such that `A * M + B` applies free-precession to the magnetization `M`.

For an in-place version, see [`freeprecess!`](@ref).

# Arguments
- `t::Real`: Duration of free-precession (ms)
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant (ms)
- `T2::Real`: Spin-spin recovery time constant (ms)
- `Δf::Real`: Off-resonance frequency (Hz)

# Examples
```jldoctest
julia> (A, B) = freeprecess(100, 1, 1000, 100, 3.75); A * Magnetization(1, 0, 0) + B
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404048
```
"""
function freeprecess(t::Real, M0::Real, T1::Real, T2::Real, Δf::Real)

    A = FreePrecessionMatrix()
    B = Magnetization()
    freeprecess!(A, B, t, M0, T1, T2, Δf)
    return (A, B)

end

"""
    freeprecess!(A, B, t, M0, T1, T2, Δf)
    freeprecess!(A, B, spin, t, [nothing])
    freeprecess!(A, B, spinmc, t, [workspace])
    freeprecess!(A, B, spin, t, grad, [nothing])
    freeprecess!(A, B, spinmc, t, grad, [workspace])

Simulate free-precession, overwriting `A` and `B` (in-place version of
[`freeprecess`](@ref)).
"""
function freeprecess!(A, B, t, M0, T1, T2, Δf)

    A.E1 = exp(-t / T1)
    A.E2 = exp(-t / T2)
    A.θ = 2π * Δf * t / 1000

    B.x = 0
    B.y = 0
    B.z = M0 * (1 - A.E1)

    return nothing

end

"""
    freeprecess(spin, t, [nothing])
    freeprecess(spinmc, t, [workspace])
    freeprecess(spin, t, grad, [nothing])
    freeprecess(spinmc, t, grad, [workspace])

Simulate free-precession for the given spin for time `t` ms, optionally in the
presence of a B0 gradient. Returns `(A, B)` such that `A * M + B` applies
free-precession to the magnetization `M`.

For `SpinMC` objects, `workspace isa BlochMcConnellWorkspace`. Pass in `nothing`
instead to use an approximate matrix exponential to solve the Bloch-McConnell
equation.

For an in-place version, see [`freeprecess!`](@ref).

# Examples
```jldoctest
julia> s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 3.75)
Spin{Float64}:
 M = Magnetization(1.0, 0.0, 0.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 3.75 Hz
 pos = Position(0.0, 0.0, 0.0) cm

julia> (A, B) = freeprecess(s, 100); A * s.M + B
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404048

julia> s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 0, Position(0, 0, 3.75))
Spin{Float64}:
 M = Magnetization(1.0, 0.0, 0.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 0.0 Hz
 pos = Position(0.0, 0.0, 3.75) cm

julia> (A, B) = freeprecess(s, 100, Gradient(0, 0, 1/GAMBAR)); A * s.M + B
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404048
```
"""
function freeprecess(spin::Spin, t, ::Nothing = nothing)

    A = FreePrecessionMatrix{eltype(spin)}()
    B = Magnetization{eltype(spin)}()
    freeprecess!(A, B, spin, t)
    return (A, B)

end

function freeprecess(spin::SpinMC{T,N}, t, workspace::Union{Nothing,BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin)) where {T,N}

    A = BlochMcConnellMatrix{T}(N)
    B = MagnetizationMC{T}(N)
    freeprecess!(A, B, spin, t, workspace)
    return (A, B)

end

function freeprecess!(
    A::FreePrecessionMatrix,
    B::Magnetization,
    spin::Spin,
    t::Real,
    workspace::Nothing = nothing
)

    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf)

end

function freeprecess!(
    A::BlochMcConnellMatrix,
    B::MagnetizationMC,
    spin::SpinMC,
    t::Real,
    workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin)
)

    expm!(A, workspace, spin, t)
    subtractmul!(B, I, A, spin.Meq)
    return nothing

end

function freeprecess(spin::Spin, t::Real, grad::Gradient, ::Nothing = nothing)

    A = FreePrecessionMatrix{eltype(spin)}()
    B = Magnetization{eltype(spin)}()
    freeprecess!(A, B, spin, t, grad)
    return (A, B)

end

function freeprecess(spin::SpinMC{T,N}, t::Real, grad::Gradient, workspace = BlochMcConnellWorkspace(spin)) where {T,N}

    A = BlochMcConnellMatrix{T}(N)
    B = MagnetizationMC{T}(N)
    freeprecess!(A, B, spin, t, grad, workspace)
    return (A, B)

end

function freeprecess!(
    A::FreePrecessionMatrix,
    B::Magnetization,
    spin::Spin,
    t::Real,
    grad::Gradient,
    workspace::Nothing = nothing
)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

function freeprecess!(
    A::BlochMcConnellMatrix,
    B::MagnetizationMC,
    spin::SpinMC,
    t::Real,
    grad::Gradient,
    workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin)
)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    expm!(A, workspace, spin, t, gradfreq)
    subtractmul!(B, I, A, spin.Meq)
    return nothing

end

struct FreePrecessionWorkspace{T1,T2,T3}
    Af::T1
    Bf::T2
    tmpA::T1
    tmpB::T2
    bm_workspace::T3
end

function FreePrecessionWorkspace(
    spin::AbstractSpin,
    bm_workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    FreePrecessionWorkspace(typeof(spin), bm_workspace)

end

function FreePrecessionWorkspace(
    spin::Union{Type{Spin{T}},Type{SpinMC{T,N}}},
    bm_workspace = spin <: Spin ? nothing : BlochMcConnellWorkspace(spin)
) where{T,N}

    if spin <: Spin
        Af = FreePrecessionMatrix{T}()
        Bf = Magnetization{T}()
        tmpA = FreePrecessionMatrix{T}()
        tmpB = Magnetization{T}()
    else
        Af = BlochMcConnellMatrix{T}(N)
        Bf = MagnetizationMC{T}(N)
        tmpA = BlochMcConnellMatrix{T}(N)
        tmpB = MagnetizationMC{T}(N)
    end
    FreePrecessionWorkspace(Af, Bf, tmpA, tmpB, bm_workspace)

end

function freeprecess(spin::Spin, t, grads, workspace = FreePrecessionWorkspace(spin))

    A = FreePrecessionMatrix{eltype(spin)}()
    B = Magnetization{eltype(spin)}()
    freeprecess!(A, B, spin, t, grads, workspace)
    return (A, B)

end

function freeprecess(spin::SpinMC{T,N}, t, grads, workspace = FreePrecessionWorkspace(spin)) where {T,N}

    A = BlochMcConnellMatrix{T}(N)
    B = MagnetizationMC{T}(N)
    freeprecess!(A, B, spin, t, grads, workspace)
    return (A, B)

end

function freeprecess!(
    A,
    B,
    spin,
    t,
    grads, # Collection of Gradients
    workspace::FreePrecessionWorkspace = FreePrecessionWorkspace(spin)
)

    Δt = t / length(grads)
    make_identity!(workspace.tmpA)
    fill!(workspace.tmpB, 0)
    for grad in grads
        freeprecess!(workspace.Af, workspace.Bf, spin, Δt, grad, workspace.bm_workspace)
        combine!(A, B, workspace.tmpA, workspace.tmpB, workspace.Af, workspace.Bf)
        copyto!(workspace.tmpA, A)
        copyto!(workspace.tmpB, B)
    end
    return nothing

end

# Exact matrix exponential
function expm!(expAt, workspace, spin, t, gradfreq = 0)

    index = 0
    for i = 1:spin.N, j = 1:spin.N

        if i == j

            r_out = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms

            A = workspace.A.A[i]
            A.R1 = -1 / spin.T1[i] - r_out # 1/ms
            A.R2 = -1 / spin.T2[i] - r_out # 1/ms
            A.Δω = 2π * (spin.Δf[i] + gradfreq) / 1000 # rad/ms

        else

            index += 1
            workspace.A.E[index].r = spin.r[i][j] # 1/ms

        end

    end

    mul!(workspace.A, t)
    expm!(expAt, workspace.A, workspace.expmworkspace)
    return nothing

end

# Approximate matrix exponential
# See page 2 of http://doi.org/10.1137/0714065
# (unnumbered equation with o(||E||) term)
# Referring to A and E from the paper:
# - A is block diagonal where each block along the diagonal is the traditional
#   Bloch matrix for the corresponding compartment (except the 1/T1 and 1/T2
#   terms also include exchange out of the compartment)
# - E is an "inverse" block diagonal matrix (i.e., the diagonal blocks are 0)
#   where the off-diagonal blocks are diagonal matrices with the exchange terms
# With this choice of A, e^At has an analytical solution. The code in this
# function follows the analytical solution of the aforementioned equation after
# manually taking e^At and doing the matrix-matrix products and integrating
# over s.
function expm!(expAt, ::Nothing, spin, t, gradfreq = 0)

    for j = 1:spin.N, i = 1:spin.N

        A = expAt.A[i][j]

        if i == j

            r_out = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms
            E1 = exp(-t * (1 / spin.T1[i] + r_out))
            E2 = exp(-t * (1 / spin.T2[i] + r_out))
            θ = 2π * (spin.Δf[i] + gradfreq) * t / 1000 # rad
            (s, c) = sincos(θ)
            E2c = E2 * c
            E2s = E2 * s

            A.a11 = E2c
            A.a21 = -E2s
            A.a31 = 0
            A.a12 = E2s
            A.a22 = E2c
            A.a32 = 0
            A.a13 = 0
            A.a23 = 0
            A.a33 = E1

        else

            r_out_i = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms
            r_out_j = sum(spin.r[j][k] for k = 1:spin.N) # 1/ms
            R1i = 1 / spin.T1[i] + r_out_i # 1/ms
            R1j = 1 / spin.T1[j] + r_out_j # 1/ms
            R2i = 1 / spin.T2[i] + r_out_i # 1/ms
            R2j = 1 / spin.T2[j] + r_out_j # 1/ms
            rji = spin.r[j][i]

            R2ji = R2j - R2i
            R2ji² = R2ji^2
            E2i = exp(-t * R2i)
            E2ji = exp(-t * R2ji)
            Δωji = 2π * (spin.Δf[j] - spin.Δf[i]) / 1000 # rad/ms
            Δωji² = Δωji^2
            θi = 2π * (spin.Δf[i] + gradfreq) * t / 1000 # rad
            θj = 2π * (spin.Δf[j] + gradfreq) * t / 1000 # rad
            (si, ci) = sincos(θi)
            (sj, cj) = sincos(θj)
            tmpc = rji * E2i * ((ci - E2ji * cj) / R2ji + Δωji * (E2ji * sj - si) / R2ji²) / (1 + Δωji² / R2ji²)
            tmps = rji * E2i * ((si - E2ji * sj) / R2ji + Δωji * (ci - E2ji * cj) / R2ji²) / (1 + Δωji² / R2ji²)

            A.a11 = tmpc
            A.a21 = -tmps
            A.a31 = 0
            A.a12 = tmps
            A.a22 = tmpc
            A.a32 = 0
            A.a13 = 0
            A.a23 = 0
            A.a33 = rji * E2i * (1 - E2ji) / R2ji

        end

    end

    return nothing

end

"""
    freeprecess!(spin, ...)

Apply free-precession to the given spin, overwriting the spin's magnetization
vector.
"""
function freeprecess!(spin::AbstractSpin, args...)

    (A, B) = freeprecess(spin, args...)
    BtoM = copy(spin.M)
    applydynamics!(spin, BtoM, A, B)

end

"""
    combine!(A, B, A1, B1, A2, B2)
    combine!(A, A1, A2)

Combine the matrices and vectors that describe the dynamics of a spin into one
matrix and one vector, overwriting `A` and `B`. The dynamics described by `A1`
and `B1` apply first, then those described by `A2` and `B2`. In other words,
`A = A2 * A1` and `B = A2 * B1 + B2`.

# Examples
```jldoctest
julia> s = Spin(1, 1000, 100, 3.75);

julia> A = BlochMatrix(); B = Magnetization();

julia> (A1, B1) = excite(s, InstantaneousRF(π/2));

julia> (A2, B2) = freeprecess(s, 100);

julia> combine!(A, B, A1, B1, A2, B2); A * s.M + B
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404054
```
"""
function combine!(A, B, A1, B1, A2, B2)

    mul!(A, A2, A1)
    mul!(B, A2, B1)
    add!(B, B2)

end

function combine!(A, B, A1, B1, A2, ::Nothing)

    mul!(A, A2, A1)
    mul!(B, A2, B1)

end

function combine!(A, B, A1, ::Nothing, A2, B2)

    mul!(A, A2, A1)
    copyto!(B, B2)

end

function combine!(A, B, A1, B1, ::Nothing, ::Nothing)

    copyto!(A, A1)
    copyto!(B, B1)

end

function combine!(A, B, ::Nothing, ::Nothing, A2, B2)

    copyto!(A, A2)
    copyto!(B, B2)

end

combine!(A, A1, A2) = mul!(A, A2, A1)

"""
    applydynamics!(spin, BtoM, A, [B])

Apply dynamics to the given spin, overwriting the spin's magnetization vector.
`BtoM` is used to store intermediate results (and is thus overwritten).

# Examples
```jldoctest
julia> s = Spin(1, 1000, 100, 3.75); s.M
Magnetization vector with eltype Float64:
 Mx = 0.0
 My = 0.0
 Mz = 1.0

julia> BtoM = Magnetization();

julia> (A,) = excite(s, InstantaneousRF(π/2)); applydynamics!(s, BtoM, A); s.M
Magnetization vector with eltype Float64:
 Mx = 1.0
 My = 0.0
 Mz = 6.123233995736766e-17

julia> (A, B) = freeprecess(s, 100); applydynamics!(s, BtoM, A, B); s.M
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404054
```
"""
function applydynamics!(spin::AbstractSpin, BtoM, A, B)

    copyto!(BtoM, B)
    muladd!(BtoM, A, spin.M)
    copyto!(spin.M, BtoM)
    return nothing

end

function applydynamics!(spin::AbstractSpin, BtoM, A, ::Nothing = nothing)

    mul!(BtoM, A, spin.M)
    copyto!(spin.M, BtoM)
    return nothing

end
