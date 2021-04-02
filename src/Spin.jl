"""
    GAMBAR

Gyromagnetic ratio for ¹H with units Hz/G.
"""
const GAMBAR = 4258

"""
    GAMMA

Gyromagnetic ratio for ¹H with units rad/s/G.
"""
const GAMMA  = 2π * GAMBAR

const DualFloat64 = Union{Float64,<:ForwardDiff.Dual{T,Float64,N} where {T,N}}

struct Magnetization{T<:Real}
    x::T
    y::T
    z::T

    function Magnetization(x, y, z)

        args = promote(x, y, z)
        T = typeof(args[1])
        new{T}(args...)

    end
end

Base.eltype(::Magnetization{T}) where {T} = T
Base.convert(::Type{Magnetization{T}}, m::Magnetization) where {T} = Magnetization(T.((m.x, m.y, m.z))...)

struct Position{T<:Real}
    x::T
    y::T
    z::T

    function Position(x, y, z)

        args = promote(x, y, z)
        T = typeof(args[1])
        new{T}(args...)

    end
end

Base.convert(::Type{Position{T}}, p::Position) where {T} = Position(T.((p.x, p.y, p.z))...)

"""
    AbstractSpin

Abstract type for representing individual spins.
"""
abstract type AbstractSpin end

"""
    Spin([M,] M0, T1, T2, Δf[, pos]) <: AbstractSpin

Create an object that represents a single spin.

# Properties
- `M::Vector{Float64} = [0.0, 0.0, M0]`: Magnetization vector [Mx, My, Mz]
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant (ms)
- `T2::Real`: Spin-spin recovery time constant (ms)
- `Δf::Real`: Off-resonance frequency (Hz)
- `pos::Vector{<:Real} = [0, 0, 0]`: Spatial location [x, y, z] (cm)
- `signal::Complex{Float64}`: Signal produced by the spin

# Examples
```jldoctest
julia> spin = Spin([1.0, 2.0, 3.0], 1, 1000, 100, 3); spin.signal
1.0 + 2.0im
```
"""
#struct Spin{T<:Union{Float64,<:ForwardDiff.Dual}} <: AbstractSpin
#    M::Vector{T}
#    M0::T
#    T1::T
#    T2::T
#    Δf::T
#    pos::Vector{T}
#
#    # Default constructor with optional argument pos
#    Spin(M::Vector{<:Real}, M0, T1, T2, Δf, pos = [0,0,0]) = begin
#        T = promote_type(eltype(M), typeof(M0), typeof(T1), typeof(T2),
#                         typeof(Δf), eltype(pos))
#        T = T <: ForwardDiff.Dual ? float(T) : Float64
#        new{T}(M, M0, T1, T2, Δf, pos)
#    end
#
#    # If magnetization vector is not specified then use equilibrium
#    Spin(M0::Real, T1, T2, Δf, pos = [0,0,0]) =
#        Spin([0,0,M0], M0, T1, T2, Δf, pos)
#end

struct Spin{T<:DualFloat64} <: AbstractSpin
    M::Magnetization{T}
    M0::T
    T1::T
    T2::T
    Δf::T
    pos::Position{Float64}

    # Default constructor with optional argument pos
    function Spin(
        M::Magnetization,
        M0,
        T1,
        T2,
        Δf,
        pos::Position = Position(0.0, 0.0, 0.0)
    )

        args = (M0, T1, T2, Δf)
        T = promote_type(Float64, eltype(M), typeof.(args)...)
        new{T}(M, args..., pos)

    end
end

# If magnetization vector is not specified then use equilibrium
Spin(M0, T1, T2, Δf, pos::Position = Position(0.0, 0.0, 0.0)) =
    Spin(Magnetization(0, 0, M0), M0, T1, T2, Δf, pos)

"""
    SpinMC([M,] M0, frac, T1, T2, Δf, τ[, pos]) <: AbstractSpin

Create an object that represents a single spin with multiple compartments.

# Properties
- `N::Integer = length(frac)`: Number of compartments
- `Meq::Vector{Float64} = vcat([[0, 0, frac[n] * M0] for n = 1:N]...)`:
    Equilibrium magnetization vector
- `M::Vector{Float64} = Meq`: Magnetization vector
    [M1x, M1y, M1z, M2x, M2y, M2z, ...]
- `M0::Real`: Equilibrium magnetization
- `frac::Vector{<:Real}`: Volume fraction of each compartment
- `T1::Vector{<:Real}`: Spin-lattice recovery time constants (ms)
- `T2::Vector{<:Real}`: Spin-spin recovery time constants (ms)
- `Δf::Vector{<:Real}`: Off-resonance frequencies (Hz)
- `τ::Vector{<:Real}`: Residence times (inverse exchange rates) (ms)
    [τ12, τ13, ..., τ1N, τ21, τ23, ..., τ2N, ...]
- `pos::Vector{<:Real} = [0, 0, 0]`: Spatial location [x, y, z] (cm)
- `A::Matrix{Float64} = ...`: Matrix of system dynamics; see slide 22 in
    https://web.stanford.edu/class/rad229/Notes/B1-Bloch-Simulations.pdf
- `signal::Complex{Float64}`: Signal produced by the spin

# Examples
```jldoctest
julia> M = [1.0, 2.0, 3.0, 0.2, 0.1, 1.0];

julia> frac = [0.8, 0.2];

julia> τ = [Inf, Inf];

julia> spin = SpinMC(M, 1, frac, [900, 400], [80, 20], [3, 13], τ); spin.signal
1.2 + 2.1im
```
"""
#struct SpinMC{T<:Union{Float64,<:ForwardDiff.Dual}} <: AbstractSpin
#    N::Int
#    M::Vector{T} # [3N]
#    Meq::Vector{T} # [3N]
#    M0::T
#    frac::Vector{T} # [N]
#    T1::Vector{T} # [N]
#    T2::Vector{T} # [N]
#    Δf::Vector{T} # [N]
#    τ::Vector{T} # [N*(N-1)]
#    pos::Vector{T}
#    A::Matrix{T} # [3N,3N]
#
#    SpinMC(M::Vector{<:Real}, M0, frac, T1, T2, Δf, τ, pos = [0,0,0]) = begin
#        T = promote_type(eltype(M), typeof(M0), eltype(frac), eltype(T1),
#                         eltype(T2), eltype(Δf), eltype(τ), eltype(pos))
#        T = T <: ForwardDiff.Dual ? float(T) : Float64
#        N = length(frac)
#        Meq = vcat([[0, 0, frac[n] * M0] for n = 1:N]...)
#        r = zeros(T, N, N) # 1/ms
#        itmp = 1
#        for j = 1:N, i = 1:N
#            if i != j
#                r[i,j] = 1 / τ[itmp]
#                itmp += 1
#            end
#        end
#        A = zeros(T, 3N, 3N)
#        for j = 1:N, i = 1:N
#            ii = 3i-2:3i
#            jj = 3j-2:3j
#            if i == j
#                tmp = sum(r[:,i]) # 1/ms
#                r1 = -1 / T1[i] - tmp # 1/ms
#                r2 = -1 / T2[i] - tmp # 1/ms
#                Δω = 2π * Δf[i] / 1000 # rad/ms
#                A[ii,jj] = [r2 Δω 0; -Δω r2 0; 0 0 r1] # Left-handed rotation
#            else
#                A[ii,jj] = r[i,j] * Diagonal(ones(Bool, 3))
#            end
#        end
#        new{T}(N, M, Meq, M0, frac, T1, T2, Δf, τ, pos, A)
#    end
#
#    # If magnetization vector is not specified then use equilibrium
#    SpinMC(M0::Real, frac, T1, T2, Δf, τ, pos = [0,0,0]) = begin
#        M = vcat([[0, 0, frac[n] * M0] for n = 1:length(frac)]...)
#        SpinMC(M, M0, frac, T1, T2, Δf, τ, pos)
#    end
#end

struct SpinMC{T<:DualFloat64,N} <: AbstractSpin
    M::NTuple{N,Magnetization{T}}
    Meq::NTuple{N,Magnetization{T}}
    M0::T
    frac::NTuple{N,T}
    T1::NTuple{N,T}
    T2::NTuple{N,T}
    Δf::NTuple{N,T}
    r::NTuple{N,NTuple{N,T}}
    pos::Position{Float64}

    function SpinMC(
        M::NTuple{N,Magnetization},
        M0,
        frac,
        T1,
        T2,
        Δf,
        τ,
        pos::Position = Position(0, 0, 0)
    ) where {N}

        N > 1 || error(sprint(print, "SpinMC expects 2 or more compartments, got ", N))
        Meq = ntuple(i -> Magnetization(0, 0, frac[i] * M0), Val(N))
        (frac, T1, T2, Δf) = Tuple.((frac, T1, T2, Δf))
        τ = promote(τ...)
        itmp = 1
        r = ntuple(Val(N)) do j
            ntuple(Val(N)) do i
                out = 1 / τ[min(itmp, N * (N - 1))]
                if i == j
                    out = zero(out)
                else
                    itmp += 1
                end
                out
            end
        end
        T = promote_type(Float64, eltype(M[1]), eltype(Meq[1]), typeof(M0),
            typeof.(frac)..., typeof.(T1)..., typeof.(T2)..., typeof.(Δf)...,
            eltype(eltype(r)))
        new{T,N}(M, Meq, M0, frac, T1, T2, Δf, r, pos)

    end
end

SpinMC(M0::Real, frac, T1, T2, Δf, τ, pos::Position = Position(0.0, 0.0, 0.0)) =
    SpinMC(ntuple(i -> Magnetization(0, 0, frac[i] * M0), length(frac)), M0, frac, T1, T2, Δf, τ, pos)

# Override Base.getproperty to allow users to type spin.signal to compute the
# signal generated by the spin
Base.getproperty(spin::Spin, s::Symbol) = begin
    if s == :signal
        M = getfield(spin, :M)
        return complex(M[1], M[2])
    else
        return getfield(spin, s)
    end
end

Base.getproperty(spin::SpinMC, s::Symbol) = begin
    if s == :signal
        M = getfield(spin, :M)
        return complex(sum(M[1:3:end]), sum(M[2:3:end]))
    else
        return getfield(spin, s)
    end
end

struct Gradient{T<:Real}
    x::T
    y::T
    z::T

    function Gradient(x, y, z)

        args = promote(x, y, z)
        T = typeof(args[1])
        new{T}(args...)

    end
end

function gradient_frequency(grad::Gradient, pos::Position)

    return GAMBAR * (grad.x * pos.x + grad.y * pos.y + grad.z * pos.z) # Hz

end

abstract type AbstractSpoiling end

struct GradientSpoiling{T<:Real} <: AbstractSpoiling
    gradient::Gradient{T}
end

GradientSpoiling(x, y, z) = GradientSpoiling(Gradient(x, y, z))

struct RFSpoiling{T<:Real} <: AbstractSpoiling
    Δθ::T
end

struct RFandGradientSpoiling{T1<:Real,T2<:Real} <: AbstractSpoiling
    gradient::GradientSpoiling{T1}
    rf::RFSpoiling{T2}
end

RFandGradientSpoiling(grad::Gradient, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad), rf)
RFandGradientSpoiling(grad::NTuple{3,Real}, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad...), rf)
RFandGradientSpoiling(gx, gy, gz, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(gx, gy, gz), rf)
RFandGradientSpoiling(grad::Union{<:GradientSpoiling,<:Gradient,<:NTuple{3,Real}}, Δθ) = RFandGradientSpoiling(grad, RFSpoiling(Δθ))
RFandGradientSpoiling(rf::Union{<:RFSpoiling,<:Real}, grad::Union{<:GradientSpoiling,<:Gradient,<:NTuple{3,Real}}) = RFandGradientSpoiling(grad, rf)
RFandGradientSpoiling(rf::RFSpoiling, gx, gy, gz) = RFandGradientSpoiling(gx, gy, gz, rf)

spoiler_gradient(s::GradientSpoiling) = s.gradient
spoiler_gradient(s::RFandGradientSpoiling) = spoiler_gradient(s.gradient)
rfspoiling_increment(s::RFSpoiling) = s.Δθ
rfspoiling_increment(s::RFandGradientSpoiling) = rfspoiling_increment(s.rf)

struct BlochMcConnellWorkspace{T<:Real}
    A::Matrix{T}
end

"""
    freeprecess(spin, t)

Simulate free-precession for the given spin.

# Arguments
- `spin::AbstractSpin`: Spin that is free-precessing
- `t::Real`: Duration of free-precession (ms)

# Return
- `A::Matrix`: Matrix that describes relaxation and precession
- `B::Vector`: Vector that describes recovery

# Examples
```jldoctest
julia> spin = Spin([1.0, 0.0, 0.0], 1, 1000, 100, 3.75)
Spin([1.0, 0.0, 0.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, B) = freeprecess(spin, 100); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404048
```
"""
freeprecess(spin::Spin, t::Real) =
    freeprecess(t, spin.M0, spin.T1, spin.T2, spin.Δf)

function freeprecess(spin::SpinMC, t::Real)

    E = expm(t * spin.A)
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess!(A, B, spin::Spin, t, workspace::Nothing = nothing)

    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf)

end

function freeprecess!(A, B, spin::SpinMC, t, workspace)

    expm!(A, workspace, spin, t)
    mul!(B, Diagonal(ones(Bool, size(A, 1))) - A, spin.Meq)
    return nothing

end

"""
    freeprecess(spin, t, grad)

Simulate free-precession for the given spin in the presence of a gradient.

# Arguments
- `spin::AbstractSpin`: Spin that is free-precessing
- `t::Real`: Duration of free-precession (ms)
- `grad::AbstractVector{<:Real}`: Gradient amplitudes [gx, gy, gz] (G/cm)

# Return
- `A::Matrix`: Matrix that describes relaxation and precession
- `B::Vector`: Vector that describes recovery

# Examples
```jldoctest
julia> spin = Spin([1.0, 0.0, 0.0], 1, 1000, 100, 0, [0, 0, 3.75])
Spin([1.0, 0.0, 0.0], 1.0, 1000.0, 100.0, 0.0, [0.0, 0.0, 3.75])

julia> (A, B) = freeprecess(spin, 100, [0, 0, 1/GAMBAR]); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404048
```
"""
function freeprecess(spin::Spin, t::Real, grad::AbstractArray{<:Real,1})

    gradfreq = GAMBAR * sum(grad[i] * spin.pos[i] for i = 1:3) # Hz
    freeprecess(t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

# See equation (6.9) in Gopal Nataraj's PhD thesis
function freeprecess(spin::SpinMC, t::Real, grad::AbstractArray{<:Real,1})

    gradfreq = GAMMA * sum(grad[i] * spin.pos[i] for i = 1:3) / 1000 # rad/ms
    ΔA = diagm(1 => repeat([gradfreq, 0, 0], spin.N), # Left-handed rotation
              -1 => repeat([-gradfreq, 0, 0], spin.N))[1:3spin.N,1:3spin.N]
    E = expm(t * (spin.A + ΔA))
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess!(A, B, spin::Spin, t, grad::Gradient, workspace::Nothing = nothing)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

function freeprecess!(A, B, spin::SpinMC, t, grad, workspace)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    expm!(A, workspace, spin, t, gradfreq)
    mul!(B, Diagonal(ones(Bool, size(A, 1))) - A, spin.Meq)
    return nothing

end

# Exact matrix exponential
function expm!(expAt, workspace, spin, t, gradfreq = 0)

    for j = 1:spin.N, i = 1:spin.N

        if i == j

            r_out = sum(spin.r[k,i] for k = 1:spin.N) # 1/ms
            r1 = -1 / spin.T1[i] - r_out # 1/ms
            r2 = -1 / spin.T2[i] - r_out # 1/ms
            Δω = 2π * (spin.Δf[i] + gradfreq) / 1000 # rad/ms

            workspace.A[3i-2,3j-2] = r2
            workspace.A[3i-1,3j-2] = -Δω
            workspace.A[3i,  3j-2] = 0
            workspace.A[3i-2,3j-1] = Δω
            workspace.A[3i-1,3j-1] = r2
            workspace.A[3i,  3j-1] = 0
            workspace.A[3i-2,3j]   = 0
            workspace.A[3i-1,3j]   = 0
            workspace.A[3i,  3j]   = r1

        else

            workspace.A[3i-2,3j-2] = spin.r[i,j]
            workspace.A[3i-1,3j-2] = 0
            workspace.A[3i,  3j-2] = 0
            workspace.A[3i-2,3j-1] = 0
            workspace.A[3i-1,3j-1] = spin.r[i,j]
            workspace.A[3i,  3j-1] = 0
            workspace.A[3i-2,3j]   = 0
            workspace.A[3i-1,3j]   = 0
            workspace.A[3i,  3j]   = spin.r[i,j]

        end

    end

    expAt .= expm(t * workspace.A)
    return nothing

end

# Approximate matrix exponential
# See page 2 of http://doi.org/10.1137/0714065
# (unnumbered equation with o(||E||) term)
function expm!(expAt, ::Nothing, spin, t, gradfreq = 0)

    for j = 1:spin.N, i = 1:spin.N

        if i == j

            r_out = sum(spin.r[k,i] for k = 1:spin.N) # 1/ms
            E1 = exp(-t * (1 / spin.T1[i] + r_out))
            E2 = exp(-t * (1 / spin.T2[i] + r_out))
            θ = 2π * (spin.Δf[i] + gradfreq) * t / 1000 # rad
            (s, c) = sincos(θ)
            E2c = E2 * c
            E2s = E2 * s

            expAt[3i-2,3j-2] = E2c
            expAt[3i-1,3j-2] = -E2s
            expAt[3i,  3j-2] = 0
            expAt[3i-2,3j-1] = E2s
            expAt[3i-1,3j-1] = E2c
            expAt[3i,  3j-1] = 0
            expAt[3i-2,3j]   = 0
            expAt[3i-1,3j]   = 0
            expAt[3i,  3j]   = E1

        else

            r_out_i = sum(spin.r[k,i] for k = 1:spin.N) # 1/ms
            r_out_j = sum(spin.r[k,j] for k = 1:spin.N) # 1/ms
            R1i = 1 / spin.T1[i] + r_out_i # 1/ms
            R1j = 1 / spin.T1[j] + r_out_j # 1/ms
            R2i = 1 / spin.T2[i] + r_out_i # 1/ms
            R2j = 1 / spin.T2[j] + r_out_j # 1/ms
            rji = spin.r[i,j]

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

            expAt[3i-2,3j-2] = tmpc
            expAt[3i-1,3j-2] = -tmps
            expAt[3i,  3j-2] = 0
            expAt[3i-2,3j-1] = tmps
            expAt[3i-1,3j-1] = tmpc
            expAt[3i,  3j-1] = 0
            expAt[3i-2,3j]   = 0
            expAt[3i-1,3j]   = 0
            expAt[3i,  3j]   = rji * E2i * (1 - E2ji) / R2ji

        end

    end

    return nothing

end

"""
    freeprecess!(spin, ...)

Apply free-precession to the given spin.
"""
function freeprecess!(spin::AbstractSpin, args...)

    (A, B) = freeprecess(spin, args...)
    applydynamics!(spin, A, B)

end

"""
    excitation(spin, θ, α)

Simulate instantaneous excitation with flip angle `α` about an axis that makes
angle `θ` with the positive x-axis.

# Arguments
- `spin::AbstractSpin`: Spin to excite
- `θ::Real`: Orientation of the axis about which to excite (rad)
- `α::Real`: Flip angle (rad)

# Return
- `A::Matrix`: Matrix that describes the excitation
- `B::Vector = zeros(length(spin.M))`: Not used, but included because other
    methods of `excitation` return a nontrivial value here

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, _) = excitation(spin, π/4, π/2); A * spin.M
3-element Array{Float64,1}:
  0.7071067811865476
 -0.7071067811865475
  6.123233995736766e-17
```
"""
function excitation(spin::Spin, θ::Real, α::Real)

    A = rotatetheta(θ, α)
    B = zeros(length(spin.M))
    return (A, B)

end

function excitation(spin::SpinMC, θ::Real, α::Real)

    A = kron(Diagonal(ones(Bool, spin.N)), rotatetheta(θ, α))
    B = zeros(length(spin.M))
    return (A, B)

end

"""
    excitation(spin, rf, Δθ, grad, dt)

Simulate non-instantaneous excitation using the hard pulse approximation.

# Arguments
- `spin::AbstractSpin`: Spin to excite
- `rf::Vector{<:Number}`: RF waveform (G); its magnitude determines the flip
    angle and its phase determines the axis of rotation
- `Δθ::Real`: Additional RF phase (e.g., for RF spoiling) (rad)
- `grad::Union{Matrix{<:Real},Vector{<:Real}}`: Gradients to play during
    excitation (G/cm); should be a 3-vector if the gradients are constant during
    excitation, otherwise it should be a 3×(length(rf)) matrix
- `dt::Real`: Time step (ms)

# Return
- `A::Matrix`: Matrix that describes excitation and relaxation
- `B::Vector`: Vector that describes excitation and relaxation
"""
function excitation(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                    grad::AbstractArray{<:Real,2}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    A = Diagonal(ones(Bool, length(spin.M)))
    B = zeros(length(spin.M))
    for t = 1:T
        (Af, Bf) = freeprecess(spin, dt/2, grad[:,t])
        (Ae, _) = excitation(spin, θ[t], α[t])
        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf
    end
    return (A, B)

end

# Excitation with constant gradient
function excitation(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                    grad::AbstractArray{<:Real,1}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    A = Diagonal(ones(Bool, length(spin.M)))
    B = zeros(length(spin.M))
    (Af, Bf) = freeprecess(spin, dt/2, grad)
    for t = 1:T
        (Ae, _) = excitation(spin, θ[t], α[t])
        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf
    end
    return (A, B)

end

"""
    excitation!(spin, ...)

Apply excitation to the given spin.
"""
function excitation!(spin::AbstractSpin, θ::Real, α::Real)

    (A, _) = excitation(spin, θ, α)
    applydynamics!(spin, A)

end

# Use this function if using RF spoiling (because A and B need to be
# recalculated for each TR, so directly modifying the magnetization should be
# faster in this case)
function excitation!(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                     grad::AbstractArray{<:Real,2}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    for t = 1:T
        (Af, Bf) = freeprecess(spin, dt/2, grad[:,t])
        (Ae, _) = excitation(spin, θ[t], α[t])
        applydynamics!(spin, Af, Bf)
        applydynamics!(spin, Ae)
        applydynamics!(spin, Af, Bf)
    end

end

function excitation!(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                     grad::AbstractArray{<:Real,1}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    (Af, Bf) = freeprecess(spin, dt/2, grad)
    for t = 1:T
        (Ae, _) = excitation(spin, θ[t], α[t])
        applydynamics!(spin, Af, Bf)
        applydynamics!(spin, Ae)
        applydynamics!(spin, Af, Bf)
    end

end

"""
    spoil(spin)

Simulate ideal spoiling (i.e., setting the transverse component of the spin's
magnetization to 0).

# Arguments
- `spin::AbstractSpin`: Spin to spoil

# Return
- `S::Matrix`: Matrix that describes ideal spoiling

# Examples
```jldoctest
julia> spin = Spin([1, 0.4, 5], 1, 1000, 100, 0)
Spin([1.0, 0.4, 5.0], 1.0, 1000.0, 100.0, 0.0, [0.0, 0.0, 0.0])

julia> S = spoil(spin); S * spin.M
3-element Array{Float64,1}:
 0.0
 0.0
 5.0
```
"""
spoil(spin::Spin) = [0 0 0; 0 0 0; 0 0 1]
spoil(spin::SpinMC) = kron(Diagonal(ones(Bool, spin.N)), [0 0 0; 0 0 0; 0 0 1])

"""
    spoil!(spin)

Apply ideal spoiling to the given spin.
"""
function spoil!(spin::Spin)

    spin.M[1:2] .= 0
    return nothing

end

function spoil!(spin::SpinMC)

    spin.M[1:3:end] .= 0
    spin.M[2:3:end] .= 0
    return nothing

end

"""
    combine(D...)

Combine the matrices and vectors that describe the dynamics of a spin into one
matrix and one vector.

# Arguments
- `D::Tuple{<:AbstractArray{<:Real,2},<:AbstractVector{<:Real}}...`: List of
    pairs of matrices and vectors, i.e., ((A1, B1), (A2, B2), ...), where the
    A's are matrices and the B's are vectors

# Return
- `A::Matrix`: Matrix that describes the spin dynamics
- `B::Vector`: Vector that describes the spin dynamics

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> D1 = excitation(spin, 0, π/2);

julia> D2 = freeprecess(spin, 100);

julia> (A, B) = combine(D1, D2); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404054
```
"""
function combine(D::Tuple{<:AbstractArray{<:Real,2},<:AbstractArray{<:Real,1}}...)

  (A, B) = D[1]
  for i = 2:length(D)
    (Ai, Bi) = D[i]
    A = Ai * A
    B = Ai * B + Bi
  end
  return (A, B)

end

"""
    applydynamics!(spin, A[, B])

Apply dynamics to the given spin.

# Arguments
- `spin::AbstractSpin`: Spin to which to apply dynamics
- `A::Matrix`: Matrix with dynamics
- `B::Vector = zeros(length(spin.M))`: Vector with dynamics

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, _) = excitation(spin, 0, π/2); applydynamics!(spin, A)

julia> (A, B) = freeprecess(spin, 100); applydynamics!(spin, A, B)

julia> spin.M
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404054
```
"""
function applydynamics!(spin::AbstractSpin, A::AbstractArray{<:Real,2},
                        B::AbstractArray{<:Real,1})

  spin.M[:] = A * spin.M + B
  return nothing

end

function applydynamics!(spin::AbstractSpin, A::AbstractArray{<:Real,2})

  spin.M[:] = A * spin.M
  return nothing

end

function applydynamics!(spin::AbstractSpin, BtoM, A, B)

    BtoM .= B
    mul!(BtoM, A, spin.M, 1, 1) # BtoM .= A * spin.M + BtoM
    spin.M .= BtoM
    return nothing

end
