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

Base.show(io::IO, pos::Position) = print(io, "(", pos.x, ", ", pos.y, ", ", pos.z, ")")
Base.show(io::IO, ::MIME"text/plain", pos::Position{T}) where {T} =
    print(io, "Position{$T}:\n x = ", pos.x, " cm\n y = ", pos.y, " cm\n z = ", pos.z, " cm")

Base.convert(::Type{Position{T}}, p::Position) where {T} = Position(T(p.x), T(p.y), T(p.z))

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
struct Spin_old{T<:Union{Float64,<:ForwardDiff.Dual}} <: AbstractSpin
    M::Vector{T}
    M0::T
    T1::T
    T2::T
    Δf::T
    pos::Vector{T}

    # Default constructor with optional argument pos
    Spin_old(M::Vector{<:Real}, M0, T1, T2, Δf, pos = [0,0,0]) = begin
        T = promote_type(eltype(M), typeof(M0), typeof(T1), typeof(T2),
                         typeof(Δf), eltype(pos))
        T = T <: ForwardDiff.Dual ? float(T) : Float64
        new{T}(M, M0, T1, T2, Δf, pos)
    end

    # If magnetization vector is not specified then use equilibrium
    Spin_old(M0::Real, T1, T2, Δf, pos = [0,0,0]) =
        Spin_old([0,0,M0], M0, T1, T2, Δf, pos)
end

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

function Base.show(io::IO, ::MIME"text/plain", spin::Spin{T}) where {T}

    print(io, "Spin{$T}:")
    print(io, "\n M = ", spin.M)
    print(io, "\n M0 = ", spin.M0)
    print(io, "\n T1 = ", spin.T1)
    print(io, "\n T2 = ", spin.T2)
    print(io, "\n Δf = ", spin.Δf)
    print(io, "\n pos = ", spin.pos)

end

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
struct SpinMC_old{T<:Union{Float64,<:ForwardDiff.Dual}} <: AbstractSpin
    N::Int
    M::Vector{T} # [3N]
    Meq::Vector{T} # [3N]
    M0::T
    frac::Vector{T} # [N]
    T1::Vector{T} # [N]
    T2::Vector{T} # [N]
    Δf::Vector{T} # [N]
    τ::Vector{T} # [N*(N-1)]
    pos::Vector{T}
    A::Matrix{T} # [3N,3N]

    SpinMC_old(M::Vector{<:Real}, M0, frac, T1, T2, Δf, τ, pos = [0,0,0]) = begin
        T = promote_type(eltype(M), typeof(M0), eltype(frac), eltype(T1),
                         eltype(T2), eltype(Δf), eltype(τ), eltype(pos))
        T = T <: ForwardDiff.Dual ? float(T) : Float64
        N = length(frac)
        Meq = vcat([[0, 0, frac[n] * M0] for n = 1:N]...)
        r = zeros(T, N, N) # 1/ms
        itmp = 1
        for j = 1:N, i = 1:N
            if i != j
                r[i,j] = 1 / τ[itmp]
                itmp += 1
            end
        end
        A = zeros(T, 3N, 3N)
        for j = 1:N, i = 1:N
            ii = 3i-2:3i
            jj = 3j-2:3j
            if i == j
                tmp = sum(r[:,i]) # 1/ms
                r1 = -1 / T1[i] - tmp # 1/ms
                r2 = -1 / T2[i] - tmp # 1/ms
                Δω = 2π * Δf[i] / 1000 # rad/ms
                A[ii,jj] = [r2 Δω 0; -Δω r2 0; 0 0 r1] # Left-handed rotation
            else
                A[ii,jj] = r[i,j] * Diagonal(ones(Bool, 3))
            end
        end
        new{T}(N, M, Meq, M0, frac, T1, T2, Δf, τ, pos, A)
    end

    # If magnetization vector is not specified then use equilibrium
    SpinMC_old(M0::Real, frac, T1, T2, Δf, τ, pos = [0,0,0]) = begin
        M = vcat([[0, 0, frac[n] * M0] for n = 1:length(frac)]...)
        SpinMC_old(M, M0, frac, T1, T2, Δf, τ, pos)
    end
end

struct SpinMC{T<:DualFloat64,N} <: AbstractSpin
    M::MagnetizationMC{T,N}
    Meq::MagnetizationMC{T,N}
    M0::T
    frac::NTuple{N,T}
    T1::NTuple{N,T}
    T2::NTuple{N,T}
    Δf::NTuple{N,T}
    r::NTuple{N,NTuple{N,T}}
    pos::Position{Float64}

    function SpinMC(
        M::MagnetizationMC{S,N},
        M0,
        frac,
        T1,
        T2,
        Δf,
        τ,
        pos::Position = Position(0, 0, 0)
    ) where {S,N}

        N > 1 || error(sprint(print, "SpinMC expects 2 or more compartments, got ", N))
        Meq = MagnetizationMC((Magnetization(0, 0, frac[i] * M0) for i = 1:N)...)
        frac = promote(frac...)
        T1 = promote(T1...)
        T2 = promote(T2...)
        Δf = promote(Δf...)
        τ = promote(τ...)
        itmp = 1
        r = ntuple(N) do i
            ntuple(N) do j
                out = 1 / τ[min(itmp, N * (N - 1))]
                if i == j
                    out = zero(out)
                else
                    itmp += 1
                end
                out
            end
        end
        T = promote_type(Float64, eltype(M), eltype(Meq), typeof(M0),
            eltype(frac), eltype(T1), eltype(T2), eltype(Δf), eltype(eltype(r)))
        new{T,N}(M, Meq, M0, frac, T1, T2, Δf, r, pos)

    end
end

SpinMC(M::NTuple{N,Magnetization}, args...) where {N} = SpinMC(MagnetizationMC(M...), args...)
SpinMC(M::NTuple{N,NTuple{3,Real}}, args...) where {N} = SpinMC(MagnetizationMC(M...), args...)
SpinMC(M::NTuple{N,Real}, args...) where {N} = SpinMC(MagnetizationMC(M...), args...)
SpinMC(M0::Real, frac, T1, T2, Δf, τ, pos::Position = Position(0.0, 0.0, 0.0)) =
    SpinMC(MagnetizationMC((Magnetization(0, 0, frac[i] * M0) for i = 1:length(frac))...), M0, frac, T1, T2, Δf, τ, pos)

function Base.show(io::IO, spin::SpinMC{T,N}) where {T,N}

    print(io, "SpinMC{$T,$N}(", spin.M, ", ")
    print(io, spin.M0, ", ")
    print(io, spin.frac, ", ")
    print(io, spin.T1, ", ")
    print(io, spin.T2, ", ")
    print(io, spin.Δf, ", ")
    print(io, spin.r, ", ")
    print(io, spin.pos, ")")

end

function Base.show(io::IO, ::MIME"text/plain", spin::SpinMC{T,N}) where {T,N}

    print(io, "SpinMC{$T,$N}:")
    print(io, "\n M = ", spin.M)
    print(io, "\n M0 = ", spin.M0)
    print(io, "\n frac = ", spin.frac)
    print(io, "\n T1 = ", spin.T1)
    print(io, "\n T2 = ", spin.T2)
    print(io, "\n Δf = ", spin.Δf)
    print(io, "\n r = ", spin.r)
    print(io, "\n pos = ", spin.pos)

end

# Override Base.getproperty to allow users to type spin.signal to compute the
# signal generated by the spin
Base.getproperty(spin::Spin, s::Symbol) = begin
    if s == :signal
        M = getfield(spin, :M)
        return complex(M.x, M.y)
    elseif s == :N
        return 1
    else
        return getfield(spin, s)
    end
end

Base.getproperty(spin::SpinMC{T,N}, s::Symbol) where {T,N} = begin
    if s == :signal
        M = getfield(spin, :M)
        return complex(sum(M[i].x for i = 1:N), sum(M[i].y for i = 1:N))
    elseif s == :N
        return N
    else
        return getfield(spin, s)
    end
end
