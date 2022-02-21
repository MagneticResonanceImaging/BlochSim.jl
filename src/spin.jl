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

"""
    Position(x, y, z)

Create a mutable `Position` object representing a 3D location. Units are cm.

# Properties
- `x::Real`: x position
- `y::Real`: y position
- `z::Real`: z position
"""
mutable struct Position{T<:Real}
    x::T
    y::T
    z::T
end

Position(x, y, z) = Position(promote(x, y, z)...)

Base.show(io::IO, pos::Position) = print(io, "Position(", pos.x, ", ", pos.y, ", ", pos.z, ")")
Base.show(io::IO, ::MIME"text/plain", pos::Position{T}) where {T} =
    print(io, "Position{$T}:\n x = ", pos.x, " cm\n y = ", pos.y, " cm\n z = ", pos.z, " cm")

Base.convert(::Type{Position{T}}, p::Position) where {T} = Position(T(p.x), T(p.y), T(p.z))

Base.:(==)(p1::Position, p2::Position) = p1.x == p2.x && p1.y == p2.y && p1.z == p2.z

const DualFloat64 = Union{Float64,<:ForwardDiff.Dual,<:ReverseDiff.TrackedReal}

"""
    AbstractSpin

Abstract type for representing individual spins.
"""
abstract type AbstractSpin end

"""
    Spin([M], M0, T1, T2, Δf, [pos]) <: AbstractSpin

Create an object that represents a single spin.

# Properties
- `M::Magnetization = Magnetization(0, 0, M0)`: Magnetization vector
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant (ms)
- `T2::Real`: Spin-spin recovery time constant (ms)
- `Δf::Real`: Off-resonance frequency (Hz)
- `pos::Position = Position(0, 0, 0)`: Spatial location (cm)
- `N::Int = 1`: Number of compartments

# Examples
```jldoctest
julia> Spin(1, 1000, 100, 0, Position(1, 2, 3))
Spin{Float64}:
 M = Magnetization(0.0, 0.0, 1.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 0.0 Hz
 pos = Position(1.0, 2.0, 3.0) cm

julia> Spin(Magnetization(1, 2, 3), 1, 1000, 100, 0)
Spin{Float64}:
 M = Magnetization(1.0, 2.0, 3.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 0.0 Hz
 pos = Position(0.0, 0.0, 0.0) cm
```
"""
struct Spin{T<:DualFloat64} <: AbstractSpin
    M::Magnetization{T}
    M0::T
    T1::T
    T2::T
    Δf::T
    pos::Position{Float64}

    # Default constructor with optional argument pos
    function Spin{T}(
        M::Magnetization,
        M0,
        T1,
        T2,
        Δf,
        pos::Position = Position(0.0, 0.0, 0.0)
    ) where {T<:DualFloat64}

        new{T}(M, M0, T1, T2, Δf, pos)

    end
end

function Spin(M, M0, T1, T2, Δf, pos...)

    args = (M0, T1, T2, Δf)
    T = promote_type(Float64, eltype(M), typeof.(args)...)
    Spin{T}(M, args..., pos...)

end

# If magnetization vector is not specified then use equilibrium
Spin{T}(M0, T1, T2, Δf, pos::Position...) where {T} =
    Spin{T}(Magnetization(0.0, 0.0, M0), M0, T1, T2, Δf, pos...)
Spin(M0, T1, T2, Δf, pos::Position...) =
    Spin(Magnetization(0.0, 0.0, M0), M0, T1, T2, Δf, pos...)

function Base.show(io::IO, ::MIME"text/plain", spin::Spin{T}) where {T}

    print(io, "Spin{$T}:")
    print(io, "\n M = ", spin.M)
    print(io, "\n M0 = ", spin.M0)
    print(io, "\n T1 = ", spin.T1, " ms")
    print(io, "\n T2 = ", spin.T2, " ms")
    print(io, "\n Δf = ", spin.Δf, " Hz")
    print(io, "\n pos = ", spin.pos, " cm")

end

Base.eltype(::Spin{T}) where {T} = T

"""
    SpinMC([M], M0, frac, T1, T2, Δf, τ, [pos]) <: AbstractSpin

Create an object that represents a single spin with multiple compartments.

# Properties
- `M::MagnetizationMC = Meq`: Magnetization vector
- `Meq::MagnetizationMC = MagnetizationMC((0, 0, frac[1] * M0), ...)`:
  Equilibrium magnetization vector
- `M0::Real`: Equilibrium magnetization
- `frac::Tuple{<:Real}`: Water fraction of each compartment
- `T1::Tuple{<:Real}`: Spin-lattice recovery time constants (ms)
- `T2::Tuple{<:Real}`: Spin-spin recovery time constants (ms)
- `Δf::Tuple{<:Real}`: Off-resonance frequencies (Hz)
- `r::Tuple{Tuple{<:Real}}`: Exchange rates (1/ms); `r[i][j]` is the exchange
  rate from compartment `i` to compartment `j`
- `pos::Position = Position(0, 0, 0)`: Spatial location (cm)
- `N::Int = length(frac)`: Number of compartments

## Note
The `SpinMC` constructor takes `τ` (inverse exchange rate, or residence time) as
input, *not* `r`. Furthermore, `τ` is given as a `Tuple` with `N^2 - N`
elements, arranged like (τ12, τ13, ..., τ1N, τ21, τ23, ..., τ2N, ...).

# Examples
```jldoctest
julia> SpinMC(1, (0.2, 0.8), (400, 1000), (20, 100), (15, 0), (100, 25))
SpinMC{Float64,2}:
 M = MagnetizationMC((0.0, 0.0, 0.2), (0.0, 0.0, 0.8))
 M0 = 1.0
 frac = (0.2, 0.8)
 T1 = (400.0, 1000.0) ms
 T2 = (20.0, 100.0) ms
 Δf = (15.0, 0.0) Hz
 r = ((0.0, 0.01), (0.04, 0.0)) 1/ms
 pos = Position(0.0, 0.0, 0.0) cm
```
"""
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

    function SpinMC{T}(
        M::MagnetizationMC{S,N},
        M0,
        frac,
        T1,
        T2,
        Δf,
        τ,
        pos::Position = Position(0.0, 0.0, 0.0)
    ) where {T<:DualFloat64,S,N}

        N > 1 || error(sprint(print, "SpinMC expects 2 or more compartments, got ", N))
        Meq = MagnetizationMC(ntuple(i -> Magnetization(0, 0, frac[i] * M0), N)...)
        itmp = 0
        r = ntuple(N) do i
            ntuple(N) do j
                if i == j
                    out = zero(T)
                else
                    itmp += 1
                    out = 1 / τ[itmp]
                end
                out
            end
        end
        new{T,N}(M, Meq, M0, frac, T1, T2, Δf, r, pos)

    end
end

function SpinMC(M, M0, frac, T1, T2, Δf, τ, pos...)

    frac = promote(frac...)
    T1 = promote(T1...)
    T2 = promote(T2...)
    Δf = promote(Δf...)
    rtmp = promote((1 ./ τ)...)
    T = promote_type(Float64, eltype(M), typeof(M0),
        eltype(frac), eltype(T1), eltype(T2), eltype(Δf), eltype(rtmp))
    SpinMC{T}(M, M0, frac, T1, T2, Δf, τ, pos...)

end

SpinMC{T}(M::NTuple{N,Magnetization}, args...) where {T,N} = SpinMC{T}(MagnetizationMC(M...), args...)
SpinMC{T}(M::NTuple{N,NTuple{3,Real}}, args...) where {T,N} = SpinMC{T}(MagnetizationMC(M...), args...)
SpinMC{T}(M0::Real, frac, T1, T2, Δf, τ, pos::Position...) where {T} =
    SpinMC{T}(MagnetizationMC([Magnetization(0.0, 0.0, frac[i] * M0) for i = 1:length(frac)]...), M0, frac, T1, T2, Δf, τ, pos...)
SpinMC(M::NTuple{N,Magnetization}, args...) where {N} = SpinMC(MagnetizationMC(M...), args...)
SpinMC(M::NTuple{N,NTuple{3,Real}}, args...) where {N} = SpinMC(MagnetizationMC(M...), args...)
SpinMC(M0::Real, frac, T1, T2, Δf, τ, pos::Position...) =
    SpinMC(MagnetizationMC((Magnetization(0.0, 0.0, frac[i] * M0) for i = 1:length(frac))...), M0, frac, T1, T2, Δf, τ, pos...)

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
    print(io, "\n T1 = ", spin.T1, " ms")
    print(io, "\n T2 = ", spin.T2, " ms")
    print(io, "\n Δf = ", spin.Δf, " Hz")
    print(io, "\n r = ", spin.r, " 1/ms")
    print(io, "\n pos = ", spin.pos, " cm")

end

Base.eltype(::SpinMC{T,N}) where {T,N} = T

signal(spin::AbstractSpin) = signal(spin.M)

function Base.getproperty(spin::Spin, s::Symbol)

    if s === :N
        return 1
    else
        return getfield(spin, s)
    end

end

function Base.getproperty(spin::SpinMC{T,N}, s::Symbol) where {T,N}

    if s === :N
        return N
    else
        return getfield(spin, s)
    end

end
