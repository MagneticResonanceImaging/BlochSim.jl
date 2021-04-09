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

mutable struct Magnetization{T<:Real}
    x::T
    y::T
    z::T

    function Magnetization(x, y, z)

        args = promote(x, y, z)
        T = typeof(args[1])
        new{T}(args...)

    end
end

Base.show(io::IO, M::Magnetization) = print(io, "[", M.x, ", ", M.y, ", ", M.z, "]")
Base.show(io::IO, ::MIME"text/plain", M::Magnetization{T}) where {T} =
    print(io, "Magnetization vector with eltype $T:\n Mx = ", M.x, "\n My = ", M.y, "\n Mz = ", M.z)

Base.zero(::Union{Magnetization{T},Type{Magnetization{T}}}) where {T} = Magnetization(zero(T), zero(T), zero(T))

Base.eltype(::Magnetization{T}) where {T} = T
Base.convert(::Type{Magnetization{T}}, M::Magnetization) where {T} = Magnetization(T(M.x), T(M.y), T(M.z))
Base.convert(::Type{Magnetization{T}}, M::Magnetization{T}) where {T} = M

Base.Vector(M::Magnetization) = [M.x, M.y, M.z]

Base.:(==)(M1::Magnetization, M2::Magnetization) = M1.x == M2.x && M1.y == M2.y && M1.z == M2.z
Base.isapprox(M1::Magnetization, M2::Magnetization; kwargs...) =
    isapprox(Vector(M1), Vector(M2); kwargs...)

struct MagnetizationMC{T<:Real,N}
    M::NTuple{N,Magnetization{T}}

    function MagnetizationMC(M::Magnetization...)

        N = length(M)
        N > 1 || error(sprint(print, "MagnetizationMC expects 2 or more compartments, got ", N))
        T = promote_type(eltype.(M)...)
        args = convert.(Magnetization{T}, M)
        new{T,N}(args)

    end
end

MagnetizationMC(M::NTuple{3,Real}...) = MagnetizationMC((Magnetization(Mi...) for Mi in M)...)

function MagnetizationMC(M::Real...)

    N = length(M)
    N % 3 == 0 || error("must specify Mx, My, and Mz for each magnetization vector")
    MagnetizationMC((Magnetization(M[i], M[i+1], M[i+2]) for i = 1:3:N)...)

end

function Base.show(io::IO, M::MagnetizationMC{T,N}) where {T,N}

    print(io, "[", M[1].x, ", ", M[1].y, ", ", M[1].z, ", ")
    for i = 2:N-1
        print(io, M[i].x, ", ", M[i].y, ", ", M[i].z, ", ")
    end
    print(io, M[N].x, ", ", M[N].y, ", ", M[N].z, "]")

end

function Base.show(io::IO, ::MIME"text/plain", M::MagnetizationMC{T,N}) where {T,N}

    print(io, "$N-compartment Magnetization vector with eltype $T:")
    for i = 1:N
        print(io, "\n Compartment ", i, ":")
        print(io, "\n  Mx = ", M[i].x)
        print(io, "\n  My = ", M[i].y)
        print(io, "\n  Mz = ", M[i].z)
    end

end

Base.zero(::Union{MagnetizationMC{T,N},Type{MagnetizationMC{T,N}}}) where {T,N} =
    MagnetizationMC((zero(Magnetization{T}) for i = 1:N)...)

Base.eltype(::MagnetizationMC{T,N}) where {T,N} = T
Base.getindex(M::MagnetizationMC, i) = M.M[i]
Base.iterate(M::MagnetizationMC{T,N}, i = 1) where {T,N} =  i > N ? nothing : (M.M[i], i + 1)
Base.convert(::Type{MagnetizationMC{T,N}}, M::MagnetizationMC{S,N}) where {S,T,N} = MagnetizationMC((convert(Magnetization{T}, Mi) for Mi in M)...)
# This next definition of convert prevents StackOverflowErrors
Base.convert(::Type{MagnetizationMC{T,N}}, M::MagnetizationMC{T,N}) where {T,N} = M

Base.:(==)(M1::MagnetizationMC{T,N}, M2::MagnetizationMC{S,N}) where {S,T,N} = all(M1[i] == M2[i] for i = 1:N)
Base.isapprox(M1::MagnetizationMC{T,N}, M2::MagnetizationMC{S,N}; kwargs...) where {S,T,N} = all(isapprox(M1[i], M2[i]; kwargs...) for i = 1:N)

abstract type AbstractBlochMatrix{T<:Real} end

# TODO: Not sure if needed
mutable struct BlochMatrix{T<:Real} <: AbstractBlochMatrix{T}
    a11::T
    a21::T
    a31::T
    a12::T
    a22::T
    a32::T
    a13::T
    a23::T
    a33::T
end

BlochMatrix(a...) = BlochMatrix(promote(a...)...)

mutable struct BlochDynamicsMatrix{T<:Real} <: AbstractBlochMatrix{T}
    R1::T
    R2::T
    Δω::T
end

BlochDynamicsMatrix{T}() where {T} = BlochDynamicsMatrix(zero(T), zero(T), zero(T))
BlochDynamicsMatrix() = BlochDynamicsMatrix{Float64}()
BlochDynamicsMatrix(R1, R2, Δω) = BlochDynamicsMatrix(promote(R1, R2, Δω)...)

Base.convert(::Type{BlochDynamicsMatrix{T}}, A::BlochDynamicsMatrix) where {T} =
    BlochDynamicsMatrix(T(A.R1), T(A.R2), T(A.Δω))
Base.convert(::Type{BlochDynamicsMatrix{T}}, A::BlochDynamicsMatrix{T}) where {T} = A

mutable struct FreePrecessionMatrix{T<:Real} <: AbstractBlochMatrix{T}
    E1::T
    E2cosθ::T
    E2sinθ::T
end

FreePrecessionMatrix{T}() where {T} = FreePrecessionMatrix(zero(T), zero(T), zero(T))
FreePrecessionMatrix() = FreePrecessionMatrix{Float64}()
FreePrecessionMatrix(E1, E2cosθ, E2sinθ) = FreePrecessionMatrix(promote(E1, E2cosθ, E2sinθ)...)

function Base.Matrix(A::FreePrecessionMatrix{T}) where {T}

    mat = Matrix{T}(undef, 3, 3)
    mat[1,1] = A.E2cosθ
    mat[2,1] = -A.E2sinθ
    mat[3,1] = zero(T)
    mat[1,2] = A.E2sinθ
    mat[2,2] = A.E2cosθ
    mat[3,2] = zero(T)
    mat[1,3] = zero(T)
    mat[2,3] = zero(T)
    mat[3,3] = A.E1

    return mat

end

mutable struct ExchangeDynamicsMatrix{T<:Real} <: AbstractBlochMatrix{T}
    r::T
end

ExchangeDynamicsMatrix{T}() where {T} = ExchangeDynamicsMatrix(zero(T))
ExchangeDynamicsMatrix() = ExchangeDynamicsMatrix{Float64}()

Base.convert(::Type{ExchangeDynamicsMatrix{T}}, A::ExchangeDynamicsMatrix) where {T} =
    ExchangeDynamicsMatrix(T(A.r))
Base.convert(::Type{ExchangeDynamicsMatrix{T}}, A::ExchangeDynamicsMatrix{T}) where {T} = A

mutable struct BlochMcConnellDynamicsMatrix{T<:Real,N,M} <: AbstractBlochMatrix{T}
    A::NTuple{N,BlochDynamicsMatrix{T}}
    E::NTuple{M,ExchangeDynamicsMatrix{T}}

    function BlochMcConnellDynamicsMatrix(
        A::NTuple{N,BlochDynamicsMatrix{T}},
        E::NTuple{M,ExchangeDynamicsMatrix{S}}
    ) where {M,N,S,T}

        M == N * (N - 1) || error("exchange rates must be defined for each pair of compartments")
        Tnew = promote_type(T, S)
        A = convert.(BlochDynamicsMatrix{Tnew}, A)
        E = convert.(ExchangeDynamicsMatrix{Tnew}, E)
        new{Tnew,N,M}(A, E)

    end
end

function BlochMcConnellDynamicsMatrix{T}(N) where {T}

    A = ntuple(i -> BlochDynamicsMatrix{T}(), N)
    E = ntuple(i -> ExchangeDynamicsMatrix{T}(), N * (N - 1))
    BlochMcConnellDynamicsMatrix(A, E)

end

BlochMcConnellDynamicsMatrix(N) = BlochMcConnellDynamicsMatrix{Float64}(N)

function Base.Matrix(A::BlochMcConnellDynamicsMatrix{T,N,M}) where {T,N,M}

    mat = Matrix{T}(undef, 3N, 3N)
    index = 0
    for j = 1:N, i = 1:N

        if i == j

            mat[3i-2,3j-2] = A.A[i].R2
            mat[3i-1,3j-2] = -A.A[i].Δω
            mat[3i  ,3j-2] = zero(T)
            mat[3i-2,3j-1] = A.A[i].Δω
            mat[3i-1,3j-1] = A.A[i].R2
            mat[3i  ,3j-1] = zero(T)
            mat[3i-2,3j]   = zero(T)
            mat[3i-1,3j]   = zero(T)
            mat[3i  ,3j]   = A.A[i].R1

        else

            index += 1
            mat[3i-2,3j-2] = A.E[index].r
            mat[3i-1,3j-2] = zero(T)
            mat[3i  ,3j-2] = zero(T)
            mat[3i-2,3j-1] = zero(T)
            mat[3i-1,3j-1] = A.E[index].r
            mat[3i  ,3j-1] = zero(T)
            mat[3i-2,3j]   = zero(T)
            mat[3i-1,3j]   = zero(T)
            mat[3i  ,3j]   = A.E[index].r

        end

    end

    return mat

end

#Base.:+(M1::Magnetization, M2::Magnetization) = Magnetization(M1.x + M2.x, M1.y + M2.y, M1.z + M2.z)
function add!(M1::Magnetization, M2::Magnetization)

    M1.x += M2.x
    M1.y += M2.y
    M1.z += M2.z
    return nothing

end

function Base.:*(A::AbstractMatrix, M::Magnetization)

    Mx = A[1,1] * M.x + A[1,2] * M.y + A[1,3] * M.z
    My = A[2,1] * M.x + A[2,2] * M.y + A[2,3] * M.z
    Mz = A[3,1] * M.x + A[3,2] * M.y + A[3,3] * M.z
    return Magnetization(Mx, My, Mz)

end

Base.:*(A::AbstractMatrix, M::MagnetizationMC{T,N}) where {T,N} = MagnetizationMC(ntuple(N) do i
    Mtmp = view(A, 3i-2:3i, 1:3) * M[1]
    for j = 2:N
        add!(Mtmp, view(A, 3i-2:3i, 3j-2:3j) * M[j])
    end
    Mtmp
end...)

# M2 = A * M1
function LinearAlgebra.mul!(M2::Magnetization, A::AbstractMatrix, M1::Magnetization)

    M2.x = A[1,1] * M1.x + A[1,2] * M1.y + A[1,3] * M1.z
    M2.y = A[2,1] * M1.x + A[2,2] * M1.y + A[2,3] * M1.z
    M2.z = A[3,1] * M1.x + A[3,2] * M1.y + A[3,3] * M1.z
    return nothing

end

# M2 = A * M1 + M2
function muladd!(M2::Magnetization, A::AbstractMatrix, M1::Magnetization)

    M2.x += A[1,1] * M1.x + A[1,2] * M1.y + A[1,3] * M1.z
    M2.y += A[2,1] * M1.x + A[2,2] * M1.y + A[2,3] * M1.z
    M2.z += A[3,1] * M1.x + A[3,2] * M1.y + A[3,3] * M1.z
    return nothing

end

# M2 = A * M1
function LinearAlgebra.mul!(M2::MagnetizationMC{T,N}, A::AbstractMatrix, M1::MagnetizationMC{S,N}) where {S,T,N}

    for i = 1:N
        M = M2[i]
        mul!(M, view(A, 3i-2:3i, 1:3), M1[1])
        for j = 2:N
            muladd!(M, view(A, 3i-2:3i, 3j-2:3j), M1[j])
        end
    end
    return nothing

end

# M2 = A * M1 + M2
function muladd!(M2::MagnetizationMC{T,N}, A::AbstractMatrix, M1::MagnetizationMC{S,N}) where {S,T,N}

    for i = 1:N
        M = M2[i]
        for j = 1:N
            muladd!(M, view(A, 3i-2:3i, 3j-2:3j), M1[j])
        end
    end
    return nothing

end

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
Base.show(io::IO, ::MIME"text/plain", pos::Magnetization{T}) where {T} =
    print(io, "Position{$T}:\n x = ", pos.x, "\n y = ", pos.y, "\n z = ", pos.z)

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
        return complex(M[1], M[2])
    else
        return getfield(spin, s)
    end
end

Base.getproperty(spin::SpinMC{T,N}, s::Symbol) where {T,N} = begin
    if s == :signal
        M = getfield(spin, :M)
        return complex(sum(M[1:3:end]), sum(M[2:3:end]))
    elseif s == :N
        return N
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

Base.show(io::IO, grad::Gradient) = print(io, "[", grad.x, ", ", grad.y, ", ", grad.z, "]")
Base.show(io::IO, ::MIME"text/plain", grad::Gradient{T}) where {T} =
    print(io, "Gradient{$T}:\n Gx = ", grad.x, " G/cm\n Gy = ", grad.y, " G/cm\n Gz = ", grad.z, " G/cm")

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

struct BlochMcConnellWorkspace{T<:Real,N}
    A::BlochMcConnellDynamicsMatrix{T,N}

    # N is number of compartments
    BlochMcConnellWorkspace(T::Type{<:Real}, N) = new{T,N}(BlochMcConnellDynamicsMatrix{T}(N))
end

BlochMcConnellWorkspace(::SpinMC{T,N}) where{T,N} = BlochMcConnellWorkspace(T, N)

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

function freeprecess(spin::SpinMC{T,N}, t::Real) where {T,N}

    A = Array{T}(undef, 3N, 3N)
    for j = 1:N, i = 1:N
        ii = 3i-2:3i
        jj = 3j-2:3j
        if i == j
            tmp = sum(spin.r[i][k] for k = 1:N) # 1/ms
            r1 = -1 / spin.T1[i] - tmp # 1/ms
            r2 = -1 / spin.T2[i] - tmp # 1/ms
            Δω = 2π * spin.Δf[i] / 1000 # rad/ms
            A[ii,jj] = [r2 Δω 0; -Δω r2 0; 0 0 r1] # Left-handed rotation
        else
            A[ii,jj] = spin.r[j][i] * Diagonal(ones(Bool, 3))
        end
    end
    E = expm(t * A)
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess!(A, B, spin::Spin, t, workspace::Nothing = nothing)

    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf)

end

function freeprecess!(A, B, spin::SpinMC, t, workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin))

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

    gradfreq = GAMBAR * (grad[1] * spin.pos.x + grad[2] * spin.pos.y + grad[3] * spin.pos.z) # Hz
    freeprecess(t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

# See equation (6.9) in Gopal Nataraj's PhD thesis
function freeprecess(spin::SpinMC{T,N}, t::Real, grad::AbstractArray{<:Real,1}) where {T,N}

    A = Array{T}(undef, 3N, 3N)
    for j = 1:N, i = 1:N
        ii = 3i-2:3i
        jj = 3j-2:3j
        if i == j
            tmp = sum(spin.r[i][k] for k = 1:N) # 1/ms
            r1 = -1 / spin.T1[i] - tmp # 1/ms
            r2 = -1 / spin.T2[i] - tmp # 1/ms
            Δω = 2π * spin.Δf[i] / 1000 # rad/ms
            A[ii,jj] = [r2 Δω 0; -Δω r2 0; 0 0 r1] # Left-handed rotation
        else
            A[ii,jj] = spin.r[j][i] * Diagonal(ones(Bool, 3))
        end
    end
    gradfreq = GAMMA * (grad[1] * spin.pos.x + grad[2] * spin.pos.y + grad[3] * spin.pos.z) / 1000 # rad/ms
    ΔA = diagm(1 => repeat([gradfreq, 0, 0], spin.N), # Left-handed rotation
              -1 => repeat([-gradfreq, 0, 0], spin.N))[1:3spin.N,1:3spin.N]
    E = expm(t * (A + ΔA))
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess!(A, B, spin::Spin, t, grad::Gradient, workspace::Nothing = nothing)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

function freeprecess!(A, B, spin::SpinMC, t, grad::Gradient, workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin))

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    expm!(A, workspace, spin, t, gradfreq)
    mul!(B, Diagonal(ones(Bool, size(A, 1))) - A, spin.Meq)
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

    # TODO: Make expm(t, A) to work directly with AbstractBlochMatrix
    expAt .= expm(t * Matrix(workspace.A))
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
    mul!(BtoM, A, spin.M, true, true) # BtoM .= A * spin.M + BtoM
    spin.M .= BtoM
    return nothing

end
