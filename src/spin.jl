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

function Base.copyto!(dst::Magnetization, src::Magnetization)

    dst.x = src.x
    dst.y = src.y
    dst.z = src.z
    return nothing

end

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

Base.copyto!(dst::MagnetizationMC{T,N}, src::MagnetizationMC{S,N}) where {S,T,N} =
    foreach(i -> copyto!(dst[i], src[i]), 1:N)

Base.eltype(::MagnetizationMC{T,N}) where {T,N} = T
Base.getindex(M::MagnetizationMC, i) = M.M[i]
Base.iterate(M::MagnetizationMC{T,N}, i = 1) where {T,N} =  i > N ? nothing : (M.M[i], i + 1)
Base.convert(::Type{MagnetizationMC{T,N}}, M::MagnetizationMC{S,N}) where {S,T,N} = MagnetizationMC((convert(Magnetization{T}, Mi) for Mi in M)...)
# This next definition of convert prevents StackOverflowErrors
Base.convert(::Type{MagnetizationMC{T,N}}, M::MagnetizationMC{T,N}) where {T,N} = M

Base.:(==)(M1::MagnetizationMC{T,N}, M2::MagnetizationMC{S,N}) where {S,T,N} = all(M1[i] == M2[i] for i = 1:N)
Base.isapprox(M1::MagnetizationMC{T,N}, M2::MagnetizationMC{S,N}; kwargs...) where {S,T,N} = all(isapprox(M1[i], M2[i]; kwargs...) for i = 1:N)

abstract type AbstractBlochMatrix{T<:Real} end

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
BlochMatrix{T}() where {T} = BlochMatrix(zero(T), zero(T), zero(T), zero(T),
                                    zero(T), zero(T), zero(T), zero(T), zero(T))
BlochMatrix() = BlochMatrix{Float64}()

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

struct BlochMcConnellMatrix{T<:Real,N} <: AbstractBlochMatrix{T}
    A::NTuple{N,NTuple{N,BlochMatrix{T}}}
end

function BlochMcConnellMatrix{T}(N) where {T}

    A = ntuple(i -> ntuple(i -> BlochMatrix{T}(), N), N)
    BlochMcConnellMatrix(A)

end

BlochMcConnellMatrix(N) = BlochMcConnellMatrix{Float64}(N)

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

# A = A * t
function LinearAlgebra.mul!(A::BlochDynamicsMatrix, t::Real)

    A.R1 *= t
    A.R2 *= t
    A.Δω *= t
    return nothing

end

function LinearAlgebra.mul!(E::ExchangeDynamicsMatrix, t::Real)

    E.r *= t
    return nothing

end

function LinearAlgebra.mul!(A::BlochMcConnellDynamicsMatrix, t::Real)

    for A in A.A
        mul!(A, t)
    end
    for E in A.E
        mul!(E, t)
    end

end

# C = A * B
function LinearAlgebra.mul!(C::BlochMatrix, A::BlochMatrix, B::BlochMatrix)

    C.a11 = A.a11 * B.a11 + A.a12 * B.a21 + A.a13 * B.a31
    C.a21 = A.a21 * B.a11 + A.a22 * B.a21 + A.a23 * B.a31
    C.a31 = A.a31 * B.a11 + A.a32 * B.a21 + A.a33 * B.a31
    C.a12 = A.a11 * B.a12 + A.a12 * B.a22 + A.a13 * B.a32
    C.a22 = A.a21 * B.a12 + A.a22 * B.a22 + A.a23 * B.a32
    C.a32 = A.a31 * B.a12 + A.a32 * B.a22 + A.a33 * B.a32
    C.a13 = A.a11 * B.a13 + A.a12 * B.a23 + A.a13 * B.a33
    C.a23 = A.a21 * B.a13 + A.a22 * B.a23 + A.a23 * B.a33
    C.a33 = A.a31 * B.a13 + A.a32 * B.a23 + A.a33 * B.a33
    return nothing

end

# C = A * B + C
function muladd!(C::BlochMatrix, A::BlochMatrix, B::BlochMatrix)

    C.a11 += A.a11 * B.a11 + A.a12 * B.a21 + A.a13 * B.a31
    C.a21 += A.a21 * B.a11 + A.a22 * B.a21 + A.a23 * B.a31
    C.a31 += A.a31 * B.a11 + A.a32 * B.a21 + A.a33 * B.a31
    C.a12 += A.a11 * B.a12 + A.a12 * B.a22 + A.a13 * B.a32
    C.a22 += A.a21 * B.a12 + A.a22 * B.a22 + A.a23 * B.a32
    C.a32 += A.a31 * B.a12 + A.a32 * B.a22 + A.a33 * B.a32
    C.a13 += A.a11 * B.a13 + A.a12 * B.a23 + A.a13 * B.a33
    C.a23 += A.a21 * B.a13 + A.a22 * B.a23 + A.a23 * B.a33
    C.a33 += A.a31 * B.a13 + A.a32 * B.a23 + A.a33 * B.a33
    return nothing

end

function LinearAlgebra.mul!(
    C::BlochMcConnellMatrix{T1,N},
    A::BlochMcConnellMatrix{T2,N},
    B::BlochMcConnellMatrix{T3,N}
) where {T1,T2,T3,N}

    for j = 1:N, i = 1:N
        mul!(C.A[i][j], A.A[i][1], B.A[1][j])
        for k = 2:N
            muladd!(C.A[i][j], A.A[i][k], B.A[k][j])
        end
    end

end

# Vector 1-norm
function absolutesum(A::BlochDynamicsMatrix)

    return 2 * abs(A.E2cosθ) + 2 * abs(A.E2sinθ) + abs(A.E1)

end

function absolutesum(E::ExchangeDynamicsMatrix)

    return 3 * abs(E.r)

end

function absolutesum(A::BlochMcConnellDynamicsMatrix{T,N,M}) where {T,N,M}

    result = absolutesum(A.A[1])
    for i = 2:N
        result += absolutesum(A.A[i])
    end
    for E in A.E
        result += absolutesum(E)
    end
    return result

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
        return complex(M.x, M.y)
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
