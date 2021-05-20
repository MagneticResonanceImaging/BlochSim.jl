mutable struct Magnetization{T<:Real}
    x::T
    y::T
    z::T
end

Magnetization(x, y, z) = Magnetization(promote(x, y, z)...)
Magnetization{T}() where {T} = Magnetization(zero(T), zero(T), zero(T))
Magnetization() = Magnetization(0.0, 0.0, 0.0)

Base.show(io::IO, M::Magnetization) = print(io, "Magnetization(", M.x, ", ", M.y, ", ", M.z, ")")
Base.show(io::IO, ::MIME"text/plain", M::Magnetization{T}) where {T} =
    print(io, "Magnetization vector with eltype $T:\n Mx = ", M.x, "\n My = ", M.y, "\n Mz = ", M.z)

Base.copy(M::Magnetization) = Magnetization(M.x, M.y, M.z)

function Base.copyto!(dst::Magnetization, src::Magnetization)

    dst.x = src.x
    dst.y = src.y
    dst.z = src.z
    return nothing

end

Base.eltype(::Magnetization{T}) where {T} = T
Base.convert(::Type{Magnetization{T}}, M::Magnetization) where {T} = Magnetization(T(M.x), T(M.y), T(M.z))

Base.Vector(M::Magnetization) = [M.x, M.y, M.z]

function Base.copyto!(dst::Magnetization, src::AbstractVector)

    dst.x = src[1]
    dst.y = src[2]
    dst.z = src[3]
    return nothing

end

function Base.copyto!(dst::AbstractVector, src::Magnetization)

    dst[1] = src.x
    dst[2] = src.y
    dst[3] = src.z
    return nothing

end

Base.:(==)(M1::Magnetization, M2::Magnetization) = M1.x == M2.x && M1.y == M2.y && M1.z == M2.z
Base.isapprox(M1::Magnetization, M2::Magnetization; kwargs...) =
    isapprox(Vector(M1), Vector(M2); kwargs...)

signal(M::Magnetization) = complex(M.x, M.y)

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
MagnetizationMC{T}(N) where {T} = MagnetizationMC(ntuple(i -> Magnetization{T}(), N)...)
MagnetizationMC(N) = MagnetizationMC(ntuple(i -> Magnetization(), N)...)

function Base.show(io::IO, M::MagnetizationMC{T,N}) where {T,N}

    print(io, "MagnetizationMC((", M[1].x, ", ", M[1].y, ", ", M[1].z, "), (")
    for i = 2:N-1
        print(io, M[i].x, ", ", M[i].y, ", ", M[i].z, "), (")
    end
    print(io, M[N].x, ", ", M[N].y, ", ", M[N].z, "))")

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

Base.copy(M::MagnetizationMC{T,N}) where {T,N} = MagnetizationMC(ntuple(i -> copy(M[i]), N)...)
Base.copyto!(dst::MagnetizationMC{T1,N}, src::MagnetizationMC{T2,N}) where {T1,T2,N} =
    foreach(i -> copyto!(dst[i], src[i]), 1:N)
Base.copyto!(dst::MagnetizationMC{T,N}, src::AbstractVector) where {T,N} =
    foreach(i -> copyto!(dst[i], view(src, 3i-2:3i)), 1:N)
Base.copyto!(dst::AbstractVector, src::MagnetizationMC{T,N}) where {T,N} =
    foreach(i -> copyto!(view(dst, 3i-2:3i), src[i]), 1:N)

Base.eltype(::MagnetizationMC{T,N}) where {T,N} = T
Base.getindex(M::MagnetizationMC, i) = M.M[i]
Base.iterate(M::MagnetizationMC{T,N}, i = 1) where {T,N} =  i > N ? nothing : (M.M[i], i + 1)
Base.convert(::Type{MagnetizationMC{T1,N}}, M::MagnetizationMC{T2,N}) where {T1,T2,N} = MagnetizationMC((convert(Magnetization{T1}, Mi) for Mi in M)...)
# This next definition of convert prevents StackOverflowErrors
Base.convert(::Type{MagnetizationMC{T,N}}, M::MagnetizationMC{T,N}) where {T,N} = M

function Base.Vector(M::MagnetizationMC{T,N}) where {T,N}

    v = Vector{T}(undef, 3N)
    for i = 1:N
        v[3i-2] = M[i].x
        v[3i-1] = M[i].y
        v[3i]   = M[i].z
    end
    return v

end

Base.:(==)(M1::MagnetizationMC{T1,N}, M2::MagnetizationMC{T2,N}) where {T1,T2,N} = all(M1[i] == M2[i] for i = 1:N)
Base.isapprox(M1::MagnetizationMC{T1,N}, M2::MagnetizationMC{T2,N}; kwargs...) where {T1,T2,N} = all(isapprox(M1[i], M2[i]; kwargs...) for i = 1:N)

signal(M::MagnetizationMC{T,N}) where {T,N} = sum(signal(M) for M in M)
