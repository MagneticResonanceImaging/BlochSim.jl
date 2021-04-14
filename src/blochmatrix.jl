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

function LinearAlgebra.mul!(
    C::BlochMcConnellMatrix{T1,N},
    A::BlochMcConnellDynamicsMatrix{T2,N,M},
    B::BlochMcConnellDynamicsMatrix{T3,N,M}
) where {T1,T2,T3,N,M}

    # TODO: Finish
    for j = 1:N, i = 1:N
        mul!(C.A[i][j], A.A[i][1], B.A[1][j])
        for k = 2:N
            muladd!(C.A[i][j], A.A[i][k], B.A[k][j])
        end
    end

end

# Vector 1-norm
function absolutesum(A::BlochDynamicsMatrix)

    return 2 * abs(A.R2) + 2 * abs(A.Δω) + abs(A.R1)

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
