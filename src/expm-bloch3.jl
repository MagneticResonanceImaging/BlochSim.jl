#=
expm-bloch3.jl
Compute matrix exponential exp(A * t)
where A is the 3×3 Bloch matrix (for a constant RF and gradient amplitude),
using an explicit eigendecomposition
based on analytical roots of cubic characteristic polynomial.

This code was developed with the help of GPT 5.2.
=#

using LinearAlgebra: cond, norm, Diagonal, lu!, mul!, rdiv!

export expm_bloch3, expm_bloch3!

const Breal = Number # Real or ForwardDiff.Dual


#=
import LinearAlgebra: cross

"""
    cross(a::NTuple{3}, b::NTuple{3})
Return cross-product for 3-tuples;
avoids allocations.
"""
function cross(a::NTuple{3,Breal}, b::NTuple{3,Breal})
    a1, a2, a3 = a
    b1, b2, b3 = b
    return (a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1)
end
=#


"""
    Bloch3ExpmWorkspace{T<:Real}

Workspace for computing the matrix exponential
of a 3×3 Bloch matrix.
"""
struct Bloch3ExpmWorkspace{T <: Breal}
    row1::Vector{Complex{T}}
    row2::Vector{Complex{T}}
    row3::Vector{Complex{T}}
    λ::Vector{Complex{T}}
    V::Matrix{Complex{T}} # eigenvectors of A
    VD::Matrix{Complex{T}} # V * Diagonal(exp.(λ * t))
end

# constructor
Bloch3ExpmWorkspace{T}() where {T <: Breal} = Bloch3ExpmWorkspace(
    ntuple(_ -> Vector{Complex{T}}(undef, 3), 4)...,
    ntuple(_ -> Matrix{Complex{T}}(undef, 3, 3), 2)...,
)


"""
    cross!(v::Vector, a::AbstractVector, b::AbstractVector)
Compute cross-product `v = a ⊗ b` for 3-vectors;
avoids allocations by mutating input `v`.
"""
function cross!(
    v::AbstractVector{<:Number},
    a::AbstractVector{<:Number},
    b::AbstractVector{<:Number},
)
    v[1] = a[2] * b[3] - a[3] * b[2]
    v[2] = a[3] * b[1] - a[1] * b[3]
    v[3] = a[1] * b[2] - a[2] * b[1]
    return v
end


"""
    matrix_bloch3(r1, r2, w, s, c)

Return 3×3 Bloch matrix
`A = [-r2  w   s;
      -w  -r2  c;
      -s  -c  -r1]`
"""
matrix_bloch3(r1::T, r2::T, w::T, s::T, c::T) where {T <: Breal} = [
    -r2  w   s;
    -w  -r2  c;
    -s  -c  -r1]


"""
    eigvals_bloch3(r1, r2, w, s, c)

Return tuple of `Complex`-typed eigenvalues of 3×3 Bloch matrix
`A = [-r2  w   s;
      -w  -r2  c;
      -s  -c  -r1]
using Cardano/trig formulas,
with `cbrt` called only on real arguments.
"""
function eigvals_bloch3(r1::T, r2::T, w::T, s::T, c::T) where {T <: Breal}

    Tp = promote_type(T, Float32)

    # coefficients of det(λI - A) = λ^3 + a λ^2 + b λ + d = 0
    a = r1 + 2r2
    b = r2^2 + 2r1*r2 + w^2 + s^2 + c^2
    d = r2^2*r1 + r2*(s^2 + c^2) + r1*w^2

    # depressed cubic y^3 + p y + q = 0 with λ = y - a/3
    p = b - a^2 / 3
    q = 2a^3 / 27 - a*b / 3 + d

    Δ = (q/2)^2 + (p/3)^3

    #=
    Δ == 0 corresponds to corner case of real roots (-r2, -r2, -r1).
    Δ < 0 also corresponds to three real roots.
    Still, we return complex roots always for the typical cases,
    and for type stability.
    =#

    if Δ ≥ 0 # some complex roots
        sqrtΔ = sqrt(Δ)
        u = cbrt(-q/2 + sqrtΔ)
        v = cbrt(-q/2 - sqrtΔ)

        y1 = u + v
        λ1 = complex(y1 - a/3)

        re = -y1/2 - a/3
        ii = Tp(sqrt(3)/2) * (u - v)
        λ2 = complex(re,  ii)
        λ3 = complex(re, -ii)
        return (λ1, λ2, λ3)

    else # three real roots via trig
        ρ = 2*sqrt(-p/3)
        arg = (-q/2) / sqrt(-(p/3)^3)
        arg = clamp(arg, -1, 1) # clamp for numerical safety
        ϕ = acos(arg) / 3
        fun(k) = complex(ρ*cos(ϕ + Tp(2π/3)*(k-1)) - a/3)
        return ntuple(fun, Val(3))
    end
end


# helper for mismatched types
function eigvals_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc,
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal}
#   T = promote_type(T1, T2, Tw, Ts, Tc, Float32)
    r1, r2, w, s, c, _ = promote(r1, r2, w, s, c, zero(Float32))
    return eigvals_bloch3(r1, r2, w, s, c)
end


function fill3!(v::Vector, a, b, c)
    v[1] = a
    v[2] = b
    v[3] = c
end


"""
    eigvec_bloch3!(v, row1, row2, row3, r1, r2, w, s, c, λ)
Compute eigenvector `v` for eigenvalue λ by crossing two rows;
switch row pairs if needed.
Mutates `v` and workspace `row1` `row2` `row3`
"""
function eigvec_bloch3!(
    v::AbstractVector{Complex{T}},
    row1::Vector{Complex{T}},
    row2::Vector{Complex{T}},
    row3::Vector{Complex{T}},
    r1::Breal, r2::Breal, w::Breal, s::Breal, c::Breal, λ::Complex{<:Breal};
    rtol::Real = 1e-12,
) where {T <: Breal}

    fill3!(row1, -r2-λ, w,     s)
    fill3!(row2, -w,   -r2-λ,  c)
    fill3!(row3, -s,   -c,   -r1-λ)
    rmax = max(norm(row1), norm(row2), norm(row3))

    cross!(v, row1, row2)
    if norm(v) < rtol * rmax
        cross!(v, row1, row3)
    end
    if norm(v) < rtol * rmax
        cross!(v, row2, row3)
    end
    return v
end


#=
# even though this version used only small tuples, it still allocated

function eigvec_bloch3(
    r1::T, r2::T, w::T, s::T, c::T, λ::Complex{T};
    rtol::Real = 1e-12,
) where {T <: Breal}

    row1 = (-r2-λ, w,     s)
    row2 = (-w,   -r2-λ,  c)
    row3 = (-s,   -c,   -r1-λ)
    rmax = max(norm(row1), norm(row2), norm(row3))

    v = cross(row1, row2)
    if norm(v) < rtol * rmax
        v = cross(row1, row3)
    end
    if norm(v) < rtol * rmax
        v = cross(row2, row3)
    end
    return v
end


# helper for mismatched types
function eigvec_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc, λ::Complex{Tλ},
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal, Tλ <: Breal}
    T = promote_type(T1, T2, Tw, Ts, Tc, Tλ, Float32)
    return eigvec_bloch3(T(r1), T(r2), T(w), T(s), T(c), Complex{T}(λ))
end
=#


"""
    (λ, V) = eigen_bloch3!(λ, V, row1, row2, row3, r1, r2, w, s, c)
Return eigendecomposition `(λ, V)`
of 3×3 Bloch matrix
`A = [-r2  w   s;
      -w  -r2  c;
      -s  -c  -r1]`
Mutates `λ`, `V` and workspace `row1 row2 row3`, so allocation free.
"""
function eigen_bloch3!(
    λ::Vector{Complex{Tw}},
    V::Matrix{Complex{Tw}},
    row1::Vector{Complex{Tw}},
    row2::Vector{Complex{Tw}},
    row3::Vector{Complex{Tw}},
    r1::T, r2::T, w::T, s::T, c::T;
) where {Tw <: Breal, T <: Breal}

    λ[1], λ[2], λ[3] = eigvals_bloch3(r1, r2, w, s, c)

    eigvec_bloch3!((@view V[:,1]), row1, row2, row3, r1, r2, w, s, c, λ[1])
    eigvec_bloch3!((@view V[:,2]), row1, row2, row3, r1, r2, w, s, c, λ[2])
    eigvec_bloch3!((@view V[:,3]), row1, row2, row3, r1, r2, w, s, c, λ[3])

    return λ, V
end


#=
# helper for mismatched types
function eigen_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc,
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal}
    r1, r2, w, s, c, _ = promote(r1, r2, w, s, c, zero(Float32))
    return eigen_bloch3(r1, r2, w, s, c)
end
=#


"""
    expm_bloch3!(expAt, work, r1, r2, w, s, c, t)

Return exp(A*t) for the 3×3 Bloch matrix
`A = [-r2 w s; -w -r2 c; -s -c -r1].`
Uses explicit eigendecomposition
based on analytical roots of cubic characteristic polynomial.
Mutates `expAt` and `work::Bloch3ExpmWorkspace`.
Allocates very little memory; just some `lu!` overhead.
"""
function expm_bloch3!(
    expAt::Matrix{Te},
    work::Bloch3ExpmWorkspace{Tw},
    r1::T, r2::T, w::T, s::T, c::T, t::Tt;
) where {Te <: Breal, Tw <: Breal, T <: Breal, Tt <: Breal}

    eigen_bloch3!(work.λ, work.V, work.row1, work.row2, work.row3,
        r1, r2, w, s, c) # eigendecomposition of A

#=
    # If V is ill-conditioned (near defective/repeated eigenvalues), fall back
    # No! Cannot compute `cond` for dual types
    if cond(V) > 1e12 || any(isnan, V) || any(isinf, V)

        A = [-r2  w  s;
             -w -r2  c;
             -s -c  -r1]
        return exp(t*A)  # Julia's matrix exponential
    end
=#

    # compute V * D(exp(λ t)) * V^{-1}
    # using lu! and rdiv! to minimize allocations
    tmp1 = work.row1
    @. tmp1 = exp(t * work.λ)
    Dexp = Diagonal(tmp1)
    mul!(work.VD, work.V, Dexp)
    rdiv!(work.VD, lu!(work.V)) # V * D * V^-1

    # inputs are real, so result should be real up to numerical roundoff
    @. expAt = real(work.VD)
    return expAt
end


"""
    expm_bloch3(r1, r2, w, s, c, t)

Return exp(A*t) for the 3×3 Bloch matrix
`A = [-r2 w s; -w -r2 c; -s -c -r1].`
Uses explicit eigendecomposition
based on analytical roots of cubic characteristic polynomial.
"""
function expm_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc, t::Tt;
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal, Tt <: Breal}
    T = promote_type(T1, T2, Tw, Ts, Tc, Tt)
    expAt = Matrix{T}(undef, 3, 3)
    work = Bloch3ExpmWorkspace{T}()
    expm_bloch3!(expAt, work, r1, r2, w, s, c, t)
    return expAt
end
