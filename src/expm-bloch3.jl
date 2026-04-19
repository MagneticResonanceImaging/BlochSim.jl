#=
expm-bloch3.jl
Compute matrix exponential exp(A * t)
where A is the 3×3 Bloch matrix.

This code was developed with the help of GPT 5.2.
=#

using LinearAlgebra: cond, cross, norm, Diagonal

export expm_bloch3

const Breal = Number # Real or ForwardDiff.Dual


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


"""
    eigvec_bloch3(r1, r2, w, s, c, λ)
Compute eigenvector for eigenvalue λ by crossing two rows;
switch row pairs if needed.
"""
function eigvec_bloch3(
    r1::T, r2::T, w::T, s::T, c::T, λ::Complex{T};
    rtol::Real = 1e-12,
) where {T <: Breal}

    R1 = [-r2-λ, w,     s]
    R2 = [-w,   -r2-λ,  c]
    R3 = [-s,   -c,   -r1-λ]
    rmax = maximum(norm, (R1, R2, R3))

    v = cross(R1, R2)
    if norm(v) < rtol * rmax
        v = cross(R1, R3)
    end
    if norm(v) < rtol * rmax
        v = cross(R2, R3)
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


"""
    (λ, V) = eigen_bloch3(r1, r2, w, s, c)
Return eigendecomposition `(λ, V)`
of 3×3 Bloch matrix
`A = [-r2  w   s;
      -w  -r2  c;
      -s  -c  -r1]`
"""
function eigen_bloch3(r1::T, r2::T, w::T, s::T, c::T) where {T <: Breal}

    λ1, λ2, λ3 = eigvals_bloch3(r1, r2, w, s, c)

    V = hcat(
        eigvec_bloch3(r1, r2, w, s, c, λ1),
        eigvec_bloch3(r1, r2, w, s, c, λ2),
        eigvec_bloch3(r1, r2, w, s, c, λ3)
    )

    return [λ1, λ2, λ3], V
end


# helper for mismatched types
function eigen_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc,
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal}
    r1, r2, w, s, c, _ = promote(r1, r2, w, s, c, zero(Float32))
    return eigen_bloch3(r1, r2, w, s, c)
end


"""
    expm_bloch3(r1, r2, w, s, c, t)

Return exp(A*t) for the 3×3 Bloch matrix
`A = [-r2 w s; -w -r2 c; -s -c -r1].`
Uses explicit eigendecomposition
based on analytical roots of cubic characteristic polynomial.
"""
function expm_bloch3(
    r1::T1, r2::T2, w::Tw, s::Ts, c::Tc, t::Tt,
) where {T1 <: Breal, T2 <: Breal, Tw <: Breal, Ts <: Breal, Tc <: Breal, Tt <: Breal}

    (λ, V) = eigen_bloch3(r1, r2, w, s, c)

#=
    # If V is ill-conditioned (near defective/repeated eigenvalues), fall back
    # Cannot compute `cond` for dual types
    if cond(V) > 1e12 || any(isnan, V) || any(isinf, V)

        T = promote_type(T1, T2, Tw, Ts, Tc, Tt)
        A = [-r2  w  s;
             -w -r2  c;
             -s -c  -r1]
        return exp(t*A)  # Julia's matrix exponential
    end
=#

    Dexp = Diagonal(exp.(t * λ))
    tmp = V * Dexp * inv(V)
    # Inputs are real; result should be real up to numerical roundoff
    return real(tmp)
end
