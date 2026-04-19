# expm-bloch3.jl

using BlochSim: expm_bloch3
using BlochSim: matrix_bloch3, eigvals_bloch3, eigvec_bloch3, eigen_bloch3
import ForwardDiff
using ExponentialAction: expv
using LinearAlgebra: eigvals, eigvecs, eigen, Diagonal, I
using Test: @inferred, @test, @testset, @test_throws


# functions to compare eigenvalues
function compare_eigs(eig_la::Vector{<:Real}, eig_b3)
    eig_b3 = collect(eig_b3)
    sort(eig_la) ≈ sort(real(eig_b3)) && eig_b3 ≈ real.(eig_b3)
end
compare_eigs(eig_la::Vector{<:Complex}, eig_b3) =
  sort(eig_la, by=imag) ≈ sort(collect(eig_b3), by=imag)


@testset "expm bloch3" begin
    r1 = rand() * 3
    r2 = rand() * 0.1
    ϕ = rand() * 2π
    w = (rand() - 0.5) * 5
    t = rand() * 0.4
    Ω = (rand() - 0.5) * 7 # RF amplitude
    s, c = sincos(ϕ) # for speed and consistency
    s *= Ω
    c *= Ω
    A = @inferred matrix_bloch3(r1, r2, w, s, c)

    # check eigenvalues
    @inferred eigvals_bloch3(Float32.((2, 1, 1, 0, 0))...) # Δ > 0 branch
    @inferred eigvals_bloch3(Float32.((3, 1, 0, 0, 0.8))...) # Δ < 0 branch
    @inferred eigvals_bloch3(r1, r2, w, s, Float16(c)) # helper
    lam3 = @inferred eigvals_bloch3(r1, r2, w, s, c)
    eig = eigen(A)
    @test compare_eigs(eig.values, lam3)

    @inferred eigvec_bloch3(r1, r2, w, s, c, ComplexF16(0)) # helper

    # check eigendecomposition
    @inferred eigen_bloch3(r1, r2, w, s, Float16(c)) # helper
    λ3, V3 = @inferred eigen_bloch3(r1, r2, w, s, c)
    @test compare_eigs(eig.values, λ3)

    @test A ≈ V3 * Diagonal(λ3) * inv(V3)


    # check exp
    E0 = exp(A * t)
    E1 = expm_bloch3(r1, r2, w, s, c, t)
    @test E1 ≈ E0

    Ev = expv(t, A, I(3))
    @test Ev ≈ E0


    # autodiff
    x = [r1, r2, w, s, c]
    f3(x) = expm_bloch3(x..., t)
    @test f3(x) == E1
    jac3 = ForwardDiff.jacobian(f3, x)

    fv(x) = expv(t, matrix_bloch3(x...), I(3))
    jacv = ForwardDiff.jacobian(fv, x)
    @test jac3 ≈ jacv

    ft(t) = expm_bloch3(x..., only(t))
    jact = ForwardDiff.jacobian(ft, [t])
    tmp = reshape(jact, 3, 3)
    @test tmp ≈ A * ft([t]) # check d/dt exp(A*t) = A * exp(A*t)

    f1(t) = exp(A * only(t))
    @test f1(t) == E0
#   jac1 = ForwardDiff.jacobian(f1, [t]) # julia does not support this!
end
