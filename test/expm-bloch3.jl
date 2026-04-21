# expm-bloch3.jl

using BlochSim: expm_bloch3, expm_bloch3!
using BlochSim: matrix_bloch3, eigvals_bloch3, eigvec_3x3!
using BlochSim: Bloch3ExpmWorkspace, cross!, eigvec_bloch3!, eigen_bloch3!
using ExponentialAction: expv
import ForwardDiff
using LinearAlgebra: cross, eigvals, eigvecs, eigen, Diagonal, I
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
    xp = [r1, r2, w, s, c] # random parameters for testing

    # helper
    A = @inferred matrix_bloch3(xp...)


    # check eigenvalues
    x1 = (2, 1, 1, 0, 0) # Δ > 0
    x2 = (3, 1, 0, 0, 0.8) # Δ < 0
    @test 0 == @allocated eigvals_bloch3(x1...) # Δ > 0
    @test 0 == @allocated eigvals_bloch3(x2...) # Δ < 0
    @inferred eigvals_bloch3(Float32.(x1)...) # Δ > 0 branch
    @inferred eigvals_bloch3(Float32.(x2)...) # Δ < 0 branch
    @inferred eigvals_bloch3(2, 1, 0, 0, 0f0) # helper
    lam3 = @inferred eigvals_bloch3(xp...)
    eig = eigen(A)
    @test compare_eigs(eig.values, lam3)


    # check cross!
    v = Vector{ComplexF64}(undef, 3)
    a = [0, 1//2, 4f0]
    b = [1.0, 2, complex(3,4)]
    @inferred cross!(v, a, b) # warm-up too
    @test 0 == @allocated cross!(v, a, b)
    @test v ≈ cross(a, b)


    # check eigvecs
    row1 = Vector{ComplexF64}(undef, 3)
    row2 = Vector{ComplexF64}(undef, 3)
    row3 = Vector{ComplexF64}(undef, 3)
    v = Vector{ComplexF64}(undef, 3)

    @inferred eigvec_3x3!(v, row1, row2, row3)
    @test [0, 0, 1] == eigvec_3x3!(v, # code cover case 2
        ComplexF64[1, 0, 0], ComplexF64[1, 0, 0], ComplexF64[0, 1, 0])
    @test [0, 0, 0] == eigvec_3x3!(v, # code cover case 3
        ComplexF64[0, 0, 0], ComplexF64[0, 0, 0], ComplexF64[1, 0, 0])
    @test 0 == @allocated eigvec_3x3!(v, row1, row2, row3)

    λ0 = Complex(0f0)
    @inferred eigvec_bloch3!(v, row1, row2, row3, 2, 1//1, 0, 0, 0, λ0)
    eigvec_bloch3!(v, row1, row2, row3, x1..., λ0) # warm-up
    @test 0 == @allocated eigvec_bloch3!(v, row1, row2, row3, 2, 1, 1, 0, 0, λ0)
    @inferred eigvec_bloch3!(v, row1, row2, row3, r1, r2, w, s, c, λ0) # warm up
    @test 0 == @allocated eigvec_bloch3!(v, row1, row2, row3, r1, r2, w, s, c, λ0)


    # check eigendecomposition
    T = Float64
    C = Complex{T}
    work = Bloch3ExpmWorkspace{T}()
    arg = (work.λ, work.V, work.row1, work.row2, work.row3)
    xt = (2., 1., 0., 0., 0.)
    @inferred eigen_bloch3!(arg..., xt...) # warm up
    @test 0 == @allocated eigen_bloch3!(arg..., xt...)

    λ3 = Vector{C}(undef, 3)
    V3 = Matrix{C}(undef, 3, 3)
    @inferred eigen_bloch3!(λ3, V3, work.row1, work.row2, work.row3, xp...)
    @test compare_eigs(eig.values, λ3)
    @test A ≈ V3 * Diagonal(λ3) * inv(V3)


    # check exp
    expAt = Matrix{T}(undef, 3, 3)
    @inferred expm_bloch3!(expAt, work, xp..., t) # warm up
    @test 80 ≥ @allocated expm_bloch3!(expAt, work, xp..., t)

    E0 = exp(A * t)
    @test expAt ≈ E0

    E1 = expm_bloch3(xp..., t)
    @test E1 ≈ E0

    Ev = expv(t, A, I(3))
    @test Ev ≈ E0


    # 2 repeated roots
    xr2 = (2, 1, 0, 0, 0)
    @test exp(matrix_bloch3(xr2...)) ≈ expm_bloch3(xr2..., 1)
    # 3 repeated roots
    xr3 = (3, 3, 0, 0, 0)
    @test exp(matrix_bloch3(xr3...)) ≈ expm_bloch3(xr3..., 1)


    # analytical solution: relaxation
    @test expm_bloch3(r1, r2, 0, 0, 0, t) ≈
        Diagonal([exp(-r2*t), exp(-r2*t), exp(-r1*t)])

    # analytical solution: precession
    @test expm_bloch3(r1, r2, (π/2)/t, 0, 0, t) ≈
        [0 exp(-r2*t) 0; -exp(-r2*t) 0 0; 0 0 exp(-r1*t)]

    # analytical solution: tip
    expm_bloch3(0, 0, 0, π/2, 0, 1) ≈ [0 0 1; 0 1 0; -1 0 0]
    expm_bloch3(0, 0, 0, 0, π/2, 1) ≈ [1 0 0; 0 0 1; 0 -1 0]


    # autodiff
    f3(x) = expm_bloch3(x..., t)
    @test f3(xp) == E1
    jac3 = ForwardDiff.jacobian(f3, xp)

    fv(x) = expv(t, matrix_bloch3(x...), I(3))
    jacv = ForwardDiff.jacobian(fv, xp)
    @test jac3 ≈ jacv

    ft(t) = expm_bloch3(xp..., only(t))
    jact = ForwardDiff.jacobian(ft, [t])
    tmp = reshape(jact, 3, 3)
    @test tmp ≈ A * ft([t]) # check d/dt exp(A*t) = A * exp(A*t)

    f1(t) = exp(A * only(t))
    @test f1(t) == E0
#   jac1 = ForwardDiff.jacobian(f1, [t]) # julia does not support this!
end
