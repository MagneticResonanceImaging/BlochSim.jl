# test/rf-rect.jl

using BlochSim: InstantaneousRF, RectRF
using BlochSim: Gradient, Magnetization, Spin, excite, excite!
using Test: @inferred, @test, @testset

@testset "RectRF" begin

    rf = @inferred RectRF(1f0)
    rf = @inferred RectRF(1f0, 1)
    rf = @inferred RectRF(1f0, 1, π)

    α = π/3
    θ = π/5
    rf = @inferred RectRF(2f0, α, θ, Gradient(0,0,0))
    @test rf isa RectRF
    @test rf.duration isa Real
    @test rf.α isa Real
    @test rf.θ isa Real

    show(devnull, rf)
    show(devnull, "text/plain", rf)
    @test 1 == length(rf)

    rf0 = InstantaneousRF(α, θ)
    rf1 = RectRF(1e-6, α, θ)
    spin = Spin(1, 1000, 100, 7.)
    (A0, b0) = @inferred excite(spin, rf0)
    (A1, b1) = @inferred excite(spin, rf1)
    @test maximum(abs, Matrix(A0) - Matrix(A1)) ≤ 1e-7
    @test maximum(abs, #= Vector(b0) - =# Vector(b1)) ≤ 1e-9

    spin0 = deepcopy(spin)
    spin1 = deepcopy(spin)
    excite!(spin0, rf0)
#   excite!(spin1, rf1) todo
end


#=
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct


    A2 = BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, rf)

    show(devnull, rf)
    show(devnull, "text/plain", rf)
    @test A1 == A2
    return B1 == B2


    (A1, B1) = @inferred excite(spin, RF(rf, dt, Δθ, grad))
    A2 = BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, RF(rf, dt, Δθ, grad))
    @test A1 == A2
    return B1 == B2

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    grad = [Gradient(0, 0, 0) for i = 1:2]
    rf = RF(fill(exp(im * π/8), length(grad)), dt, Δθ, grad)
    excite!(s, rf)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    rf = RF(fill(exp(im * π/8), 2), dt, Δθ)
    excite!(s, rf)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)

    @test duration(rf) == 2dt
    return s.M ≈ M_correct
=#
