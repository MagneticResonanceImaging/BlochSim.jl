# test/rf-rect.jl

using BlochSim: InstantaneousRF, RectRF, RF1
using BlochSim: Gradient, Magnetization, Position, Spin, excite, excite!
using Test: @inferred, @test, @testset

@testset "RectRF" begin

    rf = @inferred RectRF(1f0)
    rf = @inferred RectRF(1f0, 1)
    rf = @inferred RectRF(1f0, 1, π)

    α = π/3
    θ = π/5

    grad = Gradient(0.1, 0.2, 0.3) # G/cm (to stress test)
    rf = @inferred RectRF(2f0, α, θ, grad)
    @test rf isa RectRF
    @test rf.duration isa Real
    @test rf.α isa Real
    @test rf.θ isa Real

    show(devnull, rf)
    show(devnull, "text/plain", rf)
    @test 1 == length(rf)

    # test with very short pulses
    rf0 = InstantaneousRF(α, θ)
    rf1 = RF1(α, 1e-8, θ, grad)
    rf3 = RectRF(1e-8, α, θ, grad)
    pos = Position(3, 2, 1) # cm
    spin = Spin(1, 1000, 100, 7., pos)
    (A0, b0) = @inferred excite(spin, rf0)
    (A1, b1) = @inferred excite(spin, rf1)
    (A3, b3) = @inferred excite(spin, rf3)
    @test maximum(abs, Matrix(A0) - Matrix(A1)) ≤ 1e-6
    @test maximum(abs, Matrix(A1) - Matrix(A3)) ≤ 1e-7
    @test maximum(abs, #= Vector(b0) - =# Vector(b1)) ≤ 1e-11


    # test with longer pulse
    tRF_ms = 2
    nsamp = 2000 # 1μs spacing
    rf1 = RF1(α, tRF_ms, θ, grad; nsamp)
    @test sum(rf1.α) ≈ α
    rf3 = RectRF(tRF_ms, α, θ, grad)
    pos = Position(3, 2, 1) # cm
    spin = Spin(1, 1000, 100, 7., pos)
    (A1, b1) = @inferred excite(spin, rf1)
    (A3, b3) = @inferred excite(spin, rf3)

    @test maximum(abs, Matrix(A1) - Matrix(A3)) ≤ 1e-6
    @test maximum(abs, Vector(b1) - Vector(b3)) ≤ 1e-9

    spin1 = deepcopy(spin)
    spin3 = deepcopy(spin)
    excite!(spin1, rf1)
    excite!(spin3, rf3)

    @test maximum(abs, Vector(spin1.M - spin3.M)) ≤ 1e-6
end
