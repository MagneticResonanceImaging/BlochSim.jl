# test/slice.jl

using BlochSim: rf_slice, duration, RF, GradientSpoiling
using Test: @inferred, @test, @testset

@testset "slice" begin
    flip(rf) = sum(rf.α .* cis.(rf.θ))

    tRF_ms = 1.01 # stress test
    rf, rephasing = @inferred rf_slice(tRF_ms)
    @test π/2 ≈ flip(rf)
    @test rf isa RF
    @test duration(rf) == tRF_ms
    @test rephasing isa GradientSpoiling
    @test rephasing.Tg == tRF_ms/2

    rf = @inferred rf_slice(Val(:built_in_rephasing), tRF_ms)
    @test π/2 ≈ flip(rf)
    @test rf isa RF
    @test duration(rf) ≈ tRF_ms * 1.5 # because of rephasing gradient

    α_rad = π/3
    rf, rephasing = @inferred rf_slice(tRF_ms; α_rad)
    @test α_rad ≈ flip(rf)
    @test rf isa RF
    @test duration(rf) == tRF_ms

    rf, rephasing = @inferred rf_slice(tRF_ms; α_rad, slice_width = Inf)
    @test α_rad ≈ flip(rf)
    @test rf isa RF
    @test duration(rf) == tRF_ms

    #=
    Could test here whether the effect of the RF pulse
    with built-in rephasing gradient has the same effect has
    the cascade of an RF pulse followed by
    rephasing = GradientSpoiling(Gradient(0, 0, -gz), duration/2)
    using `excite` and `spoil` but we leave it to the docs instead.
    =#
end
