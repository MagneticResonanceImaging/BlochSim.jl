# test/slice.jl

using BlochSim: rf_slice, duration, RF
using Test: @inferred, @test, @testset

@testset "slice" begin
    rf = @inferred rf_slice()
    @test π/2 ≈ sum(rf.α .* cis.(rf.θ))             
    @test rf isa RF
    @test duration(rf) == 1.5 # because of refocusing gradient

    α_rad = π/3
    rf = rf_slice(; α_rad, refocus = :false) # no @inferred
    @test α_rad ≈ sum(rf.α .* cis.(rf.θ))
    @test rf isa RF
    @test duration(rf) == 1

    rf = rf_slice(; α_rad, slice_width = Inf) # no @inferred
    @test α_rad ≈ sum(rf.α .* cis.(rf.θ))
    @test rf isa RF
    @test duration(rf) == 1

    #=
    Could test here whether the effect of the RF pulse
    with refocusing gradient has the same effect has
    the cascade of an RF pulse followed by
    refocus = GradientSpoiling(Gradient(0, 0, -gz), duration/2)
    using `excite` and `spoil` but we leave it to the docs instead.
    =#
end
