# test/reuse.jl

using BlochSim: Spin, reuse_spin, reuse_spins!
using Test: @inferred, @test, @testset

@testset "reuse" begin
    spin1 = Spin(2, 3., 4, 5f0)
    spin2 = @inferred reuse_spin(spin1, 7, 8, 9, 1)
    @test spin1.M === spin2.M
    @test spin1.pos === spin2.pos
    @test spin2.M0 == 7
    @test 0 == @allocated reuse_spin(spin1, 7, 8, 9, 1)

    spin3 = Spin(5, 4., 3, 1//2)
    spins = [spin1, spin3]
    @inferred reuse_spins!(spins, 7, 8, 9, 1)
    @test 0 == @allocated reuse_spins!(spins, 7, 8, 9, 1)
    @test spins[end].M0 == 7

#=
    using BenchmarkTools: @benchmark
    b1 = @benchmark reuse_spin($spin1, 7, 8, 9, 1)
    b2 = @benchmark reuse_spins!($spins, 7, 8, 9, 1)
    @test 0 == b1.allocs
    @test 0 == b2.allocs
=#
end
