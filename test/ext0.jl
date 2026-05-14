# ext0.jl test extension stubs *before* loading the extension

using BlochSim: fit_signal, optimize_multistart
using Test: @test_throws, @testset

@testset "ext0" begin
    @test_throws String fit_signal()
    @test_throws String optimize_multistart()
end
