# crb.jl

using BlochSim: crb, snr2sigma, real_imag
using LinearAlgebra: I
using Random: seed!
using Test: @inferred, @test, @testset

@testset "crb" begin
    seed!(0)
    A = randn(7, 4)
    signal(x) = A * x
    x = 1:4
    σ = 3
#   bound = @inferred crb(signal, x, σ) # impossible it seems
    bound1 = crb(signal, x, σ)
    fisher1 = A'A / σ^2
    @test bound1 ≈ inv(fisher1)

    # complex signal, real parameter vector
    A = randn(ComplexF32, 7, 4)
    signal(x) = real_imag(A * x) # stack [real; imag]
    bound2 = crb(signal, x, σ)
    As = [real(A); imag(A)]
    fisher2 = As'As / σ^2
    @test bound2 ≈ inv(fisher2)

    y = A * x
    σ = @inferred snr2sigma(30, y)
end
