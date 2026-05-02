# test/fit.jl

import Optim # extension
using BlochSim: fit_signal, optimize_multistart
using Test: @inferred, @test, @testset

@testset "optimize_multistart" begin
    costs = fill(Inf, 5)
    xmin = optimize_multistart(
        x -> sum(abs2, x .- 2), i -> [1.0i]; ntry = 3, costs)
    @test xmin == [2]
    @test isinf(costs[4])

    model(x) = x .- 2
    y = [3]
    xmin = fit_signal(model, [0.], y)
    maximum(abs, xmin - [5]) ≤ 1e-4
end
