#=
fit.jl

Functions for fitting signal models to data.
Strictly speaking these have little to do
with Bloch simulations,
but one purpose of such simulations
is to develop signal models
for quantitative MRI,
and fit those models to data.
BlochSim already depends on ForwardDiff,
so we include these helper functions.
=#

import BlochSim: fit_signal, optimize_multistart # extended here

using ADTypes: AutoForwardDiff
using Optim: optimize, LBFGS


"""
    optimize_multistart(cost, x0fun; ntry = 10, costs)

Optimize the `cost` function
`ntry` times,
with the `i`th try
initialized by `x0fun(i)`.
Return the minimizer with the lowest `cost`.

Option:
- `costs[i]` is mutated to contain the `cost` of the `i`th trial.
- `autodiff` to pass to `Optim.optimize` default `AutoForwardDiff()`
- `optim_method` optimization method used by `Optim.optimize`; default `LBFGS()`

An example initialization method is
`x0fun = i -> xtrue * (1 + 0.2 * (rand() - 0.5) / 0.5)` for # ± 20% variability

```jldoctest
julia> optimize_multistart(x -> sum(abs2, x .- 2), i -> [1.0i]; ntry = 3)
1-element Vector{Float64}:
 2.0
```
"""
function optimize_multistart(
    cost,
    x0fun;
    ntry::Integer = 10,
    costs::AbstractVector = Vector{Float64}(undef, ntry),
    autodiff = AutoForwardDiff(),
    optim_method = LBFGS(),
)

    min_cost = Inf
    xbest = nothing

    for i in 1:ntry
        x0 = x0fun(i)
        opt = optimize(cost, x0, optim_method; autodiff)
        if opt.minimum < min_cost
            min_cost = opt.minimum
            xbest = opt.minimizer
        end
        costs[i] = opt.minimum
    end

    xbest == nothing && throw("No finite cost found")

    return xbest
end


"""
    xmin = fit_signal(model::Function, x0::Array, y::Array; kwargs...)
    xmin = fit_signal(model::Function, x0fun::Function, y::Array; ntry=10, kwargs...)

Perform nonlinear LS fitting
`xmin = argmin_x ‖ model(x) - y ‖²`.

If input initial parameter vector guess `x0` is an `Array`,
then perform a single optimization.
Otherwise, use `optimize_multistart` to return the best of `ntry` runs.
Extra `kwargs` are passed to `optimize_multistart`.
"""
function fit_signal(
    model,
    x0fun,
    y::AbstractArray;
    ntry::Integer = 10,
    kwargs...,
)
    cost(x) = sum(abs2, model(x) - y) # LS cost
    return optimize_multistart(cost, x0fun; ntry, kwargs...)
end

function fit_signal(model, x0::AbstractArray, y::AbstractArray; ntry::Integer = 1, kwargs...)
    ntry > 1 && @warn("ntry > 1 pointless when x0 provided")
    return fit_signal(model, i -> x0, y; ntry, kwargs...)
end
