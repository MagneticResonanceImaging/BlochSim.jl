#=
# [Matrix exponential](@id expm)

This page illustrates using the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
to compute matrix exponentials,
i.e., `expm`,
that are needed for Bloch simulations.
=#

#srcURL


#=
### Setup

First we add the Julia packages that are need for this demo.
Change `false` to `true` in the following code block
if you are using any of the following packages for the first time.
=#

if false
    import Pkg
    Pkg.add([
        "BenchmarkTools"
        "BlochSim"
        "ForwardDiff"
        "InteractiveUtils"
        "LaTeXStrings"
        "ExponentialAction"
        "LinearAlgebra"
        "MIRTjim"
        "Plots"
        "Random"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run `Pkg.add()` in the preceding code block first, if needed.

using BenchmarkTools: @benchmark
using BlochSim: expm_bloch3, matrix_bloch3
using ExponentialAction: expv
using LinearAlgebra: I
using MIRTjim: prompt
using Plots: gui, plot, plot!, default
using Random: seed!

seed!(0)
default(titlefontsize = 10, markerstrokecolor = :auto, label="", width = 1.5)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


#=
## Benchmark `exp*`

- `ExponentialAction.expv`:
  general-purpose matrix exponential

- `BlochSim.expm_bloch3`
  specific matrix exponential for 3×3 Bloch matrix
=#

r1 = rand() * 3
r2 = rand() * 0.1
ϕ = rand() * 2π
w = (rand() - 0.5) * 5
t = rand() * 0.4
Ω = (rand() - 0.5) * 7 # RF amplitude
s, c = sincos(ϕ) # for speed and consistency
s *= Ω
c *= Ω
A = matrix_bloch3(r1, r2, w, s, c)


x = [r1, r2, w, s, c]
f3(x) = expm_bloch3(x..., t)
fv(x) = expv(t, matrix_bloch3(x...), I(3))

b3 = @benchmark f3($x) # 1.3 μs (22 allocations: 1.20 KiB)

bv = @benchmark fv($x) # 3.8 μs (147 allocations: 10.4 KiB)


include("../../../inc/reproduce.jl")
