#=
# [Rotation](@id rotation)

This page discusses the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
especially its rotation conventions.
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
        "BlochSim"
        "InteractiveUtils"
        "LaTeXStrings"
        "LinearAlgebra"
        "MIRTjim"
        "Plots"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run `Pkg.add()` in the preceding code block first, if needed.

using BlochSim: Spin, SpinMC, InstantaneousRF, RF, excite
using BlochSim: excite_bloch3
#src using LinearAlgebra: diag
using MIRTjim: prompt
using Plots: annotate!, color, default, gui, plot, plot!, scatter!, scatter3d!, text

default(titlefontsize = 10, markerstrokecolor = :auto, label="", width = 1.5,
    linewidth=2)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


#=
## Rotation convention

The RF excitation convention used in `BlochSim`
is illustrated by the following plots.
=#

round2(x) = round(x, digits = 2)

Mz0, T1_ms, T2_ms, Δf_Hz = 1, 400, 10, 9 # tissue parameters
xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tuple
α_deg = 90 # flip angle °
α_rad = deg2rad(α_deg)
spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz)


# helper
function plot_tip(θ; α=π/2, color=:blue)
    v = nothing
    tmp = [1, 1, 1]
    p = plot(
     xaxis = ("", (-1, 1), tmp, ),
     yaxis = ("", (-1, 1), tmp, ),
     zaxis = ("", (-1, 1), tmp, ),
     aspect_ratio = 1,
     framestyle = :origin,
     size = (400, 400),
     camera = (139, 30),
     # thickness_scaling = 2,
     title = "θ = $(round2(rad2deg(θ)))°",
    )
    plot!([0, -sin(θ)], [0, -cos(θ)], [0, 0], lw=3, color=:black) # B1

    for α in range(0, α, 11)
        rf = InstantaneousRF(α, θ)
        A, _ = excite(spin, rf)
        v = Matrix(A) * [0, 0, 1]
        plot!([0, v[1]], [0, v[2]], [0, v[3]], lw=2; color)
        scatter3d!([v[1]], [v[2]], [v[3]],
            markercolor = color, 
            markershape = :utriangle, # lazy arrow tip
        )
    end

    annotate!([
        (1.3, 0, 0, text("x", :left, 10, :bold)),
        (0, 1.2, 0, text("y", :left, 10, :bold)),
        (0, 0, 1.2, text("z", :bottom, 10, :bold)),
    ])
    return p
end;

ptip = plot(
 plot_tip(  0, color=:blue),
 plot_tip(π/4, color=:magenta),
 plot_tip(π/2, color=:red),
 layout = (1,3),
 size = (900, 300),
)

#
prompt()

#=

As described in the docstring for
`rotatetheta!(A, α, θ)`,
excitation here applies
"left-handed rotation by angle `α`
about an axis in the x-y plane
that makes left-handed angle `θ`
with the negative y-axis."

Comments within the function
refer to p.27 of Dwight Nishimura's
["Principles of Magnetic Resonance Imaging" (1996).](https://books.google.com/books/about/Principles_of_Magnetic_Resonance_Imaging.html?id=uz9BAQAAIAAJ)

As illustrated by black line above,
the B1 axis is
`[-sin(θ) -cos(θ) 0]`

So `α = π/2` and `θ = 0`
takes equilibrium magnetization
`[0, 0, 1]`
to
`[1, 0, 0]`.

The following excitation matrix examples illustrate.
=#

rf0 = InstantaneousRF(α_rad, 0)
A1, _ = excite(spin, rf0)
round2.(Matrix(A1))

#
round2.(Matrix(A1) * [0, 0, 1])


# ## θ=π/2

#
rf2 = InstantaneousRF(α_rad, π/2)
A2, _ = excite(spin, rf2)
round2.(Matrix(A2))

#
round2.(Matrix(A2) * [0, 0, 1])

#src include("../../../inc/reproduce.jl")
