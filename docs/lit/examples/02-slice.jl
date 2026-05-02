#=
# [Slice selective excitation](@id slice select)

This page illustrates slice-selective excitation
using the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
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
        "Random"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run `Pkg.add()` in the preceding code block first, if needed.

using BlochSim: GAMMA, Gradient, GradientSpoiling, Position, Spin, SpinMC
using BlochSim: InstantaneousRF, RF, b1_gauss
using BlochSim: excite!, spoil!
#src using LinearAlgebra: diag
using MIRTjim: prompt
using Plots: annotate!, color, default, gui, plot, plot!, scatter!, scatter3d!, text
using Random: seed!

default(titlefontsize = 10, markerstrokecolor = :auto, label="", width = 1.5,
    linewidth = 2)
seed!(0);


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


#=
### Array of spins
=#

round2(x) = round(x, digits = 2)

Mz0, T1_ms, T2_ms, Δf_Hz = 1, 400, 10, 9 # tissue parameters
xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tuple

zlist = range(-1, 1, 201) # z positions (cm)

spins = map(zlist) do z
    pos = Position(0, 0, z)
    Spin(Mz0, T1_ms, T2_ms, Δf_Hz, pos)
end;


#=
### RF pulse
=#

Δt_ms = 1e-3 # ms, so 1μs dwell
tRF_ms = 0.9
nsamp = round(Int, tRF_ms / Δt_ms)
pisinc(x) = x == 0 ? zero(x) : sin(π * x) / (π * x)

nlobe = 5
sinc5(t) = sinc(2nlobe * t / tRF_ms)

α_deg = 60 # flip angle ° (somewhat large for testing)
α_rad = deg2rad(α_deg)

t = ((0:(nsamp-1))/nsamp .- 0.5) * tRF_ms # [-tRF_ms/2, tRF_ms/2)
waveform = sinc5.(t)
waveform .*= (nsamp / sum(waveform)) * b1_gauss(α_rad, tRF_ms)
prf = plot(t, waveform,
  xaxis = ("t [ms]", (-1,1) .* (tRF_ms/2), ),
  yaxis = ("b₁(t) [Gauss]", ),
  title = "5-lobe sinc for α = $(α_deg)° and tRF = $tRF_ms ms",
)

#
prompt()


#=
### Slice-selective gradient

From Fourier analysis:
`2π/(tRF/2nlobe) = GAMMA * gz * slice_width`
so
`gz = 2nlobe*2π/tRF / GAMMA / slice_width`
=#
slice_width = 0.5 # cm
gz =  (2nlobe*2π) / (tRF_ms/1000) / GAMMA / slice_width
grad = Gradient(0, 0, gz)
rf1 = RF(waveform, Δt_ms, 0, grad)
α_total = sum(rf1.α .* cis.(rf1.θ))
@assert α_total ≈ α_rad
refocus = GradientSpoiling(Gradient(0, 0, -gz), tRF_ms/2);

#=
### Excite the spins with the RF, then refocus
=#
map(spins) do spin
    excite!(spin, rf1)
    spoil!(spin, refocus)
end;


#=
### Plot slice profile
=#
mx = map(spin -> spin.M.x, spins)
my = map(spin -> spin.M.y, spins)
mz = map(spin -> spin.M.z, spins)
mmag = @. sqrt(mx^2 + my^2)
mpha = @. atan(my, mx)

pmag = plot(
  xaxis = ("z [cm]", (-1,1), [-1, -slice_width/2, 0, slice_width/2, 1]),
  yaxis = ("", (-0.2,1), ([0, sin(α_rad)/2, sin(α_rad), 1],
    ["0", "sin($(α_deg)°)/2", "sin($(α_deg)°)", 1]) ),
  legend = :right,
)
plot!(zlist, mx, label = "Mx")
plot!(zlist, my, label = "My")
plot!(zlist, mz, label = "Mz")
plot!(zlist, mmag, label = "|Mxy|")

ppha = plot( xaxis = ("z [cm]", ), legend = :right )
plot!(zlist, mpha, label = "∠Mxy")

plot(pmag, ppha; layout = (2,1),
  plot_title = "Slice profile for 5-lobe sinc, α = $(α_deg)°",
)

#
prompt()

#=
The simulation shown above
is similar to Fig. 5 of
[Pauly et al. 1989](https://doi.org/10.1016/0022-2364(89)90265-5),
where it is noted that improved refocusing
could lead to better refocusing.
This simulation ignores gradient slew-rate constraints
that must be considered in practice.
=#

#src todo: effect on qMRI


#=
ntry = 10
rand20(x::Number) = x * (1 + 0.2 * (rand() - 0.5) / 0.5) # ± 20% variability
function do_fit(signal_ri::Function, i::Int)
    seed!(i)
    y = yb + 1σ * randn(ComplexF64, size(yb))
    cost(x) = abs2(norm(signal_ri(x) - real_imag(vec(y)))) # LS cost
    return optimize_multistart(cost, todo)
end
=#
