#=
# [bSSFP RF Duration](@id bssfp-rf-dur)

This page illustrates using the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
to calculate MRI signals
for
balanced steady-state free precession
[(bSSFP)](https://en.wikipedia.org/wiki/Steady-state_free_precession_imaging)
pulse sequences.

Specifically it examines the effects of finite RF duration.


### References

- todo

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
        "ForwardDiff"
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
#src using BlochSim: bSSFPellipse
using BlochSim: bssfp, GAMMA
#src import ForwardDiff
#src using LaTeXStrings: latexstring
#src using LinearAlgebra: Diagonal, I, cond, diag, norm
using MIRTjim: prompt
using Plots: gui, plot, plot!, default
default(titlefontsize = 10, markerstrokecolor = :auto, label="", width = 1.5)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


#=
## RF pulse duration effects for 1-pool model

Examine effects of finite RF pulse duration
for a single spin with a relatively short T2.
=#

Mz0, T1_ms, T2_ms, Δf_Hz = 1, 400, 40, 0 # tissue parameters
α_deg = 50 # flip angle
TR_ms, TE_ms, α_rad = 8, 4, deg2rad(α_deg) # scan parameters
spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz)

Δϕ_rad = range(-1, 1, 101) * π # phase-cycling factors
_bssfp(Δϕ, rf) = bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ, rf)
_bssfp(rf) = map(Δϕ -> _bssfp(Δϕ, rf), Δϕ_rad) # helper

rf0 = InstantaneousRF(α_rad)
signal0 = _bssfp(rf0) # InstantaneousRF signal


#=
Specify finite-duration (rectangular) RF pulse
- `GAMMA` has units rad/s/G
- Tip angle for constant pulse:
  `α_rad = GAMMA * b1_gauss * tRF_s`
- so `b1_gauss = α_rad / GAMMA / tRF_s`
=#
b1_gauss(α_rad, tRF_ms) = α_rad / GAMMA / (tRF_ms / 1000)


#=
### Test "nearly instantaneous" RF pulse
=#
tRF_ms = 1e-12 # super-short for first test
waveform1 = [1] * b1_gauss(α_rad, tRF_ms) # single sample i.e. instant!
rf1 = RF(waveform1, tRF_ms)
#=
signal1 = _bssfp(rf1)
@assert signal0 ≈ signal1 # should be essentially identical
@assert α_rad == rf0.α ≈ only(rf1.α)
=#


#=
### Test 2ms RF pulse
Somewhat unexpectedly (to JF),
the signal matches the instantaneous RF case,
because `excite!` for a `RF` type uses `freeprecess!`
for each sample and here there is just a single sample.
=#
tRF_ms = 2
waveform2 = [1] * b1_gauss(α_rad, tRF_ms) # single sample i.e. instant!
rf2 = RF(waveform2, tRF_ms)
#=
signal2 = _bssfp(rf2)
@assert signal0 ≈ signal2 # matches!?
@assert α_rad == rf0.α ≈ only(rf2.α)
=#


#=
Examine excitation matrices
=#
A0, B0 = excite(spin, rf0)
A1, B1 = excite(spin, rf1)
@assert B0 === nothing
@assert maximum(abs, Vector(B1)) ≤ 3e-15 # 15*eps()
@assert Matrix(A0) ≈ Matrix(A1)

A2, B2 = excite(spin, rf2)
@assert !(Matrix(A1) ≈ Matrix(A2)) # huh!?
Matrix(A1) - Matrix(A2)

@which excite(spin, rf2)
throw()

#=
A1 and A2 differ, so why are is signal0 ≈ signal2
=#


# Plot
xaxis = ("phase cycling increment Δϕ (rad)", (-π, π), ((-1:1).*π, ["-π", "0", "π"]))
prfm = plot( ; xaxis, ylabel = "bSSFP signal mag",)
prfa = plot( ; xaxis, ylabel = "bSSFP signal phase",)

plot!(prfm, Δϕ_rad, abs.(signal0), label="Instantaneous")
plot!(prfa, Δϕ_rad, angle.(signal0), label="Instantaneous")

#src plot!(Δϕ_rad, abs.(signal1), label="tRF = $tRF_ms")
nw = 1000 # approximately 1μs dwell time
for tRF_ms in [1e-2 1 2]
    waveform4 = ones(nw) * b1_gauss(α_rad, tRF_ms)
    rf4 = RF(waveform4, tRF_ms / nw)
    @assert rf0.α ≈ sum(rf4.α)
    signal4 = _bssfp(rf4)
    label = "tRF = $tRF_ms ms, nw=$nw"
    plot!(prfm, Δϕ_rad, abs.(signal4); label)
    plot!(prfa, Δϕ_rad, angle.(signal4); label)
end;

prf = plot(prfm, prfa, layout=(2,1),
 plot_title = "RF pulse duration effect, T2=$T2_ms (ms)")


# todo include("../../../inc/reproduce.jl")
