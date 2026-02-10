#=
# [bSSFP](@id 01-bssfp)

This page illustrates using the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
to calculate MRI signals
for
balanced steady-state free precession
[(bSSFP)](https://en.wikipedia.org/wiki/Steady-state_free_precession_imaging)
pulse sequences.

This demo facilitates
understanding bSSFP sequences,
multi-compartment spins,
and myelin water exchange.

This demo recreates Figure 3 from [1] and Figure 2 from [2].


### References

- [1] Hargreaves, B., Vasanawala, S., Pauly, J., & Nishimura, D. (2001).
  Characterization and reduction of the transient response
  in steady‐state MR imaging.
  [MRM 46(1), 149-158](https://doi.org/10.1002/mrm.1170).

- [2] Murthy, N., Nielsen, J., Whitaker, S., Haskell, M., Swanson, S.,
  Seiberlich, N., & Fessler, J. (2022).
  Quantifying myelin water exchange using optimized bSSFP sequences.
  [Proc. Intl. Soc. Mag. Res. Med (#2068)](https://cds.ismrm.org/protected/22MProceedings/PDFfiles/2068.html).

- [3] Hinshaw, W. S. (1976).
  Image formation by nuclear magnetic resonance: the sensitive‐point method.
  [J. of Applied Physics, 47(8), 3709-21](https://doi.org/10.1063/1.323136).

- [4] Whitaker, S. T., Nataraj, G., Nielsen, J. F., & Fessler, J. A. (2020).
  Myelin water fraction estimation using small‐tip fast recovery MRI.
  [MRM 84(4), 1977-90](https://doi.org/10.1002/mrm.28259).
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
        "LaTeXStrings"
        "LinearAlgebra"
        "MIRTjim"
        "Plots"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run `Pkg.add()` in the preceding code block first, if needed.

using BlochSim: Spin, SpinMC, InstantaneousRF, RF, excite, freeprecess
using BlochSim: bssfp, bssfp_ellipse, GAMMA
import ForwardDiff
using InteractiveUtils: versioninfo
using LaTeXStrings: latexstring
using LinearAlgebra: Diagonal, I, cond, diag, norm
using MIRTjim: prompt
using Plots: gui, plot, plot!, default
default(titlefontsize = 10, markerstrokecolor = :auto, label="", width = 1.5)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


# Define some useful helper functions.

Hz_to_kHz(Δf_Hz) = Δf_Hz * 10^(-3) # convert frequencies in Hz to kHz
kHz_to_Hz(Δf_kHz) = Δf_kHz * 10^(3); # convert frequencies in kHz to Hz


#=
The bSSFP pulse sequence in Figure 2 in [1] starts with a RF pulse,
then
- `a` is at time TE
- `b` is TR-TE later, right before next RF pulse
- `c` is immediately after the next RF pulse
- `d` is TE after that next RF pulse
We use this to generate Figure 3 in [1] in two different ways.
The RF excitation is repeated periodically and, in steady-state,
the magnetization at point *a* is the same as at point *d*.
=#


#=
## Method 1: Use matrices

Use Equations 1 and 2 and Appendix A from [1]

Calculate the steady-state value at point *d*
using the method from [1] using Equations 1 and 2 and Appendix A.
=#

"""
    bssfp_matrix(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, α_rad, θ_rf_rad=0)
)

Return steady-state magnetization signal value
at the echo time
for a bSSFP sequence
using method of
[Hargreaves et al., MRM 2001](https://doi.org/10.1002/mrm.1170).

## In tissue:
- `Mz0` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `Δf_Hz` off-resonance value (Hz) (default 0)

## In scan:
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- `α_rad` flip angle of RF pulse (radians)
- `θ_rf_rad` RF pulse phase (radians) (default 0)

## Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_matrix(
    Mz0, T1_ms, T2_ms, Δf_Hz,
    TR_ms, TE_ms, α_rad, θ_rf_rad::Number = 0,
)

    M0 = [0; 0; Mz0] # initial magnetization vector

    ## rotation matrix for RF excitation about the y-axis
    if θ_rf_rad == 0 # y-axis
        R = [cos(α_rad) 0 sin(α_rad); 0 1 0; -sin(α_rad) 0 cos(α_rad)]
    elseif θ_rf_rad == -π/2 # x-axis
        R = [1 0 0; 0 cos(α_rad) sin(α_rad); 0 -sin(α_rad) cos(α_rad)]
    else
        throw("θ_rf_rad = $θ_rf_rad not implemented")
    end

    ## free precession matrix
    Pz(angle) = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1]
    P(τ_ms) = Pz( 2π * Hz_to_kHz(Δf_Hz) * τ_ms ) # angle in radians

    ## matrices for T1 and T2 relaxation over a time τ
    C(E1, E2) = [E2 0 0; 0 E2 0; 0 0 E1]
    C(τ_ms) = C(exp(-τ_ms / T1_ms), exp(-τ_ms / T2_ms))
    D(τ_ms) = (I - C(τ_ms)) * [0 ; 0 ; Mz0]

    ## matrices for various values of τ
    P1 = P(TE_ms)
    P2 = P(TR_ms - TE_ms)
    C1 = C(TE_ms)
    C2 = C(TR_ms - TE_ms)
    d1 = D(TE_ms)
    d2 = D(TR_ms - TE_ms)

    ## matrix A and vector b for steady-state calculation
    A = P1 * C1 * R * P2 * C2
    b = P1 * C1 * R * d2 + d1

    Mss = (I - A) \ b # steady-state magnetization

    return complex(Mss[1], Mss[2]) # return the complex signal
end;



#=
## Recreate Figure 3 from [1]
And compare Method 1 above
with Method 2 (using `bssfp` in `BlochSim`)
=#

TR_ms, TE_ms = 10, 5 # scan parameters
Mz0, T1_ms, T2_ms = 1, 400, 100 # tissue parameters

num_off_res_values = 401 # vector of off-resonance values
Δf_arr_Hz = kHz_to_Hz(range(-1, 1, num_off_res_values) / TR_ms) # 2 periods

flip_ang_arr_deg = [15, 30, 60, 90]; # vector of flip angles

#src θ_rf_rad = -π/2 # x-axis rotation
θ_rf_rad = 0 # y-axis rotation

# Helper functions for broadcast:
bssfp_matrix(α_rad, Δf_Hz) = # method 1
    bssfp_matrix(
        Mz0, T1_ms, T2_ms, Δf_Hz,
        TR_ms, TE_ms, α_rad, θ_rf_rad,
    )
_bssfp(α_rad, Δf_Hz, Δϕ_rad) = # method 2
    bssfp(
        Mz0, T1_ms, T2_ms, Δf_Hz,
        TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad,
    )
_bssfp_ellipse(α_rad, Δf_Hz, Δϕ_rad) = # method 3
    bssfp_ellipse(
        Mz0, T1_ms, T2_ms, Δf_Hz,
        TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad,
    )

#=
Call `bssfp` and `bssfp_matrix` and `bssfp_ellipse`
for various flip angles and off-resonance values
and verify that the calculations match.
=#
sig_matrix = bssfp_matrix.(deg2rad.(flip_ang_arr_deg)', Δf_arr_Hz)
sig_blochsim = _bssfp.(deg2rad.(flip_ang_arr_deg)', Δf_arr_Hz, 0. #= Δϕ_rad =#)
@assert sig_matrix ≈ sig_blochsim # yes they match!

sig_ellipse = _bssfp_ellipse.(deg2rad.(flip_ang_arr_deg)', Δf_arr_Hz, 0. #= Δϕ_rad =#)
@assert sig_matrix ≈ sig_ellipse # yes they match!

# Plot 1-pool signal magnitude and phase
label = reshape(map(a -> "α = $(a)°", flip_ang_arr_deg), 1, :) # row!
p_m = plot(Δf_arr_Hz, abs.(sig_blochsim); label,
    ylabel = "Signal Magnitude")
p_p = plot(Δf_arr_Hz, angle.(sig_blochsim); label,
    xlabel = "Resonant Frequency (Hz)",
    ylabel = "Signal Phase")
pmp = plot(p_m, p_p, layout=(2,1), plot_title = "bSSFP single pool",
    plot_titlefontsize = 13)

#
prompt()


#=
Explore T1 dependence of bSSFP signal model
=#
T1_ms_arr = range(0.90, 1.1, 3) * T1_ms
α_deg = 20
sig_t1 = bssfp.(Mz0, T1_ms_arr', T2_ms, Δf_arr_Hz,
    TR_ms, TE_ms, 0. #= Δϕ_rad =#, deg2rad(α_deg))

label = reshape(map(t -> "T1 = $t ms", T1_ms_arr), 1, :) # row!
pt1_m = plot(Δf_arr_Hz, abs.(sig_t1); label,
    ylabel = "Signal Magnitude")
pt1_p = plot(Δf_arr_Hz, angle.(sig_t1); label,
    xlabel = "Resonant Frequency (Hz)",
    ylabel = "Signal Phase")
pt1 = plot(pt1_m, pt1_p, layout=(2,1), plot_title = "bSSFP single pool for T1",
    plot_titlefontsize = 13)

#
prompt()

# helper functions for CRB
real_imag(x) = [real(x); imag(x)] # stacker
snr2sigma(db::Real, yb::AbstractArray{<:Complex}) =
    10^(-db/20) * norm(yb) / sqrt(length(yb));

#=
## CRB for 1-pool bSSFP
Using arbitrary flip angles and phase cycling increments as scan "design"
=#
kappa = 1 # also estimate the B1+ factor
M0_phase = π/3 # just to make it non-trivial
x = [M0_phase, Mz0, T1_ms, T2_ms, kappa] # unknowns
Δf_Hz = -7 # known

Δϕ_rad = (0:7)/8 * 2π
α_rad = [π/7, π/3]
tmp = Iterators.product(Δϕ_rad, α_rad)
design = (
  Δϕ_rad = vec(map(x->x[1], tmp)),
  α_rad = vec(map(x->x[2], tmp)),
)

signal_c = (x) -> cis(x[1] #= M0_phase =#) *
    bssfp.(x[2:end-1]..., Δf_Hz,
        TR_ms, TE_ms, design.Δϕ_rad, design.α_rad * x[end] #= kappa =#)
signal_ri(x) = real_imag(signal_c(x))

#src tmp = signal_ri(x) # test run

dB = 40 # SNR
σ = snr2sigma(dB, signal_c(x)) # noise level

# Jacobian
jac = ForwardDiff.jacobian(signal_ri, x)
fish = jac' * jac / σ^2
cond(fish) # 9e8

# CRB and standard deviation:
crb = inv(fish)
crb_std = sqrt.(diag(crb))
cov1 = round.(crb_std ./ x ; digits=3) # coefficient of variation



#=
## RF pulse duration
The preceding calculations were all
for hypothetical instantaneous" RF pulses.

Examine effects of finite RF pulse duration
for a spin with a relatively short T2.
=#

Mz0, T1_ms, T2_ms, Δf_Hz = 1, 400, 40, 0 # tissue parameters
TR_ms, TE_ms, α_rad = 8, 4, deg2rad(50) # scan parameters

Δϕ_rad = range(-1,1,101) * π
rf0 = InstantaneousRF(α_rad)
_bssfp(Δϕ_rad, rf) = bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, rf);

#=
Specify finite-duration hard (rectangular) RF pulse
- `GAMMA` has units rad/s/G
- Tip angle for constant pulse:
  `α_rad = GAMMA * b1_gauss * tRF_s`
- so `b1_gauss = α_rad / GAMMA / tRF_s`
=#
#src tRF_ms = 1e-12 # super-short for first test
tRF_ms = 1
waveform1 = [1] * α_rad / (tRF_ms / 1000) / GAMMA # single sample i.e. instant!
rf1 = RF(waveform1, tRF_ms)
signal0 = map(Δϕ -> _bssfp(Δϕ, rf0), Δϕ_rad)
signal1 = map(Δϕ -> _bssfp(Δϕ, rf1), Δϕ_rad)
@assert signal0 ≈ signal1
@assert α_rad == rf0.α ≈ only(rf1.α)

# Plot
prfm = plot(
 xlabel = "phase cycling increment Δϕ (rad)",
 ylabel = "bSSFP signal mag",
)
prfa = plot(
 xlabel = "phase cycling increment Δϕ (rad)",
 ylabel = "bSSFP signal phase",
)
plot!(prfm, Δϕ_rad, abs.(signal0), label="Instantaneous")
plot!(prfa, Δϕ_rad, angle.(signal0), label="Instantaneous")
#src plot!(Δϕ_rad, abs.(signal1), label="tRF = $tRF_ms")
nw = 1000 # approximately 1μs dwell time
for tRF_ms in [1e-2 1 2]
    waveform2 = ones(nw) * α_rad / (tRF_ms / 1000) / GAMMA
    rf2 = RF(waveform2, tRF_ms/nw)
    @assert rf0.α ≈ sum(rf2.α)
    signal2 = map(Δϕ -> _bssfp(Δϕ, rf2), Δϕ_rad)
    label = "tRF = $tRF_ms ms, nw=$nw"
    plot!(prfm, Δϕ_rad, abs.(signal2); label)
    plot!(prfa, Δϕ_rad, angle.(signal2); label)
end;

prf = plot(prfm, prfa, layout=(2,1),
 plot_title="RF pulse duration effect, T2=$T2_ms (ms)")



#=
## Multi-compartment spins and myelin water exchange

Generate Figure 2 from [2] using `BlochSim`.
First define some useful helper functions.
These functions put the parameters in the correct format
for the multi-compartment spin object constructors.
=#

"""
- in: `f_f` fast fraction (myelin fraction)
- out: `mwf_tuple` tuple with fast and slow fractions
"""
get_mwf_tuple(f_f) = (f_f, 1-f_f)


"""
    get_τ_tuple(τ_fs_ms, f_f)
## In:
- `τ_fs_ms` residence time for exchange from myelin to non-myelin water (ms)
- `f_f` fast fraction (myelin fraction)
## Out:
- `τ_tuple_ms` tuple with fast-to-slow and slow-to-fast residence times
"""
function get_τ_tuple(τ_fs_ms, f_f)
    τ_sf_ms = (1-f_f) * τ_fs_ms / f_f
    return (τ_fs_ms, τ_sf_ms)
end


"""
## In:
- `Δϕ_rad` RF phase cycling value (radians)
- `Δf0_Hz` off-resonance value (Hz)
- `Δf_myelin_Hz` # additional off-resonance value only experienced by myelin water (Hz)
- `TR_ms` repetition time (ms)
## Out:
- `Δf_tuple_Hz` tuple with off-resonance values for fast and slow compartments

[Hinshaw, J. Appl. Phys. 1976](https://doi.org/10.1063/1.323136).
"""
function get_Δf_tuple(Δϕ_rad, Δf0_Hz, Δf_myelin_Hz, TR_ms)

    ## convert the RF phase cycling value to Hz from radians
    Δϕ_Hz = kHz_to_Hz((Δϕ_rad) / (2π*TR_ms))

    ## subtract the RF phase cycling value from the off-resonance value
    Δf_RF_Hz = Δf0_Hz - Δϕ_Hz

    ## add the myelin off-resonance for the myelin term
    Δf_myelin_RF_Hz = Δf_RF_Hz + Δf_myelin_Hz

    ## create and return tuple for the spin object constructor
    Δf_tuple_Hz = (Δf_myelin_RF_Hz, Δf_RF_Hz)
    return Δf_tuple_Hz
end;


# Define a function similar to Method 2 above,
# but for multi-compartment spin objects.

"""
    bssfp2(α_deg, TR_ms, TE_ms, spin, spin_no_rf_phase_fact)

Return steady-state magnetization signal value
at the echo time
for a phase-cycled bSSFP sequence
using BlochSim.

Ref: Murthy, N., Nielsen, J. F., Whitaker, S. T., Haskell, M. W.,
Swanson, S. D., Seiberlich, N., & Fessler, J. A. (2022).
Quantifying myelin water exchange using optimized bSSFP
sequences. In Proc. Intl. Soc. Mag. Res. Med (p. 2068). [2]

## In
- `α_deg` flip angle of RF pulse (degrees)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- 'spin' multi-compartment spin object with RF phase cycling factor
- 'spin_no_rf_phase_fact' multi-compartment spin object without RF phase cycling factor

## Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp2(
    α_deg::Number, TR_ms::Number, TE_ms::Number,
    spin::SpinMC, spin_no_rf_phase_fact::SpinMC,
)

    ## convert flip angle α from degrees to radians
    α_rad = deg2rad(α_deg)

    ## excite the spin and reshape R to be the correct dimensions for a SpinMC object
    (R,) = excite(spin, InstantaneousRF(α_rad))
    R = Matrix(R.A)
    R = kron(I(2),R)

    ## precession/relaxation of the spin for TR
    (PC_TR_A, PC_TR_B) = freeprecess(spin, TR_ms)

    ## precession/relaxation of the spin for TE
    ## assume receiver modulates signal and uses the receiver phase as the RF phase
    (PC_TE_A, PE_TE_B) = freeprecess(spin_no_rf_phase_fact, TE_ms)

    ## calculate A matrix and b vector
    A = Matrix(PC_TR_A) * R
    b = Vector(PC_TR_B)

    Mss = (I - A) \ b # steady-state just before tip down
    M = R * Mss # magnetization after tip-down

    ## steady-state magnetization at the echo time
    M = Matrix(PC_TE_A) * M + Vector(PE_TE_B)

    return complex(M[1]+M[4], M[2]+M[5]) # return the complex signal
end;


# Intermediate helper
function bssfp2(
    Mz0::Number,
    frac::NTuple{2, Number},
    T1_ms::NTuple{2, Number},
    T2_ms::NTuple{2, Number},
    Δf_Hz::NTuple{2, Number},
    Δf_no_rf_phase_Hz::NTuple{2, Number},
    τ_ms::NTuple{2, Number},
    α_deg::Number, TR_ms::Number, TE_ms::Number,
)

    ## create spin (with and without RF phase-cycling factor)
    spin = SpinMC(Mz0, frac, T1_ms, T2_ms, Δf_Hz, τ_ms)
global spin1 = spin
    spin_no_rf_phase = SpinMC(Mz0, frac, T1_ms, T2_ms, Δf_no_rf_phase_Hz, τ_ms)

    return bssfp2(α_deg, TR_ms, TE_ms, spin, spin_no_rf_phase)
end;


"""
    bssfp2(...)
Version with scalar arguments (convenient for autodiff)
"""
function bssfp2(
    M0_phase::Number, # radians
    Mz0::Number,
    f_f::Number,
    T1_f_ms::Number,
    T1_s_ms::Number,
    T2_f_ms::Number,
    T2_s_ms::Number,
    τ_fs_ms::Number,
    Δff_Hz::Number, # fast component frequency shift
    ## system parameter (sometimes known):
    Δf0_Hz::Number, # B0 off resonance
    ## scan parameters (always known):
    α_deg::Number, TR_ms::Number, TE_ms::Number,
    Δϕ_deg::Number, # RF phase cycling value (degrees)
)

    τ_tuple_ms = get_τ_tuple(τ_fs_ms, f_f) # fast-to-slow and slow-to-fast residence times

    ## tuple of values incorporating off-resonance and RF phase cycling for both compartments
    Δϕ_rad = deg2rad(Δϕ_deg)
    Δf_tuple_Hz = get_Δf_tuple(Δϕ_rad, Δf0_Hz, Δff_Hz, TR_ms)
    Δf_tuple_Hz_no_rf_phase = get_Δf_tuple(0, Δf0_Hz, Δff_Hz, TR_ms)

    return cis(M0_phase) * bssfp2(
     Mz0,
     (f_f, 1 - f_f),
     (T1_f_ms, T1_s_ms),
     (T2_f_ms, T2_s_ms),
     Δf_tuple_Hz,
     Δf_tuple_Hz_no_rf_phase,
     τ_tuple_ms,
     α_deg, TR_ms, TE_ms,
    )
end;


#=
Define variables to be used in the following plots.
Values taken from [2] and [4].
=#

Mz0 = 0.77 # initial condition for longitudinal magnetization (constant)
Δf_myelin_Hz = 5.0 # frequency shift of myelin water
f_f = 0.15; # myelin water fraction (MWF), fast fraction

## T1 and T2 values
T1_f_ms = 400 # T1 for fast-relaxing, myelin water compartment
T1_s_ms = 832 # T1 for slow-relaxing, non-myelin water compartment

T2_f_ms = 20 # T2 for fast-relaxing, myelin water compartment
T2_s_ms = 80 # T2 for slow-relaxing, non-myelin water compartment

TR_ms, TE_ms = 20, 4; # scan parameters (note: TE = TR/2 would be more logical)


#=
Generate plots similar to Figure 3 from [1]
but with three different RF phase cycling factor values:
(0, 90, and 180 degrees).

For this example, choose one exchange rate:
=#
τ_fs = 50.0 # this will be varied in the next plot

flip_ang_arr_deg = [10, 40] # flip angles

Δϕ_arr_deg = [0, 90, 180] # RF phase cycling value (degrees)

## vector of off-resonance values
num_samples = 401 # number of samples (resonant frequencies)
Δf_arr_Hz = kHz_to_Hz(range(-1, 1, num_samples) / TR_ms);


# Broadcast via `map` using helper functions
bssfp_mc(Δf0_Hz::Number, Δϕ_deg::Number, α_deg::Number, τ_fs::Number) =
    bssfp2(0 #= phase =#, Mz0, f_f, T1_f_ms, T1_s_ms, T2_f_ms, T2_s_ms, τ_fs,
        Δf_myelin_Hz, Δf0_Hz, α_deg, TR_ms, TE_ms, Δϕ_deg)

tmp = Iterators.product(Δf_arr_Hz, Δϕ_arr_deg, flip_ang_arr_deg)
bssfp_mc(t3::NTuple{3, Any}) = bssfp_mc(t3..., τ_fs) # (Δf_Hz, Δϕ_deg, α_deg)
signal_mc = map(bssfp_mc, tmp);

# Plot 2-pool signals
pmcm = plot(ylabel = "Signal Magnitude")
pmcp = plot( ylabel = "Signal Phase", xlabel = "Resonant Frequency (Hz)")
for i in 1:size(signal_mc,3), j in 1:size(signal_mc,2)
    label2 = "α = $(flip_ang_arr_deg[i])°, Δϕ = $(Δϕ_arr_deg[j])°"
    tmp2 = signal_mc[:,j,i]
    plot!(pmcm, Δf_arr_Hz, abs.(tmp2))
    plot!(pmcp, Δf_arr_Hz, angle.(tmp2); label = label2)
end

p2 = plot(pmcm, pmcp, layout = (2,1), titlefontsize = 12,
    plot_title = "bSSFP 2-pool magnitude and phase")

#
prompt()


#=
## Recreate Figure 2 (magnitude plot) from [2]
and also plot the phase.
=#

Δf_Hz = 0.0 # set off-resonance to zero for this plot
tau_arr_ms = [250, 150, 50] # array of exchange values
tau_arr_marker = [:circle, :star5, :utriangle]

Δϕ_design_deg = ( # designed RF phase cycling increments
 [-176.4, -159.5, -142.1, -124.4, -107.6, -90.54, -73.62, -56.13, -39.41, -22.52, -5.272, 11.63, 28.93, 45.76, 63.08, 79.91, 96.97, 113.9, 131.3, 148.5, 166.1],
 [-168.8, -150.3, -130.1, -111.5, -93.19, -74.18, -54.68 , -37.15, -18.01, 1.342, 18.82, 38.64, 57.88, 76.48, 95.2, 113.3, 133.3, 153.1, 172.1],
)

α_arr_deg = [10.0, 40.0] # flip angles for plot
α_design_deg = [
  fill(α_arr_deg[1], length(Δϕ_design_deg[1]));
  fill(α_arr_deg[2], length(Δϕ_design_deg[2]))
]
design = (Δϕ_deg = vcat(Δϕ_design_deg...), α_deg = α_design_deg)
num_scans = length(design.Δϕ_deg) # number of different scans = 40

tmp = (τ_fs) -> (Δϕ_deg, α_deg) -> bssfp_mc(Δf_Hz, Δϕ_deg, α_deg, τ_fs)
bssfp_signal(τ_fs) = map(splat(tmp(τ_fs)), zip(design...))

signal = bssfp_signal.(tau_arr_ms)

## Plot
scan_idx = 1:num_scans
p_m = plot(title="Signal Magnitude vs. Scan Index", ylabel = "Signal Magnitude")
p_p = plot(title="Signal Phase vs. Scan Index", ylabel = "Signal Phase")
for j = 1:length(signal) # iterate over exchange values
    markershape = tau_arr_marker[j]
    local label = latexstring("\$τ_{\\mathrm{fs}}\$ = $τ_fs ms")
    plot!(p_m, scan_idx, abs.(signal[j]), linewidth=0; markershape, label)
    plot!(p_p, scan_idx, angle.(signal[j]), linewidth=0; markershape, label)
end
p3 = plot(p_m, p_p, layout = (2,1), xlabel = "Scan Index")

#
prompt()


#=
## Cramer-Rao Bound
for the designed scan
=#
kappa = 1 # also estimate the B1+ factor
M0_phase = π/3 # just to make it non-trivial
x = [M0_phase, Mz0, f_f, T1_f_ms, T1_s_ms, T2_f_ms, T2_s_ms, τ_fs, Δf_myelin_Hz, kappa] # unknowns
signal_c = (x) -> bssfp.(x[1:end-1]..., Δf_Hz, x[end]*design.α_deg, TR_ms, TE_ms, design.Δϕ_deg)
signal_ri(x) = real_imag(signal_c(x));
#src tmp = signal_ri(x) # test run

# Noise level
dB = 40 # SNR
σ = snr2sigma(dB, signal_c(x))

# Jacobian
jac = ForwardDiff.jacobian(signal_ri, x)
fish = jac' * jac / σ^2
cond(fish) # 5e10

# The condition number depends on units, so remove units:
D = Diagonal(1 ./ sqrt.(diag(fish)))
tmp = D * fish * D
cond(tmp) # 1e6

# CRB and standard deviation:
crb = inv(fish)
crb_std = sqrt.(diag(crb))

# Coefficient of variation
round.(crb_std ./ x ; digits=2)


include("../../../inc/reproduce.jl")
