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
        "LaTeXStrings"
        "LinearAlgebra"
        "MIRTjim"
        "Plots"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run `Pkg.add()` in the preceding code block first, if needed.

using BlochSim: Spin, SpinMC, InstantaneousRF, excite, freeprecess
using InteractiveUtils: versioninfo
using LaTeXStrings: latexstring
using LinearAlgebra: I
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
    bssfp_matrix(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_kHz=0)

Return steady-state magnetization signal value
at the echo time
for a bSSFP sequence
using method of
[Hargreaves et al., MRM 2001](https://doi.org/10.1002/mrm.1170).

## In
- `α_deg` flip angle of RF pulse (degrees)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- `mo` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `Δf_Hz` off-resonance value (Hz) (default 0)

## Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_matrix(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_Hz=0)

    M0 = [0; 0; mo] # initial magnetization vector

    ## rotation matrix for RF excitation about the x-axis
    α_rad = deg2rad(α_deg) # convert flip angle α from degrees to radians
    R = [1 0 0; 0 cos(α_rad) sin(α_rad); 0 -sin(α_rad) cos(α_rad)]

    ## free precession matrix
    Pz(angle) = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1]
    P(τ_ms) = Pz( 2π * Hz_to_kHz(Δf_Hz) * τ_ms ) # angle in radians

    ## matrices for T1 and T2 relaxation over a time τ
    C(E1, E2) = [E2 0 0; 0 E2 0; 0 0 E1]
    C(τ_ms) = C(exp(-τ_ms / T1_ms), exp(-τ_ms / T2_ms))
    D(τ_ms) = (I - C(τ_ms)) * [0 ; 0 ; mo]

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


# ## Method 2: Use BlochSim

"""
    bssfp_blochsim(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_kHz=0)
    bssfp_blochsim(α_deg, TR_ms, TE_ms, spin)

Return steady-state magnetization signal value
at the echo time
for a bSSFP sequence
using BlochSim.
See [Hargreaves et al., MRM 2001](https://doi.org/10.1002/mrm.1170).

## In
- `α_deg` flip angle of RF pulse (degrees)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- `mo` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `Δf_Hz` off-resonance value (Hz)

## Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_blochsim(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_Hz=0)
    spin = Spin(mo, T1_ms, T2_ms, Δf_Hz) # create a spin
    return bssfp_blochsim(α_deg, TR_ms, TE_ms, spin)
end;


function bssfp_blochsim(α_deg, TR_ms, TE_ms, spin::Spin)
    α_rad = deg2rad(α_deg) # convert flip angle α from degrees to radians

    #=
    Matrix to excite the spin
    include RF phase for instantaneous RF because above code flips over x axis
    and blochsim flips over -y axis and want to make them consistent
    =#
    (R,) = excite(spin, InstantaneousRF(α_rad, -π/2))

    ## Matrices for precession/relaxation for various time period values
    (PC1_A, PC1_B) = freeprecess(spin, TE_ms)
    (PC2_A, PC2_B) = freeprecess(spin, TR_ms-TE_ms)
    (PC_TR_A, PC_TR_B) = freeprecess(spin, TR_ms)

    ## calculate the A matrix and b vector
    A = Matrix(PC1_A * R.A * PC2_A)
    b = Matrix(PC1_A * R.A) * Vector(PC2_B) + Vector(PC1_B)

    ## calculate the steady-state magnetization at the echo time
    Mss = (I - A) \ b

    return complex(Mss[1], Mss[2]) # return the complex signal
end;


# ## Recreate Figure 3 from [1] using Methods 1 and 2

TR_ms, TE_ms = 10, 5 # scan parameters
mo, T1_ms, T2_ms = 1, 400, 100 # tissue parameters

num_off_res_values = 401 # vector of off-resonance values
Δf_arr_Hz = kHz_to_Hz(range(-1, 1, num_off_res_values) / TR_ms) # 2 periods

flip_ang_arr_deg = [15, 30, 60, 90]; # vector of flip angles

# Helper functions for broadcast:
bssfp_matrix(α_deg, Δf_Hz) =
    bssfp_matrix(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_Hz)
bssfp_blochsim(α_deg, Δf_Hz) =
    bssfp_blochsim(α_deg, TR_ms, TE_ms, mo, T1_ms, T2_ms, Δf_Hz);

#=
Call `bssfp_matrix` and `bssfp_blochsim`
for various flip angles and off-resonance values
and verify that the calculations match.
=#
sig_matrix = bssfp_matrix.(flip_ang_arr_deg', Δf_arr_Hz)
sig_blochsim = bssfp_blochsim.(flip_ang_arr_deg', Δf_arr_Hz)
@assert sig_matrix ≈ sig_blochsim # yes they match!

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
## In:
- `τ_fs_ms` residence time for exchange from myelin to non-myelin water (ms)
- `f_f` fast fraction (myelin fraction)
## Out:
- `τ_tuple_ms` tuple with fast-to-slow and slow-to-fast residence times
"""
function get_τ_tuple(τ_fs_ms, f_f)
    τ_sf_ms = (1-f_f) * τ_fs_ms / f_f
    τ_tuple_ms = (τ_fs_ms, τ_sf_ms)
    return τ_tuple_ms
end


"""
## In:
- `ΔΦ_rad` RF phase cycling value (radians)
- `Δf_Hz` off-resonance value (Hz)
- `Δf_myelin_Hz` # additional off-resonance value only experienced by myelin water (Hz)
- `TR_ms` repetition time (ms)
## Out:
- `Δf_tuple_Hz` tuple with off-resonance values for fast and slow compartments

[Hinshaw, J. Appl. Phys. 1976](https://doi.org/10.1063/1.323136).
"""
function get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)

    ## convert the RF phase cycling value to Hz from radians
    ΔΦ_Hz = kHz_to_Hz((ΔΦ_rad)/(2π*TR_ms))

    ## subtract the RF phase cycling value from the off-resonance value
    Δf_RF_Hz = Δf_Hz - ΔΦ_Hz

    ## add the myelin off-resonance for the myelin term
    Δf_myelin_RF_Hz = Δf_RF_Hz + Δf_myelin_Hz

    ## create and return tuple for the spin object constructor
    Δf_tuple_Hz = (Δf_myelin_RF_Hz, Δf_RF_Hz)
    return Δf_tuple_Hz
end;


# Define a function similar to Method 2 above,
# but for multi-compartment spin objects.

"""
    bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)

Return steady-state magnetization signal value
at the echo time
for a bSSFP sequence
using BlochSim.

Ref: Murthy, N., Nielsen, J. F., Whitaker, S. T., Haskell, M. W.,
Swanson, S. D., Seiberlich, N., & Fessler, J. A. (2022).
Quantifying myelin water exchange using optimized bSSFP
sequences. In Proc. Intl. Soc. Mag. Res. Med (p. 2068). [2]

## In
- `α_deg` flip angle of RF pulse (degrees)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- 'spin_mc' multi-compartment spin object with RF phase cycling factor
- 'spin_mc_no_rf_phase_fact' multi-compartment spin object without RF phase cycling factor

## Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)

    ## convert flip angle α from degrees to radians
    α_rad = deg2rad(α_deg)

    ## excite the spin and reshape R to be the correct dimensions for a SpinMC object
    (R,) = excite(spin_mc, InstantaneousRF(α_rad))
    R = Matrix(R.A)
    R = kron(I(2),R)

    ## precession/relaxation of the spin for TR
    (PC_TR_A, PC_TR_B) = freeprecess(spin_mc, TR_ms)

    ## precession/relaxation of the spin for TE
    ## assume receiver modulates signal and uses the receiver phase as the RF phase
    (PC_TE_A, PE_TE_B) = freeprecess(spin_mc_no_rf_phase_fact, TE_ms)

    ## calculate A matrix and b vector
    A = Matrix(PC_TR_A) * R
    b = Vector(PC_TR_B)

    Mss = (I - A) \ b # steady-state just before tip down
    M = R * Mss # magnetization after tip-down

    ## steady-state magnetization at the echo time
    M = Matrix(PC_TE_A) * M + Vector(PE_TE_B)

    return complex(M[1]+M[4], M[2]+M[5]) # return the complex signal
end;


#=
Define variables to be used in the following plots.
Values taken from [2] and [4].
=#

mo = 0.77 # initial condition for longitudinal magnetization (constant)
Δf_myelin_Hz = 5.0 # frequency shift of myelin water
f_f = 0.15; # myelin water fraction (MWF), fast fraction

## T1 and T2 values
T1_f_ms = 400 # T1 for fast-relaxing, myelin water compartment
T1_s_ms = 832 # T1 for slow-relaxing, non-myelin water compartment
T1_ms_tuple = (T1_f_ms, T1_s_ms)

T2_f_ms = 20 # T2 for fast-relaxing, myelin water compartment
T2_s_ms = 80 # T2 for slow-relaxing, non-myelin water compartment
T2_ms_tuple = (T2_f_ms, T2_s_ms)

TR_ms, TE_ms = 20, 4; # scan parameters

mwf_tuple = get_mwf_tuple(f_f) # tuple with fast and slow relaxing fractions


#=
Generate plots similar to Figure 3 from [1]
but with three different RF phase cycling factor values:
(0, 90, and 180 degrees).

For this example, choose one exchange rate:
=#
τ_fs = 50.0 # this will be varied in the next plot
## tuple with fast-to-slow and slow-to-fast residence times
τ_tuple_ms = get_τ_tuple(τ_fs, f_f)

flip_ang_arr_deg = [10, 40] # flip angles

ΔΦ_arr_deg = [0, 90, 180] # RF phase cycling value (degrees)

## vector of off-resonance values
num_samples = 401 # number of samples (resonant frequencies)
Δf_arr_Hz = kHz_to_Hz(range(-1, 1, num_samples) / TR_ms)

# Helper function for broadcast
function bssfp_blochsim_MC(Δf_Hz, ΔΦ_deg, α_deg)
    ## convert inputted RF phase cycling angle from degrees to radians
    ΔΦ_rad = deg2rad(ΔΦ_deg)

    ## get tuple of values incorporating off-resonance and RF phase cycling for both compartments
    Δf_tuple_Hz = get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)
    Δf_tuple_Hz_no_rf_phase_fact = get_Δf_tuple(0, Δf_Hz, Δf_myelin_Hz, TR_ms)

    ## create spin (with and without RF phase-cycling factor)
    spin_mc = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz, τ_tuple_ms)
    spin_mc_no_rf_phase_fact = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz_no_rf_phase_fact, τ_tuple_ms)

    return bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)
end
bssfp_blochsim_MC(t::Tuple) = bssfp_blochsim_MC(t...);

tmp = Iterators.product(Δf_arr_Hz, ΔΦ_arr_deg, flip_ang_arr_deg)
signal_MC = map(bssfp_blochsim_MC, tmp);

# Plot 2-pool signals
pmcm = plot(ylabel = "Signal Magnitude")
pmcp = plot( ylabel = "Signal Phase", xlabel = "Resonant Frequency (Hz)")
for i in 1:size(signal_MC,3), j in 1:size(signal_MC,2)
    label2 = "α = $(flip_ang_arr_deg[i])°, ΔΦ = $(ΔΦ_arr_deg[j])°"
    tmp2 = signal_MC[:,j,i]
    plot!(pmcm, Δf_arr_Hz, abs.(tmp2))
    plot!(pmcp, Δf_arr_Hz, angle.(tmp2); label = label2)
end

p2 = plot(pmcm, pmcp, layout = (2,1), titlefontsize = 12,
    plot_title = "bSSFP 2-pool magnitude and phase")

#
prompt()


#=
# Recreate Figure 2 (magnitude plot) from [2] and also add the phase plot.
=#

## Initialize the plot:
p_m = plot(title="Signal Magnitude vs. Scan Index", ylabel = "Signal Magnitude")
p_p = plot(title="Signal Phase vs. Scan Index", ylabel = "Signal Phase")

num_scans = 40 # number of different scans
scan_idx = 1:num_scans

flip_ang_arr_deg = [10.0, 40.0] # flip angles for plot
num_flip_angles = length(flip_ang_arr_deg)

Δf_Hz = 0.0 # set off-resonance to zero

tau_arr_ms = [250, 150, 50] # array of exchange values
tau_arr_marker = [:circle, :star5, :utriangle]
num_taus = length(tau_arr_ms)

sig_arr = zeros(num_scans,num_taus) # arrays to store results
sig_arr_phase = zeros(num_scans,num_taus)

ΔΦ_design_deg = ( # designed RF phase cycling increments
 [-176.4, -159.5, -142.1, -124.4, -107.6, -90.54, -73.62, -56.13, -39.41, -22.52, -5.272, 11.63, 28.93, 45.76, 63.08, 79.91, 96.97, 113.9, 131.3, 148.5, 166.1],
 [-168.8, -150.3, -130.1, -111.5, -93.19, -74.18, -54.68 , -37.15, -18.01, 1.342, 18.82, 38.64, 57.88, 76.48, 95.2, 113.3, 133.3, 153.1, 172.1],
)

curr_scan = 1

for j in 1:num_taus # iterate over exchange values
    local τ_fs = tau_arr_ms[j]
    local τ_tuple_ms = get_τ_tuple(tau_arr_ms[j], f_f)
    tau_marker = tau_arr_marker[j]

    for k in 1:num_flip_angles # iterate over flip angles
        α_deg = flip_ang_arr_deg[k]

        ## different RF phases for different flip angles - from Figure 1 in [2]
        local ΔΦ_arr_deg = ΔΦ_design_deg[k]

        for i in 1:length(ΔΦ_arr_deg) # iterate over RF phases
            ΔΦ_deg = ΔΦ_arr_deg[i]

            ## convert RF phase cycling angle from degrees to radians
            ΔΦ_rad = deg2rad(ΔΦ_deg)

            ## tuple of values incorporating off-resonance and RF phase cycling for both compartments
            Δf_tuple_Hz = get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)
            Δf_tuple_Hz_no_rf_phase_fact = get_Δf_tuple(0, Δf_Hz, Δf_myelin_Hz, TR_ms)

            ## create a spin (with and without RF phase-cycling factor)
            spin_mc = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz, τ_tuple_ms)
            spin_mc_no_rf_phase_fact = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz_no_rf_phase_fact, τ_tuple_ms)

            ## run the bSSFP blochsim and add to result array
            signal = bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)
            sig_arr[curr_scan,j] = abs(signal)
            sig_arr_phase[curr_scan,j] = angle(signal)

            global curr_scan += 1
        end
    end

    global curr_scan = 1

    plot!(p_m, scan_idx, sig_arr[:,j], linewidth=0, markershape=tau_marker,
        label = latexstring("\$τ_{\\mathrm{fs}}\$ = $τ_fs ms"))
    plot!(p_p, scan_idx, sig_arr_phase[:,j], linewidth=0, markershape=tau_marker,
        label = latexstring("\$τ_{\\mathrm{fs}}\$ = $τ_fs ms"))
end

## plot results and label axes
p3 = plot(p_m, p_p, layout = (2,1), xlabel = "Scan Index")

#
prompt()

## gui(); throw(); # xx

include("../../../inc/reproduce.jl")
