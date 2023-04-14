#=
# [bSSFP](@id 01-bssfp)

This page illustrates using the Julia package
[`BlochSim`](https://github.com/StevenWhitaker/BlochSim.jl)
to calculate MRI signals
for balanced steady-state free precession (bSSFP) pulse sequences.

This demo facilitates
understanding bSSFP sequences,
multi-compartment spins,
and myelin water exchange.

This demo recreates Figure 3 from [1] and Figure 2 from [2].
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
        "LinearAlgebra"
        "Plots"
    ])
end


# Tell this Julia session to use the following packages for this example.
# Run Pkg.add() in the preceding code block first, if needed.

using BlochSim:
using InteractiveUtils: versioninfo
using LinearAlgebra:
using MIRTjim: prompt
using Plots: plot


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() || prompt(:draw);


# Define some useful helper functions.

# convert angles in degrees to radians
function deg_to_rad(θ_deg)
    θ_rad = θ_deg*pi/180
    return θ_rad
end;

# convert frequencies in Hz to kHz
function Hz_to_kHz(Δf_Hz)
    Δf_kHz = Δf_Hz*10^(-3)
    return Δf_kHz
end;

# convert frequencies in kHz to Hz
function kHz_to_Hz(Δf_kHz)
    Δf_Hz = Δf_kHz*10^(3)
    return Δf_Hz
end;


#=
This is the bSSFP pulse sequence from Figure 2 in [1]
that we will use to generate Figure 3 in [1] in two different ways.
The RF excitation is repeated periodically and, in steady-state,
the magnetization at point *a* is the same as at point *d*.

todo <img src="hargreaves_fig2.png" width="500"/>
=#


#=
## Method 1

Use Equations 1 and 2 and Appendix A from [1]

Calculate the steady-state value at point *d*
using the method from [1] using Equations 1 and 2 and Appendix A.
=#

"""
    bssfp_matrix(α_deg, Δf_kHz, mo, T1_ms, T2_ms, TR_ms, TE_ms)

Return steady-state magentization signal value for a bssfp sequence

Ref: Hargreaves, B. A., Vasanawala, S. S., Pauly, J. M., & Nishimura,
D. G. (2001). Characterization and reduction of the transient response
in steady‐state MR imaging. Magnetic Resonance in Medicine: An Official
Journal of the International Society for Magnetic Resonance in Medicine,
46(1), 149-158. [1]

# In
- `α_deg` flip angle of RF pulse (degrees)
- `Δf_Hz` off-resonance value (Hz)
- `mo` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)


# Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_matrix(α_deg, Δf_Hz, mo, T1_ms, T2_ms, TR_ms, TE_ms)

    # convert off-resonance value to kHz
    Δf_kHz = Hz_to_kHz(Δf_Hz)

    # initial magnetization vector
    M0 = [0; 0; mo]

    # convert inputted flip angle α from degrees to radians
    α_rad = deg_to_rad(α_deg)

    # create rotation matrix for RF excitation about the x-axis
    R = [1 0 0; 0 cos(α_rad) sin(α_rad); 0 -sin(α_rad) cos(α_rad)]

    # create free precession matrix
    P(τ_ms) = [cos(2*pi*Δf_kHz*τ_ms) sin(2*pi*Δf_kHz*τ_ms) 0 ; -sin(2*pi*Δf_kHz*τ_ms) cos(2*pi*Δf_kHz*τ_ms) 0 ; 0 0 1]

    # create matrices for T1 and T2 relaxation over a time τ
    C(τ_ms) = [exp(-τ_ms/T2_ms) 0 0 ; 0 exp(-τ_ms/T2_ms) 0 ; 0 0 exp(-τ_ms/T1_ms)]
    D(τ_ms) = (I - C(τ_ms))*[0 ; 0 ; mo]

    # create matrices for various values of τ
    P1 = P(TE_ms)
    P2 = P(TR_ms - TE_ms)
    C1 = C(TE_ms)
    C2 = C(TR_ms - TE_ms)
    D1 = D(TE_ms)
    D2 = D(TR_ms - TE_ms)

    # formulate matrices A and B for steady-state calculation
    A = P1*C1*R*P2*C2
    B = P1*C1*R*D2 + D1

    # calculate the steady-state magnetization
    Mss = (I - A)\B

    # return the complex signal
    signal = complex(Mss[1], Mss[2])
    return signal

end


# ## Method 2: Use Julia Package BlochSim [3]

"""
    bssfp_blochsim(α_deg, Δf_kHz, mo, T1_ms, T2_ms, TR_ms, TE_ms)
    bssfp_blochsim(α_deg, TR_ms, TE_ms, spin)

Return steady-state magentization signal value for a bSSFP sequence using [3]

Ref: Hargreaves, B. A., Vasanawala, S. S., Pauly, J. M., & Nishimura,
D. G. (2001). Characterization and reduction of the transient response
in steady‐state MR imaging. Magnetic Resonance in Medicine: An Official
Journal of the International Society for Magnetic Resonance in Medicine,
46(1), 149-158. [1]

# In
- `α_deg` flip angle of RF pulse (degrees)
- `Δf_Hz` off-resonance value (Hz)
- `mo` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)


# Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_blochsim(α_deg, Δf_Hz, mo, T1_ms, T2_ms, TR_ms, TE_ms)
    spin = Spin(mo,T1_ms,T2_ms,Δf_Hz) # create a spin
    return bssfp_blochsim(α_deg, TR_ms, TE_ms, spin)
end


function bssfp_blochsim(α_deg, TR_ms, TE_ms, spin::Spin)

    # convert inputted flip angle α from degrees to radians
    α_rad = deg_to_rad(α_deg)

    # excite the spin
    # included RF phase for instantaneous RF because above code flips over x axis
    # and blochsim flips over -y axis and want to make them consistent
    (R,) = excite(spin, InstantaneousRF(α_rad, -pi/2))
    R = Matrix(R.A)

    # put spin through precession/relaxation for various time period values
    (PC1_A,PC1_B) = freeprecess(spin, TE_ms)
    (PC2_A,PC2_B) = freeprecess(spin, TR_ms-TE_ms)
    (PC_TR_A,PC_TR_B) = freeprecess(spin, TR_ms)

    # calculate the A and B matrices
    A = Matrix(PC1_A)*R*Matrix(PC2_A)
    B = Matrix(PC1_A)*R*Vector(PC2_B)+Vector(PC1_B)

    # calculate the steady-state magnetization at the echo time
    Mss = (I - A) \ B

    # return the complex signal
    signal = complex(Mss[1], Mss[2])
    return signal

end


# ## Recreate Figure 3 from [1] using Methods 1 and 2

# set parameter values
mo = 1
T1_ms = 400
T2_ms = 100
TR_ms = 10
TE_ms = 5

# array of off-resonance values
num_off_res_values = 100
Δf_arr_kHz = range(-1/TR_ms, 1/TR_ms, num_off_res_values)

# array of flip angles
flip_ang_arr_deg = [15 30 60 90]
num_flip_angles = length(flip_ang_arr_deg)

# array to store calculated results for both plots (methods 1 and 2)
num_plots = 2
sig_arr = zeros(num_flip_angles,num_off_res_values,num_plots)

# initialize the plot
lay = @layout [a ; b]
p_m = plot(title="Steady-State Signal Magnitude vs. Resonant Frequency: Matrix Version", titlefontsize=10)
p_b = plot(title="Steady-State Signal Magnitude vs. Resonant Frequency: BlochSim Version", titlefontsize=10)

#=
call bssfp_matrix and bssfp_blochsim
for various flip angles and off-resonance values
iterate over flip angles
=#
for i = 1:num_flip_angles
    α_deg = flip_ang_arr_deg[i]
    α_str = string(α_deg)

    # iterate over off-resonance values
    for j = 1:num_off_res_values
        Δf_kHz = Δf_arr_kHz[j]

        # convert from kHz to Hz before input into function
        Δf_Hz = kHz_to_Hz(Δf_kHz)

        # call both implementations (methods 1 and 2) of bssfp signal model
        signal_matrix = bssfp_matrix(α_deg, Δf_Hz, mo, T1_ms, T2_ms, TR_ms, TE_ms)
        signal_blochsim = bssfp_blochsim(α_deg, Δf_Hz, mo, T1_ms, T2_ms, TR_ms, TE_ms)

        # add results for methods 1 and 2 to results array
        sig_arr[i,j,1] = abs(signal_matrix)
        sig_arr[i,j,2] = abs(signal_blochsim)
    end

    # plot results for current flip angle
    plot!(p_m, Δf_arr_kHz, sig_arr[i,:,1],label="α="*α_str)
    plot!(p_b, Δf_arr_kHz, sig_arr[i,:,2],label="α="*α_str)

end

# plot results and label axes
p = plot(p_m,p_b, layout = lay)
xlabel!(p, "Resonant Frequency (kHz)")
ylabel!(p, "Signal Magnitude")


# ## Multi-Compartment Spins and Myelin Water Exchange

#=
We will generate Figure 2 from [2] using [3].
First define some useful helper functions.
These functions put the parameters in the correct format
for the multi-compartment spin object constructors.
=#

# - in: `f_f` fast fraction (myelin fraction)
# - out: `mwf_tuple` tuple with fast and slow fractions
function get_mwf_tuple(f_f)
    mwf_tuple = (f_f, 1-f_f)
    return mwf_tuple
end


# in:
# - `τ_fs_ms` residence time for exchange from myelin to non-myelin water (ms)
# - `f_f` fast fraction (myelin fraction)
# out: `τ_tuple_ms` tuple with fast-to-slow and slow-to-fast residence times
function get_τ_tuple(τ_fs_ms,f_f)
    τ_sf_ms = (1-f_f)*τ_fs_ms/f_f
    τ_tuple_ms = (τ_fs_ms,τ_sf_ms)
    return τ_tuple_ms
end


# in:
# `T1_f_ms` T1 value for fast-relaxing (myelin) compartment (ms)
# `T1_s_ms` T1 value for slow-relaxing (non-myelin) compartment (ms)
# out: `T1_tuple_ms` tuple with both T1 values
function get_T1_tuple(T1_f_ms, T1_s_ms)
    T1_tuple_ms = (T1_f_ms, T1_s_ms)
    return T1_tuple_ms
end


# in:
# `T2_f_ms` T2 value for fast-relaxing (myelin) compartment (ms)
# `T2_s_ms` T2 value for slow-relaxing (non-myelin) compartment (ms)
# out: `T2_tuple_ms` tuple with both T2 values
function get_T2_tuple(T2_f_ms, T2_s_ms)
    T2_tuple_ms = (T2_f_ms, T2_s_ms)
    return T2_tuple_ms
end


# in:
# `ΔΦ_rad` RF phase cycling value (radians)
# `Δf_Hz` off-resonance value (Hz)
# `Δf_myelin_Hz` # additional off-resonance value only experienced by myelin water (Hz)
# `TR_ms` repetition time (ms)
# out: `Δf_tuple_Hz` tuple with off-resonance values for fast and slow compartments

# Ref: W. S. Hinshaw, "Image formation by nuclear magnetic resonance:
# The sensitive-point method", J. of Appl. Phys., 1976. [4]

function get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)

    # convert the RF phase cycling value to Hz from radians
    ΔΦ_Hz = kHz_to_Hz((ΔΦ_rad)/(2*pi*TR_ms))

    # subtract the RF phase cycling value from the off-resonance value
    Δf_RF_Hz = Δf_Hz - ΔΦ_Hz

    # add the myelin off-resonance for the myelin term
    Δf_myelin_RF_Hz = Δf_RF_Hz + Δf_myelin_Hz

    # create and return tuple for the spin object constructor
    Δf_tuple_Hz = (Δf_myelin_RF_Hz, Δf_RF_Hz)
    return Δf_tuple_Hz
end


# Define a function similar to Method 2 above,
# but for multi-compartment spin objects.

"""
bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact )

Return steady-state magentization signal value for a bssfp sequence using [3]

Ref: Murthy, N., Nielsen, J. F., Whitaker, S. T., Haskell, M. W.,
Swanson, S. D., Seiberlich, N., & Fessler, J. A. (2022).
Quantifying myelin water exchange using optimized bSSFP
sequences. In Proc. Intl. Soc. Mag. Res. Med (p. 2068). [2]

# In
- `α_deg` flip angle of RF pulse (degrees)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- 'spin_mc' multi-compartment spin object with RF phase cycling factor
- 'spin_mc_no_rf_phase_fact' multi-compartment spin object without RF phase cycling factor

# Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)

    # convert inputted flip angle α from degrees to radians
    α_rad = deg_to_rad(α_deg)

    # excite the spin and reshape R to be the correct dimensions for a SpinMC object
    (R,) = excite(spin_mc, InstantaneousRF(α_rad))
    R = Matrix(R.A)
    R = kron(I(2),R)

    # precession/relaxation of the spin for TR
    (PC_TR_A,PC_TR_B) = freeprecess(spin_mc, TR_ms)

    # precession/relaxation of the spin for TE
    # assume receiver modulates signal and uses the receiver phase as the RF phase
    (PC_TE_A,PE_TE_B) = freeprecess(spin_mc_no_rf_phase_fact, TE_ms)

    # calculate the A and B matrices
    A = Matrix(PC_TR_A) * R
    B = Vector(PC_TR_B)

    # steady-state just before tip down
    Mss = (I - A)\B

    # magnetization after tip-down
    M = R * Mss

    # # calculate the steady-state magnetization at the echo time
    M = Matrix(PC_TE_A) * M + Vector(PE_TE_B)

    # return the complex signal
    signal = (complex(M[1]+M[4], M[2]+M[5]))
    return signal

end


# Define variables to be used in the following plots.
# Values taken from [2] and [5].

# initial condition for magnetization in the z-direction (constant)
mo = 0.77

# T1 and T2 values
T1_f_ms = 400 # T1 for fast-relaxing, myelin water compartment
T1_s_ms = 832 #T1 for slow-relaxing, non-myelin water compartment
T1_ms_tuple = get_T1_tuple(T1_f_ms,T1_s_ms)

T2_f_ms = 20 # T2 for fast-relaxing, myelin water compartment
T2_s_ms = 80 #T2 for slow-relaxing, non-myelin water compartment
T2_ms_tuple = get_T2_tuple(T2_f_ms,T2_s_ms)

# additional off-resonance that is only experienced by myelin water
Δf_myelin_Hz = 5.0

# myelin water fraction (MWF), fast fraction
f_f = 0.15
# tuple with fast and slow relaxing fractions
mwf_tuple = get_mwf_tuple(f_f)

# TR and TE
TR_ms = 20.0
TE_ms = 4.0


# Generate plots similar to Figure 3 from [1]
# but with three different RF phase cycling factor values:
# (0, 90, and 180 degrees).

# for this example, choose one exchange rate
τ_fs = 50.0 # this will be varied in the next plot

# tuple with fast-to-slow and slow-to-fast residence times
τ_tuple_ms = get_τ_tuple(τ_fs, f_f)

# number of samples (resonant frequencies)
num_samples = 100

# flips angles for example
flip_ang_arr_deg = [10 40]
num_flip_angles = length(flip_ang_arr_deg)

# RF phase cycling value (degrees)
ΔΦ_arr_deg = [0 90 180]
Δϕ_arr_marker = [:circle :star5 :utriangle]
num_phases = length(ΔΦ_arr_deg)

# array to store result data
sig_arr = zeros(num_flip_angles,num_phases,num_samples)

# array with off-resonance values
Δf_arr_kHz = range(-1/TR_ms, 1/TR_ms, num_samples)

p = plot(title="Steady-State Signal Magnitude vs. Resonant Frequency", titlefontsize=12)

# iterate over flip angles
for i = 1:num_flip_angles
    α_deg = flip_ang_arr_deg[i]
    α_str = string(α_deg)

    # iterate over RF phases
    for j = 1:num_phases
        ΔΦ_deg = ΔΦ_arr_deg[j]
        ΔΦ_str = string(ΔΦ_deg)

        # iterate over resonant frequencies
        for k = 1:num_samples
            Δf_kHz = Δf_arr_kHz[k]

            # convert off-resonance from kHz to Hz before input into function
            Δf_Hz = kHz_to_Hz(Δf_kHz)

            # convert inputted RF phase cycling angle from degrees to radians
            ΔΦ_rad = deg_to_rad(ΔΦ_deg)

            # get tuple of values incorporating off-resonance and RF phase cycling for both compartments
            Δf_tuple_Hz = get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)
            Δf_tuple_Hz_no_rf_phase_fact = get_Δf_tuple(0, Δf_Hz, Δf_myelin_Hz, TR_ms)

            # create a spin (with and without RF phase-cycling factor)
            spin_mc = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz, τ_tuple_ms)
            spin_mc_no_rf_phase_fact = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz_no_rf_phase_fact, τ_tuple_ms)

            # run the bssfp blochsim and add to result array
            signal = bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)
            sig_arr[i,j,k] = abs(signal)
        end

        plot!(p, Δf_arr_kHz, sig_arr[i,j,:],label="α="*α_str*"°, ΔΦ="*ΔΦ_str*"°")

    end

end

xlabel!("Resonant Frequency (kHz)")
ylabel!("Signal Magnitude")
p

# Recreate Figure 2 from [2] (magnitude plot) and also add the phase plot.

# initialize the plot
lay = @layout [a ; b]
p_m = plot(title="Signal Magnitude vs. Scan Index", titlefontsize=10)
p_p = plot(title="Signal Phase vs. Scan Index", titlefontsize=10)

# number of different scans
num_scans = 40
scan_idx = range(1,num_scans,num_scans)

# flip angles for plot
flip_ang_arr_deg = [10.0 40.0]
num_flip_angles = length(flip_ang_arr_deg)

# set off-resonance to zero
Δf_Hz = 0.0

# array of exchange values
tau_arr_ms = [250, 150, 50]
tau_arr_marker = [:circle :star5 :utriangle]
num_taus = length(tau_arr_ms)

# initialize arrays to store results
sig_arr = zeros(num_scans,num_taus)
sig_arr_phase = zeros(num_scans,num_taus)

global curr_scan = 1

# iterate over exchange values
for j = 1:num_taus
    τ_fs = tau_arr_ms[j]
    τ_tuple_ms = get_τ_tuple(tau_arr_ms[j], f_f)
    tau_str = string(τ_fs)
    tau_marker = tau_arr_marker[j]

    # iterate over flip angles
    for k = 1:num_flip_angles
        α_deg = flip_ang_arr_deg[k]

        # different RF phases for different flip angles - from Figure 1 in [2]
        if k == 1
            ΔΦ_arr_deg = [-176.4, -159.5, -142.1, -124.4, -107.6, -90.54, -73.62, -56.13, -39.41, -22.52, -5.272, 11.63, 28.93, 45.76, 63.08, 79.91, 96.97, 113.9, 131.3, 148.5, 166.1]
            itr_max = length(ΔΦ_arr_deg)

        else
            ΔΦ_arr_deg = [-168.8, -150.3, -130.1, -111.5, -93.19, -74.18, -54.68 , -37.15, -18.01, 1.342, 18.82, 38.64, 57.88, 76.48, 95.2, 113.3, 133.3, 153.1, 172.1]
            itr_max = length(ΔΦ_arr_deg)

        end

        # iterate over RF phases
        for i = 1:itr_max
            ΔΦ_deg = ΔΦ_arr_deg[i]

            # convert RF phase cycling angle from degrees to radians
            ΔΦ_rad = deg_to_rad(ΔΦ_deg)

            # get tuple of values incorporating off-resonance and RF phase cycling for both compartments
            Δf_tuple_Hz = get_Δf_tuple(ΔΦ_rad, Δf_Hz, Δf_myelin_Hz, TR_ms)
            Δf_tuple_Hz_no_rf_phase_fact = get_Δf_tuple(0, Δf_Hz, Δf_myelin_Hz, TR_ms)

            # create a spin (with and without RF phase-cycling factor)
            spin_mc = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz, τ_tuple_ms)
            spin_mc_no_rf_phase_fact = SpinMC(mo, mwf_tuple, T1_ms_tuple, T2_ms_tuple, Δf_tuple_Hz_no_rf_phase_fact, τ_tuple_ms)

            # run the bssfp blochsim and add to result array
            signal = bssfp_blochsim_MC(α_deg, TR_ms, TE_ms, spin_mc, spin_mc_no_rf_phase_fact)
            sig_arr[curr_scan,j] = abs(signal)
            sig_arr_phase[curr_scan,j] = angle(signal)

            global curr_scan += 1
        end
    end

    global curr_scan = 1

    plot!(p_m, scan_idx,sig_arr[:,j], linewidth=0, markershape=tau_marker, label="τ_fs="*tau_str)
    plot!(p_p, scan_idx,sig_arr_phase[:,j], linewidth=0, markershape=tau_marker, label="τ_fs="*tau_str)
end

# plot results and label axes
p = plot(p_m,p_p, layout = lay)
plot!(size=(600,800))
xlabel!(p, "Scan Index")
ylabel!(p_m, "Signal Magnitude")
ylabel!(p_p, "Signal Phase")
p


#=
### References

- [1] Hargreaves, B. A., Vasanawala, S. S., Pauly, J. M., & Nishimura, D. G. (2001). Characterization and reduction of the transient response in steady‐state MR imaging. Magnetic Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine, 46(1), 149-158.

- [2] Murthy, N., Nielsen, J. F., Whitaker, S. T., Haskell, M. W., Swanson, S. D., Seiberlich, N., & Fessler, J. A. (2022). Quantifying myelin water exchange using optimized bSSFP sequences. In Proc. Intl. Soc. Mag. Res. Med (p. 2068).

- [3] https://github.com/StevenWhitaker/BlochSim.jl

- [4] Hinshaw, W. S. (1976). Image formation by nuclear magnetic resonance: the sensitive‐point method. Journal of Applied Physics, 47(8), 3709-3721.

- [5] Whitaker, S. T., Nataraj, G., Nielsen, J. F., & Fessler, J. A. (2020). Myelin water fraction estimation using small‐tip fast recovery MRI. Magnetic resonance in medicine, 84(4), 1977-1990.
=#



include("../../../inc/reproduce.jl")
