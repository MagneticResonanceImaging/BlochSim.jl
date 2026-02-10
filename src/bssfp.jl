#=
bssfp.jl
Balanced steady-state free precession (bSSFP) signal
for an isochromat (1-pool) and for a 2-pool model with exchange.
=#

export bssfp, bssfp_ellipse, BSSFPTuple1

using BlochSim: Position, Spin, SpinMC, excite, freeprecess, duration
using BlochSim: AbstractRF, InstantaneousRF
using LinearAlgebra: I


"""
    bssfp_ellipse(Mz0, T1_ms, T2_ms, Δf_Hz,
       TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad=0)

Elliptical signal model for bSSFP.
This is the analytical solution to the 1-pool bSSFP signal.

Xiang et al. MRM 2014;
https://doi.org/10.1002/mrm.25098

Keskin et al. IEEE T-MI 2022;
https://doi.org/10.1109/TMI.2021.3102852
"""
function bssfp_ellipse(
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::Number, Δϕ_rad::Number,
    α_rad::Number, θ_rf_rad::Number = 0,
)
    s, c = sincos(α_rad)
    E1 = exp(-TR_ms / T1_ms)
    E2 = exp(-TR_ms / T2_ms)
    denom = 1 - E1 * c - E2^2 * (E1 - c)
    M = Mz0 * (1 - E1) * s / denom
    a = E2
    b = E2 * (1 - E1) * (1 + c) / denom
    Δω_kHz = 2π * Δf_Hz / 1000
    θ = Δω_kHz * TR_ms - Δϕ_rad # per eqn. (2) of keskin-22-cef
    post_rf = M * (1 - a * cis(θ)) / (1 - b * cos(θ)) # [6] of xiang-14-bar
    signal = post_rf * exp(-TE_ms / T2_ms) * cis(-Δω_kHz * TE_ms) # echo
    return cis(-θ_rf_rad) * signal # RF phase
end


"""
    bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad=0)
    bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, rf, [pos])
    bssfp(spin, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, rf)

Return steady-state magnetization signal value
at the echo time
for a bSSFP sequence
using `BlochSim`.

This generalizes
[Hargreaves et al., MRM 2001](https://doi.org/10.1002/mrm.1170)
by accounting for finite RF duration.

# In (tissue parameters):
- `Mz0` initial condition for magnetization in the z-direction (constant)
- `T1_ms` MRI tissue parameter for T1 relaxation (ms)
- `T2_ms` MRI tissue parameter for T2 relaxation (ms)
- `Δf_Hz` off-resonance value (Hz)
# In (scan parameters):
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms), measured from *middle* of RF pulse
- `Δϕ_rad` RF phase cycling increment (radians)
- `α_rad` flip angle of RF pulse (radians)
- `θ_rf_rad` phase of RF pulse (radians)
Or, instead of `α_rad` and `θ_rf_rad`, provide:
- `rf::AbstractRF`, e.g., `InstantaneousRF(α_rad, θ_rf_rad)`

# Option:
- `pos::Position = Position(0.0, 0.0, 0.0)` option if `rf` specified

# Out
- `signal` steady-state transverse magnetization (as a complex number)
"""
function bssfp(
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::Number, Δϕ_rad::Number,
    α_rad::Number, θ_rf_rad::Number = 0,
)
    rf = InstantaneousRF(α_rad, θ_rf_rad)
    return bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, rf)
end

function bssfp(
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::Number, Δϕ_rad::Number,
    rf::AbstractRF,
    pos::Position = Position(0.0, 0.0, 0.0),
)
    tRF_ms = duration(rf)
    tRF_ms/2 < TE_ms < TR_ms - tRF_ms/2 ||
        throw("bad TE=$TE_ms for TR=$TR_ms and tRF=$tRF_ms")
    spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz, pos)
    return bssfp(spin, TR_ms, TE_ms, Δϕ_rad, rf)
end


"""
    function bssfp(spin, TR_ms, TE_ms, rf::AbstractRF)
Classic version with no phase cycling increment,
for InstantaneousRF only.
"""
function bssfp(spin, TR_ms::Number, TE_ms::Number, rf::AbstractRF)

    (R,) = excite(spin, rf) # matrix for spin excitation
    rf isa InstantaneousRF || throw("unsupported")

    #=
    Matrices for precession/relaxation for various time period values
    The bSSFP pulse sequence in Figure 2 in [1] starts with a RF pulse, then
    - `a` is at time TE
    - `b` is TR-TE later, right before next RF pulse
    - `c` is immediately after the next RF pulse
    - `d` is TE after that next RF pulse
    =#
    (PC1_A, PC1_B) = freeprecess(spin, TE_ms)
    (PC2_A, PC2_B) = freeprecess(spin, TR_ms - TE_ms)
    (PC_TR_A, PC_TR_B) = freeprecess(spin, TR_ms)

    # calculate the A matrix and b vector
    A = Matrix(PC1_A * R.A * PC2_A)
    b = Matrix(PC1_A * R.A) * Vector(PC2_B) + Vector(PC1_B)

    Mss = (I - A) \ b # steady-state magnetization at the echo time
    return complex(Mss[1], Mss[2]) # complex signal
end


"""
    function bssfp(spin, TR_ms, TE_ms, Δϕ_rad, rf::AbstractRF)
Signal accounting for  phase cycling increment `Δϕ_rad`,
allowing for finite duration `rf` pulse.
"""
function bssfp(spin,
    TR_ms::Number, TE_ms::Number, Δϕ_rad::Number, rf::AbstractRF,
)

    (A0, d0) = excite(spin, rf) # matrix for spin excitation

    tRF_ms = duration(rf)
    (A1, d1) = freeprecess(spin, TR_ms - tRF_ms)
    Rz = FreePrecessionMatrix(1, 1, -Δϕ_rad) # phase cycling
    b = isnothing(d0) ? Vector(A0 * d1) : Vector(A0 * d1 + d0)
    A = Matrix(A0 * A1 * Rz)
    Mss = (I - A) \ b # steady-state magnetization immediately after RF

    # account for free precession from end of RF to TE:
    t_free_ms = TE_ms - tRF_ms / 2
#   Efree = FreePrecessionMatrix(I, exp(-t_free_ms/spin.T2), spin.Δf_Hz)
    return complex(Mss[1], Mss[2]) * # complex signal
        exp(-t_free_ms / spin.T2) * # T2 decay
        cis(-2π/1000*spin.Δf*t_free_ms) # off resonance
end


BSSFPTuple1 = (:Mz0, :T1_ms, :T2_ms, :Δf_Hz)

bssfp(xt::NamedTuple{BSSFPTuple1}, args...) = bssfp(xt..., args...)


#=
Multi-compartment spins and myelin water exchange
These helper functions put the parameters in the correct format
for the multi-compartment spin object constructors.
=#

"""
    get_τ_tuple(τ_fs_ms, f_f)
# In:
- `τ_fs_ms` residence time for exchange from myelin to non-myelin water (ms)
- `f_f` fast fraction (myelin fraction)
# Out:
- `τ_tuple_ms` tuple with fast-to-slow and slow-to-fast residence times
"""
function get_τ_tuple(τ_fs_ms, f_f)
    τ_sf_ms = (1-f_f) * τ_fs_ms / f_f
    return (τ_fs_ms, τ_sf_ms)
end


"""
    get_Δf_tuple(Δϕ_rad, Δf0_Hz, Δf_myelin_Hz, TR_ms)

Account for phase cycling increments as an effective frequency shift

# In:
- `Δϕ_rad` RF phase cycling value (radians)
- `Δf0_Hz` off-resonance value (Hz)
- `Δf_myelin_Hz` # additional off-resonance value only experienced by myelin water (Hz)
- `TR_ms` repetition time (ms)

# Out:
- `Δf_tuple_Hz` tuple with off-resonance values for fast and slow compartments

[Hinshaw, J. Appl. Phys. 1976](https://doi.org/10.1063/1.323136).
"""
function get_Δf_tuple(Δϕ_rad, Δf0_Hz, Δf_myelin_Hz, TR_ms)

    # convert the RF phase cycling value to Hz from radians
    Δϕ_Hz = 1000 * (Δϕ_rad) / (2π*TR_ms)

    # subtract the RF phase cycling value from the off-resonance value
    Δf_RF_Hz = Δf0_Hz - Δϕ_Hz

    # add the myelin off-resonance for the myelin term
    Δf_myelin_RF_Hz = Δf_RF_Hz + Δf_myelin_Hz

    # create and return tuple for the spin object constructor
    Δf_tuple_Hz = (Δf_myelin_RF_Hz, Δf_RF_Hz)
    return Δf_tuple_Hz
end


"""
    bssfp(spin, spin_no_rf_phase_fact, TR_ms, TE_ms, α_rad, θ_rf_rad=0)
    bssfp(spin, spin_no_rf_phase_fact, TR_ms, TE_ms, rf)

Return steady-state magnetization signal value
at the echo time
for a phase-cycled bSSFP sequence
using BlochSim.

Ref: Murthy, N., Nielsen, J. F., Whitaker, S. T., Haskell, M. W.,
Swanson, S. D., Seiberlich, N., & Fessler, J. A. (2022).
Quantifying myelin water exchange using optimized bSSFP
sequences. In Proc. Intl. Soc. Mag. Res. Med (p. 2068). [2]

# In (tissue)
- 'spin' multi-compartment spin object with RF phase cycling factor
- 'spin_no_rf_phase_fact' multi-compartment spin object without RF phase cycling factor
# In (scan)
- `TR_ms` repetition time (ms)
- `TE_ms` echo time (ms)
- `α_rad` flip angle of RF pulse (radians)
- `θ_rf_rad` phase angle of RF pulse (radians) [default 0]
Or instead provide
- `rf::AbstractRF`

# Out
- `signal` steady-state magnetization (as a complex number)
"""
function bssfp(
    spin::SpinMC, spin_no_rf_phase_fact::SpinMC,
    TR_ms::Number, TE_ms::Number, α_rad::Number, θ_rf_rad::Number=0.,
)
    rf = InstantaneousRF(α_rad, θ_rf_rad)
    return bssfp(spin, spin_no_rf_phase_fact, TR_ms, TE_ms, rf)
end


function bssfp(
    spin::SpinMC, spin_no_rf_phase_fact::SpinMC,
    TR_ms::Number, TE_ms::Number, rf::AbstractRF,
)

    rf isa InstantaneousRF ||
        throw("todo: relaxation and duraction effects?  need to derive for MC!")

    # excite the spin and reshape R to be the correct dimensions for a SpinMC object
    (R,) = excite(spin, rf)
    R = Matrix(R.A)
    R = kron(I(2), R)

    # precession/relaxation of the spin for TR
    (PC_TR_A, PC_TR_B) = freeprecess(spin, TR_ms)

    # precession/relaxation of the spin for TE
    # assume receiver modulates signal and uses the receiver phase as the RF phase
    (PC_TE_A, PE_TE_B) = freeprecess(spin_no_rf_phase_fact, TE_ms)

    # calculate A matrix and b vector
    A = Matrix(PC_TR_A) * R
    b = Vector(PC_TR_B)

    Mss = (I - A) \ b # steady-state just before tip down
    M = R * Mss # magnetization after tip-down

    # steady-state magnetization at the echo time
    M = Matrix(PC_TE_A) * M + Vector(PE_TE_B)

    return complex(M[1]+M[4], M[2]+M[5]) # return the complex signal
end


# Intermediate helper
function bssfp(
    Mz0::Number,
    frac::NTuple{2, Number},
    T1_ms::NTuple{2, Number},
    T2_ms::NTuple{2, Number},
    Δf_Hz::NTuple{2, Number},
    Δf_no_rf_phase_Hz::NTuple{2, Number},
    τ_ms::NTuple{2, Number},
    scan_args...,
)

    # create spin (with and without RF phase-cycling factor)
    spin = SpinMC(Mz0, frac, T1_ms, T2_ms, Δf_Hz, τ_ms)
    spin_no_rf_phase = SpinMC(Mz0, frac, T1_ms, T2_ms, Δf_no_rf_phase_Hz, τ_ms)

    return bssfp(spin, spin_no_rf_phase, scan_args...)
end


"""
    bssfp(...)
Version with scalar arguments (convenient for autodiff)
"""
function bssfp(
    M0_phase::Number, # radians
    Mz0::Number,
    f_f::Number,
    T1_f_ms::Number,
    T1_s_ms::Number,
    T2_f_ms::Number,
    T2_s_ms::Number,
    τ_fs_ms::Number,
    Δff_Hz::Number, # fast component frequency shift
    # system parameter (sometimes known):
    Δf0_Hz::Number, # B0 off resonance
    # scan parameters (always known):
    Δϕ_rad::Number, # RF phase cycling value (radians)
    TR_ms::Number,
    args..., # remaining scan arguments
)

    τ_tuple_ms = get_τ_tuple(τ_fs_ms, f_f) # fast-to-slow and slow-to-fast residence times

    # tuple of values incorporating off-resonance and RF phase cycling for both compartments
    Δf_tuple_Hz = get_Δf_tuple(Δϕ_rad, Δf0_Hz, Δff_Hz, TR_ms)
    Δf_tuple_Hz_no_rf_phase = get_Δf_tuple(0, Δf0_Hz, Δff_Hz, TR_ms)

    return cis(M0_phase) * bssfp(
        Mz0,
        (f_f, 1 - f_f), # fast, slow fractions
        (T1_f_ms, T1_s_ms),
        (T2_f_ms, T2_s_ms),
        Δf_tuple_Hz,
        Δf_tuple_Hz_no_rf_phase,
        τ_tuple_ms,
        TR_ms,
        args...,
    )
end
