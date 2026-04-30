#=
bssfp1.jl
Balanced steady-state free precession (bSSFP) signal
for an isochromat (1-pool).
=#

export bssfp, bSSFPtuple1
export bSSFPbloch, bSSFPbloch3, bSSFPellipse

using BlochSim: Position, Spin, excite, freeprecess, duration
using BlochSim: AbstractRF, InstantaneousRF, RF1
using LinearAlgebra: I



"""
    bSSFPmode

A type used to control how bSSFP signal is calculated.
- `bSSFPbloch` use `BlochSim` matrix computations (default)
- `bSSFPbloch3` use analytical 3×3 matrix computations
- `bSSFPellipse` use ellipse model for 1-pool
"""
struct bSSFPmode{T} end
const bSSFPbloch = bSSFPmode{:Bloch}()
const bSSFPbloch3 = bSSFPmode{:Bloch3}()
const bSSFPellipse = bSSFPmode{:Ellipse}()


TE_ms_type = Any # Union{Number, Val{:midTR}, Val{:postRF}}

"""
    _TE_ms(TE_ms, TR_ms, [rf])

Helper to convert TE symbols to a number
- `Val{:midTR}` for the usual TR/2 choice
- `Val{:postRF}` immediately after RF pulse ends
"""
_TE_ms(TE_ms::Number, TR_ms::Number, args...) = TE_ms
_TE_ms(::Val{:midTR}, TR_ms::Number, args...) = TR_ms / 2
_TE_ms(::Val{:postRF}, TR_ms::Number, tRF_ms::Number) = tRF_ms / 2 # COV_EXCL_LINE
_TE_ms(::Val{:postRF}, TR_ms::Number, rf::AbstractRF) = duration(rf) / 2
_TE_ms(::Val{:postRF}, TR_ms::Number, args...) = 0 # COV_EXCL_LINE


"""
    bssfp(bSSFPellipse, Mz0, T1_ms, T2_ms, Δf_Hz,
       TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad=0)

Elliptical signal model for bSSFP.
This is the analytical solution to the 1-pool bSSFP signal
with instantaneous RF.

Xiang et al. MRM 2014;
https://doi.org/10.1002/mrm.25098

Keskin et al. IEEE T-MI 2022;
https://doi.org/10.1109/TMI.2021.3102852
"""
function bssfp(::bSSFPmode{:Ellipse},
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::TE_ms_type, Δϕ_rad::Number,
    α_rad::Number, θ_rf_rad::Number = 0,
)
    TE_ms = _TE_ms(TE_ms, TR_ms) # handle Val
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
    bssfp(spin, Δf_Hz, TR_ms, Val(:midTR) or Val(:postRF), Δϕ_rad, rf)
    `Val(:midTR)` for TE = TR/2;  Val(:postRF) for TE = tRF/2

Return steady-state magnetization signal value
at the echo time
for a single-pool tissue
for a bSSFP sequence
using `BlochSim`.

This generalizes
[Hargreaves et al., MRM 2001](https://doi.org/10.1002/mrm.1170)
by accounting for possibly finite RF duration.

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
bssfp(::bSSFPmode{:Bloch}, args...) = bssfp(args...)

function bssfp(
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::TE_ms_type, Δϕ_rad::Number,
    α_rad::Number, θ_rf_rad::Number = 0,
)
    rf = InstantaneousRF(α_rad, θ_rf_rad)
    return bssfp(Mz0, T1_ms, T2_ms, Δf_Hz, TR_ms, TE_ms, Δϕ_rad, rf)
end

function bssfp(
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::TE_ms_type, Δϕ_rad::Number,
    rf::AbstractRF,
    pos::Position = Position(0.0, 0.0, 0.0),
)
    TE_ms = _TE_ms(TE_ms, TR_ms, rf) # handle Val
    tRF_ms = duration(rf)
    tRF_ms/2 ≤ TE_ms < TR_ms - tRF_ms/2 ||
        throw("bad TE=$TE_ms for TR=$TR_ms and tRF=$tRF_ms")
    spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz, pos)
    return bssfp(spin, TR_ms, TE_ms, Δϕ_rad, rf)
end


"""
    function bssfp(spin, TR_ms, TE_ms, rf::AbstractRF)
Classic version with no phase cycling increment,
for `InstantaneousRF` only.
"""
function bssfp(spin::Spin, TR_ms::Number, TE_ms::TE_ms_type, rf::AbstractRF)

    TE_ms = _TE_ms(TE_ms, TR_ms, rf) # handle Val
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
Signal accounting for phase cycling increment `Δϕ_rad`,
allowing for finite duration `rf` pulse.
"""
function bssfp(spin::Spin,
    TR_ms::Number, TE_ms::TE_ms_type, Δϕ_rad::Number, rf::AbstractRF,
)

    TE_ms = _TE_ms(TE_ms, TR_ms, rf) # handle Val
    (A1, b1) = excite(spin, rf) # matrix for spin excitation

    tRF_ms = duration(rf)
    (A0, d0) = freeprecess(spin, TR_ms - tRF_ms)
    Rz = FreePrecessionMatrix(1, 1, -Δϕ_rad) # phase cycling
    b = isnothing(b1) ? Vector(A1 * d0) : Vector(A1 * d0 + b1)
    A = Matrix(A1 * A0 * Rz)
    Mss = (I - A) \ b # steady-state magnetization immediately after RF

    # account for free precession from end of RF to TE:
    t_free_ms = TE_ms - tRF_ms / 2
#   Efree = FreePrecessionMatrix(I, exp(-t_free_ms/spin.T2), spin.Δf_Hz)
    return complex(Mss[1], Mss[2]) * # complex signal
        exp(-t_free_ms / spin.T2) * # T2 decay
        cis(-2π/1000*spin.Δf*t_free_ms) # off resonance
end


"""
    bSSFPtuple1
Named tuple tissue parameter keys:
`(:Mz0, :T1_ms, :T2_ms, :Δf_Hz)`
"""
bSSFPtuple1 = (:Mz0, :T1_ms, :T2_ms, :Δf_Hz)

# helpers
bssfp(xt::NamedTuple{bSSFPtuple1}, args...) = bssfp(xt..., args...)
bssfp(::bSSFPmode{:Ellipse}, xt::NamedTuple{bSSFPtuple1}, args...) =
    bssfp(bSSFPellipse, xt..., args...)


"""
     bssfp(bSSFPbloch3, tRF_ms, args...)

1-pool version for a constant RF of duration `tRF_ms`.
(Account for gradient effects in `Δf_Hz`.)
"""
function bssfp(::bSSFPmode{:Bloch3},
    tRF_ms::Number, # rectangular waveform only for now
    Mz0::Number, T1_ms::Number, T2_ms::Number, Δf_Hz::Number,
    TR_ms::Number, TE_ms::TE_ms_type, Δϕ_rad::Number,
    α_rad::Number, θ_rf_rad::Number = 0,
)

    TE_ms = _TE_ms(TE_ms, TR_ms, tRF_ms) # handle Val
    spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz)
    rf = RF1(α_rad, tRF_ms, θ_rf_rad)
    (A1, b1) = excite_bloch3(spin, rf)
    (A0, d0) = freeprecess(spin, TR_ms - tRF_ms)

    Rz = FreePrecessionMatrix(1, 1, -Δϕ_rad) # phase cycling
    b = A1 * Vector(d0) + b1
    A = A1 * Matrix(A0) * Matrix(Rz)
    Mss = (I - A) \ b # steady-state magnetization immediately after RF

    # account for free precession from end of RF to TE:
    t_free_ms = TE_ms - tRF_ms / 2
    return complex(Mss[1], Mss[2]) * # complex signal
        exp(-t_free_ms / T2_ms) * # T2 decay
        cis(-2π/1000*Δf_Hz*t_free_ms) # off resonance
end
