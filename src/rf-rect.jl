# rf-rect

using BlochSim: AbstractRF, GAMMA, Gradient
export RectRF



# Helpers for making rectangular RF pulses for use with analytical Bloch solvers

"""
    b1 = b1_gauss(α_rad, tRF_ms)

Return finite-duration (rectangular) RF pulse amplitude
# In
- `α_rad` tip angle (radians)
- `tRF_ms` pulse length (ms)

# Notes:
- `GAMMA` has units rad/s/G
- Tip angle for constant pulse:
  `α_rad = GAMMA * b1_gauss * tRF_s`
- so `b1_gauss = α_rad / GAMMA / tRF_s`
"""
b1_gauss(α_rad, tRF_ms) = α_rad / GAMMA / (tRF_ms / 1000)


"""
    rf = RF1(α_rad, tRF_ms, θ = 0)
RF "rectangular" pulse
of duration `tRF_ms`
for flip angle `α_rad`
and phase `θ` (in radians),
represented by a single-sample "waveform".
"""
function RF1(α_rad, tRF_ms, θ_rad = 0, args...; kwargs...)
    waveform = [cis(θ_rad)] * b1_gauss(α_rad, tRF_ms) # single sample "waveform"
    return RF(waveform, tRF_ms, args...; kwargs...)
end


"""
    RectRF(duration_ms, [α = π/2], [θ = 0], [grad = zeros(3)]) <: AbstractRF

Represent a rectangular RF pulse
of duration `duration_ms`
for flip angle `α`
and phase `θ`
and constant gradient 3-vector `grad`.

- `grad` B0 gradient that is turned on during the RF pulse
  (defaults to `Gradient(0, 0, 0)`, i.e., turned off).
- `grad`: Gradient applied during the RF pulse
  - `::Gradient`: Constant gradient
  - `::Vector{<:Gradient}`: Time-varying gradient
"""
struct RectRF{Ta <: Real, Td <: Real, G <: Gradient} <: AbstractRF
    duration::Td
    α::Ta
    θ::Ta
    grad::G

    function RectRF(
        duration::Td,
        α::Real = π/2,
        θ::Real = zero(α),
        grad::Gradient = Gradient(0, 0, 0),
    ) where {Td <: Real}

        grad == Gradient(0, 0, 0) || throw("todo: grad unsupported")
        Ta = promote_type(eltype(α), eltype(θ))
        new{Ta, Td, typeof(grad)}(duration, α, θ, grad)
    end
end


# helpers
Base.show(io::IO, rf::RectRF) =
    print(io, "RF(", rf.duration, rf.α, ", ", rf.θ, ", ", rf.grad, ")")

function Base.show(io::IO, ::MIME"text/plain", rf::RectRF{Ta,Td,G}) where {Ta,Td,G}

    print(io, "RectRF{$Ta,$Td,$G}:")
    print(io, "\n duration = ", rf.duration, " ms")
    print(io, "\n α = ", rf.α, " rad")
    print(io, "\n θ = ", rf.θ, " rad")
    print(io, "\n grad = ", rf.grad, " G/cm")

end

Base.length(rf::RectRF) = 1


"""
    duration(rf)

Return the duration (ms) of the RF pulse.
"""
duration(rf::RectRF) = rf.duration ## COV_EXCL_LINE


#=
"""
    ExcitationWorkspace
Struct for in-place excitation operations.
"""
struct ExcitationWorkspace{T1,T2,T3,T4,T5}
    Af::T1
    Bf::T2
    Ae::T3
    tmpA1::T4
    tmpA2::T4
    tmpB1::T2
    tmpB2::T2
    bm_workspace::T5
end

function ExcitationWorkspace(
    spin::AbstractSpin,
    bm_workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    ExcitationWorkspace(typeof(spin), bm_workspace)

end

function ExcitationWorkspace(
    spin::Union{Type{Spin{T}},Type{SpinMC{T,N}}},
    bm_workspace = spin <: Spin ? nothing : BlochMcConnellWorkspace(spin)
) where {T,N}

    Ae = ExcitationMatrix{T}()
    if spin <: Spin
        Af = FreePrecessionMatrix{T}()
        Bf = Magnetization{T}()
        tmpA1 = BlochMatrix{T}()
        tmpA2 = BlochMatrix{T}()
        tmpB1 = Magnetization{T}()
        tmpB2 = Magnetization{T}()
    else
        Af = BlochMcConnellMatrix{T}(N)
        Bf = MagnetizationMC{T}(N)
        tmpA1 = BlochMcConnellMatrix{T}(N)
        tmpA2 = BlochMcConnellMatrix{T}(N)
        tmpB1 = MagnetizationMC{T}(N)
        tmpB2 = MagnetizationMC{T}(N)
    end
    ExcitationWorkspace(Af, Bf, Ae, tmpA1, tmpA2, tmpB1, tmpB2, bm_workspace)

end

#todo    excite(spin, rf::RectRF, [workspace])


"""
    excite!(A, B, spin, rf::RectRF, [workspace])
"""
function excite(spin::Spin, rf::RectRF,
    workspace = RectRFExcitationWorkspace(spin),
)

    T = eltype(spin)
    A = BlochMatrix{T}()
    B = Magnetization{T}()
    excite!(A, B, spin, rf, workspace)
    return (A, B)
end


"""
    excite!(A, B, spin, rf::RectRF, [workspace])
"""
function excite!(
    A::Union{<:BlochMatrix,<:BlochMcConnellMatrix},
    B::Union{<:Magnetization,<:MagnetizationMC},
    spin::AbstractSpin,
    rf::RectRF{T,G},
    workspace::ExcitationWorkspace = ExcitationWorkspace(spin)
) where {T,G}

    w = workspace # shorthand
    freeprecess!(w.Af, w.Bf, spin, Δtby2, rf.grad, w.bm_workspace)
    excite!(w.Ae, spin, InstantaneousRF(rf.α[1], rf.θ[1] + rf.Δθ[])) # R₁

    combine!(w.tmpA1, w.tmpB1, w.Af, w.Bf, w.Ae, nothing) # F₁^½ R₁
    combine!(A, B, w.tmpA1, w.tmpB1, w.Af, w.Bf) # F₁^½ R₁ F₁^½
    return (A, B)
end
=#
