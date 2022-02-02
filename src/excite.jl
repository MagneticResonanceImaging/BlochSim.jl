"""
    AbstractRF

Abstract type for representing radiofrequency (RF) pulses.
"""
abstract type AbstractRF end

"""
    InstantaneousRF(α, θ = 0) <: AbstractRF

Represents an idealized instantaneous RF pulse with flip angle `α` and phase
`θ`.
"""
struct InstantaneousRF{T<:Real} <: AbstractRF
    α::T
    θ::T
end

InstantaneousRF(α, θ = zero(α)) = InstantaneousRF(promote(α, θ)...)

Base.show(io::IO, rf::InstantaneousRF) = print(io, "InstantaneousRF(", rf.α, ", ", rf.θ, ")")
Base.show(io::IO, ::MIME"text/plain", rf::InstantaneousRF{T}) where {T} =
    print(io, "Instantaneous RF pulse with eltype $T:\n α = ", rf.α, " rad\n θ = ", rf.θ, " rad")

"""
    RF(waveform, Δt, [Δθ], [grad]) <: AbstractRF

Represents an RF pulse with the given (possibly complex-valued) waveform (G) and
time step `Δt` (ms). `Δθ` is additional phase added to the waveform (defaults to
`0`), and `grad` is the B0 gradient that is turned on during the RF pulse
(defaults to `Gradient(0, 0, 0)`, i.e., turned off).

# Properties
- `α::Vector{<:Real}`: Instantaneous flip angles (rad) at each time point;
  computed from the magnitude of `waveform`
- `θ::Vector{<:Real}`: Instantaneous phase (rad) at each time point; computed
  from the phase of `waveform`
- `Δt::Real`: Time step (ms)
- `Δθ_initial::Real`: Phase added to `θ` before any phase-cycling increment has
  been applied
- `Δθ::Ref{<:Real}`: Phase to be added to `θ`; can be updated to simulate
  phase-cycling/RF spoiling
- `grad`: Gradient applied during the RF pulse
  - `::Gradient`: Constant gradient
  - `::Vector{<:Gradient}`: Time-varying gradient
"""
struct RF{T<:Real,G<:Union{<:Gradient,<:AbstractVector{<:Gradient}}} <: AbstractRF
    α::Vector{T}
    θ::Vector{T}
    Δt::Float64
    Δθ_initial::T
    Δθ::Ref{T} # Type Ref to enable RF-spoiling, which requires updating Δθ
    grad::G

    function RF(
        α::AbstractVector{<:Real},
        θ::AbstractVector{<:Real},
        Δt::Real,
        Δθ::Real,
        grad::Union{<:Gradient,<:AbstractVector{<:Gradient}}
    )

        length(α) == length(θ) || error("α and θ must have the same number of elements")
        grad isa AbstractVector && (length(grad) == length(α) ||
            error("grad is a vector but has a different number of elements than α"))
        T = promote_type(eltype(α), eltype(θ), typeof(Δθ))
        new{T,typeof(grad)}(α, θ, Δt, Δθ, Δθ, grad)

    end
end

# waveform in Gauss, Δt in ms
RF(waveform, Δt, Δθ, grad) = RF(GAMMA .* abs.(waveform) .* (Δt / 1000), angle.(waveform), Δt, Δθ, grad)
RF(waveform, Δt, Δθ::Real) = RF(waveform, Δt, Δθ, Gradient(0, 0, 0))
RF(waveform, Δt, grad) = RF(waveform, Δt, 0, grad)
RF(waveform, Δt) = RF(waveform, Δt, 0, Gradient(0, 0, 0))

Base.show(io::IO, rf::RF) = print(io, "RF(", rf.α, ", ", rf.θ, ", ", rf.Δt, ", ", rf.Δθ_initial, ", ", rf.grad, ")")

function Base.show(io::IO, ::MIME"text/plain", rf::RF{T,G}) where {T,G}

    print(io, "RF{$T,$G}:")
    print(io, "\n α = ", rf.α, " rad")
    print(io, "\n θ = ", rf.θ, " rad")
    print(io, "\n Δt = ", rf.Δt, " ms")
    print(io, "\n Δθ (initial) = ", rf.Δθ_initial, " rad")
    print(io, "\n Δθ (current) = ", rf.Δθ[], " rad")
    print(io, "\n grad = ", rf.grad, " G/cm")

end

Base.length(rf::RF) = length(rf.α)

"""
    duration(rf)

Return the duration (ms) of the RF pulse.
"""
duration(::InstantaneousRF) = 0
duration(rf::RF) = length(rf) * rf.Δt

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

"""
    rotatetheta!(A, θ, α)

Simulate left-handed rotation by angle `α` about an axis in the x-y plane that
makes angle `θ` with the negative y-axis, overwriting `A`.
"""
function rotatetheta!(A, θ, α)

    (sinθ, cosθ) = sincos(θ)
    (sinα, cosα) = sincos(α)

    A.a11 = sinθ^2 + cosα * cosθ^2
    A.a21 = sinθ * cosθ - cosα * sinθ * cosθ
    A.a31 = -sinα * cosθ
    A.a12 = sinθ * cosθ - cosα * sinθ * cosθ
    A.a22 = cosθ^2 + cosα * sinθ^2
    A.a32 = sinα * sinθ
    A.a13 = sinα * cosθ
    A.a23 = -sinα * sinθ
    A.a33 = cosα

    return nothing

end

"""
    excite(spin, rf::InstantaneousRF)
    excite(spin, rf::RF)

Simulate excitation for the given spin. Returns `(A, B)` such that `A * M + B`
applies excitation to the magnetization `M`. If `isnothing(B)` (as is the case
for `InstantaneousRF`s), then `A * M` applies excitation to `M`.

For an in-place version, see [`excite!`](@ref).

# Examples
```jldoctest
julia> s = Spin(1, 1000, 100, 3.75); s.M
Magnetization vector with eltype Float64:
 Mx = 0.0
 My = 0.0
 Mz = 1.0

julia> (A,) = excite(s, InstantaneousRF(π/2, π/4)); A * s.M
Magnetization vector with eltype Float64:
 Mx = 0.7071067811865476
 My = -0.7071067811865475
 Mz = 6.123233995736766e-17
```
"""
function excite(spin::AbstractSpin, rf::InstantaneousRF)

    (sinθ, cosθ) = sincos(rf.θ)
    (sinα, cosα) = sincos(rf.α)

    A = ExcitationMatrix(BlochMatrix(
        sinθ^2 + cosα * cosθ^2,
        sinθ * cosθ - cosα * sinθ * cosθ,
        -sinα * cosθ,
        sinθ * cosθ - cosα * sinθ * cosθ,
        cosθ^2 + cosα * sinθ^2,
        sinα * sinθ,
        sinα * cosθ,
        -sinα * sinθ,
        cosα
    ))

    return (A, nothing)

end

"""
    excite!(A, [nothing], spin, rf::InstantaneousRF, [nothing])
    excite!(A, B, spin, rf::RF, [workspace])

Simulate excitation, overwriting `A` and `B` (in-place version of
[`excite`](@ref)).

For `RF` objects, `workspace isa ExcitationWorkspace`. For `SpinMC` objects, use
`workspace = ExcitationWorkspace(spin, nothing)` to use an approximate matrix
exponential to solve the Bloch-McConnell equation.
"""
function excite!(A::ExcitationMatrix, spin::AbstractSpin, rf::InstantaneousRF)

    rotatetheta!(A.A, rf.θ, rf.α)

end

function excite!(A::ExcitationMatrix, ::Nothing, spin::AbstractSpin, rf::InstantaneousRF, ::Nothing = nothing)

    excite!(A, spin, rf)

end

function excite(spin::AbstractSpin, rf::RF{T,G}) where {T,G}

    Δtby2 = rf.Δt / 2
    if G <: AbstractVector
        (Af, Bf) = freeprecess(spin, Δtby2, rf.grad[1])
    else
        (Af, Bf) = freeprecess(spin, Δtby2, rf.grad)
    end
    (Ae,) = excite(spin, InstantaneousRF(rf.α[1], rf.θ[1] + rf.Δθ[]))

    A = Af * Ae * Af
    B = Af * (Ae * Bf) + Bf

    for t = 2:length(rf)

        if G <: AbstractVector
            (Af, Bf) = freeprecess(spin, Δtby2, rf.grad[t])
        end
        (Ae,) = excite(spin, InstantaneousRF(rf.α[t], rf.θ[t] + rf.Δθ[]))

        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf

    end

    return (A, B)

end

function excite!(
    A::Union{<:BlochMatrix,<:BlochMcConnellMatrix},
    B::Union{<:Magnetization,<:MagnetizationMC},
    spin::AbstractSpin,
    rf::RF{T,G},
    workspace::ExcitationWorkspace = ExcitationWorkspace(spin)
) where {T,G}

    Δtby2 = rf.Δt / 2
    if G <: AbstractVector
        freeprecess!(workspace.Af, workspace.Bf, spin, Δtby2, rf.grad[1], workspace.bm_workspace)
    else
        freeprecess!(workspace.Af, workspace.Bf, spin, Δtby2, rf.grad, workspace.bm_workspace)
    end
    excite!(workspace.Ae, spin, InstantaneousRF(rf.α[1], rf.θ[1] + rf.Δθ[]))

    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Af, workspace.Bf, workspace.Ae, nothing)
    combine!(A, B, workspace.tmpA1, workspace.tmpB1, workspace.Af, workspace.Bf)

    for t = 2:length(rf)

        if G <: AbstractVector
            freeprecess!(workspace.Af, workspace.Bf, spin, Δtby2, rf.grad[t], workspace.bm_workspace)
        end
        excite!(workspace.Ae, spin, InstantaneousRF(rf.α[t], rf.θ[t] + rf.Δθ[]))

        combine!(workspace.tmpA1, workspace.tmpB1, A, B, workspace.Af, workspace.Bf)
        combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Ae, nothing)
        combine!(A, B, workspace.tmpA2, workspace.tmpB2, workspace.Af, workspace.Bf)

    end

end

"""
    excite!(spin, ...)

Apply excitation to the given spin, overwriting the spin's magnetization vector.
"""
function excite!(spin::AbstractSpin, args...)

    (A, B) = excite(spin, args...)
    BtoM = copy(spin.M)
    applydynamics!(spin, BtoM, A, B)

end
