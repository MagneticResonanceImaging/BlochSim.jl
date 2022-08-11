"""
    Gradient(x, y, z)

Create a `Gradient` object representing x, y, and z B0 gradients. Units are
G/cm.

# Properties
- `x::Real`: x gradient
- `y::Real`: y gradient
- `z::Real`: z gradient
"""
struct Gradient{T<:Real}
    x::T
    y::T
    z::T
end

Gradient(x, y, z) = Gradient(promote(x, y, z)...)

Base.show(io::IO, grad::Gradient) = print(io, "Gradient(", grad.x, ", ", grad.y, ", ", grad.z, ")")
Base.show(io::IO, ::MIME"text/plain", grad::Gradient{T}) where {T} =
    print(io, "Gradient{$T}:\n Gx = ", grad.x, " G/cm\n Gy = ", grad.y, " G/cm\n Gz = ", grad.z, " G/cm")

"""
    gradient_frequency(grad, pos)

Compute the off-resonance frequency in Hz induced by the given B0 gradient
`grad` at position `pos`.
"""
function gradient_frequency(grad::Gradient, pos::Position)

    return GAMBAR * (grad.x * pos.x + grad.y * pos.y + grad.z * pos.z) # Hz

end

"""
    AbstractSpoiling

Abstract type for representing spoiling.
"""
abstract type AbstractSpoiling end

"""
    IdealSpoiling() <: AbstractSpoiling

Represents ideal spoiling, i.e., setting the transverse (x and y) components of
a spin's magnetization to 0.
"""
struct IdealSpoiling <: AbstractSpoiling end

"""
    GradientSpoiling(grad, Tg) <: AbstractSpoiling
    GradientSpoiling(gx, gy, gz, Tg)

Represents gradient spoiling, e.g., applying a gradient
`grad = Gradient(gx, gy, gz)` for time `Tg` (ms).
`grad` can be a `Gradient`
(or `gx`, `gy`, and `gz` can be scalars),
representing a constant gradient,
or `grad` can be a collection of `Gradient`s
(or `gx`, `gy`, and `gz` can be collections of values),
representing a gradient waveform
with a constant time step.
"""
struct GradientSpoiling{T} <: AbstractSpoiling
    gradient::T
    Tg::Float64

    function GradientSpoiling(gradient::T, Tg::Real) where {T}

        T <: Gradient || eltype(T) <: Gradient ||
            error("gradient must be a Gradient or a collection of Gradients")
        new{T}(gradient, Tg)

    end
end

GradientSpoiling(gx::Real, gy::Real, gz::Real, Tg) = GradientSpoiling(Gradient(gx, gy, gz), Tg)
GradientSpoiling(gx, gy, gz, Tg) = GradientSpoiling(map(Gradient, gx, gy, gz), Tg)

Base.show(io::IO, s::GradientSpoiling) = print(io, "GradientSpoiling(", s.gradient, ", ", s.Tg, ")")

function Base.show(io::IO, ::MIME"text/plain", s::GradientSpoiling{T}) where {T}

    print(io, "GradientSpoiling{$T}:")
    print(io, "\n gradient = ", s.gradient, " G/cm")
    print(io, "\n Tg = ", s.Tg, " ms")

end

# I almost removed RFSpoiling but decided against it because it can be used to
# implement, e.g., phase cycling for bSSFP
"""
    RFSpoiling(Δθ = 117°) <: AbstractSpoiling

Represents RF spoiling, i.e., quadratically incrementing the phase of the RF
pulses from TR to TR.
"""
struct RFSpoiling{T<:Real} <: AbstractSpoiling
    Δθ::T
end

RFSpoiling() = RFSpoiling(deg2rad(117))

Base.show(io::IO, s::RFSpoiling) = print(io, "RFSpoiling(", s.Δθ, ")")
Base.show(io::IO, ::MIME"text/plain", s::RFSpoiling{T}) where {T} =
    print(io, "RFSpoiling{$T}:\n Δθ = ", s.Δθ, " rad")

"""
    RFandGradientSpoiling(grad_spoiling, rf_spoiling) <: AbstractSpoiling

Represents both RF and gradient spoiling.
"""
struct RFandGradientSpoiling{T1<:Real,T2<:Real} <: AbstractSpoiling
    gradient::GradientSpoiling{T1}
    rf::RFSpoiling{T2}
end

RFandGradientSpoiling(grad::Gradient, Tg::Real, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad, Tg), rf)
RFandGradientSpoiling(grad::NTuple{3,Real}, Tg::Real, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad..., Tg), rf)
RFandGradientSpoiling(gx, gy, gz, Tg, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(gx, gy, gz, Tg), rf)
RFandGradientSpoiling(grad::GradientSpoiling, Δθ) = RFandGradientSpoiling(grad, RFSpoiling(Δθ))
RFandGradientSpoiling(grad::Union{<:Gradient,<:NTuple{3,Real}}, Tg, Δθ) = RFandGradientSpoiling(grad, Tg, RFSpoiling(Δθ))
RFandGradientSpoiling(grad::GradientSpoiling) = RFandGradientSpoiling(grad, RFSpoiling())
RFandGradientSpoiling(grad::Union{<:Gradient,<:NTuple{3,Real}}, Tg) = RFandGradientSpoiling(grad, Tg, RFSpoiling())
RFandGradientSpoiling(rf::Union{<:RFSpoiling,<:Real}, grad::GradientSpoiling) = RFandGradientSpoiling(grad, rf)
RFandGradientSpoiling(rf::Union{<:RFSpoiling,<:Real}, grad::Union{<:Gradient,<:NTuple{3,Real}}, Tg) = RFandGradientSpoiling(grad, Tg, rf)
RFandGradientSpoiling(rf::RFSpoiling, gx, gy, gz, Tg) = RFandGradientSpoiling(gx, gy, gz, Tg, rf)

Base.show(io::IO, s::RFandGradientSpoiling) = print(io, "RFandGradientSpoiling(", s.gradient, ", ", s.rf, ")")

function Base.show(io::IO, ::MIME"text/plain", s::RFandGradientSpoiling{T1,T2}) where {T1,T2}

    print(io, "RFandGradientSpoiling{$T1,$T2}:")
    print(io, "\n gradient = GradientSpoiling{$T1}:")
    print(io, "\n  gradient = ", s.gradient.gradient, " G/cm")
    print(io, "\n  Tg = ", s.gradient.Tg, " ms")
    print(io, "\n rf = RFSpoiling{$T2}:")
    print(io, "\n  Δθ = ", s.rf.Δθ, " rad")

end

"""
    spoiler_gradient(spoiling)

Get the `Gradient` object used for gradient spoiling.
"""
spoiler_gradient(s::GradientSpoiling) = s.gradient
spoiler_gradient(s::RFandGradientSpoiling) = spoiler_gradient(s.gradient)

"""
    spoiler_gradient_duration(spoiling)

Return the duration of the spoiler gradient (ms).
"""
spoiler_gradient_duration(::AbstractSpoiling) = 0
spoiler_gradient_duration(s::GradientSpoiling) = s.Tg
spoiler_gradient_duration(s::RFandGradientSpoiling) = spoiler_gradient_duration(s.gradient)

"""
    rfspoiling_increment(spoiling)

Return the quadratic phase increment used for RF spoiling.
"""
rfspoiling_increment(::AbstractSpoiling) = 0
rfspoiling_increment(s::RFSpoiling) = s.Δθ
rfspoiling_increment(s::RFandGradientSpoiling) = rfspoiling_increment(s.rf)

"""
    spoil(spin, spoiling, [nothing])
    spoil(spinmc, spoiling, [workspace])

Simulate gradient or ideal spoiling for the given spin. Returns `(A, B)`, such
that `A * M + B` applies spoiling to the magnetization `M`. If `B` is `nothing`
(as is the case for `IdealSpoiling`), then `A * M` applies spoiling, and if both
`A` and `B` are `nothing` (as is the case for `RFSpoiling`) then there is no
spoiling.

For `SpinMC` objects and for `GradientSpoiling` and `RFandGradientSpoiling`,
`workspace isa BlochMcConnellWorkspace`. Pass in `nothing` instead to use an
approximate matrix exponential to solve the Bloch-McConnell equation.

## Note
This function only simulates gradient or ideal spoiling, *not* RF spoiling. RF
spoiling must be implemented by updating the phase of the RF pulse(s) in your
sequence from TR to TR.

For an in-place version, see [`spoil!`](@ref).
```
"""
spoil(::AbstractSpin, ::IdealSpoiling = IdealSpoiling(), ::Any = nothing) = (idealspoiling, nothing)
spoil(::AbstractSpin, ::RFSpoiling, ::Any = nothing) = (nothing, nothing)

function spoil(
    spin::AbstractSpin,
    spoiling::AbstractSpoiling,
    workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    return freeprecess(spin, spoiler_gradient_duration(spoiling), spoiler_gradient(spoiling), workspace)

end

"""
    spoil!(spin)

Apply ideal spoiling to the given spin.
"""
function spoil!(spin::Spin{T}) where {T}

    spin.M.x = zero(T)
    spin.M.y = zero(T)
    return nothing

end

function spoil!(spin::SpinMC{T}) where {T}

    for Mc in spin.M
        Mc.x = zero(T)
        Mc.y = zero(T)
    end

end

applydynamics!(spin::AbstractSpin, ::IdealSpoilingMatrix) = spoil!(spin)

"""
    spoil!(A, B, spin, spoiling, [nothing])
    spoil!(A, B, spinmc, spoiling, [workspace])

Simulate gradient or ideal spoiling, overwriting `A` and `B` (in-place version
of [`spoil`](@ref)).
"""
spoil!(A::IdealSpoilingMatrix, ::Nothing, ::AbstractSpin, ::IdealSpoiling, ::Any = nothing) = nothing
spoil!(::Nothing, ::Nothing, ::AbstractSpin, ::RFSpoiling, ::Any = nothing) = nothing

function spoil!(
    A,
    B,
    spin::AbstractSpin,
    spoiling::AbstractSpoiling,
    workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    freeprecess!(A, B, spin, spoiler_gradient_duration(spoiling), spoiler_gradient(spoiling), workspace)

end
