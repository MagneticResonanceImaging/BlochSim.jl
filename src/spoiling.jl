struct Gradient{T<:Real}
    x::T
    y::T
    z::T
end

Gradient(x, y, z) = Gradient(promote(x, y, z)...)

Base.show(io::IO, grad::Gradient) = print(io, "[", grad.x, ", ", grad.y, ", ", grad.z, "]")
Base.show(io::IO, ::MIME"text/plain", grad::Gradient{T}) where {T} =
    print(io, "Gradient{$T}:\n Gx = ", grad.x, " G/cm\n Gy = ", grad.y, " G/cm\n Gz = ", grad.z, " G/cm")

function gradient_frequency(grad::Gradient, pos::Position)

    return GAMBAR * (grad.x * pos.x + grad.y * pos.y + grad.z * pos.z) # Hz

end

abstract type AbstractSpoiling end

struct IdealSpoiling <: AbstractSpoiling end

struct GradientSpoiling{T<:Real} <: AbstractSpoiling
    gradient::Gradient{T}
    Tg::Float64

    function GradientSpoiling(gradient::Gradient{T}, Tg::Real) where {T}

        new{T}(gradient, Tg)

    end
end

GradientSpoiling(x, y, z, Tg) = GradientSpoiling(Gradient(x, y, z), Tg)

struct RFSpoiling{T<:Real} <: AbstractSpoiling
    Δθ::T
end

RFSpoiling() = RFSpoiling(deg2rad(117))

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

spoiler_gradient(s::GradientSpoiling) = s.gradient
spoiler_gradient(s::RFandGradientSpoiling) = spoiler_gradient(s.gradient)
spoiler_gradient_duration(::AbstractSpoiling) = 0
spoiler_gradient_duration(s::GradientSpoiling) = s.Tg
spoiler_gradient_duration(s::RFandGradientSpoiling) = spoiler_gradient_duration(s.gradient)
rfspoiling_increment(::AbstractSpoiling) = 0
rfspoiling_increment(s::RFSpoiling) = s.Δθ
rfspoiling_increment(s::RFandGradientSpoiling) = rfspoiling_increment(s.rf)

"""
    spoil(spin)

Simulate ideal spoiling (i.e., setting the transverse component of the spin's
magnetization to 0).

# Arguments
- `spin::AbstractSpin`: Spin to spoil

# Return
- `S::Matrix`: Matrix that describes ideal spoiling

# Examples
```jldoctest
julia> spin = Spin([1, 0.4, 5], 1, 1000, 100, 0)
Spin([1.0, 0.4, 5.0], 1.0, 1000.0, 100.0, 0.0, [0.0, 0.0, 0.0])

julia> S = spoil(spin); S * spin.M
3-element Array{Float64,1}:
 0.0
 0.0
 5.0
```
"""
spoil(::Any) = [0 0 0; 0 0 0; 0 0 1]

spoil(::AbstractSpin) = idealspoiling

"""
    spoil!(spin)

Apply ideal spoiling to the given spin.
"""
spoil!(spin::AbstractSpin) = mul!(spin.M, idealspoiling)

applydynamics!(spin::AbstractSpin, ::IdealSpoilingMatrix) = spoil!(spin)

spoil!(A::IdealSpoilingMatrix, ::Nothing, spin::AbstractSpin, ::IdealSpoiling, ::Any = nothing) = nothing
spoil!(::Nothing, ::Nothing, spin::AbstractSpin, ::RFSpoiling, ::Any = nothing) = nothing

function spoil!(
    A,
    B,
    spin::AbstractSpin,
    spoiling::AbstractSpoiling,
    workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    freeprecess!(A, B, spin, spoiler_gradient_duration(spoiling), spoiler_gradient(spoiling), workspace)

end
