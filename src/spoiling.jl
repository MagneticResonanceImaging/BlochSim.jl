struct Gradient{T<:Real}
    x::T
    y::T
    z::T

    function Gradient(x, y, z)

        args = promote(x, y, z)
        T = typeof(args[1])
        new{T}(args...)

    end
end

Base.show(io::IO, grad::Gradient) = print(io, "[", grad.x, ", ", grad.y, ", ", grad.z, "]")
Base.show(io::IO, ::MIME"text/plain", grad::Gradient{T}) where {T} =
    print(io, "Gradient{$T}:\n Gx = ", grad.x, " G/cm\n Gy = ", grad.y, " G/cm\n Gz = ", grad.z, " G/cm")

function gradient_frequency(grad::Gradient, pos::Position)

    return GAMBAR * (grad.x * pos.x + grad.y * pos.y + grad.z * pos.z) # Hz

end

abstract type AbstractSpoiling end

struct GradientSpoiling{T<:Real} <: AbstractSpoiling
    gradient::Gradient{T}
end

GradientSpoiling(x, y, z) = GradientSpoiling(Gradient(x, y, z))

struct RFSpoiling{T<:Real} <: AbstractSpoiling
    Δθ::T
end

struct RFandGradientSpoiling{T1<:Real,T2<:Real} <: AbstractSpoiling
    gradient::GradientSpoiling{T1}
    rf::RFSpoiling{T2}
end

RFandGradientSpoiling(grad::Gradient, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad), rf)
RFandGradientSpoiling(grad::NTuple{3,Real}, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(grad...), rf)
RFandGradientSpoiling(gx, gy, gz, rf::RFSpoiling) = RFandGradientSpoiling(GradientSpoiling(gx, gy, gz), rf)
RFandGradientSpoiling(grad::Union{<:GradientSpoiling,<:Gradient,<:NTuple{3,Real}}, Δθ) = RFandGradientSpoiling(grad, RFSpoiling(Δθ))
RFandGradientSpoiling(rf::Union{<:RFSpoiling,<:Real}, grad::Union{<:GradientSpoiling,<:Gradient,<:NTuple{3,Real}}) = RFandGradientSpoiling(grad, rf)
RFandGradientSpoiling(rf::RFSpoiling, gx, gy, gz) = RFandGradientSpoiling(gx, gy, gz, rf)

spoiler_gradient(s::GradientSpoiling) = s.gradient
spoiler_gradient(s::RFandGradientSpoiling) = spoiler_gradient(s.gradient)
rfspoiling_increment(s::RFSpoiling) = s.Δθ
rfspoiling_increment(s::RFandGradientSpoiling) = rfspoiling_increment(s.rf)
