struct ExcitationMatrix{T<:Real}
    A::BlochMatrix{T}
end

ExcitationMatrix{T}() where {T} = ExcitationMatrix(BlochMatrix{T}())
ExcitationMatrix() = ExcitationMatrix(BlochMatrix())

abstract type AbstractRF end

struct InstantaneousRF{T<:Real} <: AbstractRF
    α::T
    θ::T
end

InstantaneousRF(α, θ = zero(α)) = InstantaneousRF(promote(α, θ)...)

struct RF{T<:Real} <: AbstractRF
    α::Vector{T}
    θ::Vector{T}
    Δt::Float64
end

# waveform in Gauss, Δt in ms
RF(waveform, Δt) = RF(GAMMA .* abs.(waveform) .* (Δt / 1000), angle.(waveform), Δt)

Base.length(rf::RF) = length(rf.α)

abstract type AbstractExcitationWorkspace end

struct ExcitationWorkspace{T<:Real} <: AbstractExcitationWorkspace
    Af::FreePrecessionMatrix{T}
    Bf::Magnetization{T}
    Ae::ExcitationMatrix{T}
    tmpA1::BlochMatrix{T}
    tmpA2::BlochMatrix{T}
    tmpB1::Magnetization{T}
    tmpB2::Magnetization{T}
    freeprecess_workspace::Nothing
end

ExcitationWorkspace(::Spin{T}) where {T} = ExcitationWorkspace(
                                                    FreePrecessionMatrix{T}(),
                                                    Magnetization{T}(),
                                                    ExcitationMatrix{T}(),
                                                    BlochMatrix{T}(),
                                                    BlochMatrix{T}(),
                                                    Magnetization{T}(),
                                                    Magnetization{T}(),
                                                    nothing)

struct ExcitationWorkspaceMC{T<:Real,N} <: AbstractExcitationWorkspace
    Af::BlochMcConnellMatrix{T,N}
    Bf::MagnetizationMC{T,N}
    Ae::ExcitationMatrix{T}
    tmpA1::BlochMcConnellMatrix{T,N}
    tmpA2::BlochMcConnellMatrix{T,N}
    tmpB1::MagnetizationMC{T,N}
    tmpB2::MagnetizationMC{T,N}
    freeprecess_workspace::BlochMcConnellWorkspace{T,N}
end

function ExcitationWorkspace(
    spin::SpinMC{T,N},
    freeprecess_workspace::BlochMcConnellWorkspace = BlochMcConnellWorkspace(spin)
) where {T,N}

    ExcitationWorkspaceMC(BlochMcConnellMatrix{T}(N),
                          MagnetizationMC{T}(N),
                          ExcitationMatrix{T}(),
                          BlochMcConnellMatrix{T}(N),
                          BlochMcConnellMatrix{T}(N),
                          MagnetizationMC{T}(N),
                          MagnetizationMC{T}(N),
                          freeprecess_workspace)

end

struct ExcitationWorkspaceMCApproximate{T<:Real,N} <: AbstractExcitationWorkspace
    Af::BlochMcConnellMatrix{T,N}
    Bf::MagnetizationMC{T,N}
    Ae::ExcitationMatrix{T}
    tmpA1::BlochMcConnellMatrix{T,N}
    tmpA2::BlochMcConnellMatrix{T,N}
    tmpB1::MagnetizationMC{T,N}
    tmpB2::MagnetizationMC{T,N}
    freeprecess_workspace::Nothing
end

function ExcitationWorkspace(
    spin::SpinMC{T,N},
    ::Nothing
) where {T,N}

    ExcitationWorkspaceMC(BlochMcConnellMatrix{T}(N),
                          MagnetizationMC{T}(N),
                          ExcitationMatrix{T}(),
                          BlochMcConnellMatrix{T}(N),
                          BlochMcConnellMatrix{T}(N),
                          MagnetizationMC{T}(N),
                          MagnetizationMC{T}(N),
                          nothing)

end

"""
    excitation(spin, θ, α)

Simulate instantaneous excitation with flip angle `α` about an axis that makes
angle `θ` with the positive x-axis.

# Arguments
- `spin::AbstractSpin`: Spin to excite
- `θ::Real`: Orientation of the axis about which to excite (rad)
- `α::Real`: Flip angle (rad)

# Return
- `A::Matrix`: Matrix that describes the excitation
- `B::Vector = zeros(length(spin.M))`: Not used, but included because other
    methods of `excitation` return a nontrivial value here

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, _) = excitation(spin, π/4, π/2); A * spin.M
3-element Array{Float64,1}:
  0.7071067811865476
 -0.7071067811865475
  6.123233995736766e-17
```
"""
function excitation(spin::Spin, θ::Real, α::Real)

    A = rotatetheta(θ, α)
    B = zeros(length(spin.M))
    return (A, B)

end

function excitation(spin::SpinMC, θ::Real, α::Real)

    A = kron(Diagonal(ones(Bool, spin.N)), rotatetheta(θ, α))
    B = zeros(length(spin.M))
    return (A, B)

end

function excite!(A::ExcitationMatrix, spin::Spin, rf::InstantaneousRF)

    rotatetheta!(A.A, rf.θ, rf.α)

end

"""
    excitation(spin, rf, Δθ, grad, dt)

Simulate non-instantaneous excitation using the hard pulse approximation.

# Arguments
- `spin::AbstractSpin`: Spin to excite
- `rf::Vector{<:Number}`: RF waveform (G); its magnitude determines the flip
    angle and its phase determines the axis of rotation
- `Δθ::Real`: Additional RF phase (e.g., for RF spoiling) (rad)
- `grad::Union{Matrix{<:Real},Vector{<:Real}}`: Gradients to play during
    excitation (G/cm); should be a 3-vector if the gradients are constant during
    excitation, otherwise it should be a 3×(length(rf)) matrix
- `dt::Real`: Time step (ms)

# Return
- `A::Matrix`: Matrix that describes excitation and relaxation
- `B::Vector`: Vector that describes excitation and relaxation
"""
function excitation(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                    grad::AbstractArray{<:Real,2}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    A = Diagonal(ones(Bool, length(spin.M)))
    B = zeros(length(spin.M))
    for t = 1:T
        (Af, Bf) = freeprecess(spin, dt/2, grad[:,t])
        (Ae, _) = excitation(spin, θ[t], α[t])
        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf
    end
    return (A, B)

end

# Excitation with constant gradient
function excitation(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                    grad::AbstractArray{<:Real,1}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    A = Diagonal(ones(Bool, length(spin.M)))
    B = zeros(length(spin.M))
    (Af, Bf) = freeprecess(spin, dt/2, grad)
    for t = 1:T
        (Ae, _) = excitation(spin, θ[t], α[t])
        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf
    end
    return (A, B)

end

function excite!(
    A::Union{<:BlochMatrix,<:BlochMcConnellMatrix},
    B::Union{<:Magnetization,<:MagnetizationMC},
    spin::AbstractSpin,
    rf::RF,
    Δθ::Real,
    grad::Union{<:AbstractVector{<:Gradient},<:Gradient},
    workspace::AbstractExcitationWorkspace = ExcitationWorkspace(spin)
)

    if grad isa AbstractVector
        freeprecess!(workspace.Af, workspace.Bf, spin, rf.Δt / 2, grad[1], workspace.freeprecess_workspace)
    else
        freeprecess!(workspace.Af, workspace.Bf, spin, rf.Δt / 2, grad, workspace.freeprecess_workspace)
    end
    rotatetheta!(workspace.Ae.A, rf.θ[1] + Δθ, rf.α[1])

    mul!(workspace.tmpA2, workspace.Ae, workspace.Af)
    mul!(A, workspace.Af, workspace.tmpA2)

    mul!(workspace.tmpB2, workspace.Ae, workspace.Bf)
    mul!(B, workspace.Af, workspace.tmpB2)
    add!(B, workspace.Bf)

    T = length(rf)
    for t = 2:T

        if grad isa AbstractVector
            freeprecess!(workspace.Af, workspace.Bf, spin, rf.Δt / 2, grad[t], workspace.freeprecess_workspace)
        end
        rotatetheta!(workspace.Ae.A, rf.θ[t] + Δθ, rf.α[t])

        mul!(workspace.tmpA1, workspace.Af, A)
        mul!(workspace.tmpA2, workspace.Ae, workspace.tmpA1)
        mul!(A, workspace.Af, workspace.tmpA2)

        mul!(workspace.tmpB1, workspace.Af, B)
        add!(workspace.tmpB1, workspace.Bf)
        mul!(workspace.tmpB2, workspace.Ae, workspace.tmpB1)
        mul!(B, workspace.Af, workspace.tmpB2)
        add!(B, workspace.Bf)

    end

    return nothing

end

"""
    excitation!(spin, ...)

Apply excitation to the given spin.
"""
function excitation!(spin::AbstractSpin, θ::Real, α::Real)

    (A, _) = excitation(spin, θ, α)
    applydynamics!(spin, A)

end

# Use this function if using RF spoiling (because A and B need to be
# recalculated for each TR, so directly modifying the magnetization should be
# faster in this case)
function excitation!(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                     grad::AbstractArray{<:Real,2}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    for t = 1:T
        (Af, Bf) = freeprecess(spin, dt/2, grad[:,t])
        (Ae, _) = excitation(spin, θ[t], α[t])
        applydynamics!(spin, Af, Bf)
        applydynamics!(spin, Ae)
        applydynamics!(spin, Af, Bf)
    end

end

function excitation!(spin::AbstractSpin, rf::AbstractArray{<:Number,1}, Δθ::Real,
                     grad::AbstractArray{<:Real,1}, dt::Real)

    T = length(rf)
    α = GAMMA * abs.(rf) * dt/1000 # Flip angle in rad
    θ = angle.(rf) .+ Δθ # RF phase in rad
    (Af, Bf) = freeprecess(spin, dt/2, grad)
    for t = 1:T
        (Ae, _) = excitation(spin, θ[t], α[t])
        applydynamics!(spin, Af, Bf)
        applydynamics!(spin, Ae)
        applydynamics!(spin, Af, Bf)
    end

end
