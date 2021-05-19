abstract type AbstractRF end

struct InstantaneousRF{T<:Real} <: AbstractRF
    α::T
    θ::T
end

InstantaneousRF(α, θ = zero(α)) = InstantaneousRF(promote(α, θ)...)

struct RF{T<:Real,G<:Union{<:Gradient,<:AbstractVector{<:Gradient}}} <: AbstractRF
    α::Vector{T}
    θ::Vector{T}
    Δt::Float64
    Δθ_initial::T
    grad::G
    Δθ::Ref{T} # Type Ref to enable RF-spoiling, which requires updating Δθ

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
        new{T,typeof(grad)}(α, θ, Δt, Δθ, grad, Δθ)

    end
end

# waveform in Gauss, Δt in ms
RF(waveform, Δt, Δθ, grad) = RF(GAMMA .* abs.(waveform) .* (Δt / 1000), angle.(waveform), Δt, Δθ, grad)
RF(waveform, Δt, Δθ::Real) = RF(waveform, Δt, Δθ, Gradient(0, 0, 0))
RF(waveform, Δt, grad) = RF(waveform, Δt, 0, grad)
RF(waveform, Δt) = RF(waveform, Δt, 0, Gradient(0, 0, 0))

Base.length(rf::RF) = length(rf.α)

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
    rotatetheta(θ, α)

Simulate left-handed rotation about an axis in the x-y plane that makes angle
`θ` with the negative y-axis.

# Arguments
- `θ::Real`: Orientation of the axis about which to rotate (rad)
- `α::Real`: Rotation angle (rad)

# Return
- `R::Matrix`: 3×3 matrix that describes rotation by angle `α` about an axis in
    the x-y plane that makes angle `θ` with the negative y-axis

## Note
`rotatetheta(θ, α) == rotatez(θ) * rotatey(-α) * rotatez(-θ)`

# Examples
```jldoctest
julia> R = rotatetheta(π/4, π/2); R * [0, 0, 1]
3-element Array{Float64,1}:
  0.7071067811865476
 -0.7071067811865475
  6.123233995736766e-17
```
"""
function rotatetheta(θ::Real, α::Real)

    A = BlochMatrix()
    rotatetheta!(A, θ, α)
    return A

end

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
function excitation(spin::Spin_old, θ::Real, α::Real)

    A = rotatetheta(θ, α)
    B = zeros(3)
    return (A, B)

end

function excitation(spin::SpinMC_old, θ::Real, α::Real)

    A = kron(Diagonal(ones(Bool, spin.N)), rotatetheta(θ, α))
    B = zeros(3spin.N)
    return (A, B)

end

function excite(spin::AbstractSpin, rf::InstantaneousRF)

    A = ExcitationMatrix{eltype(spin)}()
    excite!(A, spin, rf)
    return (A, nothing)

end

function excite!(A::ExcitationMatrix, spin::AbstractSpin, rf::InstantaneousRF)

    rotatetheta!(A.A, rf.θ, rf.α)

end

function excite!(A::ExcitationMatrix, ::Nothing, spin::AbstractSpin, rf::InstantaneousRF, ::Nothing = nothing)

    excite!(A, spin, rf)

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
    A = Diagonal(ones(Bool, 3spin.N))
    B = zeros(3spin.N)
    for t = 1:T
        (Af, Bf) = freeprecess(spin, dt/2, grad[:,t])
        Bf = Vector(Bf)
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
    A = Diagonal(ones(Bool, 3spin.N))
    B = zeros(3spin.N)
    (Af, Bf) = freeprecess(spin, dt/2, grad)
    Bf = Vector(Bf)
    for t = 1:T
        (Ae, _) = excitation(spin, θ[t], α[t])
        A = Af * Ae * Af * A
        B = Af * (Ae * (Af * B + Bf)) + Bf
    end
    return (A, B)

end

function excite(spin::AbstractSpin, rf::RF)

    T = eltype(spin)
    if spin isa Spin
        A = BlochMatrix{T}()
        B = Magnetization{T}()
    else
        N = spin.N
        A = BlochMcConnellMatrix{T}(N)
        B = MagnetizationMC{T}(N)
    end
    excite!(A, B, spin, rf)
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
