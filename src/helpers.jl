"""
    relax(t, M0, T1, T2)

Simulate relaxation.

# Arguments
- `t::Real`: Duration of relaxation
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant
- `T2::Real`: Spin-spin recovery time constant

## Note
Inputs `t`, `T1`, and `T2` should all have the same units, typically ms.

# Return
- `A::Matrix`: 3×3 matrix that describes relaxation
- `B::Vector`: 3-vector that describes recovery

# Examples
```jldoctest
julia> (A, B) = relax(100, 1, 1000, 100); A * [1, 0, 0] + B
3-element Array{Float64,1}:
 0.36787944117144233
 0.0
 0.09516258196404048
```
"""
function relax(t::Real, M0::Real, T1::Real, T2::Real)

    E1 = exp(-t/T1)
    E2 = exp(-t/T2)
    A = [E2 0 0; 0 E2 0; 0 0 E1]
    B = [0, 0, M0 * (1 - E1)]
    return (A, B)

end

"""
    rotatex(α)

Simulate left-handed rotation about the x-axis.

# Arguments
- `α::Real`: Rotation angle (rad)

# Return
- `R::Matrix`: 3×3 matrix that describes rotation about the x-axis by angle `α`

# Examples
```jldoctest
julia> R = rotatex(π/2); R * [0, 0, 1]
3-element Array{Float64,1}:
 0.0
 1.0
 6.123233995736766e-17
```
"""
function rotatex(α::Real)

    return [1 0 0; 0 cos(α) sin(α); 0 -sin(α) cos(α)]

end

"""
    rotatey(α)

Simulate left-handed rotation about the y-axis.

# Arguments
- `α::Real`: Rotation angle (rad)

# Return
- `R::Matrix`: 3×3 matrix that describes rotation about the y-axis by angle `α`

# Examples
```jldoctest
julia> R = rotatey(π/2); R * [0, 0, 1]
3-element Array{Float64,1}:
 -1.0
  0.0
  6.123233995736766e-17
```
"""
function rotatey(α::Real)

    return [cos(α) 0 -sin(α); 0 1 0; sin(α) 0 cos(α)]

end

"""
    rotatez(ϕ)

Simulate left-handed rotation about the z-axis.

# Arguments
- `ϕ::Real`: Rotation angle (rad)

# Return
- `R::Matrix`: 3×3 matrix that describes rotation about the z-axis by angle `ϕ`

# Examples
```jldoctest
julia> R = rotatez(π/2); R * [0, 1, 0]
3-element Array{Float64,1}:
 1.0
 6.123233995736766e-17
 0.0
```
"""
function rotatez(ϕ::Real)

    return [cos(ϕ) sin(ϕ) 0; -sin(ϕ) cos(ϕ) 0; 0 0 1]

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

    return [ sin(θ)^2+cos(α)*cos(θ)^2           sin(θ)*cos(θ)-cos(α)*sin(θ)*cos(θ)  sin(α)*cos(θ);
             sin(θ)*cos(θ)-cos(α)*sin(θ)*cos(θ) cos(θ)^2+cos(α)*sin(θ)^2           -sin(α)*sin(θ);
            -sin(α)*cos(θ)                      sin(α)*sin(θ)                       cos(α)]

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
    freeprecess(t, M0, T1, T2, Δf)

Simulate free-precession, i.e., relaxation and off-resonance precession.

# Arguments
- `t::Real`: Duration of free-precession (ms)
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant (ms)
- `T2::Real`: Spin-spin recovery time constant (ms)
- `Δf::Real`: Off-resonance frequency (Hz)

# Return
- `A::Matrix`: 3×3 matrix that describes relaxation and precession
- `B::Vector`: 3-vector that describes recovery

# Examples
```jldoctest
julia> (A, B) = freeprecess(100, 1, 1000, 100, 3.75); A * [1, 0, 0] + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404048
```
"""
function freeprecess(t::Real, M0::Real, T1::Real, T2::Real, Δf::Real)

    (A, B) = relax(t, M0, T1, T2)
    C = rotatez(2π * Δf * t/1000)
    return (C * A, B)

end

function freeprecess!(A, B, t, M0, T1, T2, Δf)

    E2 = exp(-t / T2)
    θ = 2π * Δf * t / 1000
    (s, c) = sincos(θ)

    A.E1 = exp(-t / T1)
    A.E2cosθ = E2 * c
    A.E2sinθ = E2 * s

    B.x = 0
    B.y = 0
    B.z = M0 * (1 - A.E1)

    return nothing

end
