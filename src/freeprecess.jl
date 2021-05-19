struct BlochMcConnellWorkspace{T<:Real,N}
    A::BlochMcConnellDynamicsMatrix{T,N}
    expmworkspace::MatrixExponentialWorkspace{T,N}

    # N is number of compartments
    BlochMcConnellWorkspace(T::Type{<:Real}, N) =
        new{T,N}(BlochMcConnellDynamicsMatrix{T}(N),
                 MatrixExponentialWorkspace{T}(N))
end

BlochMcConnellWorkspace(::Type{SpinMC{T,N}}) where {T,N} = BlochMcConnellWorkspace(T, N)
BlochMcConnellWorkspace(spin::SpinMC) = BlochMcConnellWorkspace(typeof(spin))

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

    A = FreePrecessionMatrix()
    B = Magnetization()
    freeprecess!(A, B, t, M0, T1, T2, Δf)
    return (A, B)

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

"""
    freeprecess(spin, t)

Simulate free-precession for the given spin.

# Arguments
- `spin::AbstractSpin`: Spin that is free-precessing
- `t::Real`: Duration of free-precession (ms)

# Return
- `A::Matrix`: Matrix that describes relaxation and precession
- `B::Vector`: Vector that describes recovery

# Examples
```jldoctest
julia> spin = Spin([1.0, 0.0, 0.0], 1, 1000, 100, 3.75)
Spin([1.0, 0.0, 0.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, B) = freeprecess(spin, 100); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404048
```
"""
freeprecess_old(spin::Spin_old, t::Real) =
    freeprecess(t, spin.M0, spin.T1, spin.T2, spin.Δf)

function freeprecess_old(spin::SpinMC_old, t::Real)

    E = expm(t * spin.A)
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess(spin::Spin, t)

    A = FreePrecessionMatrix{eltype(spin)}()
    B = Magnetization{eltype(spin)}()
    freeprecess!(A, B, spin, t)
    return (A, B)

end

function freeprecess(spin::SpinMC{T,N}, t, workspace = BlochMcConnellWorkspace(spin)) where {T,N}

    A = BlochMcConnellMatrix{T}(N)
    B = MagnetizationMC{T}(N)
    freeprecess!(A, B, spin, t, workspace)
    return (A, B)

end

function freeprecess!(
    A::FreePrecessionMatrix,
    B::Magnetization,
    spin::Spin,
    t::Real,
    workspace::Nothing = nothing
)

    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf)

end

function freeprecess!(
    A::BlochMcConnellMatrix,
    B::MagnetizationMC,
    spin::SpinMC,
    t::Real,
    workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin)
)

    expm!(A, workspace, spin, t)
    subtractmul!(B, I, A, spin.Meq)
    return nothing

end

"""
    freeprecess(spin, t, grad)

Simulate free-precession for the given spin in the presence of a gradient.

# Arguments
- `spin::AbstractSpin`: Spin that is free-precessing
- `t::Real`: Duration of free-precession (ms)
- `grad::AbstractVector{<:Real}`: Gradient amplitudes [gx, gy, gz] (G/cm)

# Return
- `A::Matrix`: Matrix that describes relaxation and precession
- `B::Vector`: Vector that describes recovery

# Examples
```jldoctest
julia> spin = Spin([1.0, 0.0, 0.0], 1, 1000, 100, 0, [0, 0, 3.75])
Spin([1.0, 0.0, 0.0], 1.0, 1000.0, 100.0, 0.0, [0.0, 0.0, 3.75])

julia> (A, B) = freeprecess(spin, 100, [0, 0, 1/GAMBAR]); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404048
```
"""
function freeprecess(spin::Spin, t::Real, grad::AbstractArray{<:Real,1})

    gradfreq = GAMBAR * (grad[1] * spin.pos.x + grad[2] * spin.pos.y + grad[3] * spin.pos.z) # Hz
    freeprecess(t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

# See equation (6.9) in Gopal Nataraj's PhD thesis
function freeprecess(spin::SpinMC{T,N}, t::Real, grad::AbstractArray{<:Real,1}) where {T,N}

    A = Array{T}(undef, 3N, 3N)
    for j = 1:N, i = 1:N
        ii = 3i-2:3i
        jj = 3j-2:3j
        if i == j
            tmp = sum(spin.r[i][k] for k = 1:N) # 1/ms
            r1 = -1 / spin.T1[i] - tmp # 1/ms
            r2 = -1 / spin.T2[i] - tmp # 1/ms
            Δω = 2π * spin.Δf[i] / 1000 # rad/ms
            A[ii,jj] = [r2 Δω 0; -Δω r2 0; 0 0 r1] # Left-handed rotation
        else
            A[ii,jj] = spin.r[j][i] * Diagonal(ones(Bool, 3))
        end
    end
    gradfreq = GAMMA * (grad[1] * spin.pos.x + grad[2] * spin.pos.y + grad[3] * spin.pos.z) / 1000 # rad/ms
    ΔA = diagm(1 => repeat([gradfreq, 0, 0], spin.N), # Left-handed rotation
              -1 => repeat([-gradfreq, 0, 0], spin.N))[1:3spin.N,1:3spin.N]
    E = expm(t * (A + ΔA))
    B = (Diagonal(ones(Bool, size(E, 1))) - E) * spin.Meq
    return (E, B)

end

function freeprecess!(
    A::FreePrecessionMatrix,
    B::Magnetization,
    spin::Spin,
    t::Real,
    grad::Gradient,
    workspace::Nothing = nothing
)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    freeprecess!(A, B, t, spin.M0, spin.T1, spin.T2, spin.Δf + gradfreq)

end

function freeprecess!(
    A::BlochMcConnellMatrix,
    B::MagnetizationMC,
    spin::SpinMC,
    t::Real,
    grad::Gradient,
    workspace::Union{Nothing,<:BlochMcConnellWorkspace} = BlochMcConnellWorkspace(spin)
)

    gradfreq = gradient_frequency(grad, spin.pos) # Hz
    expm!(A, workspace, spin, t, gradfreq)
    subtractmul!(B, I, A, spin.Meq)
    return nothing

end

# Exact matrix exponential
function expm!(expAt, workspace, spin, t, gradfreq = 0)

    index = 0
    for i = 1:spin.N, j = 1:spin.N

        if i == j

            r_out = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms

            A = workspace.A.A[i]
            A.R1 = -1 / spin.T1[i] - r_out # 1/ms
            A.R2 = -1 / spin.T2[i] - r_out # 1/ms
            A.Δω = 2π * (spin.Δf[i] + gradfreq) / 1000 # rad/ms

        else

            index += 1
            workspace.A.E[index].r = spin.r[i][j] # 1/ms

        end

    end

    mul!(workspace.A, t)
    expm!(expAt, workspace.A, workspace.expmworkspace)
    return nothing

end

# Approximate matrix exponential
# See page 2 of http://doi.org/10.1137/0714065
# (unnumbered equation with o(||E||) term)
function expm!(expAt, ::Nothing, spin, t, gradfreq = 0)

    for j = 1:spin.N, i = 1:spin.N

        A = expAt.A[i][j]

        if i == j

            r_out = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms
            E1 = exp(-t * (1 / spin.T1[i] + r_out))
            E2 = exp(-t * (1 / spin.T2[i] + r_out))
            θ = 2π * (spin.Δf[i] + gradfreq) * t / 1000 # rad
            (s, c) = sincos(θ)
            E2c = E2 * c
            E2s = E2 * s

            A.a11 = E2c
            A.a21 = -E2s
            A.a31 = 0
            A.a12 = E2s
            A.a22 = E2c
            A.a32 = 0
            A.a13 = 0
            A.a23 = 0
            A.a33 = E1

        else

            r_out_i = sum(spin.r[i][k] for k = 1:spin.N) # 1/ms
            r_out_j = sum(spin.r[j][k] for k = 1:spin.N) # 1/ms
            R1i = 1 / spin.T1[i] + r_out_i # 1/ms
            R1j = 1 / spin.T1[j] + r_out_j # 1/ms
            R2i = 1 / spin.T2[i] + r_out_i # 1/ms
            R2j = 1 / spin.T2[j] + r_out_j # 1/ms
            rji = spin.r[j][i]

            R2ji = R2j - R2i
            R2ji² = R2ji^2
            E2i = exp(-t * R2i)
            E2ji = exp(-t * R2ji)
            Δωji = 2π * (spin.Δf[j] - spin.Δf[i]) / 1000 # rad/ms
            Δωji² = Δωji^2
            θi = 2π * (spin.Δf[i] + gradfreq) * t / 1000 # rad
            θj = 2π * (spin.Δf[j] + gradfreq) * t / 1000 # rad
            (si, ci) = sincos(θi)
            (sj, cj) = sincos(θj)
            tmpc = rji * E2i * ((ci - E2ji * cj) / R2ji + Δωji * (E2ji * sj - si) / R2ji²) / (1 + Δωji² / R2ji²)
            tmps = rji * E2i * ((si - E2ji * sj) / R2ji + Δωji * (ci - E2ji * cj) / R2ji²) / (1 + Δωji² / R2ji²)

            A.a11 = tmpc
            A.a21 = -tmps
            A.a31 = 0
            A.a12 = tmps
            A.a22 = tmpc
            A.a32 = 0
            A.a13 = 0
            A.a23 = 0
            A.a33 = rji * E2i * (1 - E2ji) / R2ji

        end

    end

    return nothing

end

"""
    freeprecess!(spin, ...)

Apply free-precession to the given spin.
"""
function freeprecess!(spin::AbstractSpin, args...)

    (A, B) = freeprecess(spin, args...)
    applydynamics!(spin, A, B)

end

"""
    combine(D...)

Combine the matrices and vectors that describe the dynamics of a spin into one
matrix and one vector.

# Arguments
- `D::Tuple{<:AbstractArray{<:Real,2},<:AbstractVector{<:Real}}...`: List of
    pairs of matrices and vectors, i.e., ((A1, B1), (A2, B2), ...), where the
    A's are matrices and the B's are vectors

# Return
- `A::Matrix`: Matrix that describes the spin dynamics
- `B::Vector`: Vector that describes the spin dynamics

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> D1 = excitation(spin, 0, π/2);

julia> D2 = freeprecess(spin, 100);

julia> (A, B) = combine(D1, D2); A * spin.M + B
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404054
```
"""
function combine(D::Tuple{<:AbstractArray{<:Real,2},<:AbstractArray{<:Real,1}}...)

  (A, B) = D[1]
  for i = 2:length(D)
    (Ai, Bi) = D[i]
    A = Ai * A
    B = Ai * B + Bi
  end
  return (A, B)

end

function combine!(A, B, A1, B1, A2, B2)

    mul!(A, A2, A1)
    mul!(B, A2, B1)
    add!(B, B2)

end

function combine!(A, B, A1, B1, A2, ::Nothing)

    mul!(A, A2, A1)
    mul!(B, A2, B1)

end

function combine!(A, B, A1, ::Nothing, A2, B2)

    mul!(A, A2, A1)
    copyto!(B, B2)

end

function combine!(A, B, A1, B1, ::Nothing, ::Nothing)

    copyto!(A, A1)
    copyto!(B, B1)

end

function combine!(A, B, ::Nothing, ::Nothing, A2, B2)

    copyto!(A, A2)
    copyto!(B, B2)

end

combine!(A, A2, A1) = mul!(A, A2, A1)

"""
    applydynamics!(spin, A[, B])

Apply dynamics to the given spin.

# Arguments
- `spin::AbstractSpin`: Spin to which to apply dynamics
- `A::Matrix`: Matrix with dynamics
- `B::Vector = zeros(length(spin.M))`: Vector with dynamics

# Examples
```jldoctest
julia> spin = Spin(1, 1000, 100, 3.75)
Spin([0.0, 0.0, 1.0], 1.0, 1000.0, 100.0, 3.75, [0.0, 0.0, 0.0])

julia> (A, _) = excitation(spin, 0, π/2); applydynamics!(spin, A)

julia> (A, B) = freeprecess(spin, 100); applydynamics!(spin, A, B)

julia> spin.M
3-element Array{Float64,1}:
 -0.2601300475114444
 -0.2601300475114445
  0.09516258196404054
```
"""
function applydynamics!(spin::AbstractSpin, A::AbstractArray{<:Real,2},
                        B::AbstractArray{<:Real,1})

    tmp = A * spin.M
    add!(tmp, B)
    copyto!(spin.M, tmp)
    return nothing

end

function applydynamics!(spin::AbstractSpin, A::AbstractArray{<:Real,2})

    copyto!(spin.M, A * spin.M)
    return nothing

end

function applydynamics!(spin::AbstractSpin, BtoM, A, B)

    copyto!(BtoM, B)
    muladd!(BtoM, A, spin.M) # BtoM .= A * spin.M + BtoM
    copyto!(spin.M, BtoM)
    return nothing

end

function applydynamics!(spin::AbstractSpin, BtoM, A, ::Nothing = nothing)

    mul!(BtoM, A, spin.M)
    copyto!(spin.M, BtoM)
    return nothing

end
