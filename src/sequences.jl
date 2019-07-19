"""
    mese!(spin, TR, TE, nechoes; αex, αinv)

Simulate the steady-state signal acquired from a multi-echo spin echo (MESE)
sequence, assuming instantaneous excitations and ideal spoiling.

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `nechoes::Integer`: Number of spin echoes to acquire
- `αex::Real = π/2`: Flip angle for the initial excitation (rad)
- `αinv::Real = π`: Flip angle for the refocussing (inversion) pulse (rad)

# Return
- `signal::Vector{ComplexF64}`: Steady-state signal at each of the `nechoes`
    echo times
"""
function mese!(spin::AbstractSpin, TR::Real, TE::Real, nechoes::Integer;
               αex::Real = π/2, αinv::Real = π)

    TR >= TE * nechoes ||
        throw(ArgumentError("TR must be greater than or equal to TE * nechoes"))

    # Precompute spin dynamics
    Dex = excitation(spin, 0, αex)
    Dinv = excitation(spin, -π/2, αinv)
    Dte = freeprecess(spin, TE/2)
    Dtr = freeprecess(spin, TR - TE * nechoes)
    Decho = combine(Dte, Dinv, Dte)
    Dspoil = (spoil(spin), zeros(length(spin.M)))

    # Calculate steady-state magnetization immediately following excitation
    D = Decho
    for e = 2:nechoes
        D = combine(D, Decho)
    end
    (A, B) = combine(D, Dtr, Dspoil, Dex)
    spin.M[:] = (I - A) \ B

    # Calculate steady-state signal at each echo
    signal = zeros(ComplexF64, nechoes)
    for e = 1:nechoes
        applydynamics!(spin, Decho...)
        signal[e] = spin.signal
    end

    return signal

end

"""
    mese(spin, TR, TE, nechoes; αex, αinv, nspins, ncycles)

Simulate the steady-state signal acquired from a multi-echo spin echo (MESE)
sequence, assuming instantaneous excitations and nonideal gradient spoiling.

The output is a vector of `nechoes` signal values that are the complex mean of
`nspins` copies of `spin` at slightly offset positions to induce `ncycles`
cycles of phase across the spins.

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `nechoes::Integer`: Number of spin echoes to acquire
- `αex::Real = π/2`: Flip angle for the initial excitation (rad)
- `αinv::Real = π`: Flip angle for the refocussing (inversion) pulse (rad)
- `nspins::Integer = 10`: Number of intra-voxel spins to simulate
- `ncycles::Real = 1`: Number of cycles of phase across a voxel induced by each
    spoiler gradient (typically integer-valued)

# Return
- `signal::Vector{ComplexF64}`: Ensemble steady-state signal at each of the
    `nechoes` echo times
"""
function mese(spin::AbstractSpin, TR::Real, TE::Real, nechoes::Integer;
              αex::Real = π/2, αinv::Real = π, nspins::Integer = 10,
              ncycles::Real = 1)

    # Pick an arbitrary gradient strength and duration and compute the spatial
    # locations that will provide a uniform ncycles of phase
    gradz = 0.3 # G/cm
    Tg = TE / 4 # ms
    zmax = ncycles * 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:nspins)/nspins  * zmax

    # Make nspins copies of the provided spin at different spatial locations
    spins = varyz(spin, z)

    # Compute the MESE signal for each spin
    signals = map(spin -> mese!(spin, TR, TE, nechoes, [0,0,gradz], Tg,
                            αex = αex, αinv = αinv), spins) # [nspins][nechoes]

    # Find the average signal for each echo
    signal = map(s -> sum(s) / nechoes, signals)

    return signal

end

"""
    varyz(spin, z)

Copy `spin`, offsetting the z position of each copy by elements of z.
"""
function varyz(spin::Spin, z)

    map(z -> Spin(spin.M0, spin.T1, spin.T2, spin.Δf, spin.pos + [0,0,z]), z)

end

function varyz(spin::SpinMC, z)

    map(z -> SpinMC(spin.M0, spin.frac, spin.T1, spin.T2, spin.Δf, spin.τ,
                    spin.pos + [0,0,z]), z)

end

"""
    mese!(spin, TR, TE, nechoes, grad, Tg; αex, αinv)

Simulate the steady-state signal acquired from a multi-echo spin echo (MESE)
sequence, assuming instantaneous excitations and nonideal gradient spoiling.

The spoiler gradient is played before and after each refocusing pulse. There is
no RF spoiling. The phase of the excitation pulse is not reversed after each TR.

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `nechoes::Integer`: Number of spin echoes to acquire
- `grad::AbstractArray{<:Real}`: Spoiler gradient amplitudes [gx, gy, gz] (G/cm)
- `Tg::Real`: Spoiler gradient duration (ms)
- `αex::Real = π/2`: Flip angle for the initial excitation (rad)
- `αinv::Real = π`: Flip angle for the refocussing (inversion) pulse (rad)

# Return
- `signal::Vector{ComplexF64}`: Steady-state signal at each of the `nechoes`
    echo times
"""
function mese!(spin::AbstractSpin, TR::Real, TE::Real, nechoes::Integer,
               grad::AbstractArray{<:Real}, Tg::Real; αex::Real = π/2,
               αinv::Real = π)

    TR >= TE * nechoes ||
        throw(ArgumentError("TR must be greater than or equal to TE * nechoes"))
    Tg <= TE / 2 ||
        throw(ArgumentError("Tg must be less than or equal to TE / 2"))

    # Precompute spin dynamics
    Dex = excitation(spin, 0, αex)
    Dinv = excitation(spin, -π/2, αinv)
    Dg = freeprecess(spin, Tg, grad)
    Dte = freeprecess(spin, TE/2 - Tg)
    Dtr = freeprecess(spin, TR - TE * nechoes)
    Decho = combine(Dte, Dg, Dinv, Dg, Dte)

    # Calculate steady-state magnetization immediately following excitation
    D = Decho
    for e = 2:nechoes
        D = combine(D, Decho)
    end
    (A, B) = combine(D, Dtr, Dex)
    spin.M[:] = (I - A) \ B

    # Calculate steady-state signal at each echo
    signal = zeros(ComplexF64, nechoes)
    for e = 1:nechoes
        applydynamics!(spin, Decho...)
        signal[e] = spin.signal
    end

    return signal

end

"""
    spgr!(spin, TR, TE, α)

Simulate the steady-state signal acquired from a spoiled gradient echo (SPGR)
sequence, assuming instantaneous excitations and ideal spoiling.

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `α::Real`: Flip angle (rad)

# Return
- `signal::ComplexF64`: Steady-state signal
"""
function spgr!(spin::AbstractSpin, TR::Real, TE::Real, α::Real)

    TR >= TE || throw(ArgumentError("TR must be greater than or equal to TE"))

    # Precompute spin dynamics
    Dex = excitation(spin, 0, α)
    Dte = freeprecess(spin, TE)
    Dtr = freeprecess(spin, TR - TE)
    Dspoil = (spoil(spin), zeros(length(spin.M)))

    # Calculate steady-state magnetization immediately following excitation
    (A, B) = combine(Dte, Dtr, Dspoil, Dex)
    spin.M[:] = (I - A) \ B

    # Calculate steady-state signal at echo time
    applydynamics!(spin, Dte...)

    return spin.signal

end

"""
    spgr!(spin, TR, TE, α, grad, Tg; Δθinc, nTR)

Simulate the steady-state signal acquired from a spoiled gradient echo (SPGR)
sequence, assuming instantaneous excitations and nonideal spoiling.

The steady-state signal is a true steady-state when `Δθinc` is 0, otherwise
`nTR` TR's are simulated to produce a pseudo steady-state. (Note, it is only a
pseudo steady-state when many spins are averaged together.)

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `α::Real`: Flip angle (rad)
- `grad::AbstractArray{<:Real}`: Spoiler gradient amplitudes [gx, gy, gz] (G/cm)
- `Tg::Real`: Spoiler gradient duration (ms)
- `Δθinc::Real = deg2rad(117)`: RF phase increment for RF spoiling (rad)
- `nTR::Integer = 100`: Number of TR's to simulate when `Δθinc` is not 0

# Return
- `signal::ComplexF64`: Steady-state signal
"""
function spgr!(spin::AbstractSpin, TR::Real, TE::Real, α::Real,
               grad::AbstractArray{<:Real}, Tg::Real;
               Δθinc::Real = deg2rad(117), nTR::Integer = 100)

    # I would prefer to split this into two functions, but it is not possible
    # (without having multiple names like spgr_rfspoil) because I want Δθinc and
    # nTR to be optional arguments, but doing so would case the RF-spoiled case
    # and the non-RF-spoiled case to have the same method signature. And I like
    # the simplicity of having just one spgr function, so I don't want to add,
    # e.g., spgr_rfspoil.

    TR >= TE + Tg ||
        throw(ArgumentError("TR must be greater than or equal to TE + Tg"))

    # Precompute spin dynamics
    Dte = freeprecess(spin, TE)
    Dtr = freeprecess(spin, TR - TE - Tg)
    Dtg = freeprecess(spin, Tg, grad)
    D = combine(Dte, Dtr, Dtg)

    if Δθinc == 0

        # Precompute spin dynamics
        Dex = excitation(spin, 0, α)

        # Calculate steady-state magnetization immediately following excitation
        (A, B) = combine(D, Dex)
        spin.M[:] = (I - A) \ B

        # Calculate steady-state signal at echo time
        applydynamics!(spin, Dte...)

        return spin.signal

    else

        # Initialize RF spoiling parameters
        θ = 0
        Δθ = Δθinc

        # Simulate nTR TR's
        for rep = 1:nTR

            excitation!(spin, θ, α)
            applydynamics!(spin, D...)
            θ += Δθ
            Δθ += Δθinc

        end

        # Calculate signal at echo time
        excitation!(spin, θ, α)
        applydynamics!(spin, Dte...)

        return spin.signal * exp(im * θ)

    end

end
