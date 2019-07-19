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
function mese!(spin::AbstractSpin, TR, TE, nechoes; αex = π/2, αinv = π)

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
function spgr!(spin::AbstractSpin, TR, TE, α)

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
function spgr!(spin::AbstractSpin, TR, TE, α, grad, Tg; Δθinc = deg2rad(117), nTR = 100)

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
