"""
    mese!(spin, TR, TE, nechoes; αex, θex, αinv, θinv)

Simulate the steady-state signal acquired from a multi-echo spin echo (MESE)
sequence, assuming instantaneous excitations and ideal spoiling.

# Arguments
- `spin::AbstractSpin`: Spin whose signal to acquire
- `TR::Real`: Repetition time of the sequence (ms)
- `TE::Real`: Echo time of the sequence (ms)
- `nechoes::Integer`: Number of spin echoes to acquire
- `αex::Real = π/2`: Flip angle for the initial excitation (rad)
- `θex::Real = 0`: Phase of the initial excitation (rad)
- `αinv::Real = π`: Flip angle for the refocussing (inversion) pulse (rad)
- `θinv::Real = -π/2`: Phase of the refocussing (inversion) pulse (rad)

# Return
- `signal::Vector{Float64}`: Steady-state signal at each of the `nechoes` echo
    times
"""
function mese!(spin::AbstractSpin, TR, TE, nechoes; αex = π/2, θex = 0, αinv = π, θinv = -π/2)

    if TR < TE * nechoes
        error("TR < TE * nechoes")
    end

    # Precompute spin dynamics
    Dex = excitation(spin, θex, αex)
    Dinv = excitation(spin, θinv, αinv)
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
