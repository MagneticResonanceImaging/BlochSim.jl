using BlochSim, Test, MAT

# ------------------------------------------------------------------------------
# Begin tests for comparing to Brian Hargreaves' MATLAB code
# ------------------------------------------------------------------------------

function testA5b()

    answer = matread("matlabtestdata/testA5b.mat")

    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 10 # Hz
    M = [1.0, 0, 0]
    dt = 1 # ms
    time = 0:dt:1000 # ms
    spin = Spin(M, 1, T1, T2, Δf, [0,0,0])
    mag = zeros(3,length(time))
    mag[:,1] = M
    (A, B) = freeprecess(spin, dt)
    for t = 2:length(time)
      applydynamics!(spin, A, B)
      mag[:,t] = spin.M
    end

    return mag ≈ answer["M"]

end

function testF1a()

    answer = matread("matlabtestdata/testF1a.mat")

    t = 2.3 # ms
    α = π/4 # rad
    θ = 0 # rad
    Δf = -500:5:500 # Hz
    T1 = 600 # ms
    T2 = 100 # ms
    spins = map(Δf -> Spin([0,0,1.0], 1, T1, T2, Δf, [0,0,0]), Δf)
    for f = 1:length(Δf)
      (Ae, _) = excitation(spins[f], θ, α)
      (Af, Bf) = freeprecess(spins[f], t)
      applydynamics!(spins[f], Ae)
      applydynamics!(spins[f], Af, Bf)
      applydynamics!(spins[f], Ae)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF1b()

    answer = matread("matlabtestdata/testF1b.mat")

    t = 2.3 # ms
    α = π/4 # rad
    θ = 0 # rad
    Δf = 0 # Hz
    T1 = 600 # ms
    T2 = 100 # ms
    grad = [0.1, 0, 0] # G/cm
    xpos = -2:0.01:2 # cm
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    for i = 1:length(xpos)
      (Ae, _) = excitation(spins[i], θ, α)
      (Af, Bf) = freeprecess(spins[i], t, grad)
      applydynamics!(spins[i], Ae)
      applydynamics!(spins[i], Af, Bf)
      applydynamics!(spins[i], Ae)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF1c()

    answer = matread("matlabtestdata/testF1c.mat")

    t = 2.3 # ms
    α = π/4 # rad
    θ = 0 # rad
    Δf = 0 # Hz
    T1 = 600 # ms
    T2 = 100 # ms
    grad = [0.1, 0, 0] # G/cm
    xpos = -2:0.01:2 # cm
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    for i = 1:length(xpos)
      (Ae, _) = excitation(spins[i], θ, α)
      (Af, Bf) = freeprecess(spins[i], t, grad)
      (Ar, Br) = freeprecess(spins[i], t, -grad/2)
      applydynamics!(spins[i], Ae)
      applydynamics!(spins[i], Af, Bf)
      applydynamics!(spins[i], Ae)
      applydynamics!(spins[i], Ar, Br)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF2c()

    answer = matread("matlabtestdata/testF2c.mat")

    dt = 0.04 # ms
    t = 0:dt:6 # ms
    rf = 0.05 * sinc.(t .- 3) # G
    Δf = -1000:20:1000 # Hz
    spins = map(Δf -> Spin([0,0,1.0], 1, 600, 100, Δf, [0,0,0]), Δf)
    map(spin -> excitation!(spin, rf, 0, zeros(3), dt), spins)
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3a()

    answer = matread("matlabtestdata/testF3a.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -2:0.01:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1, 0, 0]
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    for i = 1:length(xpos)
      (Ae, Be) = excitation(spins[i], rf, 0, grad, dt)
      applydynamics!(spins[i], Ae, Be)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3c()

    answer = matread("matlabtestdata/testF3c.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    for i = 1:length(xpos)
        (Ae, Be) = excitation(spins[i], [rf; 0rf], 0, [grad -0.52grad], dt)
        applydynamics!(spins[i], Ae, Be)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3d()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    map(spin -> excitation!(spin, [rf; 0rf], 0, [grad -0.52grad], dt), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3f()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    map(spin -> excitation!(spin, [rf; 0rf], 0, [grad -0.52grad], dt), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testA5bMC()

    answer = matread("matlabtestdata/testA5b.mat")

    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 10 # Hz
    M = [1.0, 0, 0]
    dt = 1 # ms
    time = 0:dt:1000 # ms
    spin = SpinMC(M, 1, [1], [T1], [T2], [Δf], Vector{Int}(), [0,0,0])
    mag = zeros(3,length(time))
    mag[:,1] = M
    for t = 2:length(time)
      freeprecess!(spin, dt)
      mag[:,t] = spin.M
    end

    return mag ≈ answer["M"]

end

function testF1cMC()

    answer = matread("matlabtestdata/testF1c.mat")

    t = 2.3 # ms
    α = π/4 # rad
    θ = 0 # rad
    Δf = 0 # Hz
    T1 = 600 # ms
    T2 = 100 # ms
    grad = [0.1, 0, 0] # G/cm
    xpos = -2:0.01:2 # cm
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    excitation!.(spins, θ, α)
    map(spin -> freeprecess!(spin, t, grad), spins)
    excitation!.(spins, θ, α)
    map(spin -> freeprecess!(spin, t, -grad/2), spins)
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF2cMC()

    answer = matread("matlabtestdata/testF2c.mat")

    dt = 0.04 # ms
    t = 0:dt:6 # ms
    rf = 0.05 * sinc.(t .- 3) # G
    Δf = -1000:20:1000 # Hz
    spins = map(Δf -> SpinMC([0,0,1.0], 1, [1], [600], [100], [Δf], Vector{Int}(), [0,0,0]), Δf)
    map(spin -> excitation!(spin, rf, 0, zeros(3,length(rf)), dt), spins)
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3cMC()

    answer = matread("matlabtestdata/testF3c.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, [rf; 0rf], 0, [grad -0.52grad], dt), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3dMC()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, [rf; 0rf], 0, [grad -0.52grad], dt), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3fMC()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, [rf; 0rf], 0, [grad -0.52grad], dt), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

# ------------------------------------------------------------------------------
# End tests for comparing to Brian Hargreaves' MATLAB code
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin single compartment tests
# ------------------------------------------------------------------------------

function Spin1()

    s = Spin([1, 2, 3], 1, 1000, 100, 0)
    return s.M == [1, 2, 3] &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == [0, 0, 0] &&
           s.signal == 1 + 2im

end

function Spin2()

    s = Spin([1, 2, 3], 1.5, 1000, 100, 1.3, [0.5, 0.2, 1])
    return s.M == [1, 2, 3] &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == [0.5, 0.2, 1] &&
           s.signal == 1 + 2im

end

function Spin3()

    s = Spin(1, 1000, 100, 0)
    return s.M == [0, 0, 1] &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == [0, 0, 0] &&
           s.signal == 0

end

function Spin4()

    s = Spin(1.5, 1000, 100, 1.3, [0.5, 0.2, 1])
    return s.M == [0, 0, 1.5] &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == [0.5, 0.2, 1] &&
           s.signal == 0

end

function freeprecess1()

    s = Spin([1, 0, 0], 1, 1000, 100, 3.75)
    (A, B) = freeprecess(s, 100)
    M = A * s.M + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return M ≈ M_correct

end

function freeprecess2()

    s = Spin([1, 0, 0], 1, 1000, 100, 0, [0, 0, 3.75])
    (A, B) = freeprecess(s, 100, [0, 0, 1/GAMBAR])
    M = A * s.M + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return M ≈ M_correct

end

function freeprecess3()

    s = Spin([1, 0, 0], 1, 1000, 100, 0, [0, 0, 3.75])
    freeprecess!(s, 100, [0, 0, 1/GAMBAR])
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return s.M ≈ M_correct

end

function excitation1()

    s = Spin(1, 1000, 100, 3.75)
    (A, B) = excitation(s, π/4, π/2)
    M = A * s.M + B
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return M ≈ M_correct && B == [0, 0, 0]

end

function excitation2()

    s = Spin(1, Inf, Inf, 0)
    rf = fill(exp(im * π/8), 2)
    Δθ = π/8
    grad = zeros(3, 2)
    dt = 250π / GAMMA
    (A, B) = excitation(s, rf, Δθ, grad, dt)
    M = A * s.M + B
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return M ≈ M_correct && B == [0, 0, 0]

end

function excitation3()

    s = Spin(1, Inf, Inf, 0)
    rf = fill(exp(im * π/8), 2)
    Δθ = π/8
    grad = zeros(3)
    dt = 250π / GAMMA
    (A, B) = excitation(s, rf, Δθ, grad, dt)
    M = A * s.M + B
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return M ≈ M_correct && B == [0, 0, 0]

end

function excitation4()

    s = Spin(1, 1000, 100, 3.75)
    excitation!(s, π/4, π/2)
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return s.M ≈ M_correct

end

function excitation5()

    s = Spin(1, Inf, Inf, 0)
    rf = fill(exp(im * π/8), 2)
    Δθ = π/8
    grad = zeros(3, 2)
    dt = 250π / GAMMA
    excitation!(s, rf, Δθ, grad, dt)
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return s.M ≈ M_correct

end

function excitation6()

    s = Spin(1, Inf, Inf, 0)
    rf = fill(exp(im * π/8), 2)
    Δθ = π/8
    grad = zeros(3)
    dt = 250π / GAMMA
    excitation!(s, rf, Δθ, grad, dt)
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return s.M ≈ M_correct

end

function spoil1()

    s = Spin([1, 0.4, 5], 1, 1000, 100, 0)
    S = spoil(s)
    M = S * s.M
    M_correct = [0, 0, 5]
    return M ≈ M_correct && S == [0 0 0; 0 0 0; 0 0 1]

end

function spoil2()

    s = Spin([1, 0.4, 5], 1, 1000, 100, 0)
    spoil!(s)
    M_correct = [0, 0, 5]
    return s.M ≈ M_correct

end

function combine1()

    s = Spin(1, 1000, 100, 3.75)
    D1 = excitation(s, 0, π/2)
    D2 = freeprecess(s, 100)
    (A, B) = combine(D1, D2)
    M = A * s.M + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404054]
    return M ≈ M_correct

end

function applydynamics1()

    s = Spin(1, 1000, 100, 3.75)
    (A,) = excitation(s, 0, π/2)
    applydynamics!(s, A)
    (A, B) = freeprecess(s, 100)
    applydynamics!(s, A, B)
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404054]
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End single compartment tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin multicompartment tests
# ------------------------------------------------------------------------------

function SpinMC1()

    s = SpinMC([1, 2, 3], 1, [1], [1000], [100], [0], [])
    return s.N == 1 &&
           s.M == [1, 2, 3] &&
           s.Meq == [0, 0, 1] &&
           s.M0 == 1 &&
           s.frac == [1] &&
           s.T1 == [1000] &&
           s.T2 == [100] &&
           s.Δf == [0] &&
           isempty(s.τ) &&
           s.pos == [0, 0, 0] &&
           s.A == [-1/s.T2[1] 0 0; 0 -1/s.T2[1] 0; 0 0 -1/s.T1[1]] &&
           s.signal == 1 + 2im

end

function SpinMC2()

    s = SpinMC([1, 2, 3, 4, 5, 6], 2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], [-0.3, -1, 0])
    A11 = [-1/s.T2[1]-1/s.τ[1] 2π*s.Δf[1]/1000 0; -2π*s.Δf[1]/1000 -1/s.T2[1]-1/s.τ[1] 0; 0 0 -1/s.T1[1]-1/s.τ[1]]
    A12 = [1/s.τ[2] 0 0; 0 1/s.τ[2] 0; 0 0 1/s.τ[2]]
    A21 = [1/s.τ[1] 0 0; 0 1/s.τ[1] 0; 0 0 1/s.τ[1]]
    A22 = [-1/s.T2[2]-1/s.τ[2] 2π*s.Δf[2]/1000 0; -2π*s.Δf[2]/1000 -1/s.T2[2]-1/s.τ[2] 0; 0 0 -1/s.T1[2]-1/s.τ[2]]
    return s.N == 2 &&
           s.M == [1, 2, 3, 4, 5, 6] &&
           s.Meq == [0, 0, 0.4, 0, 0, 1.6] &&
           s.M0 == 2 &&
           s.frac == [0.2, 0.8] &&
           s.T1 == [400, 1000] &&
           s.T2 == [20, 100] &&
           s.Δf == [15, 0] &&
           s.τ == [20, 40] &&
           s.pos == [-0.3, -1, 0] &&
           s.A == [A11 A12; A21 A22] &&
           s.signal == 5 + 7im

end

function SpinMC3()

    s = SpinMC(1, [1], [1000], [100], [0], [])
    return s.N == 1 &&
           s.M == [0, 0, 1] &&
           s.Meq == [0, 0, 1] &&
           s.M0 == 1 &&
           s.frac == [1] &&
           s.T1 == [1000] &&
           s.T2 == [100] &&
           s.Δf == [0] &&
           isempty(s.τ) &&
           s.pos == [0, 0, 0] &&
           s.A == [-1/s.T2[1] 0 0; 0 -1/s.T2[1] 0; 0 0 -1/s.T1[1]] &&
           s.signal == 0

end

function SpinMC4()

    s = SpinMC(2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], [-0.3, -1, 0])
    A11 = [-1/s.T2[1]-1/s.τ[1] 2π*s.Δf[1]/1000 0; -2π*s.Δf[1]/1000 -1/s.T2[1]-1/s.τ[1] 0; 0 0 -1/s.T1[1]-1/s.τ[1]]
    A12 = [1/s.τ[2] 0 0; 0 1/s.τ[2] 0; 0 0 1/s.τ[2]]
    A21 = [1/s.τ[1] 0 0; 0 1/s.τ[1] 0; 0 0 1/s.τ[1]]
    A22 = [-1/s.T2[2]-1/s.τ[2] 2π*s.Δf[2]/1000 0; -2π*s.Δf[2]/1000 -1/s.T2[2]-1/s.τ[2] 0; 0 0 -1/s.T1[2]-1/s.τ[2]]
    return s.N == 2 &&
           s.M == [0, 0, 0.4, 0, 0, 1.6] &&
           s.Meq == [0, 0, 0.4, 0, 0, 1.6] &&
           s.M0 == 2 &&
           s.frac == [0.2, 0.8] &&
           s.T1 == [400, 1000] &&
           s.T2 == [20, 100] &&
           s.Δf == [15, 0] &&
           s.τ == [20, 40] &&
           s.pos == [-0.3, -1, 0] &&
           s.A == [A11 A12; A21 A22] &&
           s.signal == 0

end

function freeprecessMC1()

    s = SpinMC([1, 0, 0], 1, [1], [1000], [100], [3.75], [])
    (A, B) = freeprecess(s, 100)
    M = A * s.M + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return M ≈ M_correct

end

function freeprecessMC2()

    s = SpinMC([1, 0, 0], 1, [1], [1000], [100], [0], [], [0, 0, 3.75])
    (A, B) = freeprecess(s, 100, [0, 0, 1/GAMBAR])
    M = A * s.M + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return M ≈ M_correct

end

function excitationMC1()

    s = SpinMC(1.5, [1/3, 2/3], [400, 1000], [20, 100], [3.75, 3.75], [20, 40])
    (A, B) = excitation(s, π/4, π/2)
    M = A * s.M + B
    M_correct = [sqrt(2)/4, -sqrt(2)/4, 0, sqrt(2)/2, -sqrt(2)/2, 0]
    return M ≈ M_correct && B == [0, 0, 0, 0, 0, 0]

end

function spoilMC1()

    s = SpinMC([1, 0.4, 5, 0.2, 10, 0.2], 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    S = spoil(s)
    M = S * s.M
    M_correct = [0, 0, 5, 0, 0, 0.2]
    return M ≈ M_correct && S == [[0 0 0; 0 0 0; 0 0 1] zeros(3,3); zeros(3,3) [0 0 0; 0 0 0; 0 0 1]]

end

function spoilMC2()

    s = SpinMC([1, 0.4, 5, 0.2, 10, 0.2], 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    spoil!(s)
    M_correct = [0, 0, 5, 0, 0, 0.2]
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End multicompartment tests
# ------------------------------------------------------------------------------

@testset "Spin" begin

    @testset "Compare to MATLAB" begin

        @testset "Single Compartment" begin

            @test testA5b()
            @test testF1a()
            @test testF1b()
            @test testF1c()
            @test testF2c()
            @test testF3a()
            @test testF3c()
            @test testF3d()
            @test testF3f()

        end

        @testset "Multicompartment" begin

            @test testA5bMC()
            @test testF1cMC()
            @test testF2cMC()
            @test testF3cMC()
            @test testF3dMC()
            @test testF3fMC()

        end

    end

    @testset "Single Compartment" begin

        @test Spin1()
        @test Spin2()
        @test Spin3()
        @test Spin4()
        @test freeprecess1()
        @test freeprecess2()
        @test freeprecess3()
        @test excitation1()
        @test excitation2()
        @test excitation3()
        @test excitation4()
        @test excitation5()
        @test excitation6()
        @test spoil1()
        @test spoil2()
        @test combine1()
        @test applydynamics1()

    end

    @testset "Multicompartment" begin

        @test SpinMC1()
        @test SpinMC2()
        @test SpinMC3()
        @test SpinMC4()
        @test freeprecessMC1()
        @test freeprecessMC2()
        @test excitationMC1()
        @test spoilMC1()
        @test spoilMC2()

    end

end
