# ------------------------------------------------------------------------------
# Begin tests for comparing to Brian Hargreaves' MATLAB code
# ------------------------------------------------------------------------------

function testA5b()

    answer = matread("matlabtestdata/testA5b.mat")

    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 10 # Hz
    M = Magnetization(1.0, 0, 0)
    dt = 1 # ms
    time = 0:dt:1000 # ms
    spin = Spin(M, 1, T1, T2, Δf, Position(0,0,0))
    mag = zeros(3,length(time))
    mag[:,1] = Vector(M)
    A = FreePrecessionMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    freeprecess!(A, B, spin, dt)
    for t = 2:length(time)
        applydynamics!(spin, BtoM, A, B)
        mag[:,t] = Vector(spin.M)
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
    spins = map(Δf -> Spin(1, T1, T2, Δf), Δf)
    rf = InstantaneousRF(α, θ)
    Ae = ExcitationMatrix()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    for f = 1:length(Δf)
        excite!(Ae, spins[f], rf)
        freeprecess!(Af, Bf, spins[f], t)
        applydynamics!(spins[f], BtoM, Ae)
        applydynamics!(spins[f], BtoM, Af, Bf)
        applydynamics!(spins[f], BtoM, Ae)
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
    grad = Gradient(0.1, 0, 0) # G/cm
    xpos = -2:0.01:2 # cm
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = InstantaneousRF(α, θ)
    Ae = ExcitationMatrix()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(Ae, spins[i], rf)
        freeprecess!(Af, Bf, spins[i], t, grad)
        applydynamics!(spins[i], BtoM, Ae)
        applydynamics!(spins[i], BtoM, Af, Bf)
        applydynamics!(spins[i], BtoM, Ae)
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
    grad = Gradient(0.1, 0, 0) # G/cm
    grad2 = Gradient(-0.05, 0, 0) # G/cm
    xpos = -2:0.01:2 # cm
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = InstantaneousRF(α, θ)
    Ae = ExcitationMatrix()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    Ar = FreePrecessionMatrix()
    Br = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(Ae, spins[i], rf)
        freeprecess!(Af, Bf, spins[i], t, grad)
        freeprecess!(Ar, Br, spins[i], t, grad2)
        applydynamics!(spins[i], BtoM, Ae)
        applydynamics!(spins[i], BtoM, Af, Bf)
        applydynamics!(spins[i], BtoM, Ae)
        applydynamics!(spins[i], BtoM, Ar, Br)
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
    spins = map(Δf -> Spin(1, 600, 100, Δf), Δf)
    rf = RF(rf, dt)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    grad = Gradient(0, 0, 0)
    map(spins) do spin
        excite!(A, B, spin, rf, 0, grad)
        applydynamics!(spin, BtoM, A, B)
    end
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
    grad = Gradient(0.1, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = RF(rf, dt)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(A, B, spins[i], rf, 0, grad)
        applydynamics!(spins[i], BtoM, A, B)
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
    grad = [Gradient(0.1, 0, 0) for i = 1:T]
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = RF(rf, dt)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(Ae, Be, spins[i], rf, 0, grad)
        freeprecess!(Af, Bf, spins[i], t[end] - dt, grad2)
        applydynamics!(spins[i], BtoM, Ae, Be)
        applydynamics!(spins[i], BtoM, Af, Bf)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M.z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3d()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [Gradient(0.1, 0, 0) for i = 1:T]
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = RF(rf, dt)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    map(spins) do spin
        excite!(Ae, Be, spin, rf, 0, grad)
        freeprecess!(Af, Bf, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, Ae, Be)
        applydynamics!(spin, BtoM, Af, Bf)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M.z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3f()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    grad = [Gradient(0.1, 0, 0) for i = 1:T]
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> Spin(1, T1, T2, Δf, pos), pos)
    rf = RF(rf, dt)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    map(spins) do spin
        excite!(Ae, Be, spin, rf, 0, grad)
        freeprecess!(Af, Bf, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, Ae, Be)
        applydynamics!(spin, BtoM, Af, Bf)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M.z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testA5bMC()

    answer = matread("matlabtestdata/testA5b.mat")

    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 10 # Hz
    M = MagnetizationMC((1.0, 0, 0), (0, 0, 0))
    dt = 1 # ms
    time = 0:dt:1000 # ms
    spin = SpinMC(M, 1, [1, 0], [T1, T1], [T2, T2], [Δf, Δf], [Inf, Inf])
    A = BlochMcConnellMatrix(2)
    B = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    mag = zeros(3,length(time))
    mag[:,1] = Vector(M[1])
    for t = 2:length(time)
        freeprecess!(A, B, spin, dt)
        applydynamics!(spin, BtoM, A, B)
        mag[:,t] = Vector(spin.M[1])
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
    grad = Gradient(0.1, 0, 0) # G/cm
    grad2 = Gradient(-0.05, 0, 0) # G/cm
    xpos = -2:0.01:2 # cm
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> SpinMC(1, [1, 0], [T1, T1], [T2, T2], [Δf, Δf], [Inf, Inf], pos), pos)
    rf = InstantaneousRF(α)
    Ae = ExcitationMatrix()
    Af = BlochMcConnellMatrix(2)
    Bf = MagnetizationMC(2)
    Ar = BlochMcConnellMatrix(2)
    Br = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(Ae, spin, rf)
        freeprecess!(Af, Bf, spin, t, grad)
        freeprecess!(Ar, Br, spin, t, grad2)
        applydynamics!(spin, BtoM, Ae)
        applydynamics!(spin, BtoM, Af, Bf)
        applydynamics!(spin, BtoM, Ae)
        applydynamics!(spin, BtoM, Ar, Br)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF2cMC()

    answer = matread("matlabtestdata/testF2c.mat")

    dt = 0.04 # ms
    t = 0:dt:6 # ms
    rf = 0.05 * sinc.(t .- 3) # G
    Δf = -1000:20:1000 # Hz
    grad = [Gradient(0, 0, 0) for r in rf]
    spins = map(Δf -> SpinMC(1, [1, 0], [600, 10], [100, 10], [Δf, 0], [Inf, Inf]), Δf)
    rf = RF(rf, dt)
    A = BlochMcConnellMatrix(2)
    B = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(A, B, spin, rf, 0, grad)
        applydynamics!(spin, BtoM, A, B)
    end
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
    rf = RF(rf, dt)
    rf0 = RF(zeros(length(rf)), dt)
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> SpinMC(1, [1, 0], [T1, T1], [T2, T2], [Δf, Δf], [Inf, Inf], pos), pos)
    A1 = BlochMcConnellMatrix(2)
    B1 = MagnetizationMC(2)
    A2 = BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(A1, B1, spin, rf, 0, grad1)
        excite!(A2, B2, spin, rf0, 0, grad2)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[1].z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3dMC()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    rf = RF(rf, dt)
    rf0 = RF(zeros(length(rf)), dt)
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> SpinMC(1, [1, 0], [T1, T1], [T2, T2], [Δf, Δf], [Inf, Inf], pos), pos)
    A1 = BlochMcConnellMatrix(2)
    B1 = MagnetizationMC(2)
    A2 = BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(A1, B1, spin, rf, 0, grad1)
        excite!(A2, B2, spin, rf0, 0, grad2)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[1].z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3fMC()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    rf = RF(rf, dt)
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [Position(x, 0, 0) for x in xpos]
    spins = map(pos -> SpinMC(1, [1, 0], [T1, T1], [T2, T2], [Δf, Δf], [Inf, Inf], pos), pos)
    A1 = BlochMcConnellMatrix(2)
    B1 = MagnetizationMC(2)
    A2 = BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(A1, B1, spin, rf, 0, grad1)
        freeprecess!(A2, B2, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[1].z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

# ------------------------------------------------------------------------------
# End tests for comparing to Brian Hargreaves' MATLAB code
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin single compartment tests
# ------------------------------------------------------------------------------

function Spin1()

    s = Spin(Magnetization(1, 2, 3), 1, 1000, 100, 0)
    return s.M == Magnetization(1, 2, 3) &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == Position(0, 0, 0) &&
           s.signal == 1 + 2im

end

function Spin2()

    s = Spin(Magnetization(1, 2, 3), 1.5, 1000, 100, 1.3, Position(0.5, 0.2, 1))
    return s.M == Magnetization(1, 2, 3) &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == Position(0.5, 0.2, 1) &&
           s.signal == 1 + 2im

end

function Spin3()

    s = Spin(1, 1000, 100, 0)
    return s.M == Magnetization(0, 0, 1) &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == Position(0, 0, 0) &&
           s.signal == 0

end

function Spin4()

    s = Spin(1.5, 1000, 100, 1.3, Position(0.5, 0.2, 1))
    return s.M == Magnetization(0, 0, 1.5) &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == Position(0.5, 0.2, 1) &&
           s.signal == 0

end

function freeprecess1()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 3.75)
    A = FreePrecessionMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    freeprecess!(A, B, s, 100)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function freeprecess2()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 0, Position(0, 0, 3.75))
    A = FreePrecessionMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    freeprecess!(A, B, s, 100, Gradient(0, 0, 1/GAMBAR))
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function excitation1()

    s = Spin(1, 1000, 100, 3.75)
    A = ExcitationMatrix()
    BtoM = Magnetization()
    rf = InstantaneousRF(π/2, π/4)
    excite!(A, s, rf)
    applydynamics!(s, BtoM, A)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct

end

function excitation2()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    rf = RF(fill(exp(im * π/8), 2), dt)
    Δθ = π/8
    grad = [Gradient(0, 0, 0) for i = 1:length(rf)]
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    excite!(A, B, s, rf, Δθ, grad)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct && B == Magnetization(0, 0, 0)

end

function excitation3()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    rf = RF(fill(exp(im * π/8), 2), dt)
    Δθ = π/8
    grad = Gradient(0, 0, 0)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    excite!(A, B, s, rf, Δθ, grad)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct && B == Magnetization(0, 0, 0)

end

function spoil1()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct && S === BlochSim.IdealSpoilingMatrix()

end

function spoil2()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    spoil!(s)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct

end

function applydynamics1()

    s = Spin(1, 1000, 100, 3.75)
    A = ExcitationMatrix()
    BtoM = Magnetization()
    excite!(A, s, InstantaneousRF(π/2))
    applydynamics!(s, BtoM, A)
    A = FreePrecessionMatrix()
    B = Magnetization()
    freeprecess!(A, B, s, 100)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404054)
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End single compartment tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin multicompartment tests
# ------------------------------------------------------------------------------

function SpinMC1()

    s = SpinMC(MagnetizationMC((1, 2, 3), (4, 5, 6)), 2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], Position(-0.3, -1, 0))
    return s.N == 2 &&
           s.M == MagnetizationMC((1, 2, 3), (4, 5, 6)) &&
           s.Meq == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.M0 == 2 &&
           s.frac == (0.2, 0.8) &&
           s.T1 == (400, 1000) &&
           s.T2 == (20, 100) &&
           s.Δf == (15, 0) &&
           s.r == ((0, 1 / 20), (1 / 40, 0)) &&
           s.pos == Position(-0.3, -1, 0) &&
           s.signal == 5 + 7im

end

function SpinMC2()

    s = SpinMC(2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], Position(-0.3, -1, 0))
    return s.N == 2 &&
           s.M == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.Meq == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.M0 == 2 &&
           s.frac == (0.2, 0.8) &&
           s.T1 == (400, 1000) &&
           s.T2 == (20, 100) &&
           s.Δf == (15, 0) &&
           s.r == ((0, 1 / 20), (1 / 40, 0)) &&
           s.pos == Position(-0.3, -1, 0) &&
           s.signal == 0

end

function excitationMC1()

    s = SpinMC(1.5, [1/3, 2/3], [400, 1000], [20, 100], [3.75, 3.75], [20, 40])
    A = ExcitationMatrix()
    BtoM = MagnetizationMC(2)
    excite!(A, s, InstantaneousRF(π/2, π/4))
    applydynamics!(s, BtoM, A)
    M_correct = MagnetizationMC((sqrt(2)/4, -sqrt(2)/4, 0), (sqrt(2)/2, -sqrt(2)/2, 0))
    return s.M ≈ M_correct

end

function spoilMC1()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    return s.M ≈ M_correct && S === BlochSim.IdealSpoilingMatrix()

end

function spoilMC2()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    spoil!(s)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End multicompartment tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin automatic differentiation tests
# ------------------------------------------------------------------------------

function autodiff1()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = FreePrecessionMatrix{T}()
        Bf = Magnetization{T}()
        BtoM = Magnetization{T}()
        s = Spin(1, T1, T2, 10)
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 10)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 0.0009048374180359595]
    return grad ≈ correct

end

function autodiff2()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 10)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 0.0007660512555728833]
    return grad ≈ correct

end

function autodiff3()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 2)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 0.00016657611260161996]
    return grad ≈ correct

end

function autodiff4()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 1)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 8.414639502122038e-5]
    return grad ≈ correct

end

function autodiff5()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 0.1)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 8.491495820718459e-6]
    return grad ≈ correct

end

function autodiff6()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 0.01)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(s.signal)
    end
    correct = [0.0, 8.499149957624547e-7]
    return grad ≈ correct

end

# ------------------------------------------------------------------------------
# End automatic differentiation tests
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
        @test excitation1()
        @test excitation2()
        @test excitation3()
        @test spoil1()
        @test spoil2()
        @test applydynamics1()

    end

    @testset "Multicompartment" begin

        @test SpinMC1()
        @test SpinMC2()
        @test excitationMC1()
        @test spoilMC1()
        @test spoilMC2()

    end

    @testset "Automatic Differentiation" begin

        @test autodiff1()
        @test autodiff2()
        @test autodiff3()
        @test autodiff4()
        @test autodiff5()
        @test autodiff6()

    end

end
