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
    sig = map(signal, spins)

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
    sig = map(signal, spins)

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
    sig = map(signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF2c()

    answer = matread("matlabtestdata/testF2c.mat")

    dt = 0.04 # ms
    t = 0:dt:6 # ms
    rf = 0.05 * sinc.(t .- 3) # G
    Δf = -1000:20:1000 # Hz
    spins = map(Δf -> Spin(1, 600, 100, Δf), Δf)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    rf = RF(rf, dt)
    map(spins) do spin
        excite!(A, B, spin, rf)
        applydynamics!(spin, BtoM, A, B)
    end
    sig = map(signal, spins)

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
    rf = RF(rf, dt, grad)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(A, B, spins[i], rf)
        applydynamics!(spins[i], BtoM, A, B)
    end
    sig = map(signal, spins)

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
    rf = RF(rf, dt, grad)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    for i = 1:length(xpos)
        excite!(Ae, Be, spins[i], rf)
        freeprecess!(Af, Bf, spins[i], t[end] - dt, grad2)
        applydynamics!(spins[i], BtoM, Ae, Be)
        applydynamics!(spins[i], BtoM, Af, Bf)
    end
    sig = map(signal, spins)
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
    rf = RF(rf, dt, grad)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    map(spins) do spin
        excite!(Ae, Be, spin, rf)
        freeprecess!(Af, Bf, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, Ae, Be)
        applydynamics!(spin, BtoM, Af, Bf)
    end
    sig = map(signal, spins)
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
    rf = RF(rf, dt, grad)
    Ae = BlochMatrix()
    Be = Magnetization()
    Af = FreePrecessionMatrix()
    Bf = Magnetization()
    BtoM = Magnetization()
    map(spins) do spin
        excite!(Ae, Be, spin, rf)
        freeprecess!(Af, Bf, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, Ae, Be)
        applydynamics!(spin, BtoM, Af, Bf)
    end
    sig = map(signal, spins)
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
    sig = map(signal, spins)

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
    rf = RF(rf, dt, grad)
    A = BlochMcConnellMatrix(2)
    B = MagnetizationMC(2)
    BtoM = MagnetizationMC(2)
    map(spins) do spin
        excite!(A, B, spin, rf)
        applydynamics!(spin, BtoM, A, B)
    end
    sig = map(signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3cMC()

    answer = matread("matlabtestdata/testF3c.mat")

    dt = 0.05 # ms
    t = 0.1:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    rf = RF(rf, dt, grad1)
    rf0 = RF(zeros(length(rf)), dt, grad2)
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
        excite!(A1, B1, spin, rf)
        excite!(A2, B2, spin, rf0)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(signal, spins)
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
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    rf = RF(rf, dt, grad1)
    rf0 = RF(zeros(length(rf)), dt, grad2)
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
        excite!(A1, B1, spin, rf)
        excite!(A2, B2, spin, rf0)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(signal, spins)
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
    grad1 = Gradient(0.1, 0, 0)
    grad2 = Gradient(-0.052, 0, 0)
    rf = RF(rf, dt, grad1)
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
        excite!(A1, B1, spin, rf)
        freeprecess!(A2, B2, spin, t[end] - dt, grad2)
        applydynamics!(spin, BtoM, A1, B1)
        applydynamics!(spin, BtoM, A2, B2)
    end
    sig = map(signal, spins)
    Mz = map(spin -> spin.M[1].z, spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testB2c()

    answer = matread("matlabtestdata/testB2c.mat")

    spin = Spin(1, 600, 100, 0)
    mese! = MESEBlochSim(1000, 50, 1)
    mag = mese!(spin)[1]

    return signal(mag) ≈ answer["sig"]

end

function testB2d()

    answer = matread("matlabtestdata/testB2d.mat")

    spin = Spin(1, 600, 100, 0)
    mese! = MESEBlochSim(1000, 50, 8)
    mag = mese!(spin)

    return map(signal, mag) ≈ vec(answer["sig"])

end

function testB2dMC()

    answer = matread("matlabtestdata/testB2d.mat")

    spin = SpinMC(1, [0, 1], [Inf, 600], [Inf, 100], [0, 0], [Inf, Inf])
    mese! = MESEBlochSim(1000, 50, 8)
    mag = mese!(spin)

    return map(signal, mag) ≈ vec(answer["sig"])

end

function testB3a()

    answer = matread("matlabtestdata/testB3a.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    z = π/2 / (GAMMA * gradz * Tg/1000) # cm
    spin = Spin(1, 600, 100, 0, Position(0,0,z))
    spgr! = SPGRBlochSim(10, 2, π/3, GradientSpoiling(0, 0, gradz, Tg))
    spgr!(spin)

    return Vector(spin.M) ≈ vec(answer["M"])

end

function testB3b()

    answer = matread("matlabtestdata/testB3b.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 4π / (GAMMA * gradz * Tg/1000) # cm
    z = ((1:100)/100 .- 0.5) * zmax
    spins = map(z -> Spin(1, 600, 100, 0, Position(0,0,z)), z)
    spgr! = SPGRBlochSim(10, 2, π/3, GradientSpoiling(0, 0, gradz, Tg))
    foreach(spgr!, spins)
    M = mean(spin.M for spin in spins)

    return Vector(M) ≈ vec(answer["M"])

end

function testB3c()

    answer = matread("matlabtestdata/testB3c.mat")

    spin = Spin(1, 600, 100, 0)
    spgr! = SPGRBlochSim(10, 2, π/3)
    spgr!(spin)

    return Vector(spin.M) ≈ vec(answer["M"])

end

function testB5a()

    answer = matread("matlabtestdata/testB5a.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:100)/100 * zmax
    spins = map(z -> Spin(1, 600, 100, 0, Position(0,0,z)), z)
    spgr! = SPGRBlochSim(10, 2, π/6, RFandGradientSpoiling((0, 0, gradz), Tg), Val(100))
    foreach(spgr!, spins)
    sig = mean(signal(spin.M) for spin in spins)

    return sig ≈ answer["sig"]

end

function testB5b()

    answer = matread("matlabtestdata/testB5b.mat")

    α = LinRange(0, π/2, 51) # rad
    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:100)/100 * zmax
    srf = Vector{ComplexF64}(undef, 51)
    for i = 1:51
        spins = map(z -> Spin(1, 600, 100, 0, Position(0,0,z)), z)
        spgr! = SPGRBlochSim(10, 2, α[i], RFandGradientSpoiling((0, 0, gradz), Tg), Val(100))
        foreach(spgr!, spins)
        srf[i] = mean(signal(spin.M) for spin in spins)
    end
    spin = Spin(1, 600, 100, 0)
    sideal = map(α) do α
        spgr! = SPGRBlochSim(10, 2, α)
        spgr!(spin)
        signal(spin)
    end

    return srf ≈ vec(answer["sig1"]) && sideal ≈ vec(answer["sig2"])

end

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

    @testset "Sequences" begin
    
        @testset "MESE" begin
    
            @test testB2c()
            @test testB2d()
            @test testB2dMC()
    
        end
    
        @testset "SPGR" begin
    
    
            @test testB3a()
            @test testB3b()
            @test testB3c()
            @test testB5a()
            @test testB5b()
    
        end
    
    end

end
