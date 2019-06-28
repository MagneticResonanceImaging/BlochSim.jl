using BlochSim, Test, MAT

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
    spins = map(Δf -> Spin([0,0,1.0], 1, Inf, Inf, Δf, [0,0,0]), Δf)
    map(spin -> excitation!(spin, rf, 0, zeros(3,length(rf)), dt), spins)
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
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    # map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    for i = 1:length(xpos)
      (Ae, Be) = excitation(spins[i], rf, 0, grad, dt)
      applydynamics!(spins[i], Ae, Be)
    end
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3c()

    answer = matread("matlabtestdata/testF3c.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    # map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    for i = 1:length(xpos)
      (Ae, Be) = excitation(spins[i], rf, 0, grad[:,1], dt)
      applydynamics!(spins[i], Ae, Be)
    end
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3d()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3f()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> Spin([0,0,1.0], 1, T1, T2, Δf, pos), pos)
    map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
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
    spins = map(Δf -> SpinMC([0,0,1.0], 1, [1], [Inf], [Inf], [Δf], Vector{Int}(), [0,0,0]), Δf)
    map(spin -> excitation!(spin, rf, 0, zeros(3,length(rf)), dt), spins)
    sig = map(spin -> spin.signal, spins)

    return sig ≈ vec(answer["sig"])

end

function testF3cMC()

    answer = matread("matlabtestdata/testF3c.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3dMC()

    answer = matread("matlabtestdata/testF3d.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -2:0.1:2 # cm
    rf = 0.05 * sinc.(t .- 3) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 100 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

function testF3fMC()

    answer = matread("matlabtestdata/testF3f.mat")

    dt = 0.1 # ms
    t = 0:dt:6 # ms
    T = length(t)
    xpos = -5:0.1:5 # cm
    rf = 0.05 * sinc.(t .- 3) .* (exp.(im * 2π * 900 * t/1000) + exp.(-im * 2π * 900 * t/1000)) # G
    grad = [0.1ones(T)'; zeros(T)'; zeros(T)']
    T1 = 600 # ms
    T2 = 100 # ms
    Δf = 0 # Hz
    pos = [[x, 0, 0] for x in xpos]
    spins = map(pos -> SpinMC([0,0,1.0], 1, [1], [T1], [T2], [Δf], Vector{Int}(), pos), pos)
    map(spin -> excitation!(spin, rf, 0, grad, dt), spins)
    map(spin -> freeprecess!(spin, t[end] - t[1], -0.52grad[:,1]), spins)
    sig = map(spin -> spin.signal, spins)
    Mz = map(spin -> spin.M[3], spins)

    return sig ≈ vec(answer["sig"]) && Mz ≈ vec(answer["mz"])

end

@testset "Spin" begin

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
