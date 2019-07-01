using BlochSim, Test, MAT

function testB2c()

    answer = matread("matlabtestdata/testB2c.mat")

    spin = Spin(1, 600, 100, 0)
    sig = mese!(spin, 1000, 50, 1)[1]

    return sig ≈ answer["sig"]

end

function testB2d()

    answer = matread("matlabtestdata/testB2d.mat")

    spin = Spin(1, 600, 100, 0)
    sig = mese!(spin, 1000, 50, 8)

    return sig ≈ vec(answer["sig"])

end

function testB3a()

    answer = matread("matlabtestdata/testB3a.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    z = π/2 / (GAMMA * gradz * Tg/1000) # cm
    spin = Spin(1, 600, 100, 0, [0,0,z])
    spgr!(spin, 10, 2, π/3, [0,0,gradz], Tg; Δθinc = 0)

    return spin.M ≈ vec(answer["M"])

end

function testB3b()

    answer = matread("matlabtestdata/testB3b.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 4π / (GAMMA * gradz * Tg/1000) # cm
    z = ((1:100)/100 .- 0.5) * zmax
    spins = map(z -> Spin(1, 600, 100, 0, [0,0,z]), z)
    map(spin -> spgr!(spin, 10, 2, π/3, [0,0,gradz], Tg; Δθinc = 0), spins)
    M = zeros(ComplexF64, 3)
    for spin in spins
        M .+= spin.M
    end
    M ./= 100

    return M ≈ vec(answer["M"])

end

function testB3c()

    answer = matread("matlabtestdata/testB3c.mat")

    spin = Spin(1, 600, 100, 0)
    spgr!(spin, 10, 2, π/3)

    return spin.M ≈ vec(answer["M"])

end

function testB5a()

    answer = matread("matlabtestdata/testB5a.mat")

    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:100)/100 * zmax
    spins = map(z -> Spin(1, 600, 100, 0, [0,0,z]), z)
    s = map(spin -> spgr!(spin, 10, 2, π/6, [0,0,gradz], Tg, nTR = 99), spins)
    sig = sum(s) / 100

    return sig ≈ answer["sig"]

end

function testB5b()

    answer = matread("matlabtestdata/testB5b.mat")

    α = LinRange(0, π/2, 51) # rad
    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:100)/100 * zmax
    srf = zeros(ComplexF64, 51)
    for i = 1:51
        spins = map(z -> Spin(1, 600, 100, 0, [0,0,z]), z)
        s = map(spin -> spgr!(spin, 10, 2, α[i], [0,0,gradz], Tg, nTR = 99), spins)
        srf[i] = sum(s) / 100
    end
    spin = Spin(1, 600, 100, 0)
    sideal = map(α -> spgr!(spin, 10, 2, α), α)

    return srf ≈ vec(answer["sig1"]) && sideal ≈ vec(answer["sig2"])

end

@testset "Sequences" begin

    @testset "MESE" begin

        @test testB2c()
        @test testB2d()

    end

    @testset "SPGR" begin


        @test testB3a()
        @test testB3b()
        @test testB3c()
        @test testB5a()
        @test testB5b()

    end

end
