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

    α = π/3 # rad
    TE = 2 # ms
    TR = 10 # ms
    gradz = 0.3 # G/cm
    Tg = 3 # ms
    zmax = 4π / (GAMMA * gradz * Tg/1000) # cm
    z = ((1:100)/100 .- 0.5) * zmax
    spins = map(z -> Spin(1, 600, 100, 0, [0,0,z]), z)
    map(spin -> spgr!(spin, 10, 2, π/3, [0,0,gradz], Tg; Δθinc = 0), spins)
    M = zeros(ComplexF64, 100)
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

@testset "Sequences" begin

  @testset "MESE" begin

    @test testB2c()
    @test testB2d()
    @test testB3a()
    @test testB3b()
    @test testB3c()

  end

end
