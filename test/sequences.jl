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

@testset "Sequences" begin

  @testset "MESE" begin

    @test testB2c()
    @test testB2d()

  end

end
