function GradientSpoiling1()

    GradientSpoiling(0, 1, 0, 1)
    GradientSpoiling(0.0, 1, 0f0, 2//3)
    return true

end

function GradientSpoiling2()

    grad = Gradient(0, 0, 0)
    spoil = GradientSpoiling(grad, 3.0)
    return spoiler_gradient(spoil) === grad

end

function RFSpoiling1()

    RFSpoiling(deg2rad(117))
    return true

end

function RFSpoiling2()

    Δθ = 1
    spoil = RFSpoiling(Δθ)
    return rfspoiling_increment(spoil) == Δθ

end

function RFandGradientSpoiling1()

    RFandGradientSpoiling(GradientSpoiling(0, 0, 0, 3), RFSpoiling(deg2rad(117)))
    RFandGradientSpoiling(Gradient(0, 0, 0), 3, RFSpoiling(deg2rad(117)))
    RFandGradientSpoiling((0, 0, 0), 3, RFSpoiling(deg2rad(117)))
    RFandGradientSpoiling(0, 0, 0, 3, RFSpoiling(deg2rad(117)))
    RFandGradientSpoiling(GradientSpoiling(0, 0, 0, 3), deg2rad(117))
    RFandGradientSpoiling(Gradient(0, 0, 0), 3, deg2rad(117))
    RFandGradientSpoiling((0, 0, 0), 3, deg2rad(117))
    RFandGradientSpoiling(RFSpoiling(deg2rad(117)), GradientSpoiling(0, 0, 0, 3))
    RFandGradientSpoiling(RFSpoiling(deg2rad(117)), Gradient(0, 0, 0), 3)
    RFandGradientSpoiling(RFSpoiling(deg2rad(117)), (0, 0, 0), 3)
    RFandGradientSpoiling(RFSpoiling(deg2rad(117)), 0, 0, 0, 3)
    RFandGradientSpoiling(deg2rad(117), GradientSpoiling(0, 0, 0, 3))
    RFandGradientSpoiling(deg2rad(117), Gradient(0, 0, 0), 3)
    RFandGradientSpoiling(deg2rad(117), (0, 0, 0), 3)
    return true

end

function RFandGradientSpoiling2()

    grad = Gradient(1.0, 2.0, 1//2)
    Δθ = 2f0
    spoil = RFandGradientSpoiling(grad, 1.0, Δθ)
    @test spoiler_gradient(spoil) === grad
    return rfspoiling_increment(spoil) == Δθ

end

function spoil1()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = Magnetization(0, 0, 5)
    @test s.M ≈ M_correct
    return S === BlochSim.IdealSpoilingMatrix()

end

function spoil2()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    spoil!(s)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct

end

function spoilmc1()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    @test s.M ≈ M_correct
    return S === BlochSim.IdealSpoilingMatrix()

end

function spoilmc2()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    spoil!(s)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    return s.M ≈ M_correct

end

@testset "AbstractSpoiling" begin

    @testset "GradientSpoiling" begin

        @test GradientSpoiling1()
        @test GradientSpoiling2()

    end

    @testset "RFSpoiling" begin

        @test RFSpoiling1()
        @test RFSpoiling2()

    end

    @testset "RFandGradientSpoiling" begin

        @test RFandGradientSpoiling1()
        @test RFandGradientSpoiling2()

    end

end

@testset "Spoiling" begin

    @test spoil1()
    @test spoil2()
    @test spoilmc1()
    @test spoilmc2()

end
