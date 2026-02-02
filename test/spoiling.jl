using Test: @inferred, @test, @testset

function GradientSpoiling1()

    @inferred GradientSpoiling(0, 1, 0, 1)
    g = @inferred GradientSpoiling(0.0, 1, 0f0, 2//3)

    show(devnull, g.gradient)
    show(devnull, "text/plain", g.gradient)
    show(devnull, g)
    show(devnull, "text/plain", g)
    return rfspoiling_increment(g) == 0

end

function GradientSpoiling2()

    grad = @inferred Gradient(0, 0, 0)
    spoil = @inferred GradientSpoiling(grad, 3.0)
    return @inferred spoiler_gradient(spoil) === grad

end

function RFSpoiling1()

    s = @inferred RFSpoiling(deg2rad(117))

    show(devnull, s)
    show(devnull, "text/plain", s)
    return true

end

function RFSpoiling2()

    Δθ = 1
    spoil = @inferred RFSpoiling(Δθ)
    return @inferred rfspoiling_increment(spoil) == Δθ

end

function RFandGradientSpoiling1()

    @inferred RFandGradientSpoiling(GradientSpoiling(0, 0, 0, 3), RFSpoiling(deg2rad(117)))
    s = @inferred RFandGradientSpoiling(RFSpoiling(), GradientSpoiling(0, 0, 0, 3))

    show(devnull, s)
    show(devnull, "text/plain", s)
    return true

end

function RFandGradientSpoiling2()

    grad = Gradient(1.0, 2.0, 1//2)
    Δθ = 2f0
    spoil = @inferred RFandGradientSpoiling(GradientSpoiling(grad, 1.0), RFSpoiling(Δθ))
    @test @inferred spoiler_gradient(spoil) === grad
    return @inferred rfspoiling_increment(spoil) == Δθ

end

function spoil1()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    (S,) = @inferred spoil(s)
    applydynamics!(s, S)
    M_correct = Magnetization(0, 0, 5)
    @test s.M ≈ M_correct
    return S === @inferred BlochSim.IdealSpoilingMatrix()

end

function spoil2()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    spoil!(s)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct

end

function spoil3()

    s = Spin(1, 1000, 100, 0)
    spoil!(nothing, nothing, s, RFSpoiling())
    (A,) = spoil(s, RFSpoiling())

    return A === nothing

end

function spoil4()

    s = Spin(1, 1000, 100, 0)
    (A1, B1) = @inferred spoil(s, GradientSpoiling(1, 1, 1, 10))
    (A2, B2) = freeprecess(s, 10)

    @test A1 == A2
    return B1 == B2

end

function spoilmc1()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    (S,) = @inferred spoil(s)
    applydynamics!(s, S)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    @test s.M ≈ M_correct
    return S === @inferred BlochSim.IdealSpoilingMatrix()

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
    @test spoil3()
    @test spoil4()
    @test spoilmc1()
    @test spoilmc2()

end
