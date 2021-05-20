function Spin1()

    Spin{Float64}(1, 1000, 100, 0)
    Spin{Float64}(1, Int32(1000), 100, 0)
    Spin{Float64}(1f0, 1000, 100, 0, Position(0, 1, 2))
    Spin{Float64}(1, Int32(1000), 100, 0, Position(0.0, -1.0, -1))
    Spin{Float64}(1f0, 1000f0, 100f0, 0f0)
    Spin{Float64}(Magnetization(1, 1, Int32(1)), 1, 1000, 100, 0)
    Spin(1, 1000, 100, 0)
    Spin(1, Int32(1000), 100, 0)
    Spin(1f0, 1000, 100, 0, Position(0, 1, 2))
    Spin(1, Int32(1000), 100, 0, Position(0.0, -1.0, -1))
    Spin(1f0, 1000f0, 100f0, 0f0)
    Spin(Magnetization(1, 1, Int32(1)), 1, 1000, 100, 0)
    return true

end

function SpinMC1()

    SpinMC{Float64}(1, (0.8, 0.2), (1000, 100), (100, 10), (0, 0), (100, 25))
    SpinMC{Float64}(1, (0.8, 0.2), (1000f0, 100), (100, Int16(10)), (0, 0), (100, 25))
    SpinMC{Float64}((Magnetization(0, 0, 1.0), Magnetization(0, 0, 1)), 1, (1, 1),
           (100, 100), (10.0, 10), (0, 0.0), (1f0, 1.0))
    SpinMC(1, (0.8, 0.2), (1000, 100), (100, 10), (0, 0), (100, 25))
    SpinMC(1, [0.8, 0.2], (1000f0, 100), (100, Int16(10)), (0, 0), (100, 25))
    SpinMC((Magnetization(0, 0, 1.0), Magnetization(0, 0, 1)), 1, (1, 1),
           [100, 100], (10.0, 10), [0, 0.0], (1f0, 1.0))
    return true

end

function SpinMCError1()

    SpinMC((Magnetization(0,0,1),), 1, (1,), (1000,), (100,), (0,), ())

end

function SpinMCError2()

    SpinMC(1, 1, 1000, 100, 0, ())

end

function SpinMCError3()

    SpinMC(1, (1, 1), 1000, 100, 0, (2, 3))

end

@testset "AbstractSpin" begin

    @testset "Spin" begin

        @test Spin1()

    end

    @testset "SpinMC" begin

        @test SpinMC1()
        @test_throws ErrorException SpinMCError1()
        @test_throws ErrorException SpinMCError2()
        @test_throws MethodError SpinMCError3()

    end

end
