function compare_to_exp1()

    A = randn(200, 200)
    correct = exp(A)
    result = BlochSim.expm(A)
    return result ≈ correct

end

function compare_to_exp2()

    T11 = 1
    T12 = 0.4
    T21 = 0.08
    T22 = 0.02
    τ12 = 0.05
    τ21 = 0.1
    B1 = [-1/T21-1/τ12 2π 0.1; -2π -1/T21-1/τ12 0.1; -0.1 -0.1 -1/T11-1/τ12]
    B2 = [-1/T22-1/τ21 2π 0.1; -2π -1/T22-1/τ21 0.1; -0.1 -0.1 -1/T12-1/τ21]
    E12 = [1/τ12 0 0; 0 1/τ12 0; 0 0 1/τ12]
    E21 = [1/τ21 0 0; 0 1/τ21 0; 0 0 1/τ21]
    A = [B1 E21; E12 B2]
    correct = exp(A)
    result = BlochSim.expm(A)
    return result ≈ correct

end

function test_derivative_frexp1()

    f = x -> 2x * BlochSim.frexp1(x)
    x = LinRange(-5, 5, 1000)
    df_correct = 2BlochSim.frexp1.(x) .+ 2x .* BlochSim.dfrexp1.(x)
    df_forwarddiff = ForwardDiff.derivative.(f, x)
    return df_forwarddiff ≈ df_correct

end

function test_derivative_frexp2()

    f = x -> 2x * BlochSim.frexp2(x)
    x = LinRange(-5, 5, 1000)
    df_correct = 2BlochSim.frexp2.(x) .+ 2x .* BlochSim.dfrexp2.(x)
    df_forwarddiff = ForwardDiff.derivative.(f, x)
    return df_forwarddiff ≈ df_correct

end

@testset "Matrix Exponential" begin

    @testset "Compare to exp" begin

        @test compare_to_exp1()
        @test compare_to_exp2()

    end

    @testset "frexp Derivatives" begin

        @test test_derivative_frexp1()
        @test test_derivative_frexp2()

    end

end
