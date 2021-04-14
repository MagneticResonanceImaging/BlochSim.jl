function PadeApproximantOfDegree1(m)

    T11 = 1
    T12 = 0.4
    T21 = 0.08
    T22 = 0.02
    τ12 = 0.05
    τ21 = 0.1

    expA = BlochMcConnellMatrix(2)
    A = BlochSim.BlochMcConnellDynamicsMatrix(2)
    A.A[1].R2 = -1 / T21 - 1 / τ12
    A.A[1].Δω = 2π
    A.A[1].R1 = -1 / T11 - 1 / τ12
    A.A[2].R2 = -1 / T22 - 1 / τ21
    A.A[2].Δω = 2π
    A.A[2].R1 = -1 / T12 - 1 / τ21
    A.E[1].r = 1 / τ12
    A.E[2].r = 1 / τ21
    workspace = BlochSim.MatrixExponentialWorkspace{Float64}(2)

    correct = BlochSim.PadeApproximantOfDegree(Matrix(A), m)
    BlochSim.PadeApproximantOfDegree!(expA, A, workspace, m)

    return Matrix(expA) ≈ correct

end

function expm1()

    T11 = 1
    T12 = 0.4
    T21 = 0.08
    T22 = 0.02
    τ12 = 0.05
    τ21 = 0.1

    expA = BlochMcConnellMatrix(2)
    A = BlochSim.BlochMcConnellDynamicsMatrix(2)
    A.A[1].R2 = -1 / T21 - 1 / τ12
    A.A[1].Δω = 2π
    A.A[1].R1 = -1 / T11 - 1 / τ12
    A.A[2].R2 = -1 / T22 - 1 / τ21
    A.A[2].Δω = 2π
    A.A[2].R1 = -1 / T12 - 1 / τ21
    A.E[1].r = 1 / τ12
    A.E[2].r = 1 / τ21

    correct = exp(Matrix(A))
    BlochSim.expm!(expA, A)

    return Matrix(expA) ≈ correct

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

    @testset "expm Accuracy" begin

        @test PadeApproximantOfDegree1(9)
        @test PadeApproximantOfDegree1(13)
        @test expm1()

    end

    @testset "frexp Derivatives" begin

        @test test_derivative_frexp1()
        @test test_derivative_frexp2()

    end

end
