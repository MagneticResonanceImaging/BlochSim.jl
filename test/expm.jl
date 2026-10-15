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

function dfrexp1()

    f = x -> 2x * BlochSim.frexp1(x)
    x = LinRange(-5, 5, 1000)
    df_correct = 2BlochSim.frexp1.(x) .+ 2x .* BlochSim.dfrexp1.(x)
    df_forwarddiff = ForwardDiff.derivative.(f, x)

    return df_forwarddiff ≈ df_correct

end

function dfrexp2()

    f = x -> 2x * BlochSim.frexp2(x)
    x = LinRange(-5, 5, 1000)
    df_correct = 2BlochSim.frexp2.(x) .+ 2x .* BlochSim.dfrexp2.(x)
    df_forwarddiff = ForwardDiff.derivative.(f, x)

    return df_forwarddiff ≈ df_correct

end

function autodiff1()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = FreePrecessionMatrix{T}()
        Bf = Magnetization{T}()
        BtoM = Magnetization{T}()
        s = Spin(1, T1, T2, 10)
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 10)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 0.0009048374180359595]

    return grad ≈ correct

end

function autodiff2()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 10)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 0.0007660512555728833]

    return grad ≈ correct

end

function autodiff3()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 2)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 0.00016657611260161996]

    return grad ≈ correct

end

function autodiff4()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 1)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 8.414639502122038e-5]

    return grad ≈ correct

end

function autodiff5()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 0.1)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 8.491495820718459e-6]

    return grad ≈ correct

end

function autodiff6()

    Ae = ExcitationMatrix()
    rf = InstantaneousRF(π/2)
    grad = ForwardDiff.gradient([1000.0, 100.0]) do x
        T1, T2 = x
        T = typeof(T1)
        Af = BlochMcConnellMatrix{T}(2)
        Bf = MagnetizationMC{T}(2)
        BtoM = MagnetizationMC{T}(2)
        s = SpinMC(1, [0.15, 0.85], [400, T1], [20, T2], [25, 10], [Inf, Inf])
        excite!(Ae, s, rf)
        applydynamics!(s, BtoM, Ae)
        freeprecess!(Af, Bf, s, 0.01)
        applydynamics!(s, BtoM, Af, Bf)
        abs.(signal(s))
    end
    correct = [0.0, 8.499149957624547e-7]

    return grad ≈ correct

end

@testset "Matrix Exponential" begin

    @testset "expm Accuracy" begin

        @test expm1()

    end

    @testset "frexp Derivatives" begin

        @test dfrexp1()
        @test dfrexp2()

    end

    @testset "AutoDiff" begin

        @test autodiff1()
        @test autodiff2()
        @test autodiff3()
        @test autodiff4()
        @test autodiff5()
        @test autodiff6()

    end

end
