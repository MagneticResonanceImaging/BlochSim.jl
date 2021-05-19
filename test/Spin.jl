# ------------------------------------------------------------------------------
# Begin single compartment tests
# ------------------------------------------------------------------------------

function Spin1()

    s = Spin(Magnetization(1, 2, 3), 1, 1000, 100, 0)
    return s.M == Magnetization(1, 2, 3) &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == Position(0, 0, 0) &&
           signal(s) == 1 + 2im

end

function Spin2()

    s = Spin(Magnetization(1, 2, 3), 1.5, 1000, 100, 1.3, Position(0.5, 0.2, 1))
    return s.M == Magnetization(1, 2, 3) &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == Position(0.5, 0.2, 1) &&
           signal(s) == 1 + 2im

end

function Spin3()

    s = Spin(1, 1000, 100, 0)
    return s.M == Magnetization(0, 0, 1) &&
           s.M0 == 1 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 0 &&
           s.pos == Position(0, 0, 0) &&
           signal(s) == 0

end

function Spin4()

    s = Spin(1.5, 1000, 100, 1.3, Position(0.5, 0.2, 1))
    return s.M == Magnetization(0, 0, 1.5) &&
           s.M0 == 1.5 &&
           s.T1 == 1000 &&
           s.T2 == 100 &&
           s.Δf == 1.3 &&
           s.pos == Position(0.5, 0.2, 1) &&
           signal(s) == 0

end

function freeprecess1()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 3.75)
    A = FreePrecessionMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    freeprecess!(A, B, s, 100)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function freeprecess2()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 0, Position(0, 0, 3.75))
    A = FreePrecessionMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    freeprecess!(A, B, s, 100, Gradient(0, 0, 1/GAMBAR))
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function excitation1()

    s = Spin(1, 1000, 100, 3.75)
    A = ExcitationMatrix()
    BtoM = Magnetization()
    rf = InstantaneousRF(π/2, π/4)
    excite!(A, s, rf)
    applydynamics!(s, BtoM, A)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct

end

function excitation2()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    grad = [Gradient(0, 0, 0) for i = 1:2]
    rf = RF(fill(exp(im * π/8), length(grad)), dt, Δθ, grad)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    excite!(A, B, s, rf)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct && B == Magnetization(0, 0, 0)

end

function excitation3()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    rf = RF(fill(exp(im * π/8), 2), dt, Δθ)
    A = BlochMatrix()
    B = Magnetization()
    BtoM = Magnetization()
    excite!(A, B, s, rf)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct && B == Magnetization(0, 0, 0)

end

function spoil1()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct && S === BlochSim.IdealSpoilingMatrix()

end

function spoil2()

    s = Spin(Magnetization(1, 0.4, 5), 1, 1000, 100, 0)
    spoil!(s)
    M_correct = Magnetization(0, 0, 5)
    return s.M ≈ M_correct

end

function applydynamics1()

    s = Spin(1, 1000, 100, 3.75)
    A = ExcitationMatrix()
    BtoM = Magnetization()
    excite!(A, s, InstantaneousRF(π/2))
    applydynamics!(s, BtoM, A)
    A = FreePrecessionMatrix()
    B = Magnetization()
    freeprecess!(A, B, s, 100)
    applydynamics!(s, BtoM, A, B)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404054)
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End single compartment tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin multicompartment tests
# ------------------------------------------------------------------------------

function SpinMC1()

    s = SpinMC(MagnetizationMC((1, 2, 3), (4, 5, 6)), 2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], Position(-0.3, -1, 0))
    return s.N == 2 &&
           s.M == MagnetizationMC((1, 2, 3), (4, 5, 6)) &&
           s.Meq == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.M0 == 2 &&
           s.frac == (0.2, 0.8) &&
           s.T1 == (400, 1000) &&
           s.T2 == (20, 100) &&
           s.Δf == (15, 0) &&
           s.r == ((0, 1 / 20), (1 / 40, 0)) &&
           s.pos == Position(-0.3, -1, 0) &&
           signal(s) == 5 + 7im

end

function SpinMC2()

    s = SpinMC(2, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40], Position(-0.3, -1, 0))
    return s.N == 2 &&
           s.M == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.Meq == MagnetizationMC((0, 0, 0.4), (0, 0, 1.6)) &&
           s.M0 == 2 &&
           s.frac == (0.2, 0.8) &&
           s.T1 == (400, 1000) &&
           s.T2 == (20, 100) &&
           s.Δf == (15, 0) &&
           s.r == ((0, 1 / 20), (1 / 40, 0)) &&
           s.pos == Position(-0.3, -1, 0) &&
           signal(s) == 0

end

function excitationMC1()

    s = SpinMC(1.5, [1/3, 2/3], [400, 1000], [20, 100], [3.75, 3.75], [20, 40])
    A = ExcitationMatrix()
    BtoM = MagnetizationMC(2)
    excite!(A, s, InstantaneousRF(π/2, π/4))
    applydynamics!(s, BtoM, A)
    M_correct = MagnetizationMC((sqrt(2)/4, -sqrt(2)/4, 0), (sqrt(2)/2, -sqrt(2)/2, 0))
    return s.M ≈ M_correct

end

function spoilMC1()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    S = spoil(s)
    applydynamics!(s, S)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    return s.M ≈ M_correct && S === BlochSim.IdealSpoilingMatrix()

end

function spoilMC2()

    s = SpinMC(MagnetizationMC((1, 0.4, 5), (0.2, 10, 0.2)), 1, [0.2, 0.8], [400, 1000], [20, 100], [15, 0], [20, 40])
    spoil!(s)
    M_correct = MagnetizationMC((0, 0, 5), (0, 0, 0.2))
    return s.M ≈ M_correct

end

# ------------------------------------------------------------------------------
# End multicompartment tests
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Begin automatic differentiation tests
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# End automatic differentiation tests
# ------------------------------------------------------------------------------

@testset "Spin" begin

    @testset "Compare to MATLAB" begin

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

    @testset "Single Compartment" begin

        @test Spin1()
        @test Spin2()
        @test Spin3()
        @test Spin4()
        @test freeprecess1()
        @test freeprecess2()
        @test excitation1()
        @test excitation2()
        @test excitation3()
        @test spoil1()
        @test spoil2()
        @test applydynamics1()

    end

    @testset "Multicompartment" begin

        @test SpinMC1()
        @test SpinMC2()
        @test excitationMC1()
        @test spoilMC1()
        @test spoilMC2()

    end

    @testset "Automatic Differentiation" begin

        @test autodiff1()
        @test autodiff2()
        @test autodiff3()
        @test autodiff4()
        @test autodiff5()
        @test autodiff6()

    end

end
