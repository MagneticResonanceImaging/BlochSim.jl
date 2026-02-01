using BlochSim: rotatetheta

function excite1()

    θ = π/4
    α = π/2
    rf = InstantaneousRF(α, θ)
    spin = Spin(1, 1000, 100, 1.25)
    (A1,) = excite(spin, rf)
    A2 = ExcitationMatrix()
    excite!(A2, spin, rf)

    show(devnull, rf)
    show(devnull, "text/plain", rf)
    return A1.A == A2.A

end

function excite2()

    s = Spin(1, 1000, 100, 3.75)
    rf = InstantaneousRF(π/2, π/4)
    excite!(s, rf)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct

end

function excite3()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = [Gradient(0, 0, z) for z = 0:0.5:6]
    dt = 0.1
    rf = RF(rf, dt, Δθ, grad)
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = excite(spin, rf)
    A2 = BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, rf)

    show(devnull, rf)
    show(devnull, "text/plain", rf)
    @test A1 == A2
    return B1 == B2

end

function excite4()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = Gradient(0, 0, 0.5)
    dt = 0.1
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = excite(spin, RF(rf, dt, Δθ, grad))
    A2 = BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, RF(rf, dt, Δθ, grad))
    @test A1 == A2
    return B1 == B2

end

function excite5()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    grad = [Gradient(0, 0, 0) for i = 1:2]
    rf = RF(fill(exp(im * π/8), length(grad)), dt, Δθ, grad)
    excite!(s, rf)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)
    return s.M ≈ M_correct

end

function excite6()

    s = Spin(1, Inf, Inf, 0)
    dt = 250π / GAMMA
    Δθ = π/8
    rf = RF(fill(exp(im * π/8), 2), dt, Δθ)
    excite!(s, rf)
    M_correct = Magnetization(sqrt(2)/2, -sqrt(2)/2, 0)

    @test duration(rf) == 2dt
    return s.M ≈ M_correct

end

function excitemc1()

    θ = π/4
    α = π/2
    spin = SpinMC(1, (0.7, 0.2, 0.1), (1000, 400, 1000), (100, 20, 0.01), (0, 15, 0), (100, 100, 25, Inf, Inf, Inf))
    (A1,) = excite(spin, InstantaneousRF(α, θ))
    A2 = ExcitationMatrix()
    excite!(A2, spin, InstantaneousRF(α, θ))
    return A1.A == A2.A

end

function excitemc2()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = [Gradient(0, 0, z) for z = 0:0.5:6]
    dt = 0.1
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (100, 25), Position(0, 0, 1))
    (A1, B1) = excite(spin, RF(rf, dt, Δθ, grad))
    A2 = BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    excite!(A2, B2, spin, RF(rf, dt, Δθ, grad))
    @test A1 == A2
    return B1 == B2

end

function excitemc3()

    s = SpinMC(1.5, [1/3, 2/3], [400, 1000], [20, 100], [3.75, 3.75], [20, 40])
    excite!(s, InstantaneousRF(π/2, π/4))
    M_correct = MagnetizationMC((sqrt(2)/4, -sqrt(2)/4, 0), (sqrt(2)/2, -sqrt(2)/2, 0))
    return s.M ≈ M_correct

end

@testset "Excitation" begin

    @test excite1()
    @test excite2()
    @test excite3()
    @test excite4()
    @test excite5()
    @test excite6()
    @test excitemc1()
    @test excitemc2()
    @test excitemc3()
    @test duration(InstantaneousRF(π/2,0)) == 0

    @test Matrix(rotatetheta()) ≈ [0 0 1; 0 1 0; -1 0 0]

end
