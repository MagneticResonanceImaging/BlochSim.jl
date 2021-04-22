function excite1()

    θ = π/4
    α = π/2
    spin = Spin(1, 1000, 100, 1.25)
    (A1,) = excite(spin, InstantaneousRF(α, θ))
    A2 = BlochSim.ExcitationMatrix()
    excite!(A2, spin, InstantaneousRF(α, θ))
    return Matrix(A1.A) == Matrix(A2.A)

end

function excite2()

    θ = π/4
    α = π/2
    spin = SpinMC(1, (0.7, 0.2, 0.1), (1000, 400, 1000), (100, 20, 0.01), (0, 15, 0), (100, 100, 25, Inf, Inf, Inf))
    (A1,) = excite(spin, InstantaneousRF(α, θ))
    A2 = BlochSim.ExcitationMatrix()
    excite!(A2, spin, InstantaneousRF(α, θ))
    return Matrix(A1.A) == Matrix(A2.A)

end

function excite3()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = [zeros(1, length(rf)); zeros(1, length(rf)); (0:0.5:6)']
    dt = 0.1
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = excite(spin, RF(rf, dt, Δθ, [Gradient(grad[:,i]...) for i = 1:length(rf)]))
    A2 = BlochSim.BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, RF(rf, dt, Δθ, [Gradient(grad[:,i]...) for i = 1:length(rf)]))
    return Matrix(A1) == Matrix(A2) && B1 == B2

end

function excite4()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = [0, 0, 0.5]
    dt = 0.1
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = excite(spin, RF(rf, dt, Δθ, Gradient(grad...)))
    A2 = BlochSim.BlochMatrix()
    B2 = Magnetization()
    excite!(A2, B2, spin, RF(rf, dt, Δθ, Gradient(grad...)))
    return Matrix(A1) == Matrix(A2) && B1 == B2

end

function excite5()

    rf = sinc.(-3:0.5:3)
    Δθ = π/3
    grad = [zeros(1, length(rf)); zeros(1, length(rf)); (0:0.5:6)']
    dt = 0.1
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (100, 25), Position(0, 0, 1))
    (A1, B1) = excite(spin, RF(rf, dt, Δθ, [Gradient(grad[:,i]...) for i = 1:length(rf)]))
    A2 = BlochSim.BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    excite!(A2, B2, spin, RF(rf, dt, Δθ, [Gradient(grad[:,i]...) for i = 1:length(rf)]))
    return Matrix(A1) == Matrix(A2) && B1 == B2

end

@testset "Excitation" begin

    @test excite1()
    @test excite2()
    @test excite3()
    @test excite4()
    @test excite5()

end
