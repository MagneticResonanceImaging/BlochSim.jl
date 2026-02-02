using Test: @inferred, @test, @testset

function freeprecess1()

    t = 100
    spin = Spin(1, 1000, 100, 1.25)
    (A1, B1) = @inferred freeprecess(spin, t)
    A2 = @inferred FreePrecessionMatrix()
    B2 = Magnetization()
    freeprecess!(A2, B2, spin, t)
    @test A1 == A2
    return B1 == B2

end

function freeprecess2()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 3.75)
    freeprecess!(s, 100)
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function freeprecess3()

    t = 100
    grad = Gradient(0, 0, 1)
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = @inferred freeprecess(spin, t, grad)
    A2 = FreePrecessionMatrix()
    B2 = Magnetization(0.0, 0.0, 0.0)
    freeprecess!(A2, B2, spin, t, grad)
    @test A1 == A2
    return B1 == B2

end

function freeprecess4()

    s = Spin(Magnetization(1, 0, 0), 1, 1000, 100, 0, Position(0, 0, 3.75))
    freeprecess!(s, 100, Gradient(0, 0, 1/GAMBAR))
    M_correct = Magnetization(-0.2601300475114444, -0.2601300475114445, 0.09516258196404048)
    return s.M ≈ M_correct

end

function freeprecess5()

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

function freeprecess6()

    s = Spin(1, 1000, 100, 3.75)
    (A1, B1) = @inferred freeprecess(s, 100)
    (A2, B2) = @inferred freeprecess(100, s.M0, s.T1, s.T2, s.Δf)

    @test A1 == A2
    return B1 == B2

end

function freeprecess7()

    t = 100
    grads = (Gradient(0, 0, 1), Gradient(0, 0, 1))
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = @inferred freeprecess(spin, t, grads)
    (A2, B2) = @inferred freeprecess(spin, t, grads[1])
    @test A1 ≈ A2
    return B1 ≈ B2

end

function freeprecessmc1()

    t = 20
    spin = SpinMC(1, (0.7, 0.2, 0.1), (1000, 400, 1000), (100, 20, 0.01), (0, 15, 0), (100, 100, 25, Inf, Inf, Inf))
    (A1, B1) = @inferred freeprecess(spin, t)
    A2 = BlochMcConnellMatrix(3)
    B2 = MagnetizationMC(3)
    freeprecess!(A2, B2, spin, t)
    @test A1 == A2
    return B1 == B2

end

function freeprecessmc2()

    t = 20
    M0 = 1
    (fa, fb) = (0.5, 0.5)
    (T1a, T2a, Δfa) = (1000, 100, 0)
    (T1b, T2b, Δfb) = (400, 20, 15)
    spina = Spin(fa * M0, T1a, T2a, Δfa)
    spinb = Spin(fb * M0, T1b, T2b, Δfb)
    spinmc = SpinMC(M0, (fa, fb), (T1a, T1b), (T2a, T2b), (Δfa, Δfb), (Inf, Inf))
    Aa = FreePrecessionMatrix()
    Ba = Magnetization()
    Ab = FreePrecessionMatrix()
    Bb = Magnetization()
    Amc = BlochMcConnellMatrix(2)
    Bmc = MagnetizationMC(2)
    freeprecess!(Aa, Ba, spina, t)
    freeprecess!(Ab, Bb, spinb, t)
    freeprecess!(Amc, Bmc, spinmc, t)
    @test Matrix(Amc) ≈ [Matrix(Aa) zeros(3, 3); zeros(3, 3) Matrix(Ab)]
    return Bmc ≈ MagnetizationMC(Ba, Bb)

end

function freeprecessmc3()

    t = 20
    grad = Gradient(0, 0, 1)
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (100, 25), Position(0, 0, 1))
    (A1, B1) = @inferred freeprecess(spin, t, grad)
    A2 = BlochMcConnellMatrix(2)
    B2 = MagnetizationMC(2)
    freeprecess!(A2, B2, spin, t, grad)
    @test A1 ≈ A2
    return B1 ≈ B2

end

function freeprecessmc4()

    s = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (Inf, Inf))
    t = 100
    (A1, B1) = @inferred freeprecess(s, t)
    (A2, B2) = @inferred freeprecess(s, t, nothing)

    @test A1 ≈ A2
    return B1 ≈ B2

end

function freeprecessmc5()

    t = 100
    grads = (Gradient(0, 0, 1), Gradient(0, 0, 1))
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (10, 100))
    (A1, B1) = @inferred freeprecess(spin, t, grads)
    (A2, B2) = @inferred freeprecess(spin, t, grads[1])
    @test A1 ≈ A2
    return B1 ≈ B2

end

function combine1()

    A = BlochMatrix()
    B = Magnetization()
    A1 = BlochMatrix(1, 2, 3, 4, 5, 6, 7, 8, 9)
    B1 = nothing
    A2 = BlochMatrix(1, 2, 3, 4, 5, 6, 7, 8, 9)
    B2 = Magnetization(1, 1, 1)
    combine!(A, B, A1, B1, A2, B2)
    Ac = A2 * A1
    Bc = B2

    @test A == Ac
    return B == Bc

end

function combine2()

    A = BlochMatrix()
    B = Magnetization()
    A1 = BlochMatrix(1, 2, 3, 4, 5, 6, 7, 8, 9)
    B1 = Magnetization(1, 1, 1)
    A2 = nothing
    B2 = nothing
    combine!(A, B, A1, B1, A2, B2)
    Ac = A1
    Bc = B1

    @test A == Ac
    return B == Bc

end

function combine3()

    A = BlochMatrix()
    B = Magnetization()
    A1 = nothing
    B1 = nothing
    A2 = BlochMatrix(1, 2, 3, 4, 5, 6, 7, 8, 9)
    B2 = Magnetization(1, 1, 1)
    combine!(A, B, A1, B1, A2, B2)
    Ac = A2
    Bc = B2

    @test A == Ac
    return B == Bc

end

function combine4()

    A = BlochMatrix()
    A1 = BlochMatrix(1, 1, 1, 1, 1, 1, 1, 1, 1)
    A2 = BlochMatrix(1, 2, 3, 4, 5, 6, 7, 8, 9)
    combine!(A, A1, A2)
    Ac = A2 * A1

    return A == Ac

end

@testset "Free Precession" begin

    @test freeprecess1()
    @test freeprecess2()
    @test freeprecess3()
    @test freeprecess4()
    @test freeprecess5()
    @test freeprecess6()
    @test freeprecess7()
    @test freeprecessmc1()
    @test freeprecessmc2()
    @test freeprecessmc3()
    @test freeprecessmc4()
    @test freeprecessmc5()

end

@testset "combine!" begin

    @test combine1()
    @test combine2()
    @test combine3()
    @test combine4()

end
