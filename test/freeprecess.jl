function freeprecess1()

    t = 100
    spin = Spin(1, 1000, 100, 1.25)
    (A1, B1) = freeprecess(spin, t)
    A2 = similar(A1)
    B2 = similar(B1)
    freeprecess!(A2, B2, spin, t)
    return A1 == A2 && B1 == B2

end

function freeprecess2()

    t = 20
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (100, 25))
    (A1, B1) = freeprecess(spin, t)
    A2 = similar(A1)
    B2 = zero(B1)
    freeprecess!(A2, B2, spin, t)
    return A1 == A2 && B1 == B2

end

function freeprecess3()

    t = 20
    M0 = 1
    (fa, fb) = (0.5, 0.5)
    (T1a, T2a, Δfa) = (1000, 100, 0)
    (T1b, T2b, Δfb) = (400, 20, 15)
    spina = Spin(fa * M0, T1a, T2a, Δfa)
    spinb = Spin(fb * M0, T1b, T2b, Δfb)
    spinmc = SpinMC(M0, (fa, fb), (T1a, T1b), (T2a, T2b), (Δfa, Δfb), (Inf, Inf))
    Aa = Array{Float64}(undef, 3, 3)
    Ba = Array{Float64}(undef, 3)
    Ab = Array{Float64}(undef, 3, 3)
    Bb = Array{Float64}(undef, 3)
    Amc = Array{Float64}(undef, 6, 6)
    Bmc = MagnetizationMC(Magnetization(0.0, 0.0, 0.0), Magnetization(0.0, 0.0, 0.0))
    freeprecess!(Aa, Ba, spina, t)
    freeprecess!(Ab, Bb, spinb, t)
    freeprecess!(Amc, Bmc, spinmc, t)
    return Amc ≈ [Aa zeros(3, 3); zeros(3, 3) Ab] && Bmc ≈ (Magnetization(Ba...), Magnetization(Bb...))

end

function freeprecess4()

    t = 100
    grad = Gradient(0, 0, 1)
    spin = Spin(1, 1000, 100, 1.25, Position(0, 0, 1))
    (A1, B1) = freeprecess(spin, t, [grad.x, grad.y, grad.z])
    A2 = similar(A1)
    B2 = similar(B1)
    freeprecess!(A2, B2, spin, t, grad)
    return A1 == A2 && B1 == B2

end

function freeprecess5()

    t = 20
    grad = Gradient(0, 0, 1)
    spin = SpinMC(1, (0.8, 0.2), (1000, 400), (100, 20), (0, 15), (100, 25), Position(0, 0, 1))
    (A1, B1) = freeprecess(spin, t, [grad.x, grad.y, grad.z])
    A2 = similar(A1)
    B2 = zero(B1)
    freeprecess!(A2, B2, spin, t, grad)
    return A1 ≈ A2 && B1 ≈ B2

end

@testset "Free Precession" begin

    @test freeprecess1()
    @test freeprecess2()
    @test freeprecess3()
    @test freeprecess4()
    @test freeprecess5()

end
