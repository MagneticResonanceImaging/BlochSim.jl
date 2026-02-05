using Test: @inferred, @test, @testset

function mese1()

    rfex = RF([π/2], [0], 0, 0, Gradient(0, 0, 0))
    rfref = RF([π], [-π/2], 0, 0, Gradient(0, 0, 0))
    rephaser1 = GradientSpoiling(0, 0, 0, 0)
    rephaser2 = GradientSpoiling([0], [0], [0], 0)
    crusher = GradientSpoiling(0, 0, 0, 0)
    spoiler = GradientSpoiling(0, 0, 0, 0)
    mese1! = @inferred MESEBlochSim(3000, 10, 32, rfex, rfref, rephaser1, crusher, spoiler)
    mese2! = @inferred MESEBlochSim(3000, 10, 32, rfex, rfref, rephaser2, crusher, nothing)
    s1 = Spin(1, 1000, 100, 0)
    s2 = SpinMC(1, (1, 0), (1000, 100), (100, 10), (0, 0), (Inf, Inf))
    M1 = mese1!(s1)
    M2 = mese1!(s2)
    M3 = mese2!(s2)

    show(devnull, mese1!)
    show(devnull, "text/plain", mese1!)
    return all(M1[i] ≈ M2[i][1] ≈ M3[i][1] for i = 1:32)

end

@testset "MESE" begin

    @test mese1()

end
