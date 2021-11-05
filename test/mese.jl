function mese1()

    rfex = RF([π/2], [0], 0, 0, Gradient(0, 0, 0))
    rfref = RF([π], [-π/2], 0, 0, Gradient(0, 0, 0))
    mese! = MESEBlochSim(3000, 10, 32, rfex, rfref)
    s1 = Spin(1, 1000, 100, 0)
    s2 = SpinMC(1, (1, 0), (1000, 100), (100, 10), (0, 0), (Inf, Inf))
    M1 = mese!(s1)
    M2 = mese!(s2)

    show(devnull, mese!)
    show(devnull, "text/plain", mese!)
    return all(M1[i] ≈ M2[i][1] for i = 1:32)

end

@testset "MESE" begin

    @test mese1()

end
