function spgr1()

    spgr! = SPGRBlochSim(10, 5, deg2rad(20), Val(1))
    s = Spin(1, 1000, 100, 0)
    spgr!(s)

    show(devnull, spgr!)
    show(devnull, "text/plain", spgr!)
    return true

end

function spgr2()

    spgr1! = SPGRBlochSim(10, 5, RF([0.1], 0.01), RFSpoiling(), Val(2), Val(true))
    s1 = Spin(1, 1000, 100, 0)
    s2 = SpinMC(1, (1, 0), (1000, 100), (100, 10), (0, 0), (Inf, Inf))
    M11 = spgr1!(s1)
    M12 = spgr1!(s2)

    spgr2! = SPGRBlochSim(10, 5, RF([0.1], 0.01), RFandGradientSpoiling(RFSpoiling(), GradientSpoiling(0, 0, 0, 0)), Val(2), Val(true))
    s1 = Spin(1, 1000, 100, 0)
    s2 = SpinMC(1, (1, 0), (1000, 100), (100, 10), (0, 0), (Inf, Inf))
    M21 = spgr2!(s1)
    M22 = spgr2!(s2)

    return all(M11[i] ≈ M12[i][1] ≈ M21[i] ≈ M22[i][1] for i = 1:2)

end

@testset "SPGR" begin

    @test spgr1()
    @test spgr2()

end
