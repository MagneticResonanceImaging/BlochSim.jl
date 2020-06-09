function relax1()

    (A, B) = BlochSim.relax(100, 1, 1000, 100)
    M = A * [1, 0, 0] + B
    M_correct = [0.36787944117144233, 0.0, 0.09516258196404048]
    return M ≈ M_correct

end

function rotatex1()

    R = BlochSim.rotatex(π/2)
    M = R * [0, 0, 1]
    M_correct = [0, 1, 0]
    return M ≈ M_correct

end

function rotatey1()

    R = BlochSim.rotatey(π/2)
    M = R * [0, 0, 1]
    M_correct = [-1, 0, 0]
    return M ≈ M_correct

end

function rotatez1()

    R = BlochSim.rotatez(π/2)
    M = R * [0, 1, 0]
    M_correct = [1, 0, 0]
    return M ≈ M_correct

end

function rotatetheta1()

    R = BlochSim.rotatetheta(π/4, π/2)
    M = R * [0, 0, 1]
    M_correct = [sqrt(2)/2, -sqrt(2)/2, 0]
    return M ≈ M_correct

end

function freeprecess1()

    (A, B) = BlochSim.freeprecess(100, 1, 1000, 100, 3.75)
    M = A * [1, 0, 0] + B
    M_correct = [-0.2601300475114444, -0.2601300475114445, 0.09516258196404048]
    return M ≈ M_correct

end

@testset "Helper Functions" begin

    @test relax1()
    @test rotatex1()
    @test rotatey1()
    @test rotatez1()
    @test rotatetheta1()
    @test freeprecess1()

end
