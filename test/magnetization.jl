function magnetization1()

    M1 = Magnetization()
    M2 = Magnetization{Int}()
    M3 = Magnetization(0, 0, 0.0)

    @test eltype(M1) == eltype(M3) != eltype(M2)
    @test M1 == M2 == M3

end

function magnetization2()

    M = Magnetization()
    show(devnull, M)
    show(devnull, "text/plain", M)

    @test true

end

function magnetization3()

    M = Magnetization(1, 2, 3)
    Mcopy = copy(M)
    Mcopyto = Magnetization()
    copyto!(Mcopyto, M)

    @test M ≈ Mcopy ≈ Mcopyto
    @test M == Mcopy == Mcopyto
    @test signal(M) == 1 + 2im

end

function magnetization4()

    M = Magnetization()
    MInt = convert(Magnetization{Int}, M)
    MF64 = convert(Magnetization{Float64}, M)

    @test M == MInt == MF64

end

function magnetization5()

    src = [0, 0, 0]
    dst = Magnetization(3, 2, 1)
    copyto!(dst, src)

    @test Vector(dst) == src

end

function magnetization6()

    src = Magnetization(3, 2, 1)
    dst = [0, 0, 0]
    copyto!(dst, src)

    @test dst == Vector(src)

end

function magnetizationmc1()

    M1 = MagnetizationMC(3)
    M2 = MagnetizationMC{Int}(3)
    M3 = MagnetizationMC((0, 0, 0), (0, 0, 0), (0, 0, 0))

    @test eltype(M1) != eltype(M2) == eltype(M3)
    @test M1 == M2 == M3

end

function magnetizationmc2()

    M = MagnetizationMC(2)
    show(devnull, M)
    show(devnull, "text/plain", M)

    @test true

end

function magnetizationmc3()

    M = MagnetizationMC((1, 2, 3), (3, 2, 1))
    Mcopy = copy(M)
    Mcopyto = MagnetizationMC(2)
    copyto!(Mcopyto, M)

    @test M ≈ Mcopy ≈ Mcopyto
    @test M == Mcopy == Mcopyto
    @test signal(M) == 4 + 4im

end

function magnetizationmc4()

    M = MagnetizationMC(3)
    MInt = convert(MagnetizationMC{Int,3}, M)
    MF64 = convert(MagnetizationMC{Float64,3}, M)

    @test M == MInt == MF64

end

function magnetizationmc5()

    src = [0, 0, 0, 0, 0, 0]
    dst = MagnetizationMC((3, 2, 1), (4, 5, 6))
    copyto!(dst, src)

    @test Vector(dst) == src

end

function magnetizationmc6()

    src = MagnetizationMC((3, 2, 1), (4, 5, 6))
    dst = [0, 0, 0, 0, 0, 0]
    copyto!(dst, src)

    @test dst == Vector(src)

end

@testset "Magnetization" begin

    magnetization1()
    magnetization2()
    magnetization3()
    magnetization4()
    magnetization5()
    magnetization6()

end

@testset "MagnetizationMC" begin

    magnetizationmc1()
    magnetizationmc2()
    magnetizationmc3()
    magnetizationmc4()
    magnetizationmc5()
    magnetizationmc6()

end
