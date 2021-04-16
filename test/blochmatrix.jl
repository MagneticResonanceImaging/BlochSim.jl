function mul1()

    M2 = Magnetization()
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    M1 = Magnetization(1, 2, 3)
    mul!(M2, A, M1)
    correct = Matrix(A) * Vector(M1)

    return Vector(M2) == correct

end

function mul2()

    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    t = 10
    correct = Matrix(A) * t
    mul!(A, t)

    return Matrix(A) == correct

end

function mul3()

    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 2
    t = 10
    correct = Matrix(A) * t
    mul!(A, t)

    return Matrix(A) == correct

end

function mul4()

    T11 = 1
    T12 = 0.4
    T21 = 0.08
    T22 = 0.02
    τ12 = 0.05
    τ21 = 0.1

    A = BlochSim.BlochMcConnellDynamicsMatrix(2)
    A.A[1].R2 = -1 / T21 - 1 / τ12
    A.A[1].Δω = 2π
    A.A[1].R1 = -1 / T11 - 1 / τ12
    A.A[2].R2 = -1 / T22 - 1 / τ21
    A.A[2].Δω = 2π
    A.A[2].R1 = -1 / T12 - 1 / τ21
    A.E[1].r = 1 / τ12
    A.E[2].r = 1 / τ21
    t = 10
    correct = Matrix(A) * t
    mul!(A, t)

    return Matrix(A) == correct

end

function mul5()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    t = 10
    mul!(C, A, t)
    correct = Matrix(A) * t

    return Matrix(C) == correct

end

function mul6()

    C = BlochMcConnellMatrix(2)
    A = BlochMcConnellMatrix(2)
    for i = 1:2, j = 1:2
        A.A[i][j].a11 = 1
        A.A[i][j].a21 = 2
        A.A[i][j].a31 = 3
        A.A[i][j].a12 = 4
        A.A[i][j].a22 = 5
        A.A[i][j].a32 = 6
        A.A[i][j].a13 = 7
        A.A[i][j].a23 = 8
        A.A[i][j].a33 = 9
    end
    t = 10
    mul!(C, A, t)
    correct = Matrix(A) * t

    return Matrix(C) == correct

end

function mul7()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul8()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul9()

    C = BlochSim.BlochMatrix()
    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 1
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul10()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    B = BlochSim.BlochDynamicsMatrix()
    B.R1 = 10
    B.R2 = 20
    B.Δω = 30
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul11()

    C = BlochSim.BlochMatrix()
    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 2
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul12()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(A) * Matrix(B)
    mul!(C, A, B)

    return Matrix(C) == correct

end

function mul13()

    C = BlochSim.BlochMatrix()
    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(B) * Matrix(A)
    mul!(C, B, A)

    return Matrix(C) == correct

end

function mul14()

    C = BlochMcConnellMatrix(2)
    A = BlochMcConnellMatrix(2)
    for i = 1:2, j = 1:2
        A.A[i][j].a11 = 1 + i + 10j
        A.A[i][j].a21 = 2 + i + 10j
        A.A[i][j].a31 = 3 + i + 10j
        A.A[i][j].a12 = 4 + i + 10j
        A.A[i][j].a22 = 5 + i + 10j
        A.A[i][j].a32 = 6 + i + 10j
        A.A[i][j].a13 = 7 + i + 10j
        A.A[i][j].a23 = 8 + i + 10j
        A.A[i][j].a33 = 9 + i + 10j
    end
    B = BlochMcConnellMatrix(2)
    for i = 1:2, j = 1:2
        B.A[i][j].a11 = 10 + i + 10j
        B.A[i][j].a21 = 20 + i + 10j
        B.A[i][j].a31 = 30 + i + 10j
        B.A[i][j].a12 = 40 + i + 10j
        B.A[i][j].a22 = 50 + i + 10j
        B.A[i][j].a32 = 60 + i + 10j
        B.A[i][j].a13 = 70 + i + 10j
        B.A[i][j].a23 = 80 + i + 10j
        B.A[i][j].a33 = 90 + i + 10j
    end
    mul!(C, A, B)
    correct = Matrix(A) * Matrix(B)

    return Matrix(C) == correct

end

function mul15()

    # Use Int because Float64 had slight floating point errors
    C = BlochMcConnellMatrix{Int}(2)
    A = BlochSim.BlochMcConnellDynamicsMatrix{Int}(2)
    A.A[1].R2 = 1
    A.A[1].Δω = 2
    A.A[1].R1 = 3
    A.A[2].R2 = 4
    A.A[2].Δω = 5
    A.A[2].R1 = 6
    A.E[1].r = 7
    A.E[2].r = 8
    B = BlochMcConnellMatrix{Int}(2)
    for i = 1:2, j = 1:2
        B.A[i][j].a11 = 10 + i + 10j
        B.A[i][j].a21 = 20 + i + 10j
        B.A[i][j].a31 = 30 + i + 10j
        B.A[i][j].a12 = 40 + i + 10j
        B.A[i][j].a22 = 50 + i + 10j
        B.A[i][j].a32 = 60 + i + 10j
        B.A[i][j].a13 = 70 + i + 10j
        B.A[i][j].a23 = 80 + i + 10j
        B.A[i][j].a33 = 90 + i + 10j
    end
    mul!(C, A, B)
    correct = Matrix(A) * Matrix(B)

    return Matrix(C) == correct

end

function mul16()

    T11 = 1
    T12 = 0.4
    T21 = 0.08
    T22 = 0.02
    τ12 = 0.05
    τ21 = 0.1

    C = BlochMcConnellMatrix(2)
    A = BlochSim.BlochMcConnellDynamicsMatrix(2)
    A.A[1].R2 = -1 / T21 - 1 / τ12
    A.A[1].Δω = 2π
    A.A[1].R1 = -1 / T11 - 1 / τ12
    A.A[2].R2 = -1 / T22 - 1 / τ21
    A.A[2].Δω = 2π
    A.A[2].R1 = -1 / T12 - 1 / τ21
    A.E[1].r = 1 / τ12
    A.E[2].r = 1 / τ21
    mul!(C, A, A)
    correct = Matrix(A) * Matrix(A)

    return Matrix(C) == correct

end

function muladd1()

    M2 = Magnetization(1, 1, 1)
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    M1 = Magnetization(1, 2, 3)
    correct = Matrix(A) * Vector(M1) + Vector(M2)
    BlochSim.muladd!(M2, A, M1)

    return Vector(M2) == correct

end

function muladd2()

    C = BlochSim.BlochMatrix()
    C.a11 = 10
    C.a21 = 20
    C.a31 = 30
    C.a12 = 40
    C.a22 = 50
    C.a32 = 60
    C.a13 = 70
    C.a23 = 80
    C.a33 = 90
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    t = 10
    correct = Matrix(A) * t + Matrix(C)
    BlochSim.muladd!(C, A, t)

    return Matrix(C) == correct

end

function muladd3()

    C = BlochMcConnellMatrix(3)
    for i = 1:3, j = 1:3
        C.A[i][j].a11 = 10
        C.A[i][j].a21 = 20
        C.A[i][j].a31 = 30
        C.A[i][j].a12 = 40
        C.A[i][j].a22 = 50
        C.A[i][j].a32 = 60
        C.A[i][j].a13 = 70
        C.A[i][j].a23 = 80
        C.A[i][j].a33 = 90
    end
    A = BlochMcConnellMatrix(3)
    for i = 1:3, j = 1:3
        A.A[i][j].a11 = 1
        A.A[i][j].a21 = 2
        A.A[i][j].a31 = 3
        A.A[i][j].a12 = 4
        A.A[i][j].a22 = 5
        A.A[i][j].a32 = 6
        A.A[i][j].a13 = 7
        A.A[i][j].a23 = 8
        A.A[i][j].a33 = 9
    end
    t = 10
    correct = Matrix(A) * t + Matrix(C)
    BlochSim.muladd!(C, A, t)

    return Matrix(C) == correct

end

function muladd4()

    C = BlochSim.BlochMatrix()
    C.a11 = 10
    C.a21 = 20
    C.a31 = 30
    C.a12 = 40
    C.a22 = 50
    C.a32 = 60
    C.a13 = 70
    C.a23 = 80
    C.a33 = 90
    t = 10
    correct = I(3) * t + Matrix(C)
    BlochSim.muladd!(C, I, t)

    return Matrix(C) == correct

end

function muladd5()

    C = BlochMcConnellMatrix(3)
    for i = 1:3, j = 1:3
        C.A[i][j].a11 = 10
        C.A[i][j].a21 = 20
        C.A[i][j].a31 = 30
        C.A[i][j].a12 = 40
        C.A[i][j].a22 = 50
        C.A[i][j].a32 = 60
        C.A[i][j].a13 = 70
        C.A[i][j].a23 = 80
        C.A[i][j].a33 = 90
    end
    t = 10
    correct = I(9) * t + Matrix(C)
    BlochSim.muladd!(C, I, t)

    return Matrix(C) == correct

end

function muladd6()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd7()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.BlochDynamicsMatrix()
    A.R2 = 1
    A.Δω = 2
    A.R1 = 3
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd8()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 1
    B = BlochSim.BlochMatrix()
    B.a11 = 10
    B.a21 = 20
    B.a31 = 30
    B.a12 = 40
    B.a22 = 50
    B.a32 = 60
    B.a13 = 70
    B.a23 = 80
    B.a33 = 90
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd9()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.BlochDynamicsMatrix()
    A.R2 = 1
    A.Δω = 2
    A.R1 = 3
    B = BlochSim.BlochDynamicsMatrix()
    B.R2 = 10
    B.Δω = 20
    B.R1 = 30
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd10()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 1
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd11()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.BlochDynamicsMatrix()
    A.R2 = 1
    A.Δω = 2
    A.R1 = 3
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(A) * Matrix(B) + Matrix(C)
    BlochSim.muladd!(C, A, B)

    return Matrix(C) == correct

end

function muladd12()

    C = BlochSim.BlochMatrix()
    C.a11 = 1
    C.a21 = 1
    C.a31 = 1
    C.a12 = 1
    C.a22 = 1
    C.a32 = 1
    C.a13 = 1
    C.a23 = 1
    C.a33 = 1
    A = BlochSim.BlochDynamicsMatrix()
    A.R2 = 1
    A.Δω = 2
    A.R1 = 3
    B = BlochSim.ExchangeDynamicsMatrix()
    B.r = 10
    correct = Matrix(B) * Matrix(A) + Matrix(C)
    BlochSim.muladd!(C, B, A)

    return Matrix(C) == correct

end

function subtractmul1()

    M2 = Magnetization()
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    M1 = Magnetization(10.0, 20.0, 30.0)
    BlochSim.subtractmul!(M2, I, A, M1)
    correct = (I - Matrix(A)) * Vector(M1)

    return Vector(M2) == correct

end

function subtractmul2()

    M2 = MagnetizationMC(2)
    A = BlochMcConnellMatrix(2)
    for i = 1:2, j = 1:2
        A.A[i][j].a11 = 1 + i + 10j
        A.A[i][j].a21 = 2 + i + 10j
        A.A[i][j].a31 = 3 + i + 10j
        A.A[i][j].a12 = 4 + i + 10j
        A.A[i][j].a22 = 5 + i + 10j
        A.A[i][j].a32 = 6 + i + 10j
        A.A[i][j].a13 = 7 + i + 10j
        A.A[i][j].a23 = 8 + i + 10j
        A.A[i][j].a33 = 9 + i + 10j
    end
    M1 = MagnetizationMC((10.0, 20.0, 30.0), (1.0, 2.0, 3.0))
    BlochSim.subtractmul!(M2, I, A, M1)
    correct = (I - Matrix(A)) * Vector(M1)

    return Vector(M2) == correct

end

function subtractmuladd1()

    M2 = Magnetization(1.0, 1.0, 1.0)
    A = BlochSim.BlochMatrix()
    A.a11 = 1
    A.a21 = 2
    A.a31 = 3
    A.a12 = 4
    A.a22 = 5
    A.a32 = 6
    A.a13 = 7
    A.a23 = 8
    A.a33 = 9
    M1 = Magnetization(10.0, 20.0, 30.0)
    BlochSim.subtractmul!(M2, I, A, M1)
    correct = (I - Matrix(A)) * Vector(M1)

    return Vector(M2) == correct

end

function getblock1()

    A = BlochSim.BlochMcConnellDynamicsMatrix(2)
    A.A[1].R2 = 1
    A.A[1].Δω = 2
    A.A[1].R1 = 3
    A.A[2].R2 = 10
    A.A[2].Δω = 20
    A.A[2].R1 = 30
    A.E[1].r = 100
    A.E[2].r = 200
    result = [Matrix(BlochSim.getblock(A, 1, 1)) Matrix(BlochSim.getblock(A, 1, 2));
              Matrix(BlochSim.getblock(A, 2, 1)) Matrix(BlochSim.getblock(A, 2, 2))]
    correct = Matrix(A)

    return result == correct

end

function getblock2()

    A = BlochMcConnellMatrix(2)
    for i = 1:2, j = 1:2
        A.A[i][j].a11 = 10 + i + 10j
        A.A[i][j].a21 = 20 + i + 10j
        A.A[i][j].a31 = 30 + i + 10j
        A.A[i][j].a12 = 40 + i + 10j
        A.A[i][j].a22 = 50 + i + 10j
        A.A[i][j].a32 = 60 + i + 10j
        A.A[i][j].a13 = 70 + i + 10j
        A.A[i][j].a23 = 80 + i + 10j
        A.A[i][j].a33 = 90 + i + 10j
    end
    result = [Matrix(BlochSim.getblock(A, 1, 1)) Matrix(BlochSim.getblock(A, 1, 2));
              Matrix(BlochSim.getblock(A, 2, 1)) Matrix(BlochSim.getblock(A, 2, 2))]
    correct = Matrix(A)

    return result == correct

end

function absolutesum1()

    A = BlochSim.BlochDynamicsMatrix()
    A.R1 = 1
    A.R2 = 2
    A.Δω = 3
    result = BlochSim.absolutesum(A)
    correct = norm(Matrix(A), 1)

    return result == correct

end

function absolutesum2()

    A = BlochSim.ExchangeDynamicsMatrix()
    A.r = 3
    result = BlochSim.absolutesum(A)
    correct = norm(Matrix(A), 1)

    return result == correct

end

function absolutesum3()

    A = BlochSim.BlochMcConnellDynamicsMatrix{Int}(2)
    A.A[1].R2 = 1
    A.A[1].Δω = 2
    A.A[1].R1 = 3
    A.A[2].R2 = 4
    A.A[2].Δω = 5
    A.A[2].R1 = 6
    A.E[1].r = 7
    A.E[2].r = 8
    result = BlochSim.absolutesum(A)
    correct = norm(Matrix(A), 1)

    return result == correct

end

@testset "Bloch Matrices" begin

    @testset "mul!" begin

        @test mul1()
        @test mul2()
        @test mul3()
        @test mul4()
        @test mul5()
        @test mul6()
        @test mul7()
        @test mul8()
        @test mul9()
        @test mul10()
        @test mul11()
        @test mul12()
        @test mul13()
        @test mul14()
        @test mul15()
        @test mul16()

    end

    @testset "muladd!" begin

        @test muladd1()
        @test muladd2()
        @test muladd3()
        @test muladd4()
        @test muladd5()
        @test muladd6()
        @test muladd7()
        @test muladd8()
        @test muladd9()
        @test muladd10()
        @test muladd11()
        @test muladd12()

    end

    @testset "subtractmul!" begin

        @test subtractmul1()
        @test subtractmul2()

    end

    @testset "subtractmuladd!" begin

        @test subtractmuladd1()

    end

    @testset "getblock" begin

        @test getblock1()
        @test getblock2()

    end

    @testset "absolutesum" begin

        @test absolutesum1()
        @test absolutesum2()
        @test absolutesum3()

    end

end
