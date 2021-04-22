struct MESEBlochSim{T1<:AbstractRF,T2<:AbstractRF}
    TR::Float64
    TE::Float64
    nechoes::Int
    rfex::T1
    rfref::T2

    function MESEBlochSim(TR, TE, nechoes, rfex::T1, rfref::T2) where {T1,T2}

        TR >= TE * nechoes || error("TR must be greater than or equal to TE * nechoes")
        # Should probably also check to make sure TE is not during RF pulse,
        # and that excitation/refocussing pulses do not overlap
        new{T1,T2}(TR, TE, nechoes, rfex, rfref)

    end
end

MESEBlochSim(TR, TE, nechoes) = MESEBlochSim(TR, TE, nechoes, InstantaneousRF(π/2), InstantaneousRF(π, -π/2))

struct MESEBlochSimWorkspace{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12}
    Aex::T1
    Bex::T2
    Aref::T3
    Bref::T4
    Ate::T5
    Bte::T6
    Atr::T5
    Btr::T6
    Aecho::T7
    Becho::T6
    tmpA1::T7
    tmpB1::T6
    tmpA2::T7
    tmpB2::T6
    mat::T8
    vec::T9
    bm_workspace::T10
    ex_workspace::T11
    ref_workspace::T12
end

function MESEBlochSimWorkspace(
    spin::AbstractSpin,
    scan::MESEBlochSim{T1,T2},
    bm_workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
) where {T1,T2}

    T = eltype(spin)
    if T1 <: InstantaneousRF
        Aex = ExcitationMatrix{T}()
        Bex = nothing
        ex_workspace = nothing
    else
        if spin isa Spin
            Aex = BlochMatrix{T}()
            Bex = Magnetization{T}()
        else
            Aex = BlochMcConnellMatrix{T}(spin.N)
            Bex = MagnetizationMC{T}(spin.N)
        end
        ex_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if T2 <: InstantaneousRF
        Aref = ExcitationMatrix{T}()
        Bref = nothing
        ref_workspace = nothing
    else
        if spin isa Spin
            Aref = BlochMatrix{T}()
            Bref = Magnetization{T}()
        else
            Aref = BlochMcConnellMatrix{T}(spin.N)
            Bref = MagnetizationMC{T}(spin.N)
        end
        ref_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if spin isa Spin
        Ate = FreePrecessionMatrix{T}()
        Bte = Magnetization{T}()
        Atr = FreePrecessionMatrix{T}()
        Btr = Magnetization{T}()
        Aecho = BlochMatrix{T}()
        Becho = Magnetization{T}()
        tmpA1 = BlochMatrix{T}()
        tmpB1 = Magnetization{T}()
        tmpA2 = BlochMatrix{T}()
        tmpB2 = Magnetization{T}()
        mat = Matrix{T}(undef, 3, 3)
        vec = Vector{T}(undef, 3)
    else
        N = spin.N
        Ate = BlochMcConnellMatrix{T}(N)
        Bte = MagnetizationMC{T}(N)
        Atr = BlochMcConnellMatrix{T}(N)
        Btr = MagnetizationMC{T}(N)
        Aecho = BlochMcConnellMatrix{T}(N)
        Becho = MagnetizationMC{T}(N)
        tmpA1 = BlochMcConnellMatrix{T}(N)
        tmpB1 = MagnetizationMC{T}(N)
        tmpA2 = BlochMcConnellMatrix{T}(N)
        tmpB2 = MagnetizationMC{T}(N)
        mat = Matrix{T}(undef, 6, 6)
        vec = Vector{T}(undef, 6)
    end
    MESEBlochSimWorkspace(Aex, Bex, Aref, Bref, Ate, Bte, Atr, Btr, Aecho,
        Becho, tmpA1, tmpB1, tmpA2, tmpB2, mat, vec, bm_workspace, ex_workspace,
        ref_workspace)

end

# This function does not correctly handle timing with non-instantaneous RF pulses
function (scan::MESEBlochSim)(spin::AbstractSpin, workspace::MESEBlochSimWorkspace = MESEBlochSimWorkspace(spin, scan))

    excite!(workspace.Aex, workspace.Bex, spin, scan.rfex, workspace.ex_workspace)
    excite!(workspace.Aref, workspace.Bref, spin, scan.rfref, workspace.ex_workspace)
    freeprecess!(workspace.Ate, workspace.Bte, spin, scan.TE / 2, workspace.bm_workspace)
    freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TR - scan.TE * scan.nechoes, workspace.bm_workspace)
    S = spoil(spin)

    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Ate, workspace.Bte, workspace.Aref, workspace.Bref)
    combine!(workspace.Aecho, workspace.Becho, workspace.tmpA1, workspace.tmpB1, workspace.Ate, workspace.Bte)

    copyto!(workspace.tmpA1, workspace.Aecho)
    copyto!(workspace.tmpB1, workspace.Becho)
    for e = 2:scan.nechoes
        combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aecho, workspace.Becho)
        copyto!(workspace.tmpA1, workspace.tmpA2)
        copyto!(workspace.tmpB1, workspace.tmpB2)
    end
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atr, workspace.Btr)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, S, nothing)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aex, workspace.Bex)
    subtract!(workspace.mat, I, workspace.tmpA2)
    copyto!(workspace.vec, workspace.tmpB2)
    F = lu!(workspace.mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    Mout = Vector{typeof(spin.M)}(undef, scan.nechoes)
    for e = 1:scan.nechoes
        applydynamics!(spin, workspace.tmpB1, workspace.Aecho, workspace.Becho)
        Mout[e] = copy(spin.M)
    end

    return Mout

end
