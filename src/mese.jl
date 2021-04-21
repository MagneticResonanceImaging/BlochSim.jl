struct MESEBlochSim
    TR
    TE
    nechoes
    αex
    αref
    config
end

struct MESEBlochSimConfig
end

struct MESEBlochSimWorkspace
    Aex
    Aref
    Ate
    Bte
    Atr
    Btr
    bm_workspace
end

function (scan::MESEBlochSim)(spin::AbstractSpin, workspace)

    excite!(workspace.Aex, spin, excitation_pulse(scan))
    excite!(workspace.Aref, spin, refocussing_pulse(scan))
    freeprecess!(workspace.Ate, workspace.Bte, spin, scan.TE / 2, workspace.bm_workspace)
    freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TR - scan.TE * scan.nechoes, workspace.bm_workspace)
    S = spoil(spin)

    combine!(workspace.Atmp1, workspace.Btmp1, workspace.Ate, workspace.Bte, workspace.Aref, nothing)
    combine!(workspace.Aecho, workspace.Becho, workspace.Atmp1, workspace.Btmp1, workspace.Ate, workspace.Bte)

    copyto!(workspace.Atmp1, workspace.Aecho)
    copyto!(workspace.Btmp1, workspace.Becho)
    for e = 2:scan.nechoes
        combine!(workspace.Atmp2, workspace.Btmp2, workspace.Atmp1, workspace.Btmp1, workspace.Aecho, workspace.Becho)
        copyto!(workspace.Atmp1, workspace.Atmp2)
        copyto!(workspace.Btmp1, workspace.Btmp2)
    end
    combine!(workspace.Atmp2, workspace.Btmp2, workspace.Atmp1, workspace.Btmp1, workspace.Atr, workspace.Btr)
    combine!(workspace.Atmp1, workspace.Btmp1, workspace.Atmp2, workspace.Btmp2, S, nothing)
    combine!(workspace.Atmp2, workspace.vec, workspace.Atmp1, workspace.Btmp1, workspace.Aex, nothing)
    subtract!(workspace.mat, I, workspace.Atmp2)
    F = lu!(mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    Mout = Vector{typeof(spin.M)}(undef, nechoes)
    for e = 1:nechoes
        applydynamics!(spin, workspace.Aecho, workspace.Becho)
        Mout[e] = copy(spin.M)
    end

    return Mout

end
