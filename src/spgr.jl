struct SPGRBlochSim{T1<:AbstractRF,T2<:AbstractSpoiling,T3}
    TR::Float64
    TE::Float64
    rf::T1
    spoiling::T2
    nTR::Val{T3}
    # TODO: Add save_transients, or maybe save_every (for saving only every specified number of TRs)

    function SPGRBlochSim(TR, TE, rf::T1, spoiling::T2, nTR::Val{T3}) where {T1<:AbstractRF,T2<:AbstractSpoiling,T3}

        Tg = spoiler_gradient_duration(spoiling)
        TR >= TE + Tg || error("TR must be greater than or equal to TE + Tg")
        (T3 isa Int && T3 >= 0) || error("nTR must be a nonnegative Int")
        T2 <: Union{<:RFSpoiling,<:RFandGradientSpoiling} && (T3 > 0 ||
            error("nTR must be positive when simulating RF spoiling"))
        # Should probably also check to make sure TE is not during RF pulse
        new{T1,T2,T3}(TR, TE, rf, spoiling)

    end
end

SPGRBlochSim(TR, TE, rf, spoiling::AbstractSpoiling) = SPGRBlochSim(TR, TE, rf, spoiling, Val(0))
SPGRBlochSim(TR, TE, rf, nTR::Val) = SPGRBlochSim(TR, TE, rf, IdealSpoiling(), nTR)
SPGRBlochSim(TR, TE, rf) = SPGRBlochSim(TR, TE, rf, IdealSpoiling(), Val(0))
SPGRBlochSim(TR, TE, α::Real, spoiling, nTR) = SPGRBlochSim(TR, TE, InstantaneousRF(α), spoiling, nTR)

struct SPGRBlochSimWorkspace{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}
    Aex::T1
    Bex::T2
    Atr::T3
    Btr::T4
    Atg::T5
    Btg::T6
    tmpA1::T7
    tmpB1::T4
    tmpA2::T7
    tmpB2::T4
    mat::T8
    vec::T9
    bm_workspace::T10
    ex_workspace::T11
end

function SPGRBlochSimWorkspace(
    spin::AbstractSpin,
    scan::SPGRBlochSim{T1,T2},
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
    if T2 <: IdealSpoiling
        Atg = idealspoiling
        Btg = nothing
    elseif spin isa Spin
        Atg = FreePrecessionMatrix{T}()
        Btg = Magnetization{T}()
    else
        Atg = BlochMcConnellMatrix{T}(spin.N)
        Btg = MagnetizationMC{T}(spin.N)
    end
    if spin isa Spin
        Atr = FreePrecessionMatrix{T}()
        Btr = Magnetization{T}()
        tmpA1 = BlochMatrix{T}()
        tmpB1 = Magnetization{T}()
        tmpA2 = BlochMatrix{T}()
        tmpB2 = Magnetization{T}()
        mat = Matrix{T}(undef, 3, 3)
        vec = Vector{T}(undef, 3)
    else
        N = spin.N
        Atr = BlochMcConnellMatrix{T}(N)
        Btr = MagnetizationMC{T}(N)
        tmpA1 = BlochMcConnellMatrix{T}(N)
        tmpB1 = MagnetizationMC{T}(N)
        tmpA2 = BlochMcConnellMatrix{T}(N)
        tmpB2 = MagnetizationMC{T}(N)
        mat = Matrix{T}(undef, 6, 6)
        vec = Vector{T}(undef, 6)
    end
    SPGRBlochSimWorkspace(Aex, Bex, Atr, Btr, Atg, Btg, tmpA1, tmpB1, tmpA2,
                          tmpB2, mat, vec, bm_workspace, ex_workspace)

end

# Case when nTR = 0
# This function does not correctly handle timing with non-instantaneous RF pulses
function (scan::SPGRBlochSim{<:AbstractRF,<:AbstractSpoiling,0})(
    spin::AbstractSpin,
    workspace::SPGRBlochSimWorkspace = SPGRBlochSimWorkspace(spin, scan)
)

    Tg = spoiler_gradient_duration(scan.spoiling)

    excite!(workspace.Aex, workspace.Bex, spin, scan.rf, workspace.ex_workspace)
    freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TR - Tg, workspace.bm_workspace)
    spoil!(workspace.Atg, workspace.Btg, spin, scan.spoiling, workspace.bm_workspace)

    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Atr, workspace.Btr, workspace.Atg, workspace.Btg)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aex, workspace.Bex)
    subtract!(workspace.mat, I, workspace.tmpA2)
    copyto!(workspace.vec, workspace.tmpB2)
    F = lu!(workspace.mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TE, workspace.bm_workspace)
    applydynamics!(spin, workspace.tmpB1, workspace.Atr, workspace.Btr)

end

# Case when nTR > 0
function (scan::SPGRBlochSim{<:AbstractRF,<:AbstractSpoiling,nTR})(
    spin::AbstractSpin,
    workspace::SPGRBlochSimWorkspace = SPGRBlochSimWorkspace(spin, scan)
) where {nTR}



end
