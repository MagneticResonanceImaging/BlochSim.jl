struct SPGRBlochSim{T1<:AbstractRF,T2<:AbstractSpoiling,nTR,save_transients}
    TR::Float64
    TE::Float64
    rf::T1
    spoiling::T2

    function SPGRBlochSim(
        TR,
        TE,
        rf::T1,
        spoiling::T2,
        nTR::Val{T3},
        save_transients::Val{T4}
    ) where {T1<:AbstractRF,T2<:AbstractSpoiling,T3,T4}

        Tg = spoiler_gradient_duration(spoiling)
        TR >= TE + Tg || error("TR must be greater than or equal to TE + Tg")
        # Should probably also check to make sure TE is not during RF pulse
        (T3 isa Int && T3 >= 0) || error("nTR must be a nonnegative Int")
        T2 <: Union{<:RFSpoiling,<:RFandGradientSpoiling} && (T3 > 0 ||
            error("nTR must be positive when simulating RF spoiling"))
        T4 isa Bool || error("save_transients must be a Bool")
        T3 == 0 && T4 &&
            @warn("save_transients is true, but nTR = 0; no transients will be saved")
        new{T1,T2,T3,T4}(TR, TE, rf, spoiling)

    end
end

SPGRBlochSim(TR, TE, rf, spoiling::AbstractSpoiling, nTR::Val) = SPGRBlochSim(TR, TE, rf, spoiling, nTR, Val(false))
SPGRBlochSim(TR, TE, rf, spoiling::AbstractSpoiling) = SPGRBlochSim(TR, TE, rf, spoiling, Val(0), Val(false))
SPGRBlochSim(TR, TE, rf, nTR::Val, save_transients::Val = Val(false)) = SPGRBlochSim(TR, TE, rf, IdealSpoiling(), nTR, save_transients)
SPGRBlochSim(TR, TE, rf) = SPGRBlochSim(TR, TE, rf, IdealSpoiling(), Val(0), Val(false))
SPGRBlochSim(TR, TE, α::Real, spoiling, nTR, save_transients) = SPGRBlochSim(TR, TE, InstantaneousRF(α), spoiling, nTR, save_transients)

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
    elseif T2 <: RFSpoiling
        Atg = nothing
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
# This function does not correctly handle timing with non-instantaneous RF pulses
function (scan::SPGRBlochSim{<:AbstractRF,T,nTR,save})(
    spin::AbstractSpin,
    workspace::SPGRBlochSimWorkspace = SPGRBlochSimWorkspace(spin, scan)
) where {T,nTR,save}

    rf = scan.rf
    rfspoiling = T <: Union{<:RFSpoiling,<:RFandGradientSpoiling}
    Tg = spoiler_gradient_duration(scan.spoiling)

    if rfspoiling
        rf isa RF && (rf.Δθ[] = rf.Δθ_initial)
        Δθinc = rfspoiling_increment(scan.spoiling)
        θ = zero(Δθinc) # For knowing how much phase to remove when recording signal
        Δθ = Δθinc
    end

    if save
        M = Vector{typeof(spin.M)}(undef, nTR)
        freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TR - scan.TE - Tg, workspace.bm_workspace)
    else
        freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TR - Tg, workspace.bm_workspace)
    end
    spoil!(workspace.Atg, workspace.Btg, spin, scan.spoiling, workspace.bm_workspace)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Atr, workspace.Btr, workspace.Atg, workspace.Btg)
    freeprecess!(workspace.Atr, workspace.Btr, spin, scan.TE, workspace.bm_workspace)

    for rep = 1:nTR-1

        excite!(workspace.Aex, workspace.Bex, spin, rf, workspace.ex_workspace)
        applydynamics!(spin, workspace.tmpB2, workspace.Aex, workspace.Bex)
        if save
            applydynamics!(spin, workspace.tmpB2, workspace.Atr, workspace.Btr)
            M[rep] = copy(spin.M)
            if rfspoiling
                if spin isa Spin
                    tmp = signal(spin.M) * exp(im * θ)
                    M[rep].x = real(tmp)
                    M[rep].y = imag(tmp)
                else
                    for i = 1:spin.N
                        tmp = signal(spin.M[i]) * exp(im * θ)
                        M[rep][i].x = real(tmp)
                        M[rep][i].y = imag(tmp)
                    end
                end
            end
        end
        applydynamics!(spin, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1)

        if rfspoiling
            if rf isa InstantaneousRF
                rf = InstantaneousRF(rf.α, rf.θ + Δθ)
            else
                rf.Δθ[] += Δθ
            end
            θ += Δθ
            Δθ += Δθinc
        end

    end

    excite!(workspace.Aex, workspace.Bex, spin, rf, workspace.ex_workspace)
    applydynamics!(spin, workspace.tmpB2, workspace.Aex, workspace.Bex)
    applydynamics!(spin, workspace.tmpB2, workspace.Atr, workspace.Btr)
    if rfspoiling
        if spin isa Spin
            tmp = signal(spin.M) * exp(im * θ)
            spin.M.x = real(tmp)
            spin.M.y = imag(tmp)
        else
            for i = 1:spin.N
                tmp = signal(spin.M[i]) * exp(im * θ)
                spin.M[i].x = real(tmp)
                spin.M[i].y = imag(tmp)
            end
        end
    end
    if save
        M[nTR] = copy(spin.M)
        return M
    end

end
