"""
    mese! = MESEBlochSim(TR, TE, nechoes, [rfex, rfref, rephaser, crusher, spoiling])
    mese!(spin, [workspace])

Simulate a multi-echo spin echo (MESE) scan on `spin`, overwriting the spin's
magnetization vector. Returns a `Vector` with the magnetization vectors at each
echo.

# Arguments
- `TR::Real`: Repetition time (ms)
- `TE::Real`: First echo time, and echo spacing (ms);
  the first echo time is measured from the middle of the excitation pulse
- `nechoes::Integer`: Number of echoes to readout
- `rfex::AbstractRF = InstantaneousRF(π/2)`: Excitation RF pulse
- `rfref::AbstractRF = InstantaneousRF(π, -π/2)`: Refocussing RF pulse
- `rephaser::Union{<:GradientSpoiling,Nothing} = nothing`: Slice-select
  excitation rephasing gradient
- `crusher::Union{<:GradientSpoiling,Nothing} = nothing`: Crusher gradient
  (placed on either side of each refocussing pulse)
- `spoiling::Union{IdealSpoiling,<:GradientSpoiling,Nothing} = IdealSpoiling()`:
  Type of spoiling to apply

`workspace isa MESEBlochSimWorkspace`.
"""
struct MESEBlochSim{T1<:AbstractRF,T2<:AbstractRF,
        T3<:Union{<:GradientSpoiling,Nothing},
        T4<:Union{<:GradientSpoiling,Nothing},
        T5<:Union{IdealSpoiling,<:GradientSpoiling,Nothing}}
    TR::Float64
    TE::Float64
    nechoes::Int
    rfex::T1
    rfref::T2
    rephaser::T3
    crusher::T4
    spoiling::T5

    # Constructor ensures sequence timing works out
    # Specifically:
    # 1. The TR must be long enough to collect all the echoes
    #    and to include spoiling at the end
    # 2. The first echo must occur after the excitation pulse, prephasing
    #    gradient, and the refocussing pulse (with its flanking crusher
    #    gradients), and the first refocussing pulse must be at TE/2
    # 3. The echo spacing must be greater than the duration of the
    #    refocussing pulse and its flanking crusher gradients
    # As written, 3. is covered by 2., but this changes if
    # a readout gradient is added to the simulation.
    function MESEBlochSim(TR, TE, nechoes, rfex::T1, rfref::T2, rephaser::T3,
            crusher::T4, spoiling::T5) where {T1,T2,T3,T4,T5}

        dur = x -> isnothing(x) ? 0.0 : spoiler_gradient_duration(x)

        TR >= TE * nechoes + dur(spoiling) + duration(rfex) / 2 ||
            error("TR must be long enough to collect all echoes and to include spoiling")
        TE / 2 >= duration(rfex) / 2 + dur(rephaser) + duration(rfref) / 2 + dur(crusher) ||
            error("first refocussing pulse must occur at TE / 2")
        new{T1,T2,T3,T4,T5}(TR, TE, nechoes, rfex, rfref, rephaser, crusher, spoiling)

    end
end

MESEBlochSim(TR, TE, nechoes) = MESEBlochSim(TR, TE, nechoes, InstantaneousRF(π/2), InstantaneousRF(π, -π/2), nothing, nothing, IdealSpoiling())

Base.show(io::IO, mese::MESEBlochSim) =
    print(io, "MESEBlochSim(", mese.TR, ", ", mese.TE, ", ", mese.nechoes, ", ", mese.rfex, ", ", mese.rfref, ", ", mese.rephaser, ", ", mese.crusher, ", ", mese.spoiling, ")")

function Base.show(io::IO, ::MIME"text/plain", mese::MESEBlochSim)

    print(io, "Multi-Echo Spin Echo (MESE) Bloch Simulation:")
    print(io, "\n TR = ", mese.TR, " ms")
    print(io, "\n TE (and echo spacing) = ", mese.TE, " ms")
    print(io, "\n nechoes = ", mese.nechoes)
    print(io, "\n rfex (excitation pulse) = ")
    show(io, "text/plain", mese.rfex)
    print(io, "\n rfref (refocussing pulses) = ")
    show(io, "text/plain", mese.rfref)
    print(io, "\n rephaser (after excitation pulse) = ")
    show(io, "text/plain", mese.rephaser)
    print(io, "\n crusher = ")
    show(io, "text/plain", mese.crusher)
    print(io, "\n spoiling = ")
    show(io, "text/plain", mese.spoiling)

end

struct MESEBlochSimWorkspace{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21}
    Aex::T1
    Bex::T2
    Aref::T3
    Bref::T4
    Are::T5
    Bre::T6
    Acrush::T7
    Bcrush::T8
    As::T9
    Bs::T10
    Ate1::T11
    Bte1::T12
    Ate::T11
    Bte::T12
    Atr::T11
    Btr::T12
    Aecho1::T13
    Becho1::T12
    Aecho::T13
    Becho::T12
    tmpA1::T13
    tmpB1::T12
    tmpA2::T13
    tmpB2::T12
    mat::T14
    vec::T15
    bm_workspace::T16
    ex_workspace::T17
    ref_workspace::T18
    re_workspace::T19
    crush_workspace::T20
    s_workspace::T21
end

function MESEBlochSimWorkspace(
    spin::AbstractSpin,
    scan::MESEBlochSim,
    bm_workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    MESEBlochSimWorkspace(typeof(spin), typeof(scan), bm_workspace)

end

function MESEBlochSimWorkspace(
    spin::Union{Type{Spin{T}},Type{SpinMC{T,N}}},
    scan::Type{MESEBlochSim{T1,T2,T3,T4,T5}},
    bm_workspace = spin <: Spin ? nothing : BlochMcConnellWorkspace(spin)
) where {T,N,T1,T2,T3,T4,T5}

    if T1 <: InstantaneousRF
        Aex = ExcitationMatrix{T}()
        Bex = nothing
        ex_workspace = nothing
    else
        if spin <: Spin
            Aex = BlochMatrix{T}()
            Bex = Magnetization{T}()
        else
            Aex = BlochMcConnellMatrix{T}(N)
            Bex = MagnetizationMC{T}(N)
        end
        ex_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if T2 <: InstantaneousRF
        Aref = ExcitationMatrix{T}()
        Bref = nothing
        ref_workspace = nothing
    else
        if spin <: Spin
            Aref = BlochMatrix{T}()
            Bref = Magnetization{T}()
        else
            Aref = BlochMcConnellMatrix{T}(N)
            Bref = MagnetizationMC{T}(N)
        end
        ref_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if T3 <: GradientSpoiling
        if spin <: Spin
            Are = FreePrecessionMatrix{T}()
            Bre = Magnetization{T}()
        else
            Are = BlochMcConnellMatrix{T}(N)
            Bre = MagnetizationMC{T}(N)
        end
        if T3 <: GradientSpoiling{<:Gradient}
            re_workspace = bm_workspace
        else
            re_workspace = FreePrecessionWorkspace(spin, bm_workspace)
        end
    else
        Are = nothing
        Bre = nothing
        re_workspace = nothing
    end
    if T4 <: GradientSpoiling
        if spin <: Spin
            Acrush = FreePrecessionMatrix{T}()
            Bcrush = Magnetization{T}()
        else
            Acrush = BlochMcConnellMatrix{T}(N)
            Bcrush = MagnetizationMC{T}(N)
        end
        if T4 <: GradientSpoiling{<:Gradient}
            crush_workspace = bm_workspace
        else
            crush_workspace = FreePrecessionWorkspace(spin, bm_workspace)
        end
    else
        Acrush = nothing
        Bcrush = nothing
        crush_workspace = nothing
    end
    if T5 <: IdealSpoiling
        As = idealspoiling
        Bs = nothing
        s_workspace = nothing
    elseif T5 <: GradientSpoiling
        if spin <: Spin
            As = FreePrecessionMatrix{T}()
            Bs = Magnetization{T}()
        else
            As = BlochMcConnellMatrix{T}(N)
            Bs = MagnetizationMC{T}(N)
        end
        if T5 <: GradientSpoiling{<:Gradient}
            s_workspace = bm_workspace
        else
            s_workspace = FreePrecessionWorkspace(spin, bm_workspace)
        end
    else
        As = nothing
        Bs = nothing
        s_workspace = nothing
    end
    if spin <: Spin
        Ate1 = FreePrecessionMatrix{T}()
        Bte1 = Magnetization{T}()
        Ate = FreePrecessionMatrix{T}()
        Bte = Magnetization{T}()
        Atr = FreePrecessionMatrix{T}()
        Btr = Magnetization{T}()
        Aecho1 = BlochMatrix{T}()
        Becho1 = Magnetization{T}()
        Aecho = BlochMatrix{T}()
        Becho = Magnetization{T}()
        tmpA1 = BlochMatrix{T}()
        tmpB1 = Magnetization{T}()
        tmpA2 = BlochMatrix{T}()
        tmpB2 = Magnetization{T}()
        mat = Matrix{T}(undef, 3, 3)
        vec = Vector{T}(undef, 3)
    else
        Ate1 = BlochMcConnellMatrix{T}(N)
        Bte1 = MagnetizationMC{T}(N)
        Ate = BlochMcConnellMatrix{T}(N)
        Bte = MagnetizationMC{T}(N)
        Atr = BlochMcConnellMatrix{T}(N)
        Btr = MagnetizationMC{T}(N)
        Aecho1 = BlochMcConnellMatrix{T}(N)
        Becho1 = MagnetizationMC{T}(N)
        Aecho = BlochMcConnellMatrix{T}(N)
        Becho = MagnetizationMC{T}(N)
        tmpA1 = BlochMcConnellMatrix{T}(N)
        tmpB1 = MagnetizationMC{T}(N)
        tmpA2 = BlochMcConnellMatrix{T}(N)
        tmpB2 = MagnetizationMC{T}(N)
        mat = Matrix{T}(undef, 3N, 3N)
        vec = Vector{T}(undef, 3N)
    end
    MESEBlochSimWorkspace(Aex, Bex, Aref, Bref, Are, Bre, Acrush, Bcrush, As,
        Bs, Ate1, Bte1, Ate, Bte, Atr, Btr, Aecho1, Becho1, Aecho, Becho, tmpA1,
        tmpB1, tmpA2, tmpB2, mat, vec, bm_workspace, ex_workspace,
        ref_workspace, re_workspace, crush_workspace, s_workspace)

end

function (scan::MESEBlochSim)(spin::AbstractSpin, workspace::MESEBlochSimWorkspace = MESEBlochSimWorkspace(spin, scan))

    dur = x -> isnothing(x) ? 0.0 : spoiler_gradient_duration(x)

    # Excitation pulse
    excite!(workspace.Aex, workspace.Bex, spin, scan.rfex, workspace.ex_workspace)

    # Rephasing gradient
    isnothing(scan.rephaser) || spoil!(workspace.Are, workspace.Bre, spin, scan.rephaser, workspace.re_workspace)

    # Time between rephasing gradient and crusher gradient
    # Compute the time such that the middle of the first refocussing pulse
    # occurs at scan.TE / 2
    t = scan.TE / 2 - duration(scan.rfref) / 2 - dur(scan.crusher) - dur(scan.rephaser) - duration(scan.rfex) / 2
    freeprecess!(workspace.Ate1, workspace.Bte1, spin, t, workspace.bm_workspace)

    # Crusher gradient
    isnothing(scan.crusher) || spoil!(workspace.Acrush, workspace.Bcrush, spin, scan.crusher, workspace.crush_workspace)

    # Refocussing pulse
    excite!(workspace.Aref, workspace.Bref, spin, scan.rfref, workspace.ref_workspace)

    # Time between crusher gradient and spin echo
    # Compute the time such that the spin echo occurs scan.TE / 2
    # after the center of the refocussing pulse
    t = scan.TE / 2 - duration(scan.rfref) / 2 - dur(scan.crusher)
    freeprecess!(workspace.Ate, workspace.Bte, spin, t, workspace.bm_workspace)

    # Time after final echo and before spoiling
    t = scan.TR - duration(scan.rfex) / 2 - dur(scan.spoiling) - scan.TE * scan.nechoes
    freeprecess!(workspace.Atr, workspace.Btr, spin, t, workspace.bm_workspace)

    # Spoiling
    isnothing(scan.spoiling) || spoil!(workspace.As, workspace.Bs, spin, scan.spoiling, workspace.s_workspace)

    # Combine dynamics that occur between echoes:
    # TE -> crusher -> refocus -> crusher -> TE
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Ate, workspace.Bte, workspace.Acrush, workspace.Bcrush)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aref, workspace.Bref)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, workspace.Acrush, workspace.Bcrush)
    combine!(workspace.Aecho, workspace.Becho, workspace.tmpA1, workspace.tmpB1, workspace.Ate, workspace.Bte)

    # Combine dynamics that occur after excitation until the first echo:
    # rephase -> wait -> crusher -> refocus -> crusher -> TE
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Are, workspace.Bre, workspace.Ate1, workspace.Bte1)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Acrush, workspace.Bcrush)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, workspace.Aref, workspace.Bref)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Acrush, workspace.Bcrush)
    combine!(workspace.Aecho1, workspace.Becho1, workspace.tmpA2, workspace.tmpB2, workspace.Ate, workspace.Bte)

    # Combine dynamics of the whole TR
    # Have excitation last so steady-state gives magnetization
    # immediately following the excitation pulse
    copyto!(workspace.tmpA1, workspace.Aecho1)
    copyto!(workspace.tmpB1, workspace.Becho1)
    for e = 2:scan.nechoes
        combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aecho, workspace.Becho)
        copyto!(workspace.tmpA1, workspace.tmpA2)
        copyto!(workspace.tmpB1, workspace.tmpB2)
    end
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atr, workspace.Btr)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, workspace.As, workspace.Bs)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Aex, workspace.Bex)

    # Compute steady-state magnetization
    subtract!(workspace.mat, I, workspace.tmpA2)
    copyto!(workspace.vec, workspace.tmpB2)
    F = lu!(workspace.mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    # Collect the multi-echo spin echo data
    Mout = Vector{typeof(spin.M)}(undef, scan.nechoes)
    applydynamics!(spin, workspace.tmpB1, workspace.Aecho1, workspace.Becho1)
    Mout[1] = copy(spin.M)
    for e = 2:scan.nechoes
        applydynamics!(spin, workspace.tmpB1, workspace.Aecho, workspace.Becho)
        Mout[e] = copy(spin.M)
    end

    return Mout

end
