#=
test/bssfp.jl
=#

using BlochSim: bssfp, BSSFPTuple1
using ForwardDiff: ForwardDiff
import ForwardDiff: derivative, gradient
using Test: @inferred, @test, @testset

real_imag(z) = [real(z); imag(z)] # stacker


"""
bSSFP signal for TE = TR/2, Δϕ = π, Δf=0
see lenz:10:lor
"""
function freeman_hill(T1, T2, TR, TE, α)
    @assert TE ≈ TR/2
    E1 = exp(-TR / T1)
    E2 = exp(-TR / T2)
    return sin(α) * (1 - E1) * sqrt(E2) /
        (1 - E1 * E2 - (E1 - E2) * cos(α))
end


# single-pool test
@testset "bssfp1" begin
    α_deg, TR_ms, TE_ms = 20, 10, 5f0 # scan parameters
    Mz0, T1_ms, T2_ms = 0.9, 400, 100 # tissue parameters
    Δf_Hz = -30.5
    xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tissue
    α_rad = deg2rad(α_deg)
    Δϕ_rad = π/3
Δϕ_rad = 0
    xs = (; TR_ms, TE_ms, Δϕ_rad, α_rad)
    tmp1 = @inferred bssfp(xt..., xs...)
    @test tmp1 isa Complex{<:AbstractFloat}

    tmp2 = @inferred bssfp(xt, xs...)
    @test tmp1 == tmp2

    fun(xt) = real_imag(bssfp(xt..., xs...))
    grad = ForwardDiff.jacobian(fun, collect(xt))
    @test grad isa Matrix{<:AbstractFloat}

    # compare to classic Freeman-Hill formula:
    tmp3 = Mz0 * freeman_hill(T1_ms, T2_ms, TR_ms, TE_ms, α_rad)
    Δϕ_rad = π
    xtf = (; Mz0, T1_ms, T2_ms, Δf_Hz = 0 - Δϕ_rad/(2π*(TR_ms/1000)))
#   xtf = (; Mz0, T1_ms, T2_ms, Δf_Hz = 0) # todo
    rf_phase_rad = π/2
    xsf = (; TR_ms, TE_ms, Δϕ_rad=0*Δϕ_rad, α_rad, rf_phase_rad) # todo
    tmp4 = @inferred bssfp(xtf, xsf...)
    @test tmp3 ≈ tmp4
end

#throw()

# two-pool test
@testset "bssfp2" begin
# todo: currently all the values must be Float for jacobian to work

    M0_phase = π/3
    Mz0 = 0.9
    f_f = 0.15 # fast fraction
    T1_f_ms = 400. # T1 for fast-relaxing, myelin water compartment
    T1_s_ms = 832. # T1 for slow-relaxing, non-myelin water compartment
    T2_f_ms = 20. # T2 for fast-relaxing, myelin water compartment
    T2_s_ms = 80. # T2 for slow-relaxing, non-myelin water compartment
    τ_fs_ms = 50. # residence time
    Δf_myelin_Hz = 5. # frequency shift of myelin water
    Δf0_Hz = -40. # from B0
    α_deg, TR_ms, TE_ms = 20, 10, 5 # scan parameters
    α_rad = deg2rad(α_deg)
    Δϕ_deg = 30 # RF phase cycling increment
    Δϕ_rad = deg2rad(Δϕ_deg)

    xt = (; M0_phase, Mz0, f_f, T1_f_ms, T1_s_ms, T2_f_ms, T2_s_ms,
        τ_fs_ms, Δf_myelin_Hz, Δf0_Hz)
    xs = (; Δϕ_rad, TR_ms, TE_ms, α_rad)
    tmp1 = @inferred bssfp(xt..., xs...)

    @test tmp1 isa Complex{<:AbstractFloat}

    fun(xt) = real_imag(bssfp(xt..., xs...))
    tmp2 = @inferred fun(xt)
    @test real_imag(tmp1) == tmp2

    grad = ForwardDiff.jacobian(fun, collect(xt))
    @test grad isa Matrix{<:AbstractFloat}
end
