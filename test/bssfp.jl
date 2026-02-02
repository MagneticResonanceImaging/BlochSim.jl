#=
test/bssfp.jl
=#

using BlochSim: bssfp, BSSFPTuple1
using ForwardDiff: ForwardDiff
import ForwardDiff: derivative, gradient
using Test: @inferred, @test, @testset

real_imag(z) = [real(z); imag(z)] # stacker

# single-pool test
@testset "bssfp1" begin
    α_deg, TR_ms, TE_ms = 20, 10, 5 # scan parameters
    Mz0, T1_ms, T2_ms = 0.9, 400, 100 # tissue parameters
    Δf_Hz = -30.5
    xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tissue
    α_rad = deg2rad(α_deg)
    xs = (; TR_ms, TE_ms, α_rad)
    tmp1 = @inferred bssfp(xt..., xs...)
    @test tmp1 isa Complex{<:AbstractFloat}

    tmp2 = @inferred bssfp(xt, xs...)
    @test tmp1 == tmp2

    fun(xt) = real_imag(bssfp(xt..., xs...))
    grad = ForwardDiff.jacobian(fun, collect(xt))
    @test grad isa Matrix{<:AbstractFloat}
end


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
    ΔΦ_deg = 30 # RF phase cycling increment
    ΔΦ_rad = deg2rad(ΔΦ_deg)

    xt = (; M0_phase, Mz0, f_f, T1_f_ms, T1_s_ms, T2_f_ms, T2_s_ms,
        τ_fs_ms, Δf_myelin_Hz, Δf0_Hz)
    xs = (; ΔΦ_rad, TR_ms, TE_ms, α_rad)
    tmp1 = @inferred bssfp(xt..., xs...)

    @test tmp1 isa Complex{<:AbstractFloat}

    fun(xt) = real_imag(bssfp(xt..., xs...))
    tmp2 = @inferred fun(xt)
    @test real_imag(tmp1) == tmp2

    grad = ForwardDiff.jacobian(fun, collect(xt))
    @test grad isa Matrix{<:AbstractFloat}
end
