#=
test/bssfp1.jl 1-pool case
=#

using BlochSim: bssfp, Spin, InstantaneousRF, duration, real_imag
using BlochSim: bSSFPbloch, bSSFPbloch3, bSSFPellipse, rf_slice
using ForwardDiff: ForwardDiff
import ForwardDiff: derivative, gradient, jacobian
using Test: @inferred, @test, @testset


"""
bSSFP signal for TE = TR/2, Δϕ = π, Δf = 0, θ_rf = 0
see lenz:10:lor
"""
function freeman_hill(T1, T2, TR, TE, α)
    @assert TE ≈ TR/2
    E1 = exp(-TR / T1)
    E2 = exp(-TR / T2)
    return sin(α) * (1 - E1) * sqrt(E2) /
        (1 - E1 * E2 - (E1 - E2) * cos(α))
end


# old vs new way for Δϕ_rad=0 (the only option the old way supports)
# todo: cut after new way works!
@testset "bssfp0" begin
    Mz0, T1_ms, T2_ms, Δf_Hz = 0.7, 400, 80, -9 # tissue
    spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz)
    TR_ms, TE_ms, α_rad, Δϕ_rad, θ_rf_rad = 10, 5f0, 0, π/5, π/7 # scan
    rf = InstantaneousRF(α_rad, θ_rf_rad)
    oldsig = @inferred bssfp(spin, TR_ms, TE_ms, rf)
    newsig = @inferred bssfp(spin, TR_ms, TE_ms, Δϕ_rad, rf)
    @test oldsig ≈ newsig
end


# single-pool test
@testset "bssfp1" begin
    Mz0, T1_ms, T2_ms, Δf_Hz = 0.7, 400, 80, -9 # tissue parameters
    xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tissue
    TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad = 20, 10, 5f0, π/3, π/5, π/7 # scan
    xs = (; TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad)

    # long argument list with no structs
    sig1 = @inferred bssfp(xt..., xs...)
    @test sig1 isa Complex{<:AbstractFloat}

    sig2 = @inferred bssfp(xt, xs...) # tuple version
    @test sig1 == sig2

    # jacobian thereof
    fun(xt) = real_imag(bssfp(xt..., xs...))
    jac = ForwardDiff.jacobian(fun, collect(xt))
    @test jac isa Matrix{<:AbstractFloat}
    @test !any(isnan, jac)

    # short vs long form for InstantaneousRF
    spin = Spin(Mz0, T1_ms, T2_ms, Δf_Hz)
    rf = InstantaneousRF(α_rad, θ_rf_rad)
    sig3 = @inferred bssfp(spin, TR_ms, TE_ms, Δϕ_rad, rf)
    @test sig1 == sig3

    # ellipse formula
    sig7 = @inferred bssfp(bSSFPellipse, xt..., xs...)
    @test isapprox(sig1, sig7; atol=1e-9)

    sig8 = @inferred bssfp(bSSFPellipse, xt, xs...)
    @test sig7 == sig8

    # helpers
    bpost = bssfp(spin, TR_ms, duration(rf)/2, 0, rf) # signal right after RF
    @test bpost == bssfp(spin, TR_ms, Val(:postRF), 0, rf)
    bmid = bssfp(spin, TR_ms, TR_ms/2, 0, rf)
    @test bmid == bssfp(spin, TR_ms, Val(:midTR), 0, rf)

    # slice-select version
    Δt_ms = 0.1 # coarse for testing
    rf, rephasing = rf_slice(1; α_rad, Δt_ms)
    sig4 = @inferred bssfp(spin, TR_ms, TE_ms, Δϕ_rad, (rephasing, rf, rephasing))

    # jacobian of slice-select version
    funs(xt) = real_imag(bssfp(Spin(xt...),
        TR_ms, TE_ms, Δϕ_rad, (rephasing, rf, rephasing)))
    @test real_imag(sig4) == funs(xt)
    jac = ForwardDiff.jacobian(funs, collect(xt))
    @test jac isa Matrix{<:AbstractFloat}
    @test !any(isnan, jac)

    # test with RF constructed (kappa fixed to 1)
    function _bssfp1(xt)
        spin1 = Spin(xt...)
        rf1, rephasing1 = rf_slice(1; α_rad = 1 * α_rad, Δt_ms) # kappa
        tmp = (rephasing1, rf1, rephasing1)
        return real_imag(bssfp(spin1, TR_ms, TE_ms, Δϕ_rad, tmp))
    end
    @test real_imag(sig4) ≈ _bssfp1(xt)
    jac1 = ForwardDiff.jacobian(_bssfp1, collect(xt))
    @test jac1 isa Matrix{<:AbstractFloat}
    @test !any(isnan, jac1)

    # test jacobian with kappa too
    function _bssfp1(xk)
        spin1 = Spin(xk[1:4]...)
        rf1, rephasing1 = rf_slice(1; α_rad = xk[5] * α_rad, Δt_ms) # kappa
        tmp = (rephasing1, rf1, rephasing1)
        return real_imag(bssfp(spin1, TR_ms, TE_ms, Δϕ_rad, tmp))
    end
    @test real_imag(sig4) ≈ _bssfp1([xt...; 1])
    kappa = 0.9
    xk = [xt...; kappa]
    jac = ForwardDiff.jacobian(_bssfp1, xk)
    @test jac isa Matrix{<:AbstractFloat}
    @test !any(isnan, jac)
end


# compare to classic Freeman-Hill formula:
@testset "freeman-hill" begin
    Mz0, T1_ms, T2_ms, Δf_Hz = 0.7, 400, 80, 0 # note Δf_Hz=0 !
    xt = (; Mz0, T1_ms, T2_ms, Δf_Hz) # tissue
    TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad = 20, 10, π, π/5, 0 # scan
    xs = (; TR_ms, TE_ms, Δϕ_rad, α_rad, θ_rf_rad)

    sig4 = Mz0 * freeman_hill(T1_ms, T2_ms, TR_ms, TE_ms, α_rad)
    sig5 = @inferred bssfp(xt, xs...)
    @test sig4 ≈ sig5

    sig6 = @inferred bssfp(bSSFPbloch, xt, xs...)
    @test sig4 ≈ sig6

    tRF_ms = 1e-3 # super-short RF pulse
    sig7 = @inferred bssfp(bSSFPbloch3, tRF_ms, xt..., xs...)
    @test isapprox(sig4, sig7, atol = 2e-6)
end
