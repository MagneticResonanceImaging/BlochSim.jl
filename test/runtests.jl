using BlochSim
using LinearAlgebra
using Test

@testset "aqua.jl" begin
    include("aqua.jl")
end

@testset "BlochSim.jl" begin
    include("matlab.jl")
    include("magnetization.jl")
    include("blochmatrix.jl")
    include("spin.jl")
    include("spoiling.jl")
    include("expm.jl")
    include("expm-bloch3.jl")
    include("freeprecess.jl")
    include("excite.jl")
    include("rf-rect.jl")
    include("bssfp1.jl")
    include("bssfp2.jl")
    include("mese.jl")
    include("spgr.jl")
    include("crb.jl")
    include("slice.jl")

    include("ext0.jl") # extension - before "fit.jl"
    include("fit.jl") # extension - should be last
end
