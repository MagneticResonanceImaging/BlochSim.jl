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
    include("freeprecess.jl")
    include("excite.jl")
    include("bssfp.jl")
    include("mese.jl")
    include("spgr.jl")
end
