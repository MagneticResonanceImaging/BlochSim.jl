using BlochSim
using ForwardDiff
using LinearAlgebra
using MAT
using Test

@testset "BlochSim.jl" begin
#    include("helpers.jl")
#    include("expm.jl")
#    include("Spin.jl")
#    include("sequences.jl")

    include("blochmatrix.jl")
    include("AbstractSpin.jl")
    include("AbstractSpoiling.jl")
    include("freeprecess.jl")
    include("expm.jl")
end
