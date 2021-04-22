using BlochSim
using ForwardDiff
using LinearAlgebra
using MAT
using Test

@testset "BlochSim.jl" begin
    include("sequences.jl")

    include("blochmatrix.jl")
    include("AbstractSpin.jl")
    include("AbstractSpoiling.jl")
    include("helpers.jl")
    include("freeprecess.jl")
    include("excite.jl")
    include("expm.jl")
    include("Spin.jl")
end
