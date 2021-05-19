using BlochSim
using ForwardDiff
using LinearAlgebra
using MAT
using Statistics
using Test

@testset "BlochSim.jl" begin
#    include("blochmatrix.jl")
#    include("AbstractSpin.jl")
#    include("AbstractSpoiling.jl")
#    include("helpers.jl")
#    include("freeprecess.jl")
#    include("excite.jl")
#    include("expm.jl")
#    include("Spin.jl")
#    include("sequences.jl")

    include("matlab.jl")
    include("magnetization.jl")
    include("blochmatrix.jl")
    include("spin.jl")
    include("spoiling.jl")
    include("expm.jl")
end
