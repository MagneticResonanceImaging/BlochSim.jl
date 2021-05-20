using BlochSim
using ForwardDiff
using LinearAlgebra
using MAT
using Statistics
using Test

@testset "BlochSim.jl" begin
    include("matlab.jl")
    include("magnetization.jl")
    include("blochmatrix.jl")
    include("spin.jl")
    include("spoiling.jl")
    include("expm.jl")
    include("freeprecess.jl")
    include("excite.jl")
end
