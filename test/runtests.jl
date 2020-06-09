using BlochSim
using ForwardDiff
using MAT
using Test

@testset "BlochSim.jl" begin
    include("helpers.jl")
    include("expm.jl")
    include("Spin.jl")
    include("sequences.jl")
end
