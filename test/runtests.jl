using BlochSim
using MAT
using Test
@testset "BlochSim.jl" begin
    include("helpers.jl")
    include("Spin.jl")
    include("sequences.jl")
end
