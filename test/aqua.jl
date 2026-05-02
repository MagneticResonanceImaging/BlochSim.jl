using BlochSim: BlochSim
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_ambiguities(BlochSim; broken = true)
    Aqua.test_unbound_args(BlochSim; broken = true)
    Aqua.test_all(BlochSim;
       ambiguities = false,
       unbound_args = false, # todo
       stale_deps = (ignore = [:ADTypes],), # used by Optim extension
    )
end
