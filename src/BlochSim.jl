module BlochSim

using LinearAlgebra: I, diagm

include("helpers.jl")
include("Spin.jl")

export GAMMA, GAMBAR
export AbstractSpin, Spin, SpinMC
export freeprecess, freeprecess!
export excitation, excitation!
export spoil, spoil!
export combine, applydynamics!

end
