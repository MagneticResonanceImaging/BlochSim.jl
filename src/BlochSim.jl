module BlochSim

using LinearAlgebra: I, diagm

include("helpers.jl")
include("Spin.jl")
include("sequences.jl")

export GAMMA, GAMBAR
export AbstractSpin, Spin, SpinMC
export freeprecess, freeprecess!
export excitation, excitation!
export spoil, spoil!
export combine, applydynamics!

export mese!, spgr!

end
