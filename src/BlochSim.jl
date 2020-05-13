module BlochSim

using ForwardDiff
using LinearAlgebra

export GAMMA
export GAMBAR
export AbstractSpin
export Spin
export SpinMC
export freeprecess
export freeprecess!
export excitation
export excitation!
export spoil
export spoil!
export applydynamics!
export mese
export mese!
export spgr!

include("helpers.jl")
include("Spin.jl")
include("sequences.jl")

end
