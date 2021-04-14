module BlochSim

using ForwardDiff
using LinearAlgebra

export GAMMA
export GAMBAR

export Magnetization
export MagnetizationMC
export Position
export AbstractSpin
export Spin
export SpinMC

export Gradient
export AbstractSpoiling
export GradientSpoiling
export RFSpoiling
export RFandGradientSpoiling
export spoiler_gradient
export rfspoiling_increment

export FreePrecessionMatrix

export BlochMcConnellWorkspace
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

include("spin.jl")
include("spoiling.jl")
include("expm.jl")
include("helpers.jl")
include("freeprecess.jl")
include("sequences.jl")

end
