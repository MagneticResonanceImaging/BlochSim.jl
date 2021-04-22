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
export signal

export Gradient
export AbstractSpoiling
export GradientSpoiling
export RFSpoiling
export RFandGradientSpoiling
export spoiler_gradient
export rfspoiling_increment
export InstantaneousRF
export RF

export BlochMatrix
export FreePrecessionMatrix
export ExcitationMatrix
export BlochMcConnellMatrix

export BlochMcConnellWorkspace
export ExcitationWorkspace
export freeprecess
export freeprecess!
export excitation
export excitation!
export excite
export excite!
export spoil
export spoil!
export applydynamics!
export MESEBlochSim
export MESEBlochSimWorkspace
export mese
export mese!
export spgr!

include("magnetization.jl")
include("blochmatrix.jl")
include("spin.jl")
include("spoiling.jl")
include("expm.jl")
include("helpers.jl")
include("freeprecess.jl")
include("excite.jl")
include("mese.jl")

end
