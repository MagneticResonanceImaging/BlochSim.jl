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
export IdealSpoiling
export GradientSpoiling
export RFSpoiling
export RFandGradientSpoiling
export idealspoiling
export spoiler_gradient
export spoiler_gradient_duration
export rfspoiling_increment

export AbstractRF
export InstantaneousRF
export RF
export duration

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

export add!
export subtract!
export muladd!
export combine!
export applydynamics!

export MESEBlochSim
export MESEBlochSimWorkspace
export SPGRBlochSim
export SPGRBlochSimWorkspace

include("magnetization.jl")
include("blochmatrix.jl")
include("spin.jl")
include("spoiling.jl")
include("expm.jl")
include("helpers.jl")
include("freeprecess.jl")
include("excite.jl")
include("mese.jl")
include("spgr.jl")

end
