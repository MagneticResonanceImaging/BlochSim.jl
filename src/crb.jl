#=
crb.jl
Functions for computing CRB.
Strictly speaking these have little to do
with Bloch simulations,
but one purpose of such simulations
is to develop signal models
and scan designs
for quantitative MRI,
and this package already depends on ForwardDiff,
so we include these helper functions.
=#

export crb

public snr2sigma
public real_imag

import ForwardDiff # jacobian


# helper functions for CRB

"""
    real_imag(y)
Return the stack `[real(y); imag(y)]`
that is needed for CRB
because `ForwardDiff`
cannot handle complex signal values.
"""
real_imag(y) = [real(y); imag(y)] # COV_EXCL_LINE


"""
    snr2sigma(snr_db, s::AbstractArray)
Convert SNR in dB to noise standard deviation `σ`
for a complex-valued signal `s`
to which complex-valued noise will be added.
"""
snr2sigma(db::Real, s::AbstractArray{<:Number}) =
     10^(-db/20) * norm(s) / sqrt(length(s))


"""
    crb(signal, x, σ)
Compute Cramer-Rao Bound (CRB)
for estimating parameter vector `x ∈ ℝᴺ`
from the noisy signal
`y = signal(x) + ϵ ∈ ℝᴹ`
for `M ≥ N`
where `ϵ` denotes additive white Gaussian noise
with standard deviation `σ`.

The output of `signal(x)` must be real!
For complex signals,
use the helper function `Blochsim.real_imag`.
"""
function crb(
    signal, x::AbstractArray{Tx}, σ::Tσ,
) where {Tx <: Number, Tσ <: Number}

    #=
    Type inference seems impossible because the compiler
    does not know the type of the output of `signal`.
    =#
    T = promote_type(Tx, Tσ, Float32)
    jac = ForwardDiff.jacobian(signal, x)              
    fish = jac' * jac / σ^2    
    crb = inv(fish)
    return T.(crb)
end
