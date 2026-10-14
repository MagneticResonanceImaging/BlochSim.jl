# reuse.jl
# Utilities for avoiding allocations with spins

export reuse_spin, reuse_spins!


"""
    reuse_spin(spin, M0, T1, T2, Δf)

Return new `Spin` struct
with given tissue parameter values,
"recycling" the memory used in `spin`
for its magnetization `M` and position `pos` vectors,
to avoid heap allocations.

Caution.
Using this leads to two spins
sharing the same `M` and `pos` vectors.
"""
function reuse_spin(spin::Spin,
     M0::Number = spin.M0,
     T1::Number = spin.T1,
     T2::Number = spin.T2,
     Δf::Number = spin.Δf,
)
     return Spin(spin.M, M0, T1, T2, Δf, spin.pos)
end


"""
    reuse_spins!(spins, Mz0, T1_ms, T2_ms, Δf_Hz)

Update tissue parameter values
for each spin in `spins`,
"recycling" the memory used
for its magnetization `M` and position `pos`
to avoid heap allocations.

This is useful when using an array of `Spin`s
to model slice profile effects.
"""
function reuse_spins!(spins::AbstractVector{<:Spin}, M0, T1, T2, Δf)
     for (i, spin) in enumerate(spins)
         spins[i] = reuse_spin(spin, M0, T1, T2, Δf)
     end
     return spins
end
