#=
slice.jl
RF and gradient waveforms for slice-selective excitation
=#

#using BlochSim: b1_gauss, GAMMA, Gradient, GradientSpoiling, RF

export rf_slice


"""
    rf_sinc(α_rad, tRF_ms, Δt_ms, shape, nlobe)

Make truncated sinc RF waveform.
"""
function rf_sinc(
    α_rad::Real,
    tRF_ms::Real,
    Δt_ms::Real,
    shape::Symbol,
    nlobe::Real,
)

    shape === :sinc || throw("unsupported shape $shape")
    nsamp = round(Int, tRF_ms / Δt_ms)
    nsamp * Δt_ms ≈ tRF_ms || @warn("Δt_ms=$Δt_ms does not divide tRF_ms=$tRF_ms")
    t = ((0:(nsamp-1))/nsamp .- 0.5) * tRF_ms # [-tRF_ms/2, tRF_ms/2)
    wave = sinc.(2nlobe * t / tRF_ms)
    wave .*= nsamp / sum(wave) # make sum to unity
    return wave * b1_gauss(α_rad, tRF_ms) # flip angle
end


"""
    gz = gz_sinc(tRF_ms, nlobe, slice_width)

Slice-selective gradient amplitude (in G/cm)

From Fourier analysis (small tip-angle approximation):
`2π / (tRF_ms/2nlobe) = GAMMA * gz * slice_width`
"""
gz_sinc(tRF_ms::Real, nlobe::Int, slice_width::Real) =
    (2nlobe * 2π) / (tRF_ms/1000) / GAMMA / slice_width


"""
    rf, rephasing = rf_slice(tRF_ms = 1 ; kwargs...)
    rf = rf_slice( Val(:built_in_rephasing), tRF_ms = 1 ; kwargs...)

Make RF pulse for slice selection
(with constant gradient along z direction),
and rephasing gradient
of length `tRF_ms/2`.

The version with `Val(:built_in_rephasing)`
pads the RF waveform with zeros
and includes the rephasing gradient;
it is just for testing!
It is not recommended because it is less efficient computationally
and its `duration(rf) = (3/2) * tRF_ms`,

Caution:
the bipolar gradient waveform used here
ignores slew-rate constraints.

# In
- `Val(:built_in_rephasing)` (see above)
- `tRF_ms` duration in ms of RF portion of the pulse; default 1 [ms]
  (excluding rephasing gradient)

# Option
- `α_rad` flip angle; default π/2 radian
- `Δt_ms` dwell time (sampling period); default 1e-3 for 1μs
- `shape` pulse shape; default `:sinc`
- `nlobe` number of sinc lobes; default 5
- `slice_width` default 1 cm
- `rephasing::Bool` include rephasing gradient? default `!isinf(slice_width)`
- `Δθ` phase passed to `RF`; default 0 radian

# Out
- `rf` a `BlochSim.RF` object
"""
function rf_slice(
    tRF_ms::Real = 1.0,
    ;
    α_rad::Real = π/2,
    Δt_ms::Real = 1e-3, # ms, so 1μs dwell
    shape::Symbol = :sinc,
    nlobe::Real = 5,
    slice_width::Real = 1, # cm
    rephasing::Bool = !isinf(slice_width),
    Δθ::Real = 0,
)

    wave = rf_sinc(α_rad, tRF_ms, Δt_ms, shape, nlobe) # RF waveform
    gz = gz_sinc(tRF_ms, nlobe, slice_width) # gradient amplitude

    # RF pulse with constant slice-selection gradient:
    rf = RF(wave, Δt_ms, Δθ, Gradient(0, 0, gz))

    rephasing = GradientSpoiling(Gradient(0, 0, -gz), tRF_ms/2)

    return rf, rephasing
end



function rf_slice(
    ::Val{:built_in_rephasing},
    tRF_ms::Real = 1.0,
    ;
    α_rad::Real = π/2,
    Δt_ms::Real = 1e-3, # ms, so 1μs dwell
    shape::Symbol = :sinc,
    nlobe::Real = 5,
    slice_width::Real = 1, # cm
    rephasing::Bool = !isinf(slice_width),
    Δθ::Real = 0,
)

    wave = rf_sinc(α_rad, tRF_ms, Δt_ms, shape, nlobe) # RF waveform
    nsamp = length(wave)
    gz = gz_sinc(tRF_ms, nlobe, slice_width) # gradient amplitude

    grad = [fill(Gradient(0, 0, gz), nsamp); # excite
            fill(Gradient(0, 0, -gz), nsamp÷2)] # rephasing
    wave = [wave; zeros(nsamp÷2)]
    rf = RF(wave, Δt_ms, Δθ, grad)

    return rf
end
