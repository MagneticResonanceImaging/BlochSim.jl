#=
slice.jl
RF and gradient waveforms for slice selection
=#

# using BlochSim: b1_gauss, GAMMA, Gradient, RF

export rf_slice


"""
    rf = rf_slice( ; kwargs...)

Make RF pulse for slice selection
(along z direction),
including refocusing gradient
of length `tRF_ms/2` by default.

Caution:
with that default,
`duration(rf) = (3/2) * tRF_ms`,

Caution:
the bipolar gradient waveform used here
ignores slew-rate constraints.

# Option
- `α_rad` flip angle; default π/2 radian
- `tRF_ms` duration in ms of RF portion of the pulse; default 1 [ms]
  (excluding refocus gradient)
- `Δt_ms` dwell time (sampling period); default 1e-3 for 1μs
- `shape` pulse shape; default `:sinc`
- `nlobe` number of sinc lobes; default 5
- `slice_width` default 1 cm (set to `Inf` for non-selective)
- `refocus::Bool` include refocusing gradient? default `!isinf(slice_width)`
- `Δθ` phase passed to `RF`; default 0 radian

# Out
- `rf` a `BlochSim.RF` object
"""
function rf_slice( ;
    α_rad::Real = π/2,
    Δt_ms::Real = 1e-3, # ms, so 1μs dwell
    tRF_ms::Real = 1.0,
    shape::Symbol = :sinc,
    nlobe::Real = 5,
    slice_width::Real = 1, # cm
    refocus::Bool = !isinf(slice_width),
    Δθ::Real = 0,
)

    # Make truncated sinc RF waveform
    shape === :sinc || throw("unsupported shape $shape")
    nsamp = round(Int, tRF_ms / Δt_ms)
    t = ((0:(nsamp-1))/nsamp .- 0.5) * tRF_ms # [-tRF_ms/2, tRF_ms/2)
    wave = sinc.(2nlobe * t / tRF_ms)
    wave .*= nsamp / sum(wave) # make sum to unity
    wave .*= b1_gauss(α_rad, tRF_ms) # flip angle

    #=
    Slice-selective gradient amplitude
    From Fourier analysis (small tip-angle approximation):
    2π / (tRF_ms/2nlobe) = GAMMA * gz * slice_width
    =#
    gz = (2nlobe * 2π) / (tRF_ms/1000) / GAMMA / slice_width

    if refocus
        grad = [fill(Gradient(0, 0, gz), nsamp); # excite
                fill(Gradient(0, 0, -gz), nsamp÷2)] # refocus
        wave = [wave; zeros(nsamp÷2)]
        rf = RF(wave, Δt_ms, Δθ, grad)

    else # no refocus gradient
        grad = Gradient(0, 0, gz)

        # RF pulse with constant slice-selection gradient:
        rf = RF(wave, Δt_ms, Δθ, grad)
    end
    return rf
end
