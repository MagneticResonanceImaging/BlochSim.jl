# BlochSim.jl

https://github.com/StevenWhitaker/BlochSim.jl

[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![action status][action-img]][action-url]
[![pkgeval status][pkgeval-img]][pkgeval-url]
[![codecov][codecov-img]][codecov-url]
[![license][license-img]][license-url]

This Julia package provides functionality
for simulating arbitrary MRI pulse sequences.
It includes support for (traditional) single-compartment Bloch simulations
(using `Spin` objects) as well as multi-compartment Bloch-McConnell simulations
(using `SpinMC` objects).

## Getting started
This package is registered in the
[General](https://github.com/JuliaRegistries/General) registry, so you can
install it at the REPL with `] add BlochSim`.

The main functionality is provided by the functions
`freeprecess`, `excite`, and `spoil`
(and their mutating variants `freeprecess!`, `excite!`, and `spoil!`).
These functions can be used to simulate a wide variety of MRI sequences.
In addition, this package provides implementations for
a multi-echo spin echo (MESE) scan (`MESEBlochSim`)
and a spoiled gradient-recalled echo (SPGR) scan (`SPGRBlochSim`).

## Examples
See the examples given in the documentation strings for how to use the provided
functions. To access the documentation for, e.g., `freeprecess`, simply type
`?freeprecess` at the Julia REPL after loading the package.

For examples of how to simulate full MRI sequences, see
[src/mese.jl](src/mese.jl) and [src/spgr.jl](src/spgr.jl) in this repo,
and [STFR.jl](https://github.com/StevenWhitaker/STFR.jl).

Additionally,
the following blog posts
cover the basics of the Bloch equations
and some examples of how to use this package.
- [Simulating MRI Physics with the Bloch Equations](https://blog.glcs.io/simulating-mri-physics-with-the-bloch-equations)
- [Mastering MRI Bloch Simulations with BlochSim.jl in Julia](https://blog.glcs.io/mastering-mri-bloch-simulations-with-blochsimjl-in-julia)
- [Simulating MRI Scans with BlochSim.jl](https://blog.glcs.io/simulating-mri-scans-with-blochsimjl)

Below are some more concrete examples of how to use this package.

```julia
julia> using BlochSim

julia> spin = Spin(1, 1000, 100, 3.75)
Spin{Float64}:
 M = Magnetization(0.0, 0.0, 1.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 3.75 Hz
 pos = Position(0.0, 0.0, 0.0) cm

julia> excite!(spin, InstantaneousRF(π/2)) # 90° excitation

julia> spin.M # Mz is not quite 0 due to numerical roundoff
Magnetization vector with eltype Float64:
 Mx = 1.0
 My = 0.0
 Mz = 6.123233995736766e-17

julia> freeprecess!(spin, 100) # Free-precess for 100 ms

julia> spin.M
Magnetization vector with eltype Float64:
 Mx = -0.2601300475114444
 My = -0.2601300475114445
 Mz = 0.09516258196404054

julia> spgr! = SPGRBlochSim(5, 2.5, deg2rad(20)) # Create an object to simulate an SPGR scan
Spoiled Gradient-Recalled Echo (SPGR) Bloch Simulation:
 TR = 5.0 ms
 TE = 2.5 ms
 rf (excitation pulse) = Instantaneous RF pulse with eltype Float64:
 α = 0.3490658503988659 rad
 θ = 0.0 rad
 spoiling = IdealSpoiling()
 steady-state

julia> spgr!(spin) # Simulate a steady-state SPGR scan applied to the given spin

julia> spin.M # Steady-state magnetization
Magnetization vector with eltype Float64:
 Mx = 0.025553542433162182
 My = -0.0015069712547712193
 Mz = 0.07442699373678281

julia> spinmc = SpinMC(1, (0.2, 0.8), (400, 1000), (20, 80), (15, 0), (100, 25))
SpinMC{Float64,2}:
 M = MagnetizationMC((0.0, 0.0, 0.2), (0.0, 0.0, 0.8))
 M0 = 1.0
 frac = (0.2, 0.8)
 T1 = (400.0, 1000.0) ms
 T2 = (20.0, 80.0) ms
 Δf = (15.0, 0.0) Hz
 r = ((0.0, 0.01), (0.04, 0.0)) 1/ms
 pos = Position(0.0, 0.0, 0.0) cm

julia> spgr!(spinmc) # The same SPGR scan can be used on multi-compartment spins

julia> spinmc.M # Steady-state magnetization
2-compartment Magnetization vector with eltype Float64:
 Compartment 1:
  Mx = -0.09359002635156467
  My = 0.02433674787041617
  Mz = -0.36973998540693054
 Compartment 2:
  Mx = 0.1541252837882581
  My = 0.00031515000730316224
  Mz = 0.5077167235922019

julia> signal(spin) # Grab the observed signal from the spin
0.025553542433162182 - 0.0015069712547712193im

julia> signal(spinmc)
0.060535257436693427 + 0.02465189787771933im
```

## Related package(s)

* [MRIgeneralizedBloch.jl](https://github.com/JakobAsslaender/MRIgeneralizedBloch.jl)
  focuses on [magnetization transfer](https://doi.org/10.1002/mrm.29071)


## Acknowledgement
This package was developed based on
[Brian Hargreaves' Bloch simulation tutorial](http://mrsrl.stanford.edu/~brian/bloch/).
All tests for this package of the form `testX0x` (like `testA5b` or `testF3d`)
are based on the corresponding section in the tutorial (see
[test/matlab.jl](test/matlab.jl)).

<!-- URLs -->
[action-img]: https://github.com/StevenWhitaker/BlochSim.jl/actions/workflows/runtests.yml/badge.svg
[action-url]: https://github.com/StevenWhitaker/BlochSim.jl/actions
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BlochSim.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BlochSim.html
[codecov-img]: https://codecov.io/gh/StevenWhitaker/BlochSim.jl/branch/main/graph/badge.svg?token=tduieBgema
[codecov-url]: https://codecov.io/gh/StevenWhitaker/BlochSim.jl
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://StevenWhitaker.github.io/BlochSim.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://StevenWhitaker.github.io/BlochSim.jl/dev
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE
