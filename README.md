# BlochSim.jl

[![Build Status](https://travis-ci.org/StevenWhitaker/BlochSim.jl.svg?branch=master)](https://travis-ci.org/StevenWhitaker/BlochSim.jl)
[![codecov](https://codecov.io/gh/StevenWhitaker/BlochSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/StevenWhitaker/BlochSim.jl)

This package provides functionality for simulating arbitrary MRI sequences.
It includes support for (traditional) single-compartment Bloch simulations
(using `Spin` objects) as well as multi-compartment Bloch-McConnell simulations
(using `SpinMC` objects).

## Getting started
This package is registered in the
[`General`](https://github.com/JuliaRegistries/General) registry, so you can
install at the REPL with `] add BlochSim`.

## Examples
See the examples given in the documentation strings for how to use the provided
functions. To access the documentation for, e.g., `freeprecess`, simply type
`?freeprecess` at the Julia REPL after loading the package.

For examples of how to simulate full MRI sequences, see
[src/sequences.jl](https://github.com/StevenWhitaker/BlochSim.jl/blob/master/src/sequences.jl).

## Acknowledgement
This package was developed based on
[Brian Hargreaves' Bloch simulation tutorial](http://mrsrl.stanford.edu/~brian/bloch/).
All tests for this package of the form `testX0x` (like `testA5b` or `testF3d`)
are based on the corresponding section in the tutorial.
