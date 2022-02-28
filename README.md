# MPIHaloArrays.jl - Work in progress


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smillerc.github.io/MPIHaloArrays.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

A high-level array type to help with halo, or ghost-cell exchanges commonly found in large-scale PDE problems. Very similar in goals and design to `MPIArrays.jl` and `ImplicitGlobalGrid.jl`.

```julia
using MPI, MPIHaloArrays

MPI.Init()

# Create the MPI topology, which here is a Cartesian 4x4 domain (using 16 cores)
topo = CartesianTopology([4,4], [true, true])
rank = MPI.Comm_rank(MPI.COMM_WORLD)

# How many halo cells in each dimension (fixed for all dimensions)
nhalo = 2
N = 200

local_data = rand(N,N) # does not include halo cells
x = MPIHaloArray{Float64}(local_data,topo,nhalo)

updatehalo!(x)

GC.gc()
MPI.Finalize()
```

## Installation

The package is registered and can be installed with `Pkg.add` as

```julia
pkg> add MPIHaloArrays
```

## Documentation

- [[**STABLE**](https://smillerc.github.io/MPIHaloArrays.jl/stable)] &mdash; **most recently tagged version of the documentation.**
- [[**DEV**](https://smillerc.github.io/MPIHaloArrays.jl/dev)] &mdash; **most recent development version of the documentation.**

## Project Status

The package is tested against Julia `v1.0` and the latest `v1.x` on Linux, macOS, and Windows.

## Contributing and Questions