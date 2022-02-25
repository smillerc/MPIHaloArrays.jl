# MPIHaloArrays.jl - Work in progress

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaLang.github.io/Example.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLang.github.io/Example.jl/dev) -->

GitHub Actions : [![Build Status](https://github.com/JuliaLang/Example.jl/workflows/CI/badge.svg)](https://github.com/JuliaLang/Example.jl/actions?query=workflow%3ACI+branch%3Amaster)

<!-- AppVeyor: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaLang/Example.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tkelman/example-jl/branch/master) -->

<!-- [![Coverage Status](https://coveralls.io/repos/JuliaLang/Example.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaLang/Example.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLang/Example.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLang/Example.jl?branch=master) -->


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